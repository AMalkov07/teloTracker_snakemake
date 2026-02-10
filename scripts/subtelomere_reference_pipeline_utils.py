"""
Utility functions for telomere extension and polishing pipeline.

Import this with: from subtelomere_reference_pipeline_utils import *
"""

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import os
import pandas as pd


def get_75th_percentile_reads(input_tsv, output_dir, output_file):
    """
    Read TSV file and get the 75th percentile read for each chr_end based on repeat_length.

    Parameters:
    input_tsv (str): Path to the input TSV file
    output_dir (str): Path to output directory
    output_file (str): Path to output file for read IDs

    Returns:
    tuple: (dict mapping chr_end to read_id, list of all selected read_ids)
    """
    print("Reading input TSV file and selecting 75th percentile reads...")

    # Read the TSV
    df = pd.read_csv(input_tsv, sep='\t')

    # Filter for telomere reads
    df_telomere_reads = df[(df['repeat_length'] >= 30) & (df['Adapter_After_Telomere'] == True)]

    # Define all chromosome ends
    chr_ends = [f'{chr_num}{side}' for chr_num in range(1, 17) for side in ['L', 'R']]

    chr_end_to_read = {}
    all_selected_read_ids = []

    for chr_end in chr_ends:
        if chr_end not in df_telomere_reads['chr_end'].values:
            print(f"Warning: {chr_end} not found in the input file.")
            continue

        df_chr_telomere_reads = df_telomere_reads[df_telomere_reads['chr_end'] == chr_end]

        if df_chr_telomere_reads.empty:
            print(f"Warning: No telomere reads found for {chr_end}.")
            continue

        # Get mode y_prime_probe_count
        mode_y_prime = df_chr_telomere_reads['y_prime_probe_count'].mode()[0]

        # Filter for mode y_prime
        df_mode_y_chr_telomere_reads = df_chr_telomere_reads[
            df_chr_telomere_reads['y_prime_probe_count'] == mode_y_prime]

        if df_mode_y_chr_telomere_reads.empty:
            print(f"Warning: No reads found with mode y_prime_probe_count for {chr_end}.")
            continue

        # Get the read ID for the read with the 75% quantile repeat_length
        sorted_df = df_mode_y_chr_telomere_reads.sort_values(by='repeat_length')
        size_of_df = len(sorted_df)
        quantile_index = int(size_of_df * 0.75)
        if quantile_index >= size_of_df:
            quantile_index = size_of_df - 1

        quantile_read_id = sorted_df.iloc[quantile_index]['read_id']
        chr_end_to_read[chr_end] = quantile_read_id
        all_selected_read_ids.append(quantile_read_id)

        print(f"{chr_end}: Selected read {quantile_read_id} (repeat_length: {sorted_df.iloc[quantile_index]['repeat_length']})")

    print(f"\nSelected {len(chr_end_to_read)} reads for {len(chr_end_to_read)} chromosome ends")

    # Save read IDs to file
    with open(output_file, 'w') as f:
        for chr_end, read_id in chr_end_to_read.items():
            f.write(f"{chr_end}\t{read_id}\n")

    # Also save just the read IDs for easy extraction
    read_ids_only_file = output_file.replace('.txt', '_only.txt')
    with open(read_ids_only_file, 'w') as f:
        for read_id in all_selected_read_ids:
            f.write(f"{read_id}\n")

    return chr_end_to_read, all_selected_read_ids


def extract_selected_reads(reads_fastq, read_ids_file, output_fastq):
    """
    Extract selected reads from FASTQ file.

    Parameters:
    reads_fastq (str): Path to input FASTQ file
    read_ids_file (str): Path to file with read IDs (one per line)
    output_fastq (str): Path to output FASTQ file
    """
    print(f"Extracting selected reads from {reads_fastq}")

    # Read the read IDs
    with open(read_ids_file, 'r') as f:
        selected_read_ids = set(line.strip() for line in f if line.strip())

    # Index the FASTQ reads by ID
    print("Indexing FASTQ file...")
    fastq_index = SeqIO.to_dict(SeqIO.parse(reads_fastq, "fastq"))

    # Extract selected reads
    matched_reads = [fastq_index[read_id] for read_id in selected_read_ids if read_id in fastq_index]

    if not matched_reads:
        raise ValueError("No FASTQ reads matched the selected read IDs!")

    print(f"Writing {len(matched_reads)} selected reads to {output_fastq}")
    with open(output_fastq, "w") as out_f:
        SeqIO.write(matched_reads, out_f, "fastq")


def get_softclip_from_cigar(read, side="start"):
    """Extract soft-clipped sequence from read based on CIGAR string"""
    if not read.cigartuples:
        return ""

    if side == "start":
        if read.cigartuples[0][0] == 4:  # 4 = soft clip
            clip_length = read.cigartuples[0][1]
            return read.query_sequence[:clip_length]
    elif side == "end":
        if read.cigartuples[-1][0] == 4:  # 4 = soft clip
            clip_length = read.cigartuples[-1][1]
            return read.query_sequence[-clip_length:]

    return ""


def extend_reference_multi(bamfile, reference, read_ids_file, output_fasta, trim):
    """
    Extend reference genome using soft-clipped bases from multiple reads

    Args:
        bamfile: Path to BAM file
        reference: Path to reference FASTA file
        read_ids_file: Path to file with read IDs (one per line)
        output_fasta: Path to output extended FASTA file
        trim: Number of bp to trim from the end of extensions to remove adapters
    """
    # Read the read IDs
    with open(read_ids_file, 'r') as f:
        read_ids = [line.strip() for line in f if line.strip()]

    samfile = pysam.AlignmentFile(bamfile, "rb")
    ref_seqs = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))

    read_ids_set = set(read_ids)
    extensions = {}  # Store extensions per reference

    print(f"\nSearching for soft-clipped bases from {len(read_ids)} reads...")
    if trim != 0:
        print(f"Will trim {trim}bp from extension ends to remove adapters")

    # Collect all extensions from all specified reads
    for read in samfile.fetch(until_eof=True):
        if read.query_name not in read_ids_set or read.is_secondary or read.is_supplementary:
            continue

        ref_name = read.reference_name

        if ref_name not in extensions:
            extensions[ref_name] = {'prefix': [], 'suffix': []}

        # Check if read has soft-clipped bases at start (5' extension)
        prefix = get_softclip_from_cigar(read, side="start")
        if prefix and len(prefix) > trim:
            trimmed_prefix = prefix[:-trim] if trim > 0 else prefix
            extensions[ref_name]['prefix'].append((read.query_name, trimmed_prefix, len(prefix)))
            #print(f"Found 5' extension from read {read.query_name} on {ref_name}: {len(prefix)}bp (trimmed to {len(trimmed_prefix)}bp)")

        # Check if read has soft-clipped bases at end (3' extension)
        suffix = get_softclip_from_cigar(read, side="end")
        if suffix and len(suffix) > trim:
            trimmed_suffix = suffix[trim:] if trim > 0 else suffix
            extensions[ref_name]['suffix'].append((read.query_name, trimmed_suffix, len(suffix)))
            #print(f"Found 3' extension from read {read.query_name} on {ref_name}: {len(suffix)}bp (trimmed to {len(trimmed_suffix)}bp)")

    samfile.close()

    # Create extended reference(s)
    extended_records = []
    total_extensions = 0

    for ref_name, ext_data in extensions.items():
        orig_ref = str(ref_seqs[ref_name].seq)

        # Use the longest prefix and suffix (after trimming) if multiple exist
        prefix = ""
        prefix_read = ""
        prefix_orig_len = 0
        if ext_data['prefix']:
            prefix_read, prefix, prefix_orig_len = max(ext_data['prefix'], key=lambda x: len(x[1]))
            print(f"Using longest 5' extension for {ref_name} from {prefix_read}: added {len(prefix)}bp")

        suffix = ""
        suffix_read = ""
        suffix_orig_len = 0
        if ext_data['suffix']:
            suffix_read, suffix, suffix_orig_len = max(ext_data['suffix'], key=lambda x: len(x[1]))
            print(f"Using longest 3' extension for {ref_name} from {suffix_read}: added {len(suffix)}bp")

        if prefix or suffix:
            new_ref = prefix + orig_ref + suffix

            # Create proper header
            new_id = f"{ref_name}_extended"
            description = ""
            if prefix:
                description += f" (added {prefix_orig_len}bp)"
            if suffix:
                description += f" (added {suffix_orig_len}bp)"

            new_record = SeqRecord(
                Seq(new_ref),
                id=new_id,
                description=description
            )

            extended_records.append(new_record)
            total_extensions += 1

    if extended_records:
        with open(output_fasta, "w") as out_fasta:
            SeqIO.write(extended_records, out_fasta, "fasta")

        print(f"\n{total_extensions} reference(s) extended and written to {output_fasta}")
        return True
    else:
        print(f"\nNo soft-clipped bases found near chromosome ends for any of the specified reads")
        # Copy original reference to output so pipeline can continue
        subprocess.run(["cp", reference, output_fasta], check=True)
        print(f"Copied original reference to {output_fasta}")
        return False
