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


def extend_reference_multi(bamfile, reference, read_ids_file, output_fasta, trim,
                           chr_arm_pairs_file=None):
    """
    Extend reference genome using soft-clipped bases from multiple reads

    Args:
        bamfile: Path to BAM file
        reference: Path to reference FASTA file
        read_ids_file: Path to file with read IDs (one per line)
        output_fasta: Path to output extended FASTA file
        trim: Number of bp to trim from the end of extensions to remove adapters
        chr_arm_pairs_file: "<chr_end>\\t<read_id>" pairs from select_reads. When
            supplied, each scaffold contributes ONLY to its own arm's side:
            an L-arm read's start-clip may become a prefix, an R-arm read's
            end-clip may become a suffix, and the opposite-side clip of each is
            ignored entirely.

    WHY ARM BINDING MATTERS
    ----------------------
    Each contig carries both arms, so extensions[contig] holds one 'prefix' and
    one 'suffix' list. Without arm awareness a clip is filed purely by which end
    of the CIGAR it sits on, and the winner is whichever is LONGEST. Arm and
    clip-side normally coincide geometrically (verified: all 32 arms map L->
    contig start, R-> contig end), but nothing enforces it, so a read clipped on
    BOTH sides lands in both lists and can win the other arm's slot.

    That is what corrupted 7172 chr16R: the chr16L scaffold is a fold-back whose
    spurious 82,350bp end-clip beat chr16R's legitimate 6,814bp clip, grafting
    chr16L's inverted subtelomere onto chr16R.

    Removing bad reads alone does not close this. Perfectly healthy R-arm reads
    carry a telomeric start-clip (1,848bp on 1R, 2,049bp on 9R) because the base
    reference is truncated before the telomeres -- larger than the legitimate
    prefix clips of chr7L (1,935bp) and chr15L (1,778bp). Binding by arm makes
    the wrong-side clip ineligible rather than merely out-competed.
    """
    # Read the read IDs
    with open(read_ids_file, 'r') as f:
        read_ids = [line.strip() for line in f if line.strip()]

    # read_id -> chr_end ("16L" / "16R"), used to bind clips to their own arm
    read_to_arm = {}
    if chr_arm_pairs_file:
        with open(chr_arm_pairs_file, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    read_to_arm[parts[1]] = parts[0]
        print(f"Arm binding ENABLED: {len(read_to_arm)} scaffold/arm pairs loaded")
    else:
        print("WARNING: no chr_arm pairs file supplied -- falling back to "
              "side-of-CIGAR filing (a read clipped on both sides can win the "
              "wrong arm's slot)")

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

        # Which side is this read allowed to contribute to?
        # L arm -> contig start -> prefix only.  R arm -> contig end -> suffix only.
        # The opposite-side clip is never examined, so it cannot compete.
        arm = read_to_arm.get(read.query_name)
        if read_to_arm and arm is None:
            print(f"  WARNING: {read.query_name} has no arm assignment, skipping")
            continue
        allow_prefix = (arm is None) or arm.upper().endswith('L')
        allow_suffix = (arm is None) or arm.upper().endswith('R')

        # Soft-clipped bases at start (5' extension) -- L arms only
        if allow_prefix:
            prefix = get_softclip_from_cigar(read, side="start")
            if prefix and len(prefix) > trim:
                trimmed_prefix = prefix[:-trim] if trim > 0 else prefix
                extensions[ref_name]['prefix'].append((read.query_name, trimmed_prefix, len(prefix)))

        # Soft-clipped bases at end (3' extension) -- R arms only
        if allow_suffix:
            suffix = get_softclip_from_cigar(read, side="end")
            if suffix and len(suffix) > trim:
                trimmed_suffix = suffix[trim:] if trim > 0 else suffix
                extensions[ref_name]['suffix'].append((read.query_name, trimmed_suffix, len(suffix)))

    samfile.close()

    # Create extended reference(s)
    extended_records = []
    total_extensions = 0

    for ref_name, ext_data in extensions.items():
        orig_ref = str(ref_seqs[ref_name].seq)

        # With arm binding there is exactly one scaffold per arm, so each side
        # must have at most one candidate. More than one means something
        # upstream is broken (duplicate rows in the pairs file, or two reads
        # assigned the same arm) -- fail loudly rather than silently picking.
        if read_to_arm:
            for side in ('prefix', 'suffix'):
                if len(ext_data[side]) > 1:
                    who = ', '.join(f'{r}({len(s)}bp)' for r, s, _ in ext_data[side])
                    raise RuntimeError(
                        f'{ref_name}: {len(ext_data[side])} {side} candidates but arm '
                        f'binding permits exactly 1 -- {who}. Check '
                        f'selected_read_chr_arm_pairs.txt for duplicate arms.')

        # Use the longest prefix and suffix (after trimming) if multiple exist
        prefix = ""
        prefix_read = ""
        prefix_orig_len = 0
        if ext_data['prefix']:
            prefix_read, prefix, prefix_orig_len = max(ext_data['prefix'], key=lambda x: len(x[1]))
            print(f"Using 5' extension for {ref_name} from {prefix_read} "
                  f"[{read_to_arm.get(prefix_read, 'unbound')}]: added {len(prefix)}bp")

        suffix = ""
        suffix_read = ""
        suffix_orig_len = 0
        if ext_data['suffix']:
            suffix_read, suffix, suffix_orig_len = max(ext_data['suffix'], key=lambda x: len(x[1]))
            print(f"Using 3' extension for {ref_name} from {suffix_read} "
                  f"[{read_to_arm.get(suffix_read, 'unbound')}]: added {len(suffix)}bp")

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
