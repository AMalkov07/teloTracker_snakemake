"""
Step 10: Run minimap2 alignment for one chr end against the day0 reference.

Produces supplementary evidence only — does NOT classify reads or gate
downstream analysis. The feature-by-feature chunk analysis in step 11
is the primary recombination detection method.

Outputs:
  {output_bam}  — sorted, indexed BAM (for IGV inspection)
  {output_tsv}  — per-read supplementary alignment info

Usage:
  python run_alignment.py \
      --reads-fasta  <per-chr-end reads FASTA> \
      --day0-ref     <day0 reference FASTA> \
      --chr-end      chr4R \
      --output-tsv   results/{base}/recombination/{base}_{chr_end}_alignment.tsv \
      --output-bam   results/{base}/recombination/{base}_{chr_end}.bam \
      --threads      4
"""

import argparse
import os
import subprocess
import sys

# Ensure scripts/ is on the import path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pysam

from recombination_utils import (
    telo_side_from_header,
    write_results_tsv,
)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description='Align reads to day0 reference (supplementary evidence)')
    p.add_argument('--reads-fasta', required=True)
    p.add_argument('--day0-ref',    required=True)
    p.add_argument('--chr-end',     required=True)
    p.add_argument('--output-tsv',  required=True)
    p.add_argument('--output-bam',  required=True)
    p.add_argument('--threads',     type=int, default=4)
    return p.parse_args()

# ---------------------------------------------------------------------------
# minimap2
# ---------------------------------------------------------------------------

def run_minimap2(reads_fasta, day0_ref, output_bam, threads):
    """Align reads to the day0 reference; produce sorted, indexed BAM."""
    print(f'  Running minimap2 ({threads} threads)...')

    minimap2_cmd = [
        'minimap2', '-a', '-x', 'map-ont',
        '--secondary=no',
        '-t', str(threads),
        day0_ref, reads_fasta,
    ]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-o', output_bam]

    with subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as mm:
        with subprocess.Popen(sort_cmd, stdin=mm.stdout, stderr=subprocess.PIPE) as st:
            mm.stdout.close()
            st.communicate()
            mm.wait()
            if mm.returncode != 0:
                print(f'  minimap2 failed: {mm.stderr.read().decode()}', file=sys.stderr)
                sys.exit(1)
            if st.returncode != 0:
                print(f'  samtools sort failed', file=sys.stderr)
                sys.exit(1)

    subprocess.run(['samtools', 'index', output_bam], check=True, capture_output=True)
    print(f'  BAM written to {output_bam}')

# ---------------------------------------------------------------------------
# Read FASTA headers (for telo_side)
# ---------------------------------------------------------------------------

def read_fasta_headers(fasta_path):
    """Return dict {read_id: full_header} from a FASTA file."""
    headers = {}
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith('>'):
                full_header = line[1:].rstrip()
                read_id = full_header.split()[0]
                headers[read_id] = full_header
    return headers

# ---------------------------------------------------------------------------
# BAM parsing — extract supplementary alignment info
# ---------------------------------------------------------------------------

def extract_supplementary_info(bam_path, chr_end, read_headers):
    """
    For each read, find supplementary alignments and record the contig(s)
    they map to. This is supplementary evidence for cross-feature reconciliation.

    Returns list of row dicts.
    """
    results = []

    bam = pysam.AlignmentFile(bam_path, 'rb')
    aln_by_read = {}
    for aln in bam.fetch(until_eof=True):
        aln_by_read.setdefault(aln.query_name, []).append(aln)
    bam.close()

    for read_id, alns in aln_by_read.items():
        full_header = read_headers.get(read_id, read_id)
        telo_side = telo_side_from_header(full_header, chr_end)

        primary = None
        supplementaries = []
        for aln in alns:
            if aln.is_secondary or aln.is_unmapped:
                continue
            if aln.is_supplementary:
                supplementaries.append(aln)
            else:
                primary = aln

        # Primary alignment info
        primary_contig = ''
        primary_mapq = 0
        if primary is not None:
            primary_contig = primary.reference_name or ''
            primary_mapq = primary.mapping_quality

        # Supplementary contigs (deduplicated, ordered by alignment length)
        supp_contigs = []
        for s in sorted(supplementaries,
                        key=lambda a: a.query_alignment_length or 0,
                        reverse=True):
            contig = s.reference_name or ''
            if contig and contig not in supp_contigs and contig != primary_contig:
                supp_contigs.append(contig)

        results.append({
            'read_id': read_id,
            'chr_end': chr_end,
            'telo_side': telo_side,
            'primary_contig': primary_contig,
            'primary_mapq': primary_mapq,
            'supplementary_contigs': ';'.join(supp_contigs),
            'n_supplementary': len(supp_contigs),
        })

    return results

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    print(f'run_alignment.py -- chr_end={args.chr_end}')

    for path in [args.reads_fasta, args.day0_ref]:
        if not os.path.exists(path):
            print(f'ERROR: missing input: {path}', file=sys.stderr)
            sys.exit(1)

    # Read FASTA headers for telo_side determination
    print('  Reading FASTA headers...')
    read_headers = read_fasta_headers(args.reads_fasta)
    n_reads = len(read_headers)
    print(f'  {n_reads} reads')

    if n_reads == 0:
        print('  No reads -- writing empty output')
        write_results_tsv([], args.output_tsv)
        return

    # Run minimap2
    os.makedirs(os.path.dirname(args.output_bam), exist_ok=True)
    run_minimap2(args.reads_fasta, args.day0_ref, args.output_bam, args.threads)

    # Extract supplementary alignment info
    print('  Extracting supplementary alignment info...')
    rows = extract_supplementary_info(args.output_bam, args.chr_end, read_headers)

    # Write output
    os.makedirs(os.path.dirname(args.output_tsv) or '.', exist_ok=True)
    write_results_tsv(rows, args.output_tsv)

    # Summary
    n_with_supp = sum(1 for r in rows if r['n_supplementary'] > 0)
    print(f'  Output: {args.output_tsv} ({len(rows)} reads, {n_with_supp} with supplementary alignments)')


if __name__ == '__main__':
    main()
