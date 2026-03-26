"""
Step 11c: Spacer recombination analysis (hybrid approach).

For reads where the minimap2 breakpoint falls in a spacer feature, uses two
complementary strategies:
  1. Whole-tail BLAST: BLAST the diverged tail against the spacer library to
     identify the source chr end (weighted 0.7).
  2. Chunk analysis: break the read's spacer region into overlapping 250 bp
     windows; BLAST each chunk to locate the exact switch point within the
     spacer (weighted 0.3).

Usage:
  python analyze_spacer_recombination.py \
      --breakpoints-tsv  results/{base}/recombination/{base}_{chr_end}_breakpoints.tsv \
      --diverged-fasta   results/{base}/recombination/{base}_{chr_end}_diverged_tails.fasta \
      --reads-fasta      results/{base}/chr_anchor_included_individual_files/{base}_blasted_{anchor}_{chr_end}_anchor_reads.fasta \
      --day0-ref         results/{day0_base}/assembly_{strain}/assembly_{strain}_dorado_reference.fasta \
      --spacer-lib-dir   references/pairings_for_spacers/{strain}_pairings/ \
      --chr-end  chr4R  --strain 6991 \
      --output-tsv  results/{base}/recombination/{base}_{chr_end}_spacer_recomb.tsv \
      --threads 4
"""

import argparse
import glob as glob_module
import os
import re
import sys
import tempfile

import pandas as pd

from recombination_utils import (
    Hypothesis,
    hypotheses_to_row_dict,
    get_first_breakpoint_feature_type,
    get_first_breakpoint_is_mid_element,
    get_first_breakpoint_pos_on_read,
    load_bed_features,
    normalize_hypotheses,
    read_fasta,
    run_blast,
    write_fasta_list,
    write_results_tsv,
    MIN_CONFIDENCE_TO_REPORT,
    SPACER_CHUNK_SIZE,
    SPACER_CHUNK_STEP,
)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description='Spacer recombination analysis')
    p.add_argument('--breakpoints-tsv', required=True)
    p.add_argument('--diverged-fasta',  required=True)
    p.add_argument('--reads-fasta',     required=True)
    p.add_argument('--day0-ref',        required=True)
    p.add_argument('--spacer-lib-dir',  required=True)
    p.add_argument('--chr-end',         required=True)
    p.add_argument('--strain',          required=True)
    p.add_argument('--output-tsv',      required=True)
    p.add_argument('--threads',         type=int, default=4)
    return p.parse_args()

# ---------------------------------------------------------------------------
# Build combined spacer library
# ---------------------------------------------------------------------------

def build_combined_spacer_fasta(spacer_lib_dir, strain, tmp_dir):
    pattern = os.path.join(spacer_lib_dir, f'{strain}_paired_*.fasta')
    fasta_files = sorted(glob_module.glob(pattern))
    if not fasta_files:
        return None
    combined_path = os.path.join(tmp_dir, f'{strain}_spacer_combined.fasta')
    with open(combined_path, 'w') as out_fh:
        for fasta_path in fasta_files:
            with open(fasta_path) as in_fh:
                for line in in_fh:
                    out_fh.write(line)
    return combined_path

# ---------------------------------------------------------------------------
# Source chr end from spacer subject name
# ---------------------------------------------------------------------------

def subject_to_chr_end(sseqid):
    m = re.search(r'(chr\d+[LR])', sseqid)
    return m.group(1) if m else sseqid

# ---------------------------------------------------------------------------
# Chunk analysis
# ---------------------------------------------------------------------------

def chunk_sequence(seq, chunk_size=SPACER_CHUNK_SIZE, step=SPACER_CHUNK_STEP):
    """Return list of (chunk_seq, start_pos) with overlapping windows."""
    chunks = []
    start = 0
    while start + chunk_size <= len(seq):
        chunks.append((seq[start:start + chunk_size], start))
        start += step
    # Include a final partial chunk if there's leftover sequence >= 50 bp
    if start < len(seq) and len(seq) - start >= 50:
        chunks.append((seq[start:], start))
    return chunks


def blast_chunks(chunks, spacer_db_fasta, tmp_dir, read_id):
    """
    Write chunks to a temp FASTA (headers = chunk_{start}) and BLAST against spacer db.
    Returns DataFrame with chunk_start column added.
    """
    if not chunks:
        return pd.DataFrame()

    chunk_fasta = os.path.join(tmp_dir, f'{read_id}_chunks.fasta')
    entries = [(f'chunk_{start}', seq) for seq, start in chunks]
    write_fasta_list(entries, chunk_fasta)

    df = run_blast(chunk_fasta, spacer_db_fasta, tmp_dir, label=f'{read_id}_chunk')
    if df.empty:
        return df

    def parse_start(qseqid):
        try:
            return int(qseqid.split('_')[1])
        except (IndexError, ValueError):
            return -1

    df['chunk_start'] = df['qseqid'].apply(parse_start)
    df['source_chr_end'] = df['sseqid'].apply(subject_to_chr_end)
    return df


def find_chunk_switch_point(chunk_blast_df):
    """
    Given per-chunk BLAST results, find where the source chr end changes.

    Returns (switch_pos, source_before, source_after, chunk_confidence).
    switch_pos=-1 if no switch found.
    """
    if chunk_blast_df.empty:
        return -1, '', '', 0.0

    # For each chunk position, pick the top-bitscore source chr end
    best = (
        chunk_blast_df.sort_values('bitscore', ascending=False)
        .groupby('chunk_start')
        .first()
        .reset_index()[['chunk_start', 'source_chr_end', 'bitscore']]
    )
    best = best.sort_values('chunk_start')

    if len(best) < 2:
        return -1, best.iloc[0]['source_chr_end'] if len(best) else '', '', 1.0

    # Find the first position where source chr end changes
    switch_pos = -1
    source_before = best.iloc[0]['source_chr_end']
    source_after = ''
    for i in range(1, len(best)):
        if best.iloc[i]['source_chr_end'] != source_before:
            switch_pos = int(best.iloc[i]['chunk_start'])
            source_after = best.iloc[i]['source_chr_end']
            break

    if switch_pos < 0:
        return -1, source_before, '', 1.0

    # Confidence = fraction of chunks agreeing with the two-segment assignment
    n_before = sum(1 for _, row in best.iterrows()
                   if row['chunk_start'] < switch_pos and row['source_chr_end'] == source_before)
    n_after = sum(1 for _, row in best.iterrows()
                  if row['chunk_start'] >= switch_pos and row['source_chr_end'] == source_after)
    n_total = len(best)
    chunk_confidence = (n_before + n_after) / n_total if n_total > 0 else 0.0

    return switch_pos, source_before, source_after, chunk_confidence

# ---------------------------------------------------------------------------
# Tail BLAST confidence
# ---------------------------------------------------------------------------

def tail_blast_confidence(blast_hits, diverged_seq_len):
    """
    Compute source chr end and confidence from BLAST of the full diverged tail.
    Returns (source_chr_end, c_tail).
    """
    if blast_hits.empty:
        return '', 0.0

    blast_hits = blast_hits.copy()
    blast_hits['source_chr_end'] = blast_hits['sseqid'].apply(subject_to_chr_end)
    best_by_source = (
        blast_hits.groupby('source_chr_end')['bitscore'].max().reset_index()
        .sort_values('bitscore', ascending=False)
    )
    top = best_by_source.iloc[0]
    top_source = top['source_chr_end']
    best_hit = blast_hits[blast_hits['bitscore'] == top['bitscore']].iloc[0]

    c_tail = (best_hit['pident'] / 100.0) * min(
        1.0, best_hit['length'] / max(diverged_seq_len, 1)
    ) ** 0.5
    return top_source, c_tail

# ---------------------------------------------------------------------------
# Hypothesis generation
# ---------------------------------------------------------------------------

def generate_spacer_hypotheses(read_id, bp_row, chunk_blast_df, tail_blast_df, diverged_seq_len):
    """Generate Hypothesis objects combining tail BLAST and chunk analysis."""
    hypotheses = []
    is_mid = get_first_breakpoint_is_mid_element(bp_row)

    # --- Tail stage ---
    tail_source, c_tail = tail_blast_confidence(tail_blast_df, diverged_seq_len)

    # --- Chunk stage ---
    switch_pos, source_before, source_after, chunk_conf = find_chunk_switch_point(chunk_blast_df)

    # If both stages have no result
    if not tail_source and switch_pos < 0:
        hypotheses.append(Hypothesis(
            h_type='ambiguous',
            description='Spacer breakpoint but no match in spacer library',
            confidence=0.2,
            ambiguous=True,
        ))
        return normalize_hypotheses(hypotheses)

    # Determine agreed source
    if tail_source and source_after and tail_source == source_after:
        # Both agree
        c_combined = c_tail * 0.7 + chunk_conf * 0.3
        ambiguous = chunk_conf < 0.6
        desc = (
            f"Spacer switch to {tail_source}"
            + (f" at ~{switch_pos} bp into spacer" if switch_pos >= 0 else '')
        )
        if is_mid:
            c_combined *= 0.7
        hypotheses.append(Hypothesis(
            h_type='spacer',
            description=desc,
            confidence=c_combined,
            source_chr_ends=[tail_source],
            ambiguous=ambiguous,
        ))

        # Compound: check for second switch in chunk analysis
        if not chunk_blast_df.empty:
            best = (
                chunk_blast_df.sort_values('bitscore', ascending=False)
                .groupby('chunk_start').first().reset_index()
                .sort_values('chunk_start')
            )
            switches = []
            prev_src = best.iloc[0]['source_chr_end'] if len(best) else ''
            for _, row in best.iloc[1:].iterrows():
                if row['source_chr_end'] != prev_src:
                    switches.append((int(row['chunk_start']), row['source_chr_end']))
                    prev_src = row['source_chr_end']
            if len(switches) >= 2:
                c1 = c_combined
                c2 = chunk_conf * 0.3
                c_compound = c1 * c2
                hypotheses.append(Hypothesis(
                    h_type='compound',
                    description=(
                        f"Compound spacer: switch to {switches[0][1]} at ~{switches[0][0]} bp, "
                        f"then to {switches[1][1]} at ~{switches[1][0]} bp"
                    ),
                    confidence=c_compound,
                    is_compound=True,
                    ambiguous=True,
                ))

    elif tail_source and source_after and tail_source != source_after:
        # Disagreement: generate separate hypotheses
        c_tail_h = c_tail * 0.7
        if is_mid:
            c_tail_h *= 0.7
        hypotheses.append(Hypothesis(
            h_type='spacer',
            description=f"Spacer switch to {tail_source} (tail BLAST evidence)",
            confidence=c_tail_h,
            source_chr_ends=[tail_source],
            ambiguous=True,
        ))
        c_chunk_h = chunk_conf * 0.3
        if is_mid:
            c_chunk_h *= 0.7
        if c_chunk_h >= MIN_CONFIDENCE_TO_REPORT:
            hypotheses.append(Hypothesis(
                h_type='spacer',
                description=(
                    f"Spacer switch to {source_after} at ~{switch_pos} bp "
                    f"(chunk analysis evidence; disagrees with tail BLAST)"
                ),
                confidence=c_chunk_h,
                source_chr_ends=[source_after],
                ambiguous=True,
            ))

    elif tail_source and switch_pos < 0:
        # Only tail BLAST result
        c = c_tail * 0.7
        if is_mid:
            c *= 0.7
        hypotheses.append(Hypothesis(
            h_type='spacer',
            description=f"Spacer switch to {tail_source} (tail BLAST only; no chunk switch point found)",
            confidence=c,
            source_chr_ends=[tail_source],
            ambiguous=True,
        ))

    elif switch_pos >= 0 and not tail_source:
        # Only chunk result
        c = chunk_conf * 0.3
        if is_mid:
            c *= 0.7
        hypotheses.append(Hypothesis(
            h_type='spacer',
            description=(
                f"Spacer switch to {source_after} at ~{switch_pos} bp "
                f"(chunk analysis only; no tail BLAST match)"
            ),
            confidence=c,
            source_chr_ends=[source_after],
            ambiguous=True,
        ))

    return normalize_hypotheses(hypotheses)

# ---------------------------------------------------------------------------
# Extract the spacer region from a read
# ---------------------------------------------------------------------------

def extract_spacer_region_from_read(read_seq, bp_row, features):
    """
    Extract the portion of the read that corresponds to the spacer region.
    Uses the breakpoint position and telo_side to determine which part of
    the read is the spacer.  Returns the subsequence (may be the full read
    if precise extraction is not possible).
    """
    bp_on_read = get_first_breakpoint_pos_on_read(bp_row)
    telo_side = bp_row.get('telo_side', 'end')

    if bp_on_read < 0:
        return read_seq

    if telo_side == 'end':
        # Spacer is anchor-proximal; breakpoint is near start of diverged region
        # Use a window around the breakpoint: 2000 bp before + 2000 bp after
        start = max(0, bp_on_read - 2000)
        end = min(len(read_seq), bp_on_read + 2000)
    else:
        start = max(0, bp_on_read - 2000)
        end = min(len(read_seq), bp_on_read + 2000)

    region = read_seq[start:end]
    return region if len(region) >= SPACER_CHUNK_SIZE else read_seq

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    print(f'analyze_spacer_recombination.py — chr_end={args.chr_end}')

    if not os.path.exists(args.breakpoints_tsv):
        print(f'ERROR: missing breakpoints TSV: {args.breakpoints_tsv}', file=sys.stderr)
        sys.exit(1)

    df_bp = pd.read_csv(args.breakpoints_tsv, sep='\t')
    if df_bp.empty:
        write_results_tsv([], args.output_tsv)
        return

    df_bp['_first_bp_feat_type'] = df_bp['breakpoint_feature_types'].astype(str).apply(lambda v: v.split(';')[0])
    spacer_reads = df_bp[df_bp['_first_bp_feat_type'] == 'spacer']
    print(f'  {len(spacer_reads)} reads with spacer breakpoints (of {len(df_bp)} total)')

    # Load full read sequences for chunk extraction
    read_seqs = {}
    if os.path.exists(args.reads_fasta):
        read_seqs = read_fasta(args.reads_fasta)

    # Load diverged tails
    diverged_seqs = {}
    if os.path.exists(args.diverged_fasta):
        diverged_seqs = read_fasta(args.diverged_fasta)

    # Load BED features for spacer boundary information
    features = load_bed_features(args.day0_bed if hasattr(args, 'day0_bed') else '', args.chr_end) \
        if hasattr(args, 'day0_bed') else []

    rows = []

    if spacer_reads.empty or not os.path.exists(args.spacer_lib_dir):
        for _, bp_row in df_bp.iterrows():
            row = {'read_id': bp_row['read_id'], 'chr_end': args.chr_end,
                   'spacer_analyzed': False, 'spacer_switch_pos': -1,
                   'spacer_source_before': '', 'spacer_source_after': '',
                   'n_hypotheses': 0}
            row.update(hypotheses_to_row_dict([]))
            rows.append(row)
        write_results_tsv(rows, args.output_tsv)
        return

    with tempfile.TemporaryDirectory() as tmp_dir:
        combined_lib = build_combined_spacer_fasta(args.spacer_lib_dir, args.strain, tmp_dir)
        if combined_lib is None:
            print(f'  WARNING: no spacer FASTA files found in {args.spacer_lib_dir}')
            for _, bp_row in df_bp.iterrows():
                row = {'read_id': bp_row['read_id'], 'chr_end': args.chr_end,
                       'spacer_analyzed': False, 'spacer_switch_pos': -1,
                       'spacer_source_before': '', 'spacer_source_after': '',
                       'n_hypotheses': 0}
                row.update(hypotheses_to_row_dict([]))
                rows.append(row)
            write_results_tsv(rows, args.output_tsv)
            return

        for _, bp_row in df_bp.iterrows():
            read_id = bp_row['read_id']
            feat_type = get_first_breakpoint_feature_type(bp_row)

            if feat_type != 'spacer':
                row = {'read_id': read_id, 'chr_end': args.chr_end,
                       'spacer_analyzed': False, 'spacer_switch_pos': -1,
                       'spacer_source_before': '', 'spacer_source_after': '',
                       'n_hypotheses': 0}
                row.update(hypotheses_to_row_dict([]))
                rows.append(row)
                continue

            # Find full read sequence
            full_seq = ''
            for header in read_seqs:
                if header.startswith(read_id) or read_id in header:
                    full_seq = read_seqs[header]
                    break

            # Find diverged tail
            diverged_seq = ''
            diverged_len = 0
            for header in diverged_seqs:
                if header.startswith(read_id) or read_id in header:
                    diverged_seq = diverged_seqs[header]
                    diverged_len = len(diverged_seq)
                    break

            # --- Tail BLAST ---
            tail_blast_df = pd.DataFrame()
            if diverged_seq:
                tail_fasta = os.path.join(tmp_dir, f'{read_id}_tail.fasta')
                write_fasta_list([(read_id, diverged_seq)], tail_fasta)
                tail_blast_df = run_blast(tail_fasta, combined_lib, tmp_dir, label=f'{read_id}_tail')

            # --- Chunk BLAST ---
            chunk_blast_df = pd.DataFrame()
            if full_seq:
                spacer_region = extract_spacer_region_from_read(full_seq, bp_row, features)
                chunks = chunk_sequence(spacer_region)
                if chunks:
                    chunk_blast_df = blast_chunks(chunks, combined_lib, tmp_dir, read_id)

            # Generate hypotheses
            hypotheses = generate_spacer_hypotheses(
                read_id, bp_row, chunk_blast_df, tail_blast_df, diverged_len
            )

            # Extract switch info for row metadata
            switch_pos, source_before, source_after, _ = find_chunk_switch_point(chunk_blast_df)

            row = {
                'read_id': read_id,
                'chr_end': args.chr_end,
                'spacer_analyzed': True,
                'spacer_switch_pos': switch_pos,
                'spacer_source_before': source_before,
                'spacer_source_after': source_after,
                'n_hypotheses': len(hypotheses),
            }
            row.update(hypotheses_to_row_dict(hypotheses))
            rows.append(row)

    os.makedirs(os.path.dirname(args.output_tsv) or '.', exist_ok=True)
    write_results_tsv(rows, args.output_tsv)
    print(f'  Output: {args.output_tsv} ({len(rows)} reads)')


if __name__ == '__main__':
    main()
