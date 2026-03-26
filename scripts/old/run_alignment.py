"""
Step 10: Run minimap2 alignment for one chr end against the day0 reference.

Classifies each read into one of four categories:
  - no_recombination:  high identity + coverage, no supplementary alignments
  - two_clean_halves:  two (or more) chained segments tiling the read
  - one_clean_half:    anchor-proximal segment is clean, but the rest is unaligned or low quality
  - no_clean_halves:   no clean alignment covering the anchor-proximal region

Produces:
  {output_tsv}  — per-read classification and segment/breakpoint TSV
  {output_bam}  — sorted, indexed BAM
  Companion diverged-tails FASTA written alongside the TSV as
      <tsv_stem>_diverged_tails.fasta

Usage:
  python run_alignment.py \
      --reads-fasta  results/{base}/chr_anchor_included_individual_files/{base}_blasted_{anchor}_{chr_end}_anchor_reads.fasta \
      --day0-ref     results/{day0_base}/assembly_{strain}/assembly_{strain}_dorado_reference.fasta \
      --day0-bed     results/{day0_base}/pretelomeric_labels/pretelomeric_regions_{strain}_simp.bed \
      --chr-end      chr4R \
      --output-tsv   results/{base}/recombination/{base}_{chr_end}_breakpoints.tsv \
      --output-bam   results/{base}/recombination/{base}_{chr_end}.bam \
      --threads      4
"""

import argparse
import os
import subprocess
import sys

import pysam

from recombination_utils import (
    classify_breakpoint,
    load_bed_features,
    telo_side_from_header,
    write_results_tsv,
    write_fasta_list,
    MIN_WHOLE_READ_COVERAGE,
    MIN_WHOLE_READ_IDENTITY,
    MIN_SEGMENT_IDENTITY,
    MAX_INTER_SEGMENT_GAP,
    MIN_MAPQ,
    MIN_DIVERGED_LEN,
)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description='Align reads to day0 reference and classify recombination')
    p.add_argument('--reads-fasta', required=True)
    p.add_argument('--day0-ref',    required=True)
    p.add_argument('--day0-bed',    required=True)
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
    print(f'  Running minimap2 ({threads} threads)…')

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
# Read FASTA sequences
# ---------------------------------------------------------------------------

def read_fasta_sequences(fasta_path):
    """Return dict {read_id: sequence} and dict {read_id: full_header}."""
    seqs = {}
    headers = {}
    current_id = None
    current_seq = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(current_seq)
                full_header = line[1:]
                current_id = full_header.split()[0]
                headers[current_id] = full_header
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = ''.join(current_seq)
    return seqs, headers

# ---------------------------------------------------------------------------
# CIGAR parsing
# ---------------------------------------------------------------------------

def cigar_aligned_and_identity(cigar_tuples, nm_tag):
    """Return (aligned_bases, identity) from pysam cigar_tuples and NM tag."""
    aligned_bases = sum(length for op, length in cigar_tuples if op in (0, 7, 8))
    if aligned_bases == 0:
        return 0, 0.0
    mismatches = nm_tag if nm_tag is not None else 0
    identity = max(0.0, 1.0 - mismatches / aligned_bases)
    return aligned_bases, identity

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def chr_end_matches_contig(chr_end, contig_name):
    """Check if a reference contig matches the expected chr end.
    chr9R → expect 'chr9' in the contig name (e.g., 'chr9_extended')."""
    chr_num = chr_end.rstrip('LR')
    return chr_num in contig_name


def _ref_pos_at_query_boundary(seg, side):
    """Get the reference coordinate at a query boundary.

    side='left'  → ref coordinate corresponding to query_start
    side='right' → ref coordinate corresponding to query_end

    Forward alignment: query_start→ref_start, query_end→ref_end
    Reverse alignment: query_start→ref_end,   query_end→ref_start
    """
    if side == 'left':
        return seg['ref_start'] if not seg['is_reverse'] else seg['ref_end']
    else:  # 'right'
        return seg['ref_end'] if not seg['is_reverse'] else seg['ref_start']


def _remove_contained_segments(segments):
    """Remove segments that are entirely contained within another segment."""
    if len(segments) <= 1:
        return segments
    # Sort by query_start ascending, then query_end descending.
    # This ensures that for same-start segments, the longest comes first.
    sorted_segs = sorted(segments, key=lambda s: (s['query_start'], -s['query_end']))
    result = [sorted_segs[0]]
    for seg in sorted_segs[1:]:
        if seg['query_end'] <= result[-1]['query_end']:
            continue  # entirely within the previous segment
        result.append(seg)
    return result

# ---------------------------------------------------------------------------
# BAM parsing → classification
# ---------------------------------------------------------------------------

def parse_bam_to_breakpoints(bam_path, chr_end, features, read_seqs, read_headers):
    """
    Parse BAM alignments, classify each read into one of four categories.
    Returns (list_of_row_dicts, list_of_(header, diverged_seq)_tuples).
    """
    results = []
    diverged_tails = []

    bam = pysam.AlignmentFile(bam_path, 'rb')
    aln_by_read = {}
    for aln in bam.fetch(until_eof=True):
        aln_by_read.setdefault(aln.query_name, []).append(aln)
    bam.close()

    for read_id, alns in aln_by_read.items():
        raw_seq = read_seqs.get(read_id, '')
        read_length = len(raw_seq)
        full_header = read_headers.get(read_id, read_id)
        telo_side = telo_side_from_header(full_header, chr_end)

        qc_flags = []
        if read_length < 500:
            qc_flags.append('too_short')

        # --- Collect non-secondary, mapped alignments ---
        segments_raw = []
        for aln in alns:
            if aln.is_secondary or aln.is_unmapped:
                continue
            ct = aln.cigartuples
            if ct is None:
                continue
            nm = aln.get_tag('NM') if aln.has_tag('NM') else 0
            aligned_bases, identity = cigar_aligned_and_identity(ct, nm)
            coverage = aligned_bases / read_length if read_length > 0 else 0.0

            segments_raw.append({
                'query_start': aln.query_alignment_start,
                'query_end': aln.query_alignment_end,
                'ref_name': aln.reference_name or '',
                'ref_start': aln.reference_start,
                'ref_end': aln.reference_end or 0,
                'identity': round(identity, 4),
                'coverage': round(coverage, 4),
                'aligned_bases': aligned_bases,
                'mapq': aln.mapping_quality,
                'is_reverse': aln.is_reverse,
            })

        # --- No alignments → no_clean_halves ---
        if not segments_raw:
            qc_flags.append('no_alignment')
            results.append(_make_row(
                read_id, chr_end, read_length, telo_side,
                alignment_class='no_clean_halves',
                chained=[], breakpoints=[], diverged_seq_length=0,
                qc_flags=qc_flags,
            ))
            continue

        # --- Sort by query_start, remove contained segments, then chain ---
        segments_raw.sort(key=lambda s: s['query_start'])
        cleaned = _remove_contained_segments(segments_raw)

        # Chain adjacent segments (left to right on the read).
        # Cross-contig segments get a generous overlap allowance because
        # minimap2 extends both alignments through similar transition zones.
        # Same-contig segments use strict gap checking.
        MIN_CHAIN_SEGMENT_SIZE = 500  # bp — ignore tiny supplementary fragments
        chained = [cleaned[0]]
        for i in range(1, len(cleaned)):
            seg = cleaned[i]
            if seg['identity'] < MIN_SEGMENT_IDENTITY:
                break
            seg_span = seg['query_end'] - seg['query_start']
            if seg_span < MIN_CHAIN_SEGMENT_SIZE:
                continue  # skip tiny fragments

            gap = seg['query_start'] - chained[-1]['query_end']

            if seg['ref_name'] != chained[-1]['ref_name']:
                # Cross-contig: allow large overlap (common at recombination junctions)
                # but limit gaps to MAX_INTER_SEGMENT_GAP
                if gap <= MAX_INTER_SEGMENT_GAP:
                    chained.append(seg)
                else:
                    break
            else:
                # Same contig: strict gap/overlap check
                if abs(gap) <= MAX_INTER_SEGMENT_GAP:
                    chained.append(seg)
                else:
                    break

        # --- Identify the anchor segment (maps to expected contig) ---
        anchor_idx = None
        for i, seg in enumerate(chained):
            if chr_end_matches_contig(chr_end, seg['ref_name']):
                anchor_idx = i
                break

        if anchor_idx is None:
            qc_flags.append('anchor_mismatch')
            # Check if any segment is clean enough for one_clean_half
            if chained[0]['identity'] >= MIN_SEGMENT_IDENTITY:
                alignment_class = 'one_clean_half'
            else:
                alignment_class = 'no_clean_halves'
                qc_flags.append('low_quality_anchor')
        else:
            anchor_seg = chained[anchor_idx]
            if anchor_seg['mapq'] < MIN_MAPQ:
                qc_flags.append('ambiguous_alignment')

            if anchor_seg['identity'] < MIN_SEGMENT_IDENTITY:
                alignment_class = 'no_clean_halves'
                qc_flags.append('low_quality_anchor')
            elif len(chained) == 1:
                if (len(cleaned) == 1
                        and anchor_seg['identity'] >= MIN_WHOLE_READ_IDENTITY
                        and anchor_seg['coverage'] >= MIN_WHOLE_READ_COVERAGE):
                    alignment_class = 'no_recombination'
                else:
                    alignment_class = 'one_clean_half'
            else:
                alignment_class = 'two_clean_halves'

        # --- Reorder chained segments: anchor first, then toward telomere ---
        # For the output, segment_1 = anchor, segment_2 = next toward telomere
        if anchor_idx is not None and anchor_idx > 0:
            # Anchor is not the first segment on the read — the telomere side
            # is to the LEFT. Reverse so anchor comes first.
            chained = list(reversed(chained))
            anchor_idx = len(chained) - 1 - anchor_idx

        # --- Build breakpoints (between consecutive chained segments) ---
        breakpoints = []
        if alignment_class == 'two_clean_halves' and len(chained) >= 2:
            for j in range(len(chained) - 1):
                prev_seg = chained[j]
                next_seg = chained[j + 1]

                # Breakpoint on the read: midpoint of overlap zone, or junction if no overlap
                overlap_start = max(prev_seg['query_start'], next_seg['query_start'])
                overlap_end = min(prev_seg['query_end'], next_seg['query_end'])
                if overlap_end > overlap_start:
                    # Segments overlap — breakpoint at midpoint of overlap
                    bp_on_read = (overlap_start + overlap_end) // 2
                else:
                    # No overlap — breakpoint at the gap
                    bp_on_read = prev_seg['query_end']

                # Breakpoint on the reference: use prev_seg's boundary facing next_seg
                if prev_seg['query_start'] < next_seg['query_start']:
                    bp_on_ref = _ref_pos_at_query_boundary(prev_seg, 'right')
                else:
                    bp_on_ref = _ref_pos_at_query_boundary(prev_seg, 'left')

                feat_name, feat_type, is_mid = classify_breakpoint(bp_on_ref, features)

                breakpoints.append({
                    'pos_on_read': bp_on_read,
                    'pos_on_ref': bp_on_ref,
                    'element': feat_name,
                    'feature_type': feat_type,
                    'is_mid_element': is_mid,
                    'post_contig': next_seg['ref_name'],
                })

        # --- Extract diverged tail ---
        diverged_seq = ''
        diverged_seq_length = 0

        if alignment_class in ('one_clean_half', 'two_clean_halves'):
            anchor_seg = chained[0] if anchor_idx is not None else chained[0]

            # Diverged tail = larger unaligned portion of the read
            left_len = anchor_seg['query_start']
            right_len = read_length - anchor_seg['query_end']

            if left_len >= right_len and left_len >= MIN_DIVERGED_LEN:
                diverged_seq = raw_seq[:anchor_seg['query_start']]
                bp_on_read = anchor_seg['query_start']
                bp_on_ref = _ref_pos_at_query_boundary(anchor_seg, 'left')
            elif right_len >= MIN_DIVERGED_LEN:
                diverged_seq = raw_seq[anchor_seg['query_end']:]
                bp_on_read = anchor_seg['query_end']
                bp_on_ref = _ref_pos_at_query_boundary(anchor_seg, 'right')
            else:
                qc_flags.append('short_tail')

            diverged_seq_length = len(diverged_seq)

            if diverged_seq and len(diverged_seq) >= MIN_DIVERGED_LEN:
                diverged_tails.append((f'{read_id} diverged_start={bp_on_read}', diverged_seq))

                # For one_clean_half, add a breakpoint at the diverged tail boundary
                if alignment_class == 'one_clean_half' and not breakpoints:
                    feat_name, feat_type, is_mid = classify_breakpoint(bp_on_ref, features)
                    breakpoints.append({
                        'pos_on_read': bp_on_read,
                        'pos_on_ref': bp_on_ref,
                        'element': feat_name,
                        'feature_type': feat_type,
                        'is_mid_element': is_mid,
                        'post_contig': '',
                    })

        results.append(_make_row(
            read_id, chr_end, read_length, telo_side,
            alignment_class=alignment_class,
            chained=chained, breakpoints=breakpoints,
            diverged_seq_length=diverged_seq_length,
            qc_flags=qc_flags,
        ))

    return results, diverged_tails

# ---------------------------------------------------------------------------
# Output row construction
# ---------------------------------------------------------------------------

def _make_row(read_id, chr_end, read_length, telo_side,
              alignment_class, chained, breakpoints,
              diverged_seq_length, qc_flags):
    """Build a flat dict for the output TSV."""
    row = {
        'read_id': read_id,
        'chr_end': chr_end,
        'read_length': read_length,
        'telo_side': telo_side,
        'alignment_class': alignment_class,
    }

    for i in range(3):
        prefix = f'segment_{i+1}_'
        if i < len(chained):
            seg = chained[i]
            row[prefix + 'contig'] = seg['ref_name']
            row[prefix + 'identity'] = seg['identity']
            row[prefix + 'coverage'] = seg['coverage']
            row[prefix + 'start'] = seg['query_start']
            row[prefix + 'end'] = seg['query_end']
        else:
            row[prefix + 'contig'] = ''
            row[prefix + 'identity'] = ''
            row[prefix + 'coverage'] = ''
            row[prefix + 'start'] = ''
            row[prefix + 'end'] = ''

    row['n_breakpoints'] = len(breakpoints)
    if breakpoints:
        row['breakpoint_elements'] = ';'.join(bp['element'] for bp in breakpoints)
        row['breakpoint_feature_types'] = ';'.join(bp['feature_type'] for bp in breakpoints)
        row['breakpoint_positions_on_read'] = ';'.join(str(bp['pos_on_read']) for bp in breakpoints)
        row['breakpoint_positions_on_ref'] = ';'.join(str(bp['pos_on_ref']) for bp in breakpoints)
        row['breakpoint_is_mid_element'] = ';'.join(str(bp['is_mid_element']) for bp in breakpoints)
        row['breakpoint_post_contigs'] = ';'.join(bp['post_contig'] for bp in breakpoints)
    else:
        row['breakpoint_elements'] = ''
        row['breakpoint_feature_types'] = ''
        row['breakpoint_positions_on_read'] = ''
        row['breakpoint_positions_on_ref'] = ''
        row['breakpoint_is_mid_element'] = ''
        row['breakpoint_post_contigs'] = ''

    row['diverged_seq_length'] = diverged_seq_length
    row['qc_flags'] = ';'.join(qc_flags)

    return row

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    print(f'run_alignment.py — chr_end={args.chr_end}')

    for path in [args.reads_fasta, args.day0_ref, args.day0_bed]:
        if not os.path.exists(path):
            print(f'ERROR: missing input: {path}', file=sys.stderr)
            sys.exit(1)

    print('  Loading BED features…')
    features = load_bed_features(args.day0_bed, args.chr_end)
    print(f'  Found {len(features)} features for {args.chr_end}')

    print('  Reading read sequences…')
    read_seqs, read_headers = read_fasta_sequences(args.reads_fasta)
    n_reads = len(read_seqs)
    print(f'  {n_reads} reads loaded')

    if n_reads == 0:
        print('  No reads — writing empty outputs and exiting')
        write_results_tsv([], args.output_tsv)
        tails_path = args.output_tsv.replace('_breakpoints.tsv', '_diverged_tails.fasta')
        open(tails_path, 'w').close()
        return

    os.makedirs(os.path.dirname(args.output_bam), exist_ok=True)
    run_minimap2(args.reads_fasta, args.day0_ref, args.output_bam, args.threads)

    print('  Parsing BAM and classifying reads…')
    rows, diverged_tails = parse_bam_to_breakpoints(
        args.output_bam, args.chr_end, features, read_seqs, read_headers,
    )

    os.makedirs(os.path.dirname(args.output_tsv) or '.', exist_ok=True)
    write_results_tsv(rows, args.output_tsv)
    print(f'  Breakpoints TSV: {args.output_tsv} ({len(rows)} reads)')

    tails_path = args.output_tsv.replace('_breakpoints.tsv', '_diverged_tails.fasta')
    write_fasta_list(diverged_tails, tails_path)
    print(f'  Diverged tails FASTA: {tails_path} ({len(diverged_tails)} sequences)')

    from collections import Counter
    class_counts = Counter(r['alignment_class'] for r in rows)
    print(f'  Classification summary:')
    for cls in ['no_recombination', 'two_clean_halves', 'one_clean_half', 'no_clean_halves']:
        print(f'    {cls}: {class_counts.get(cls, 0)}')


if __name__ == '__main__':
    main()
