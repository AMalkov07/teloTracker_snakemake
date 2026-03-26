"""
Step 11a: Y prime recombination analysis.

Runs RepeatMasker on the FULL per-chr-end reads FASTA (not just diverged tails),
then compares Y prime positions on each read against the day0 reference order
derived from the labeled BED file.  Mirrors the position-by-position logic of
the existing get_stats_of_recombination.py but adds confidence scoring,
color-group ambiguity, downstream consistency checking, and compound detection.

Usage:
  python analyze_y_prime_recombination.py \
      --breakpoints-tsv  results/{base}/recombination/{base}_{chr_end}_breakpoints.tsv \
      --reads-fasta      results/{base}/chr_anchor_included_individual_files/{base}_blasted_{anchor}_{chr_end}_anchor_reads.fasta \
      --diverged-fasta   results/{base}/recombination/{base}_{chr_end}_diverged_tails.fasta \
      --day0-bed         results/{day0_base}/pretelomeric_labels/pretelomeric_regions_{strain}_simp.bed \
      --y-prime-lib      references/extracted_yprimes_{strain}.fasta \
      --chr-end          chr4R  --strain 6991 \
      --output-tsv       results/{base}/recombination/{base}_{chr_end}_y_prime_recomb.tsv \
      --threads          4
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile

import pandas as pd

from recombination_utils import (
    Hypothesis,
    hypotheses_to_row_dict,
    is_recombination_candidate,
    get_first_breakpoint_feature_type,
    get_first_breakpoint_is_mid_element,
    get_first_breakpoint_pos_on_read,
    load_bed_features,
    normalize_hypotheses,
    write_results_tsv,
    MIN_CONFIDENCE_TO_REPORT,
)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description='Y prime recombination analysis')
    p.add_argument('--breakpoints-tsv', required=True)
    p.add_argument('--reads-fasta',     required=True)
    p.add_argument('--diverged-fasta',  required=True)
    p.add_argument('--day0-bed',        required=True)
    p.add_argument('--y-prime-lib',     required=True)
    p.add_argument('--chr-end',         required=True)
    p.add_argument('--strain',          required=True)
    p.add_argument('--output-tsv',      required=True)
    p.add_argument('--threads',         type=int, default=4)
    return p.parse_args()

# ---------------------------------------------------------------------------
# Y prime library parsing
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Read sequences (for post-Y-prime telomere check)
# ---------------------------------------------------------------------------

def load_read_sequences(fasta_path):
    """Load read sequences from FASTA, keyed by read_id (first word of header)."""
    seqs = {}
    current_id = None
    current_seq = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id is not None:
        seqs[current_id] = ''.join(current_seq)
    return seqs


# ---------------------------------------------------------------------------
# Post-Y-prime telomere validation
# ---------------------------------------------------------------------------

def _telomeric_fraction(sequence):
    """
    Estimate what fraction of a sequence is telomeric repeats.
    Yeast telomeric repeats are TG1-3 / C1-3A (irregular).
    Check if the sequence is dominated by T+G or C+A bases.
    True telomeric sequence is typically >80% T+G or >80% C+A;
    random DNA is ~50%.
    """
    seq = sequence.upper()
    if not seq:
        return 0.0
    n = len(seq)
    tg_count = sum(1 for base in seq if base in 'TG')
    ca_count = sum(1 for base in seq if base in 'CA')
    return max(tg_count / n, ca_count / n)


def check_post_y_prime_sequence(read_seq, ordered_hits, telo_side):
    """
    Check if the sequence after the last Y prime (toward the telomere) is telomeric.

    ordered_hits: anchor-first list from order_y_primes_for_read
      → last entry = telomere-proximal Y prime

    Returns (tail_length, telo_fraction, check_result):
      - tail_length: bp between last Y prime and telomere end of read
      - telo_fraction: fraction of that tail that looks telomeric (0-1)
      - check_result: 'pass' | 'fail' | 'n/a'
    """
    if not ordered_hits:
        return 0, 0.0, 'n/a'

    # Last hit in anchor-first order = telomere-proximal
    telo_hit = ordered_hits[-1]

    if telo_side == 'end':
        # Telomere at read end → tail is after the last hit
        tail_seq = read_seq[telo_hit['match_end_on_read']:]
    else:
        # Telomere at read start → tail is before the last (telomere-proximal) hit
        tail_seq = read_seq[:telo_hit['match_start_on_read']]

    tail_length = len(tail_seq)

    if tail_length < 50:
        return tail_length, 1.0, 'pass'

    telo_fraction = _telomeric_fraction(tail_seq)

    # Short tails (<200bp) with moderate telomeric content are fine
    if tail_length < 200 and telo_fraction >= 0.55:
        return tail_length, round(telo_fraction, 4), 'pass'

    if telo_fraction >= 0.65:
        return tail_length, round(telo_fraction, 4), 'pass'

    return tail_length, round(telo_fraction, 4), 'fail'


# ---------------------------------------------------------------------------
# Y prime library parsing
# ---------------------------------------------------------------------------

def parse_y_prime_header(header):
    """
    Parse a Y prime FASTA header like:
      >Y_Prime_chr2L1#Short/Solo/ID4_Green-Light
    Returns dict with keys: id, color, color_group, size, array_type, origin.
    """
    header = header.lstrip('>')
    # Split on '#' to separate name from classification
    if '#' in header:
        name_part, class_part = header.split('#', 1)
    else:
        name_part = header
        class_part = ''

    origin = name_part  # e.g. Y_Prime_chr2L1

    size = ''
    array_type = ''
    color_group = ''
    y_id = ''
    color = ''

    # class_part looks like "Short/Solo/ID4_Green-Light"
    if class_part:
        parts = class_part.split('/')
        if len(parts) >= 1:
            size = parts[0]
        if len(parts) >= 2:
            array_type = parts[1]
        if len(parts) >= 3:
            color_group = parts[2]  # e.g. "ID4_Green-Light"
            # Split ID from color
            cg_parts = color_group.split('_', 1)
            y_id = cg_parts[0]        # e.g. "ID4"
            color = cg_parts[1] if len(cg_parts) > 1 else ''

    return {
        'id': y_id,
        'color': color,
        'color_group': color_group,
        'size': size,
        'array_type': array_type,
        'origin': origin,
        'full_name': header.split()[0],
    }


def build_y_prime_color_map(y_prime_lib_fasta):
    """
    Parse the Y prime library FASTA.
    Returns:
      color_map:    {color_group: [seq_name, ...]}
      name_to_info: {seq_name: parsed_header_dict}
    """
    color_map = {}
    name_to_info = {}
    with open(y_prime_lib_fasta) as fh:
        for line in fh:
            line = line.strip()
            if not line.startswith('>'):
                continue
            info = parse_y_prime_header(line)
            seq_name = info['full_name']
            name_to_info[seq_name] = info
            cg = info['color_group']
            if cg not in color_map:
                color_map[cg] = []
            color_map[cg].append(seq_name)
    return color_map, name_to_info

# ---------------------------------------------------------------------------
# Reference Y prime order from BED
# ---------------------------------------------------------------------------

def get_reference_y_prime_order(features):
    """
    Extract Y prime features from the BED features list (for this chr_end).
    Sort by start position ascending → anchor-proximal first (Y_prime_1 has lowest coordinate).
    Returns list of {'name', 'color_group', 'id'} in anchor-first order.
    """
    yp_features = [f for f in features if f['feature_type'] == 'y_prime']
    yp_features.sort(key=lambda f: f['start'])
    result = []
    for f in yp_features:
        # Extract color_group and id from the BED feature name
        # Feature name pattern: {chr_end}_Y_Prime_{N}_{color_group}  (varies)
        # Extract the ID portion (e.g. "ID4") using a regex
        id_match = re.search(r'(ID\d+)', f['name'])
        y_id = id_match.group(1) if id_match else ''
        # Color group: everything after the last ID token if present
        cg_match = re.search(r'(ID\d+[_\w-]*)', f['name'])
        color_group = cg_match.group(1) if cg_match else f['name']
        result.append({
            'feature_name': f['name'],
            'color_group': color_group,
            'id': y_id,
        })
    return result

# ---------------------------------------------------------------------------
# RepeatMasker
# ---------------------------------------------------------------------------

def run_repeatmasker_on_full_reads(reads_fasta, y_prime_lib, tmp_dir, threads):
    """
    Run RepeatMasker on the full reads FASTA against the Y prime library.
    Returns DataFrame: read_id, y_prime_name, match_start_on_read,
                       match_end_on_read, sw_score, divergence_pct.
    """
    print('    Running RepeatMasker on full reads…')
    rm_dir = os.path.join(tmp_dir, 'repeatmasker_y_prime')
    os.makedirs(rm_dir, exist_ok=True)

    result = subprocess.run(
        [
            'RepeatMasker', reads_fasta,
            '-lib', y_prime_lib,
            '-s', '-pa', str(threads),
            '--cutoff', '1000',
            '-no_is', '-norna',
            '-dir', rm_dir,
        ],
        capture_output=True,
    )

    out_file = os.path.join(rm_dir, os.path.basename(reads_fasta) + '.out')
    if not os.path.exists(out_file):
        print('    RepeatMasker produced no output file')
        return pd.DataFrame(columns=['read_id', 'y_prime_name', 'match_start_on_read',
                                     'match_end_on_read', 'sw_score', 'divergence_pct', 'strand'])

    # Parse .out file (fixed-width / whitespace-delimited)
    rows = []
    with open(out_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('SW') or line.startswith('score'):
                continue
            parts = line.split()
            if len(parts) < 14:
                continue
            try:
                sw_score = int(parts[0])
                divergence_pct = float(parts[1])
                read_id = parts[4]
                match_start = int(parts[5])
                match_end = int(parts[6])
                strand = parts[8]
                y_prime_name = parts[9]
                rows.append({
                    'read_id': read_id,
                    'y_prime_name': y_prime_name,
                    'match_start_on_read': match_start,
                    'match_end_on_read': match_end,
                    'sw_score': sw_score,
                    'divergence_pct': divergence_pct,
                    'strand': strand,
                })
            except (ValueError, IndexError):
                continue

    df = pd.DataFrame(rows)
    print(f'    RepeatMasker: {len(df)} hits across all reads')
    return df

# ---------------------------------------------------------------------------
# Per-read Y prime ordering
# ---------------------------------------------------------------------------

def order_y_primes_for_read(rm_hits_for_read, telo_side):
    """
    Sort RepeatMasker hits by match_start_on_read.
    If telo_side=='beginning': reverse so index 0 = anchor-proximal (Y_prime_1).
    Filters out simple repeat annotations (e.g., (CA)n, (A)n) that RepeatMasker
    reports alongside real Y prime matches.
    Returns list of hit dicts in anchor-first order.
    """
    hits = rm_hits_for_read.to_dict('records')
    # Filter to real Y prime hits only — exclude simple repeats like (CA)n, (TGGAA)n, etc.
    hits = [h for h in hits if h['y_prime_name'].startswith('Y_Prime_')]
    hits.sort(key=lambda h: h['match_start_on_read'])
    if telo_side == 'beginning':
        hits = hits[::-1]
    return hits

# ---------------------------------------------------------------------------
# Position-by-position comparison (mirrors get_stats_of_recombination.py logic)
# ---------------------------------------------------------------------------

def find_first_divergence(read_y_primes, ref_y_primes, name_to_info):
    """
    Compare Y primes on read vs reference order position-by-position.

    read_y_primes: list of hit dicts (from order_y_primes_for_read), anchor-first
    ref_y_primes:  list of reference dicts (from get_reference_y_prime_order), anchor-first
    name_to_info:  {seq_name: parsed_header_dict}

    Returns (divergence_idx, actual_hit, expected_ref, downstream_consistent).
    divergence_idx=None if read matches reference exactly.
    """
    # Determine color_group for each read hit
    def hit_color_group(hit):
        info = name_to_info.get(hit['y_prime_name'], {})
        return info.get('color_group', hit['y_prime_name'])

    # Check for Y prime loss (read has fewer Y primes than reference)
    if len(read_y_primes) < len(ref_y_primes):
        divergence_idx = len(read_y_primes)
        expected_ref = ref_y_primes[divergence_idx] if divergence_idx < len(ref_y_primes) else None
        return divergence_idx, None, expected_ref, True

    # Position-by-position comparison
    for i in range(len(ref_y_primes)):
        if i >= len(read_y_primes):
            break
        read_cg = hit_color_group(read_y_primes[i])
        ref_cg = ref_y_primes[i]['color_group']
        if read_cg != ref_cg:
            # Divergence found — check downstream consistency
            downstream_consistent = True
            if i + 1 < len(read_y_primes):
                diverged_cg = read_cg
                for j in range(i + 1, len(read_y_primes)):
                    if hit_color_group(read_y_primes[j]) != diverged_cg:
                        downstream_consistent = False
                        break
            return i, read_y_primes[i], ref_y_primes[i], downstream_consistent

    return None, None, None, True  # no divergence

# ---------------------------------------------------------------------------
# Hypothesis generation
# ---------------------------------------------------------------------------

def generate_y_prime_hypotheses(
    read_id, breakpoint_info_row, divergence_idx, actual_hit, expected_ref,
    downstream_consistent, color_map, name_to_info, all_read_y_primes,
):
    """Generate Hypothesis objects for one read's Y prime analysis."""
    hypotheses = []

    recomb_flagged = is_recombination_candidate(breakpoint_info_row)
    is_mid_element = get_first_breakpoint_is_mid_element(breakpoint_info_row)
    bp_feature_type = get_first_breakpoint_feature_type(breakpoint_info_row)
    bp_pos_on_read = get_first_breakpoint_pos_on_read(breakpoint_info_row)

    # Case 1: no RepeatMasker hits on full read
    if len(all_read_y_primes) == 0:
        hypotheses.append(Hypothesis(
            h_type='ambiguous',
            description="No Y' matches found on read",
            confidence=0.3,
            ambiguous=True,
        ))
        return normalize_hypotheses(hypotheses)

    # Case 2: no divergence from reference
    if divergence_idx is None:
        if recomb_flagged:
            hypotheses.append(Hypothesis(
                h_type='no_y_prime_change',
                description="Y' array matches reference; recombination may be in spacer/X element only",
                confidence=0.8,
            ))
        else:
            hypotheses.append(Hypothesis(
                h_type='no_recombination',
                description="Y' array matches reference exactly",
                confidence=0.9,
            ))
        return normalize_hypotheses(hypotheses)

    # Case 3: Y prime loss (actual_hit is None)
    if actual_hit is None:
        exp_id = expected_ref['id'] if expected_ref else '?'
        hypotheses.append(Hypothesis(
            h_type='y_prime',
            description=f"Loss of Y' at position {divergence_idx + 1} (expected {exp_id}, found none)",
            confidence=0.7,
            breakpoint_element=expected_ref['feature_name'] if expected_ref else '',
            ambiguous=True,
        ))
        return normalize_hypotheses(hypotheses)

    # Case 4: divergence found
    y_prime_name = actual_hit['y_prime_name']
    info_at_divergence = name_to_info.get(y_prime_name, {})
    actual_cg = info_at_divergence.get('color_group', y_prime_name)
    exp_id = expected_ref['id'] if expected_ref else '?'

    # Base confidence: penalize divergence, reward high SW score
    # Use a normalized SW score relative to a typical expected value (1500)
    sw = actual_hit.get('sw_score', 500)
    div_pct = actual_hit.get('divergence_pct', 5.0)
    c_base = (1.0 - div_pct / 100.0) * min(1.0, sw / 1500.0)
    c_base = max(0.05, c_base)

    # Downstream consistency
    if downstream_consistent:
        c_base *= 1.2
    else:
        c_base *= 0.7

    # Mid-element penalty
    if is_mid_element and bp_feature_type == 'y_prime':
        c_base *= 0.7

    # Same-color confidence distribution
    same_color_seqs = color_map.get(actual_cg, [actual_cg])
    n_same_color = max(1, len(same_color_seqs))

    # Collect SW scores for same-color hits in this read to decide equal vs weighted split
    same_color_hits = [
        h for h in all_read_y_primes
        if name_to_info.get(h['y_prime_name'], {}).get('color_group', '') == actual_cg
    ]
    if len(same_color_hits) > 1:
        sw_scores = [h['sw_score'] for h in same_color_hits]
        sw_max = max(sw_scores)
        sw_min = min(sw_scores)
        use_weighted = (sw_max / max(sw_min, 1)) > 1.5
    else:
        use_weighted = False

    # Build one hypothesis per unique source chr end implied by same-color sequences
    # (color group implies which chr ends are ambiguous sources)
    # For simplicity, map each same-color sequence to its origin chr end
    source_chr_ends_seen = {}
    for seq_name in same_color_seqs:
        origin_info = name_to_info.get(seq_name, {})
        origin = origin_info.get('origin', seq_name)
        # Extract chr_end from origin (e.g. "Y_Prime_chr2L1" → "chr2L")
        m = re.search(r'(chr\d+[LR])', origin)
        origin_chr = m.group(1) if m else origin
        source_chr_ends_seen.setdefault(origin_chr, []).append(seq_name)

    if use_weighted:
        total_sw = sum(sw_scores)
        weight_for_actual = sw / max(total_sw, 1)
        c_for_actual = c_base * weight_for_actual
    else:
        c_for_actual = c_base / n_same_color

    # Primary hypothesis for the actual divergence
    source_list = list(source_chr_ends_seen.keys())
    bp_on_read = actual_hit.get('match_start_on_read', bp_pos_on_read)
    hypotheses.append(Hypothesis(
        h_type='y_prime',
        description=(
            f"Y' divergence at position {divergence_idx + 1}: "
            f"expected {exp_id}, found {actual_cg}"
        ),
        confidence=c_for_actual,
        breakpoint_element='',
        breakpoint_pos_on_read=bp_on_read,
        source_chr_ends=source_list,
        ambiguous=(n_same_color > 1),
        is_compound=False,
    ))

    # Additional hypotheses for other same-color sources
    if n_same_color > 1 and len(source_chr_ends_seen) > 1:
        remaining_sources = [s for s in source_chr_ends_seen if s not in source_list[:1]]
        c_others = c_base - c_for_actual
        if remaining_sources and c_others >= MIN_CONFIDENCE_TO_REPORT:
            hypotheses.append(Hypothesis(
                h_type='y_prime',
                description=(
                    f"Y' divergence at position {divergence_idx + 1}: "
                    f"expected {exp_id}, found {actual_cg} (alternative sources)"
                ),
                confidence=c_others,
                source_chr_ends=remaining_sources,
                ambiguous=True,
            ))

    # Compound detection: if downstream is inconsistent AND there's a second color group
    if not downstream_consistent and len(all_read_y_primes) > divergence_idx + 1:
        post_cgs = set()
        for j in range(divergence_idx + 1, len(all_read_y_primes)):
            cg_j = name_to_info.get(all_read_y_primes[j]['y_prime_name'], {}).get('color_group', '')
            if cg_j and cg_j != actual_cg:
                post_cgs.add(cg_j)
        if post_cgs:
            c_compound = c_for_actual * 0.5  # compound events get half confidence
            hypotheses.append(Hypothesis(
                h_type='compound',
                description=(
                    f"Compound: first Y' switch to {actual_cg}, "
                    f"then switch to {';'.join(post_cgs)}"
                ),
                confidence=c_compound,
                is_compound=True,
                ambiguous=True,
            ))

    return normalize_hypotheses(hypotheses)

# ---------------------------------------------------------------------------
# Cross-validate with minimap2 breakpoint
# ---------------------------------------------------------------------------

def cross_validate_with_minimap2(divergence_idx, actual_hit, bp_pos_on_read, hypotheses):
    """
    Compare where the Y prime divergence occurred on the read with where minimap2
    detected the breakpoint.  Apply a small boost or penalty to confidence.
    Returns (validation_result, updated_hypotheses).
    """
    if actual_hit is None or bp_pos_on_read < 0:
        return 'n/a', hypotheses

    y_prime_bp = actual_hit.get('match_start_on_read', -1)
    if y_prime_bp < 0:
        return 'n/a', hypotheses

    distance = abs(y_prime_bp - bp_pos_on_read)
    if distance <= 1000:
        validation = 'confirms'
        for h in hypotheses:
            h.confidence = min(1.0, h.confidence * 1.1)
    else:
        validation = 'contradicts'
        for h in hypotheses:
            h.confidence *= 0.8

    return validation, normalize_hypotheses(hypotheses)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    print(f'analyze_y_prime_recombination.py — chr_end={args.chr_end}')

    for path in [args.breakpoints_tsv, args.reads_fasta, args.day0_bed, args.y_prime_lib]:
        if not os.path.exists(path):
            print(f'ERROR: missing input: {path}', file=sys.stderr)
            sys.exit(1)

    # Load breakpoints
    df_bp = pd.read_csv(args.breakpoints_tsv, sep='\t')
    if df_bp.empty:
        print('  No reads in breakpoints TSV — writing empty output')
        write_results_tsv([], args.output_tsv)
        return

    # Load BED features and Y prime reference order
    features = load_bed_features(args.day0_bed, args.chr_end)
    ref_y_primes = get_reference_y_prime_order(features)
    print(f'  Reference Y prime order for {args.chr_end}: '
          f'{[r["id"] for r in ref_y_primes]} ({len(ref_y_primes)} Y primes)')

    # Build Y prime color map from library
    color_map, name_to_info = build_y_prime_color_map(args.y_prime_lib)
    print(f'  Y prime library: {len(name_to_info)} sequences, {len(color_map)} color groups')

    # Load read sequences for post-Y-prime telomere check
    print('  Loading read sequences for telomere validation…')
    read_seqs = load_read_sequences(args.reads_fasta)
    print(f'  {len(read_seqs)} read sequences loaded')

    # Run RepeatMasker on full reads
    with tempfile.TemporaryDirectory() as tmp_dir:
        df_rm = run_repeatmasker_on_full_reads(
            args.reads_fasta, args.y_prime_lib, tmp_dir, args.threads
        )

    # Process each read
    rows = []
    for _, bp_row in df_bp.iterrows():
        read_id = bp_row['read_id']
        telo_side = bp_row.get('telo_side', 'end')
        bp_pos_on_read = get_first_breakpoint_pos_on_read(bp_row)

        # Get all RepeatMasker hits for this read
        if df_rm.empty or 'read_id' not in df_rm.columns:
            rm_for_read = pd.DataFrame()
        else:
            rm_for_read = df_rm[df_rm['read_id'] == read_id]

        # Order Y primes anchor-first
        ordered_hits = order_y_primes_for_read(rm_for_read, telo_side)

        # Position-by-position comparison
        divergence_idx, actual_hit, expected_ref, downstream_consistent = find_first_divergence(
            ordered_hits, ref_y_primes, name_to_info
        )

        # Generate hypotheses
        hypotheses = generate_y_prime_hypotheses(
            read_id, bp_row, divergence_idx, actual_hit, expected_ref,
            downstream_consistent, color_map, name_to_info, ordered_hits,
        )

        # Cross-validate with minimap2 breakpoint
        validation_result, hypotheses = cross_validate_with_minimap2(
            divergence_idx, actual_hit, bp_pos_on_read, hypotheses
        )

        # Build output row
        row = {
            'read_id': read_id,
            'chr_end': args.chr_end,
            'y_prime_divergence_idx': divergence_idx if divergence_idx is not None else -1,
            'y_prime_expected_at_divergence': (
                expected_ref['id'] if expected_ref else ''
            ),
            'y_prime_found_at_divergence': (
                name_to_info.get(actual_hit['y_prime_name'], {}).get('color_group', '')
                if actual_hit else ''
            ),
            'y_prime_downstream_consistent': downstream_consistent,
            'y_prime_recombination_status': _recomb_status(
                divergence_idx, actual_hit, expected_ref, downstream_consistent,
                len(ordered_hits), ref_y_primes,
            ),
            'y_prime_validation_result': validation_result,
        }

        # Post-Y-prime telomere check
        read_seq = read_seqs.get(read_id, '')
        tail_len, telo_frac, telo_check = check_post_y_prime_sequence(
            read_seq, ordered_hits, telo_side
        )
        row['post_y_prime_seq_length'] = tail_len
        row['post_y_prime_telo_fraction'] = telo_frac
        row['post_y_prime_check'] = telo_check

        row['n_hypotheses'] = len(hypotheses)
        row.update(hypotheses_to_row_dict(hypotheses))
        rows.append(row)

    os.makedirs(os.path.dirname(args.output_tsv) or '.', exist_ok=True)
    write_results_tsv(rows, args.output_tsv)
    print(f'  Output: {args.output_tsv} ({len(rows)} reads)')


def _recomb_status(divergence_idx, actual_hit, expected_ref, downstream_consistent,
                   n_read_y_primes, ref_y_primes):
    """Map divergence info to a human-readable status string."""
    if divergence_idx is None:
        return 'No Change'
    if actual_hit is None:
        return "Y' Loss"
    if divergence_idx == 0:
        return "1st Y' Change"
    if not downstream_consistent:
        return 'Ambiguous'
    return "Y' Recombination"


if __name__ == '__main__':
    main()
