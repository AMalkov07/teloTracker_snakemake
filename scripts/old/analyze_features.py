"""
Step 11: Combined feature analysis — structural walk from anchor to telomere.

Replaces the three separate scripts (analyze_y_prime_recombination.py,
analyze_x_prime_recombination.py, analyze_spacer_recombination.py).

Runs batch preprocessing (BLAST) upfront, then walks each read
from anchor outward through the subtelomeric features in structural order:
  anchor → spacer → X element → Y primes → telomere

For each read, determines:
  1. feature_reach — which features the read actually covers
  2. features_in_order — are features in expected positions?
  3. Recombination analysis for each feature type
  4. Post-Y-prime telomere validation

Usage:
  python analyze_features.py \
      --breakpoints-tsv  results/{base}/recombination/{base}_{chr_end}_breakpoints.tsv \
      --reads-fasta      results/{base}/chr_anchor_included_individual_files/{base}_blasted_{anchor}_{chr_end}_anchor_reads.fasta \
      --diverged-fasta   results/{base}/recombination/{base}_{chr_end}_diverged_tails.fasta \
      --day0-bed         results/{day0}/pretelomeric_labels/pretelomeric_regions_{strain}_simp.bed \
      --y-prime-lib      references/extracted_yprimes_{strain}.fasta \
      --spacer-lib-dir   references/pairings_for_spacers/{strain}_pairings/ \
      --x-element-lib-dir references/pairings_for_x_element_ends/{strain}_pairings/ \
      --chr-end          chr4R  --strain 7302 \
      --output-tsv       results/{base}/recombination/{base}_{chr_end}_features.tsv \
      --threads          4
"""

import argparse
import os
import re
import sys
import tempfile

import pandas as pd

from recombination_utils import (
    Hypothesis,
    hypotheses_to_row_dict,
    is_recombination_candidate,
    has_structural_evidence,
    get_first_breakpoint_feature_type,
    get_first_breakpoint_element,
    get_first_breakpoint_is_mid_element,
    get_first_breakpoint_pos_on_read,
    get_post_breakpoint_contig,
    load_bed_features,
    normalize_hypotheses,
    run_blast,
    read_fasta,
    write_fasta,
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
    p = argparse.ArgumentParser(description='Combined feature analysis — structural walk')
    p.add_argument('--breakpoints-tsv', required=True)
    p.add_argument('--reads-fasta',     required=True)
    p.add_argument('--diverged-fasta',  required=True)
    p.add_argument('--day0-bed',        required=True)
    p.add_argument('--y-prime-lib',     required=True)
    p.add_argument('--spacer-lib-dir',  required=True)
    p.add_argument('--x-element-lib-dir', required=True)
    p.add_argument('--chr-end',         required=True)
    p.add_argument('--strain',          required=True)
    p.add_argument('--output-tsv',      required=True)
    p.add_argument('--threads',         type=int, default=4)
    return p.parse_args()

# ---------------------------------------------------------------------------
# Read sequence loading
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
# Feature reach — determine which BED features a read covers
# ---------------------------------------------------------------------------

FEATURE_ORDER = ['anchor', 'spacer', 'x_core', 'x_variable', 'y_prime']

def determine_feature_reach(bp_row, features):
    """
    Determine which BED features a read covers based on alignment extent.

    For no_recombination: read covers all features.
    For one_clean_half/two_clean_halves: features up to breakpoint_pos_on_ref.
    For no_clean_halves: no features reached.

    Returns (feature_reach_list, features_in_order, order_issues).
    """
    alignment_class = str(bp_row.get('alignment_class', ''))

    if alignment_class == 'no_clean_halves':
        return [], True, ''

    if not features:
        return ['anchor'], True, 'no_bed_features'

    if alignment_class == 'no_recombination':
        # Read covers all features
        reached = []
        for f in features:
            ft = f['feature_type']
            if ft == 'y_prime':
                reached.append(f['name'])
            elif ft in FEATURE_ORDER:
                reached.append(ft)
        return reached, True, ''

    # one_clean_half or two_clean_halves: check up to breakpoint
    bp_ref_str = str(bp_row.get('breakpoint_positions_on_ref', ''))
    try:
        bp_on_ref = int(bp_ref_str.split(';')[0])
    except (ValueError, IndexError):
        return ['anchor'], True, 'no_breakpoint_ref_pos'

    reached = []
    for f in features:
        if f['end'] <= bp_on_ref:
            ft = f['feature_type']
            if ft == 'y_prime':
                reached.append(f['name'])
            elif ft in FEATURE_ORDER:
                reached.append(ft)
        elif f['start'] <= bp_on_ref:
            # Breakpoint is inside this feature
            ft = f['feature_type']
            if ft == 'y_prime':
                reached.append(f['name'] + '_partial')
            elif ft in FEATURE_ORDER:
                reached.append(ft + '_partial')

    # Validate order: check features appear in expected sequence
    order_issues = []
    expected_order_types = []
    for r in reached:
        # Normalize to base type
        base = r.replace('_partial', '')
        if 'Y_Prime' in base:
            expected_order_types.append('y_prime')
        elif base in FEATURE_ORDER:
            expected_order_types.append(base)

    prev_idx = -1
    for ft in expected_order_types:
        if ft in FEATURE_ORDER:
            idx = FEATURE_ORDER.index(ft)
            if idx < prev_idx:
                order_issues.append(f'{ft}_out_of_order')
            prev_idx = max(prev_idx, idx)

    features_in_order = len(order_issues) == 0
    return reached, features_in_order, ';'.join(order_issues)

# ---------------------------------------------------------------------------
# Y prime library parsing (ported from analyze_y_prime_recombination.py)
# ---------------------------------------------------------------------------

def parse_y_prime_header(header):
    """Parse a Y prime FASTA header like >Y_Prime_chr2L1#Short/Solo/ID4_Green-Light."""
    header = header.lstrip('>')
    if '#' in header:
        name_part, class_part = header.split('#', 1)
    else:
        name_part, class_part = header, ''

    origin = name_part
    size, array_type, color_group, y_id, color = '', '', '', '', ''

    if class_part:
        parts = class_part.split('/')
        if len(parts) >= 1: size = parts[0]
        if len(parts) >= 2: array_type = parts[1]
        if len(parts) >= 3:
            color_group = parts[2]
            cg_parts = color_group.split('_', 1)
            y_id = cg_parts[0]
            color = cg_parts[1] if len(cg_parts) > 1 else ''

    return {
        'id': y_id, 'color': color, 'color_group': color_group,
        'size': size, 'array_type': array_type, 'origin': origin,
        'full_name': header.split()[0],
    }


def build_y_prime_color_map(y_prime_lib_fasta):
    """Parse Y prime library. Returns (color_map, name_to_info)."""
    color_map = {}
    name_to_info = {}
    with open(y_prime_lib_fasta) as fh:
        for line in fh:
            if not line.strip().startswith('>'):
                continue
            info = parse_y_prime_header(line.strip())
            seq_name = info['full_name']
            name_to_info[seq_name] = info
            cg = info['color_group']
            color_map.setdefault(cg, []).append(seq_name)
    return color_map, name_to_info


def get_reference_y_prime_order(features):
    """Extract Y prime features from BED, sorted anchor-first."""
    yp_features = [f for f in features if f['feature_type'] == 'y_prime']
    yp_features.sort(key=lambda f: f['start'])
    result = []
    for f in yp_features:
        id_match = re.search(r'(ID\d+)', f['name'])
        y_id = id_match.group(1) if id_match else ''
        cg_match = re.search(r'(ID\d+[_\w-]*)', f['name'])
        color_group = cg_match.group(1) if cg_match else f['name']
        result.append({
            'feature_name': f['name'], 'color_group': color_group,
            'id': y_id, 'start': f['start'], 'end': f['end'],
        })
    return result

# ---------------------------------------------------------------------------
# Y prime BLAST (batch preprocessing — replaces RepeatMasker)
# ---------------------------------------------------------------------------

MAX_MERGE_GAP = 500   # max gap between BLAST hits to merge into one Y prime
MIN_YPRIME_HIT_LEN = 500  # minimum alignment length for a Y prime hit after merging

def blast_y_primes_on_full_reads(reads_fasta, y_prime_lib, tmp_dir, threads):
    """BLAST full reads against Y prime library. Returns DataFrame of hits
    with columns: read_id, y_prime_name, match_start_on_read, match_end_on_read,
    sw_score, divergence_pct, strand."""
    print('    BLASTing full reads against Y prime library…')

    blast_df = run_blast(reads_fasta, y_prime_lib, tmp_dir, label='y_prime_full',
                         min_pident=80.0, evalue=1e-10)

    if blast_df.empty:
        print('    Y prime BLAST: 0 hits')
        return pd.DataFrame(columns=['read_id', 'y_prime_name', 'match_start_on_read',
                                     'match_end_on_read', 'sw_score', 'divergence_pct', 'strand'])

    print(f'    Y prime BLAST: {len(blast_df)} raw hits')

    # Convert to internal format and merge adjacent hits from the same Y prime
    rows = []
    for _, h in blast_df.iterrows():
        rows.append({
            'read_id': h['qseqid'],
            'y_prime_name': h['sseqid'],
            'match_start_on_read': min(int(h['qstart']), int(h['qend'])),
            'match_end_on_read': max(int(h['qstart']), int(h['qend'])),
            'sw_score': int(h['bitscore']),
            'divergence_pct': round(100.0 - h['pident'], 2),
            'strand': '+' if h['qstart'] <= h['qend'] else '-',
        })

    df = pd.DataFrame(rows)

    # Merge adjacent hits: same read + same Y prime + within MAX_MERGE_GAP
    merged = []
    for (rid, yname), grp in df.groupby(['read_id', 'y_prime_name']):
        grp = grp.sort_values('match_start_on_read')
        current = grp.iloc[0].to_dict()
        for i in range(1, len(grp)):
            row = grp.iloc[i]
            gap = row['match_start_on_read'] - current['match_end_on_read']
            if gap <= MAX_MERGE_GAP:
                # Merge: extend range, sum scores, weighted-average divergence
                total_len = (current['match_end_on_read'] - current['match_start_on_read'] +
                             row['match_end_on_read'] - row['match_start_on_read'])
                cur_len = current['match_end_on_read'] - current['match_start_on_read']
                row_len = row['match_end_on_read'] - row['match_start_on_read']
                current['match_end_on_read'] = row['match_end_on_read']
                current['sw_score'] += row['sw_score']
                if total_len > 0:
                    current['divergence_pct'] = round(
                        (current['divergence_pct'] * cur_len + row['divergence_pct'] * row_len) / total_len, 2)
            else:
                merged.append(current)
                current = row.to_dict()
        merged.append(current)

    df_merged = pd.DataFrame(merged)
    print(f'    After merging adjacent hits: {len(df_merged)} Y prime hits')

    # Filter out short hits — real Y primes are several kb
    if not df_merged.empty:
        df_merged['hit_length'] = df_merged['match_end_on_read'] - df_merged['match_start_on_read']
        df_merged = df_merged[df_merged['hit_length'] >= MIN_YPRIME_HIT_LEN].copy()
        df_merged.drop(columns='hit_length', inplace=True)
        print(f'    After filtering hits < {MIN_YPRIME_HIT_LEN}bp: {len(df_merged)} Y prime hits')

    # Deduplicate overlapping hits on the same read: multiple library sequences
    # (e.g., chr5R_Y_Prime_1, chr8L_Y_Prime_1) can match the same Y prime on
    # a read. Keep only the best-scoring hit at each position.
    if not df_merged.empty:
        deduped = []
        for rid, grp in df_merged.groupby('read_id'):
            grp = grp.sort_values('sw_score', ascending=False)
            kept = []
            for _, hit in grp.iterrows():
                hit_start = hit['match_start_on_read']
                hit_end = hit['match_end_on_read']
                # Check if this hit overlaps >50% with any already-kept hit
                overlaps = False
                for k in kept:
                    overlap_start = max(hit_start, k['match_start_on_read'])
                    overlap_end = min(hit_end, k['match_end_on_read'])
                    if overlap_end > overlap_start:
                        overlap_len = overlap_end - overlap_start
                        hit_len = hit_end - hit_start
                        if overlap_len > 0.5 * hit_len:
                            overlaps = True
                            break
                if not overlaps:
                    kept.append(hit.to_dict())
            deduped.extend(kept)
        df_merged = pd.DataFrame(deduped)
        print(f'    After deduplicating overlapping hits: {len(df_merged)} Y prime hits')

    return df_merged

# ---------------------------------------------------------------------------
# Y prime per-read analysis (ported from analyze_y_prime_recombination.py)
# ---------------------------------------------------------------------------

def order_y_primes_for_read(rm_hits_for_read, telo_side):
    """Sort and orient RepeatMasker Y prime hits for anchor-first ordering."""
    hits = rm_hits_for_read.to_dict('records')
    hits.sort(key=lambda h: h['match_start_on_read'])
    if telo_side == 'beginning':
        hits = hits[::-1]
    return hits


def find_first_divergence(read_y_primes, ref_y_primes, name_to_info):
    """Compare Y primes position-by-position. Returns (div_idx, actual_hit, expected_ref, downstream_consistent)."""
    def hit_color_group(hit):
        info = name_to_info.get(hit['y_prime_name'], {})
        return info.get('color_group', hit['y_prime_name'])

    if len(read_y_primes) < len(ref_y_primes):
        divergence_idx = len(read_y_primes)
        expected_ref = ref_y_primes[divergence_idx] if divergence_idx < len(ref_y_primes) else None
        return divergence_idx, None, expected_ref, True

    for i in range(len(ref_y_primes)):
        if i >= len(read_y_primes):
            break
        read_cg = hit_color_group(read_y_primes[i])
        ref_cg = ref_y_primes[i]['color_group']
        if read_cg != ref_cg:
            downstream_consistent = True
            if i + 1 < len(read_y_primes):
                diverged_cg = read_cg
                for j in range(i + 1, len(read_y_primes)):
                    if hit_color_group(read_y_primes[j]) != diverged_cg:
                        downstream_consistent = False
                        break
            return i, read_y_primes[i], ref_y_primes[i], downstream_consistent

    return None, None, None, True


def generate_y_prime_hypotheses(read_id, bp_row, divergence_idx, actual_hit, expected_ref,
                                downstream_consistent, color_map, name_to_info, all_read_y_primes):
    """Generate hypotheses for Y prime analysis."""
    hypotheses = []
    recomb_flagged = is_recombination_candidate(bp_row)
    is_mid_element = get_first_breakpoint_is_mid_element(bp_row)
    bp_feature_type = get_first_breakpoint_feature_type(bp_row)
    bp_pos_on_read = get_first_breakpoint_pos_on_read(bp_row)

    if len(all_read_y_primes) == 0:
        hypotheses.append(Hypothesis(h_type='ambiguous',
            description="No Y' matches found on read", confidence=0.3, ambiguous=True))
        return normalize_hypotheses(hypotheses)

    if divergence_idx is None:
        if recomb_flagged:
            hypotheses.append(Hypothesis(h_type='no_y_prime_change',
                description="Y' array matches reference; recombination may be in spacer/X element only",
                confidence=0.8))
        else:
            hypotheses.append(Hypothesis(h_type='no_recombination',
                description="Y' array matches reference exactly", confidence=0.9))
        return normalize_hypotheses(hypotheses)

    if actual_hit is None:
        exp_id = expected_ref['id'] if expected_ref else '?'
        hypotheses.append(Hypothesis(h_type='y_prime',
            description=f"Loss of Y' at position {divergence_idx + 1} (expected {exp_id}, found none)",
            confidence=0.7, breakpoint_element=expected_ref['feature_name'] if expected_ref else '',
            ambiguous=True))
        return normalize_hypotheses(hypotheses)

    # Divergence found
    y_prime_name = actual_hit['y_prime_name']
    info_at_div = name_to_info.get(y_prime_name, {})
    actual_cg = info_at_div.get('color_group', y_prime_name)
    exp_id = expected_ref['id'] if expected_ref else '?'

    sw = actual_hit.get('sw_score', 500)
    div_pct = actual_hit.get('divergence_pct', 5.0)
    c_base = (1.0 - div_pct / 100.0) * min(1.0, sw / 1500.0)
    c_base = max(0.05, c_base)

    if downstream_consistent:
        c_base *= 1.2
    else:
        c_base *= 0.7

    if is_mid_element and bp_feature_type == 'y_prime':
        c_base *= 0.7

    same_color_seqs = color_map.get(actual_cg, [actual_cg])
    n_same_color = max(1, len(same_color_seqs))
    c_for_actual = c_base / n_same_color

    source_chr_ends_seen = {}
    for seq_name in same_color_seqs:
        origin_info = name_to_info.get(seq_name, {})
        origin = origin_info.get('origin', seq_name)
        m = re.search(r'(chr\d+[LR])', origin)
        origin_chr = m.group(1) if m else origin
        source_chr_ends_seen.setdefault(origin_chr, []).append(seq_name)

    source_list = list(source_chr_ends_seen.keys())
    bp_on_read = actual_hit.get('match_start_on_read', bp_pos_on_read)
    hypotheses.append(Hypothesis(
        h_type='y_prime',
        description=f"Y' divergence at position {divergence_idx + 1}: expected {exp_id}, found {actual_cg}",
        confidence=c_for_actual, breakpoint_pos_on_read=bp_on_read,
        source_chr_ends=source_list, ambiguous=(n_same_color > 1),
    ))

    if n_same_color > 1 and len(source_chr_ends_seen) > 1:
        remaining = [s for s in source_chr_ends_seen if s not in source_list[:1]]
        c_others = c_base - c_for_actual
        if remaining and c_others >= MIN_CONFIDENCE_TO_REPORT:
            hypotheses.append(Hypothesis(
                h_type='y_prime',
                description=f"Y' divergence at position {divergence_idx + 1}: alternative sources",
                confidence=c_others, source_chr_ends=remaining, ambiguous=True,
            ))

    if not downstream_consistent and len(all_read_y_primes) > divergence_idx + 1:
        post_cgs = set()
        for j in range(divergence_idx + 1, len(all_read_y_primes)):
            cg_j = name_to_info.get(all_read_y_primes[j]['y_prime_name'], {}).get('color_group', '')
            if cg_j and cg_j != actual_cg:
                post_cgs.add(cg_j)
        if post_cgs:
            hypotheses.append(Hypothesis(
                h_type='compound',
                description=f"Compound: first Y' switch to {actual_cg}, then switch to {';'.join(post_cgs)}",
                confidence=c_for_actual * 0.5, is_compound=True, ambiguous=True,
            ))

    return normalize_hypotheses(hypotheses)


def y_prime_recomb_status(divergence_idx, actual_hit):
    """Map divergence info to a human-readable status string."""
    if divergence_idx is None:
        return 'No Change'
    if actual_hit is None:
        return "Y' Loss"
    if divergence_idx == 0:
        return "1st Y' Change"
    return "Y' Recombination"

# ---------------------------------------------------------------------------
# Post-Y-prime telomere check
# ---------------------------------------------------------------------------

def _telomeric_fraction(sequence):
    """Fraction of bases that are T+G or C+A (whichever is higher)."""
    seq = sequence.upper()
    if not seq:
        return 0.0
    n = len(seq)
    tg = sum(1 for b in seq if b in 'TG')
    ca = sum(1 for b in seq if b in 'CA')
    return max(tg / n, ca / n)


def check_post_y_prime_sequence(read_seq, ordered_hits, telo_side):
    """Check if sequence after last Y prime is telomeric."""
    if not ordered_hits:
        return 0, 0.0, 'n/a'

    telo_hit = ordered_hits[-1]
    if telo_side == 'end':
        tail_seq = read_seq[telo_hit['match_end_on_read']:]
    else:
        tail_seq = read_seq[:telo_hit['match_start_on_read']]

    tail_length = len(tail_seq)
    if tail_length < 50:
        return tail_length, 1.0, 'pass'

    telo_frac = _telomeric_fraction(tail_seq)

    if tail_length < 200 and telo_frac >= 0.55:
        return tail_length, round(telo_frac, 4), 'pass'
    if telo_frac >= 0.65:
        return tail_length, round(telo_frac, 4), 'pass'
    return tail_length, round(telo_frac, 4), 'fail'

# ---------------------------------------------------------------------------
# X element analysis (ported from analyze_x_prime_recombination.py)
# ---------------------------------------------------------------------------

def build_combined_library_fasta(lib_dir, strain, tmp_dir, label):
    """Concatenate all pairing FASTAs into one combined file."""
    combined_path = os.path.join(tmp_dir, f'{strain}_{label}_combined.fasta')
    pattern_dir = lib_dir
    if not os.path.isdir(pattern_dir):
        return None

    with open(combined_path, 'w') as out:
        for fname in sorted(os.listdir(pattern_dir)):
            if fname.endswith('.fasta'):
                fpath = os.path.join(pattern_dir, fname)
                with open(fpath) as inp:
                    out.write(inp.read())

    if os.path.getsize(combined_path) == 0:
        return None
    return combined_path


def subject_to_chr_end(sseqid):
    """Extract chr_end from a BLAST subject ID (e.g., 'chr4R_spacer_...' → 'chr4R')."""
    m = re.search(r'(chr\d+[LR])', sseqid)
    return m.group(1) if m else sseqid


def analyze_x_element(bp_row, diverged_seqs, x_blast_df):
    """Analyze X element for one read. Returns (x_found, x_source, hypotheses)."""
    read_id = bp_row['read_id']
    bp_feature = get_first_breakpoint_feature_type(bp_row)
    if bp_feature not in ('x_variable', 'x_core'):
        return False, '', []

    if x_blast_df.empty or read_id not in diverged_seqs:
        return True, '', [Hypothesis(h_type='ambiguous',
            description='X element breakpoint but no diverged sequence', confidence=0.2, ambiguous=True)]

    hits = x_blast_df[x_blast_df['qseqid'] == read_id]
    if hits.empty:
        return True, '', [Hypothesis(h_type='ambiguous',
            description='X element breakpoint but no BLAST match', confidence=0.2, ambiguous=True)]

    hits = hits.sort_values('bitscore', ascending=False)
    top_score = hits.iloc[0]['bitscore']
    ambig_hits = hits[hits['bitscore'] >= top_score * 0.9]

    sources = {}
    for _, h in ambig_hits.iterrows():
        ce = subject_to_chr_end(h['sseqid'])
        if ce not in sources or h['bitscore'] > sources[ce]['bitscore']:
            sources[ce] = h

    hypotheses = []
    div_len = len(diverged_seqs.get(read_id, ''))
    for ce, h in sources.items():
        c = (h['pident'] / 100.0) * min(1.0, h['length'] / max(div_len, 1)) ** 0.5
        c /= max(1, len(sources))
        hypotheses.append(Hypothesis(
            h_type='x_prime', description=f'X element matches {ce}',
            confidence=c, source_chr_ends=[ce], ambiguous=(len(sources) > 1),
        ))

    return True, subject_to_chr_end(hits.iloc[0]['sseqid']), normalize_hypotheses(hypotheses)

# ---------------------------------------------------------------------------
# Spacer analysis (ported from analyze_spacer_recombination.py)
# ---------------------------------------------------------------------------

def chunk_sequence(seq, chunk_size=SPACER_CHUNK_SIZE, step=SPACER_CHUNK_STEP):
    """Break sequence into overlapping chunks."""
    chunks = []
    for start in range(0, max(1, len(seq) - chunk_size + 1), step):
        chunk = seq[start:start + chunk_size]
        if len(chunk) >= chunk_size // 2:
            chunks.append((chunk, start))
    return chunks


def analyze_spacer(bp_row, diverged_seqs, spacer_blast_df, read_seqs, spacer_db, tmp_dir, features):
    """Analyze spacer for one read. Returns (spacer_found, spacer_source, hypotheses)."""
    read_id = bp_row['read_id']
    bp_feature = get_first_breakpoint_feature_type(bp_row)
    if bp_feature != 'spacer':
        return False, '', []

    hypotheses = []

    # Tail BLAST
    tail_source = ''
    tail_confidence = 0.0
    if not spacer_blast_df.empty:
        hits = spacer_blast_df[spacer_blast_df['qseqid'] == read_id]
        if not hits.empty:
            hits = hits.sort_values('bitscore', ascending=False)
            top = hits.iloc[0]
            tail_source = subject_to_chr_end(top['sseqid'])
            div_len = len(diverged_seqs.get(read_id, ''))
            tail_confidence = (top['pident'] / 100.0) * min(1.0, top['length'] / max(div_len, 1)) ** 0.5

    # Chunk analysis on spacer region
    chunk_source = ''
    chunk_confidence = 0.0
    read_seq = read_seqs.get(read_id, '')
    if read_seq and spacer_db:
        # Extract spacer region: use the aligned portion around the breakpoint
        bp_on_read = get_first_breakpoint_pos_on_read(bp_row)
        if bp_on_read > 0:
            window = 2000
            spacer_start = max(0, bp_on_read - window)
            spacer_end = min(len(read_seq), bp_on_read + window)
            spacer_seq = read_seq[spacer_start:spacer_end]
            if len(spacer_seq) >= SPACER_CHUNK_SIZE:
                chunks = chunk_sequence(spacer_seq)
                if chunks:
                    # Write chunks to temp FASTA and BLAST
                    chunk_fasta = os.path.join(tmp_dir, f'{read_id}_chunks.fasta')
                    with open(chunk_fasta, 'w') as fh:
                        for seq, start in chunks:
                            fh.write(f'>chunk_{start}\n{seq}\n')
                    chunk_hits = run_blast(chunk_fasta, spacer_db, tmp_dir,
                                          label=f'{read_id}_chunks')
                    if not chunk_hits.empty:
                        # Assign each chunk to its top-hit source
                        chunk_sources = {}
                        for _, h in chunk_hits.iterrows():
                            cid = h['qseqid']
                            if cid not in chunk_sources or h['bitscore'] > chunk_sources[cid][1]:
                                chunk_sources[cid] = (subject_to_chr_end(h['sseqid']), h['bitscore'])
                        # Find consensus source
                        from collections import Counter
                        src_counts = Counter(src for src, _ in chunk_sources.values())
                        if src_counts:
                            chunk_source = src_counts.most_common(1)[0][0]
                            chunk_confidence = src_counts[chunk_source] / max(1, len(chunk_sources))

    # Combine tail and chunk evidence
    if tail_source and chunk_source:
        if tail_source == chunk_source:
            combined_c = tail_confidence * 0.7 + chunk_confidence * 0.3
            hypotheses.append(Hypothesis(
                h_type='spacer',
                description=f'Spacer recombination from {tail_source}',
                confidence=combined_c, source_chr_ends=[tail_source],
            ))
        else:
            hypotheses.append(Hypothesis(
                h_type='spacer', description=f'Spacer tail matches {tail_source}',
                confidence=tail_confidence * 0.7, source_chr_ends=[tail_source], ambiguous=True,
            ))
            hypotheses.append(Hypothesis(
                h_type='spacer', description=f'Spacer chunks match {chunk_source}',
                confidence=chunk_confidence * 0.3, source_chr_ends=[chunk_source], ambiguous=True,
            ))
    elif tail_source:
        hypotheses.append(Hypothesis(
            h_type='spacer', description=f'Spacer recombination from {tail_source}',
            confidence=tail_confidence * 0.7, source_chr_ends=[tail_source],
        ))
    elif chunk_source:
        hypotheses.append(Hypothesis(
            h_type='spacer', description=f'Spacer chunks match {chunk_source}',
            confidence=chunk_confidence * 0.3, source_chr_ends=[chunk_source], ambiguous=True,
        ))
    else:
        hypotheses.append(Hypothesis(
            h_type='ambiguous', description='Spacer breakpoint but no match in library',
            confidence=0.2, ambiguous=True,
        ))

    source = tail_source or chunk_source
    return True, source, normalize_hypotheses(hypotheses)

# ---------------------------------------------------------------------------
# Unified hypothesis generation (replaces step 12 aggregation)
# ---------------------------------------------------------------------------

TIER_STRUCTURAL = 0.5
TIER_YPRIME = 0.35
TIER_ELEMENT = 0.15


def build_unified_hypotheses(bp_row, y_prime_hyps, spacer_hyps, x_element_hyps):
    """
    Combine structural evidence (step 10) + Y prime + spacer/X element
    into unified hypotheses.
    """
    alignment_class = str(bp_row.get('alignment_class', ''))

    if alignment_class == 'no_recombination':
        return [Hypothesis(h_type='no_recombination',
            description='No recombination detected', confidence=0.95)]

    if alignment_class == 'no_clean_halves':
        return [Hypothesis(h_type='ambiguous',
            description='No clean alignment — cannot assess recombination',
            confidence=0.1, ambiguous=True)]

    hypotheses = []

    # Structural evidence (two_clean_halves only)
    structural_source = ''
    if alignment_class == 'two_clean_halves':
        structural_source = get_post_breakpoint_contig(bp_row)
        if structural_source:
            # Extract chr_end from contig name
            m = re.search(r'chr(\d+)', structural_source)
            if m:
                structural_source_ce = f'chr{m.group(1)}'
            else:
                structural_source_ce = structural_source

    # Get top hypothesis from each element analysis
    top_yp = y_prime_hyps[0] if y_prime_hyps else None
    top_sp = spacer_hyps[0] if spacer_hyps else None
    top_xe = x_element_hyps[0] if x_element_hyps else None

    # For two_clean_halves: structural + Y prime + element
    if alignment_class == 'two_clean_halves' and structural_source:
        c = TIER_STRUCTURAL
        desc_parts = [f'Structural: {structural_source}']

        if top_yp and top_yp.h_type == 'y_prime':
            c += top_yp.confidence * TIER_YPRIME
            desc_parts.append(f"Y prime: {top_yp.description}")
            # Consistency boost
            if top_yp.source_chr_ends and any(structural_source_ce in s for s in top_yp.source_chr_ends):
                c *= 1.3

        if top_sp and top_sp.h_type == 'spacer':
            c += top_sp.confidence * TIER_ELEMENT
            desc_parts.append(f"Spacer: {top_sp.description}")
        elif top_xe and top_xe.h_type == 'x_prime':
            c += top_xe.confidence * TIER_ELEMENT
            desc_parts.append(f"X element: {top_xe.description}")

        source = [structural_source_ce] if structural_source else []
        hypotheses.append(Hypothesis(
            h_type='recombination', description='; '.join(desc_parts),
            confidence=min(1.0, c), source_chr_ends=source,
        ))

    # For one_clean_half: Y prime + element (no structural)
    elif alignment_class == 'one_clean_half':
        if top_yp and top_yp.h_type in ('y_prime', 'no_y_prime_change'):
            c = top_yp.confidence * 0.7
            hypotheses.append(Hypothesis(
                h_type=top_yp.h_type, description=top_yp.description,
                confidence=c, source_chr_ends=top_yp.source_chr_ends,
                ambiguous=top_yp.ambiguous,
            ))

        if top_sp and top_sp.h_type == 'spacer':
            c = top_sp.confidence * 0.3
            hypotheses.append(Hypothesis(
                h_type='spacer', description=top_sp.description,
                confidence=c, source_chr_ends=top_sp.source_chr_ends,
                ambiguous=top_sp.ambiguous,
            ))
        elif top_xe and top_xe.h_type == 'x_prime':
            c = top_xe.confidence * 0.3
            hypotheses.append(Hypothesis(
                h_type='x_prime', description=top_xe.description,
                confidence=c, source_chr_ends=top_xe.source_chr_ends,
                ambiguous=top_xe.ambiguous,
            ))

        if not hypotheses:
            hypotheses.append(Hypothesis(
                h_type='ambiguous', description='One clean half but no source identified',
                confidence=0.2, ambiguous=True,
            ))

    return normalize_hypotheses(hypotheses)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    print(f'analyze_features.py — chr_end={args.chr_end}')

    for path in [args.breakpoints_tsv, args.reads_fasta, args.day0_bed, args.y_prime_lib]:
        if not os.path.exists(path):
            print(f'ERROR: missing input: {path}', file=sys.stderr)
            sys.exit(1)

    # Load breakpoints from step 10
    df_bp = pd.read_csv(args.breakpoints_tsv, sep='\t')
    if df_bp.empty:
        print('  No reads — writing empty output')
        write_results_tsv([], args.output_tsv)
        return

    # Load BED features
    features = load_bed_features(args.day0_bed, args.chr_end)
    ref_y_primes = get_reference_y_prime_order(features)
    print(f'  BED: {len(features)} features, {len(ref_y_primes)} Y primes for {args.chr_end}')

    # Load read sequences
    print('  Loading read sequences…')
    read_seqs = load_read_sequences(args.reads_fasta)
    print(f'  {len(read_seqs)} reads loaded')

    # Load diverged tails
    diverged_seqs = {}
    if os.path.exists(args.diverged_fasta) and os.path.getsize(args.diverged_fasta) > 0:
        raw = read_fasta(args.diverged_fasta)
        for header, seq in raw.items():
            rid = header.split()[0]
            diverged_seqs[rid] = seq

    # Build Y prime color map
    color_map, name_to_info = build_y_prime_color_map(args.y_prime_lib)
    print(f'  Y prime library: {len(name_to_info)} sequences, {len(color_map)} color groups')

    # === Batch preprocessing ===
    with tempfile.TemporaryDirectory() as tmp_dir:
        # 1. BLAST full reads against Y prime library
        df_rm = blast_y_primes_on_full_reads(
            args.reads_fasta, args.y_prime_lib, tmp_dir, args.threads)

        # 2. BLAST diverged tails against spacer library
        spacer_db = build_combined_library_fasta(
            args.spacer_lib_dir, args.strain, tmp_dir, 'spacer')
        spacer_blast_df = pd.DataFrame()
        if spacer_db and diverged_seqs:
            print('    BLASTing diverged tails against spacer library…')
            spacer_blast_df = run_blast(
                args.diverged_fasta, spacer_db, tmp_dir, label='spacer_tail')
            print(f'    Spacer BLAST: {len(spacer_blast_df)} hits')

        # 3. BLAST diverged tails against X element library
        x_element_db = build_combined_library_fasta(
            args.x_element_lib_dir, args.strain, tmp_dir, 'x_element')
        x_blast_df = pd.DataFrame()
        if x_element_db and diverged_seqs:
            print('    BLASTing diverged tails against X element library…')
            x_blast_df = run_blast(
                args.diverged_fasta, x_element_db, tmp_dir, label='x_element_tail')
            print(f'    X element BLAST: {len(x_blast_df)} hits')

        # === Per-read structural walk ===
        print('  Walking features per read…')
        rows = []
        for _, bp_row in df_bp.iterrows():
            read_id = bp_row['read_id']
            telo_side = bp_row.get('telo_side', 'end')
            alignment_class = str(bp_row.get('alignment_class', ''))

            # 1. Determine feature reach
            feature_reach, features_in_order, order_issues = determine_feature_reach(bp_row, features)

            # 2. Y prime analysis (if read has Y prime hits)
            rm_for_read = pd.DataFrame()
            if not df_rm.empty and 'read_id' in df_rm.columns:
                rm_for_read = df_rm[df_rm['read_id'] == read_id]

            ordered_hits = order_y_primes_for_read(rm_for_read, telo_side)
            y_prime_count = len(ordered_hits)

            divergence_idx, actual_hit, expected_ref, downstream_consistent = find_first_divergence(
                ordered_hits, ref_y_primes, name_to_info)

            y_prime_hyps = generate_y_prime_hypotheses(
                read_id, bp_row, divergence_idx, actual_hit, expected_ref,
                downstream_consistent, color_map, name_to_info, ordered_hits)

            # 3. Post-Y-prime telomere check
            read_seq = read_seqs.get(read_id, '')
            tail_len, telo_frac, telo_check = check_post_y_prime_sequence(
                read_seq, ordered_hits, telo_side)

            # 4. Spacer analysis (BLAST-based, only for reads with spacer breakpoints)
            spacer_found, spacer_source, spacer_hyps = analyze_spacer(
                bp_row, diverged_seqs, spacer_blast_df, read_seqs,
                spacer_db, tmp_dir, features)

            # 5. X element analysis (BLAST-based, only for reads with X element breakpoints)
            x_found, x_source, x_hyps = analyze_x_element(
                bp_row, diverged_seqs, x_blast_df)

            # 5b. Override found flags based on alignment coverage —
            # if the aligned region covers a feature's BED coordinates,
            # mark it as found even if BLAST wasn't run for that feature.
            if any('spacer' in f for f in feature_reach):
                spacer_found = True
            if any('x_core' in f or 'x_variable' in f for f in feature_reach):
                x_found = True

            # 6. Build unified hypotheses
            unified_hyps = build_unified_hypotheses(bp_row, y_prime_hyps, spacer_hyps, x_hyps)

            # Safety net: ensure at least one hypothesis exists
            if not unified_hyps:
                unified_hyps = [Hypothesis(
                    h_type='ambiguous',
                    description='Unable to determine recombination status',
                    confidence=0.1, ambiguous=True,
                )]

            # Build output row
            row = {
                'read_id': read_id,
                'chr_end': args.chr_end,
                'read_length': bp_row.get('read_length', 0),
                'telo_side': telo_side,
                'alignment_class': alignment_class,

                # Structural walk
                'feature_reach': ','.join(feature_reach),
                'features_in_order': features_in_order,
                'feature_order_issues': order_issues,

                # Spacer
                'spacer_found': spacer_found,
                'spacer_source': spacer_source,

                # X element
                'x_element_found': x_found,
                'x_element_source': x_source,

                # Y prime
                'y_prime_count_on_read': y_prime_count,
                'y_prime_divergence_idx': divergence_idx if divergence_idx is not None else -1,
                'y_prime_expected_at_divergence': expected_ref['id'] if expected_ref else '',
                'y_prime_found_at_divergence': (
                    name_to_info.get(actual_hit['y_prime_name'], {}).get('color_group', '')
                    if actual_hit else ''),
                'y_prime_downstream_consistent': downstream_consistent,
                'y_prime_recombination_status': y_prime_recomb_status(divergence_idx, actual_hit),

                # Post-Y-prime telomere check
                'post_y_prime_seq_length': tail_len,
                'post_y_prime_telo_fraction': telo_frac,
                'post_y_prime_check': telo_check,

                # Unified hypothesis
                'n_hypotheses': len(unified_hyps),
            }
            row.update(hypotheses_to_row_dict(unified_hyps))
            rows.append(row)

    # Write output
    os.makedirs(os.path.dirname(args.output_tsv) or '.', exist_ok=True)
    write_results_tsv(rows, args.output_tsv)

    # Summary
    from collections import Counter
    class_counts = Counter(r['alignment_class'] for r in rows)
    yp_status_counts = Counter(r['y_prime_recombination_status'] for r in rows)
    print(f'  Output: {args.output_tsv} ({len(rows)} reads)')
    print(f'  Alignment classes: {dict(class_counts)}')
    print(f'  Y prime status: {dict(yp_status_counts)}')


if __name__ == '__main__':
    main()
