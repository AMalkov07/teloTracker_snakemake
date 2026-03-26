"""
Step 11: Combined feature analysis — chunk-based recombination detection.

For each read:
  1. Spacer chunk analysis:  250bp chunks BLASTed against combined spacer library
  2. X element chunk analysis: 250bp chunks BLASTed against combined X element library
  3. Y prime analysis: BLAST full read against Y prime library, position-by-position comparison
  4. Cross-feature reconciliation + confidence scoring

All analyses run on EVERY read — no gating by alignment classification or Y prime status.

Usage:
  python analyze_features.py \
      --reads-fasta       <per-chr-end reads FASTA> \
      --alignment-tsv     <step 10 supplementary alignment TSV> \
      --day0-bed          <day0 BED features file> \
      --y-prime-lib       references/extracted_yprimes_{strain}.fasta \
      --spacer-lib-dir    references/pairings_for_spacers/{strain}_pairings/ \
      --x-element-lib-dir references/pairings_for_x_element_ends/{strain}_pairings/ \
      --chr-end           chr4R  --strain 7302 \
      --output-tsv        results/{base}/recombination/{base}_{chr_end}_features.tsv \
      --threads           4
"""

import argparse
import os
import re
import sys
import tempfile

# Ensure scripts/ is on the import path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

from recombination_utils import (
    load_bed_features,
    telo_side_from_header,
    run_blast,
    write_results_tsv,
    BLAST_MIN_PIDENT,
    SPACER_CHUNK_SIZE,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CHUNK_SIZE = SPACER_CHUNK_SIZE  # 250bp, non-overlapping

# Feature distinctiveness weights (Factor 1 of confidence scoring)
DISTINCTIVENESS = {
    'spacer': 0.9,      # spacers are generally unique between chr ends
    'y_prime': 0.6,     # distinct between IDs but ID assignment can be unreliable
    'x_element': 0.4,   # often similar across chr ends
}

# Confidence penalties
COMPLEXITY_PENALTY = 0.3       # multiplied when cross-feature results are inconsistent
PARTIAL_AGREEMENT_FACTOR = 0.7 # when features partially agree

MIN_YPRIME_HIT_LEN = 500      # minimum BLAST hit length for Y prime
MAX_YPRIME_MERGE_GAP = 500    # max gap to merge adjacent Y prime hits

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description='Combined feature analysis — chunk-based')
    p.add_argument('--reads-fasta',    required=True)
    p.add_argument('--alignment-tsv',  required=True, help='Step 10 supplementary alignment TSV')
    p.add_argument('--day0-bed',       required=True)
    p.add_argument('--day0-ref',       default='', help='Day0 reference FASTA (for spacer quick check)')
    p.add_argument('--y-prime-lib',    required=True)
    p.add_argument('--spacer-lib-dir',     required=True, help='Directory with spacer pairing FASTAs')
    p.add_argument('--x-element-lib-dir',  required=True, help='Directory with X element pairing FASTAs')
    p.add_argument('--chr-end',        required=True)
    p.add_argument('--strain',         required=True)
    p.add_argument('--output-tsv',     required=True)
    p.add_argument('--threads',        type=int, default=4)
    return p.parse_args()

# ---------------------------------------------------------------------------
# FASTA loading
# ---------------------------------------------------------------------------

def load_reads(fasta_path):
    """Return (dict {read_id: sequence}, dict {read_id: full_header})."""
    seqs, headers = {}, {}
    cur_id, cur_seq = None, []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if cur_id is not None:
                    seqs[cur_id] = ''.join(cur_seq)
                hdr = line[1:]
                cur_id = hdr.split()[0]
                headers[cur_id] = hdr
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_id is not None:
        seqs[cur_id] = ''.join(cur_seq)
    return seqs, headers

# ---------------------------------------------------------------------------
# Chunking + batch BLAST
# ---------------------------------------------------------------------------

def build_combined_library(lib_dir, tmp_dir, label):
    """Concatenate all FASTA files in a directory into one deduplicated combined library.
    Removes duplicate sequences (same header) that appear across multiple pairing files."""
    if not os.path.isdir(lib_dir):
        return None
    combined_path = os.path.join(tmp_dir, f'{label}_combined.fasta')
    seen_headers = set()
    with open(combined_path, 'w') as out:
        for fname in sorted(os.listdir(lib_dir)):
            if not fname.endswith('.fasta'):
                continue
            fpath = os.path.join(lib_dir, fname)
            writing = False
            with open(fpath) as inp:
                for line in inp:
                    if line.startswith('>'):
                        header = line.strip()
                        if header in seen_headers:
                            writing = False
                        else:
                            seen_headers.add(header)
                            writing = True
                            out.write(line)
                    elif writing:
                        out.write(line)
    if os.path.getsize(combined_path) == 0:
        return None
    print(f'  {label} combined library: {len(seen_headers)} unique sequences (deduplicated)')
    return combined_path


def chunk_reads(read_seqs, chunk_size=CHUNK_SIZE):
    """Chunk all reads into non-overlapping pieces.
    Returns list of (chunk_id, sequence) where chunk_id = '{read_id}__chunk_{i}'."""
    chunks = []
    for read_id, seq in read_seqs.items():
        for i in range(0, len(seq), chunk_size):
            chunk_seq = seq[i:i + chunk_size]
            if len(chunk_seq) >= chunk_size // 2:  # keep chunks at least half-size
                chunk_id = f'{read_id}__chunk_{i}'
                chunks.append((chunk_id, chunk_seq))
    return chunks


BLAST_BATCH_SIZE = 2000  # max chunks per BLAST call to avoid OOM

def batch_blast_chunks(chunks, library_fasta, tmp_dir, label):
    """Write chunks to temp FASTA, BLAST against library in batches, return DataFrame."""
    if not chunks or not os.path.exists(library_fasta) or os.path.getsize(library_fasta) == 0:
        return pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'bitscore', 'length'])

    all_dfs = []
    for batch_idx in range(0, len(chunks), BLAST_BATCH_SIZE):
        batch = chunks[batch_idx:batch_idx + BLAST_BATCH_SIZE]
        batch_label = f'{label}_batch{batch_idx}'
        chunk_fasta = os.path.join(tmp_dir, f'{batch_label}_chunks.fasta')
        with open(chunk_fasta, 'w') as fh:
            for chunk_id, seq in batch:
                fh.write(f'>{chunk_id}\n{seq}\n')

        df = run_blast(chunk_fasta, library_fasta, tmp_dir, label=batch_label)
        if not df.empty:
            all_dfs.append(df)

    if all_dfs:
        return pd.concat(all_dfs, ignore_index=True)
    return pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'bitscore', 'length'])


def parse_chunk_results(blast_df):
    """Parse batch BLAST results into per-read, per-chunk best hits.

    Returns dict: {read_id: [(chunk_start, info_dict), ...]} sorted by chunk_start.
    """
    if blast_df.empty:
        return {}

    # First pass: collect best hit per (read, chunk)
    best_hits = {}
    second_hits = {}
    for _, row in blast_df.iterrows():
        qid = row['qseqid']
        parts = qid.rsplit('__chunk_', 1)
        if len(parts) != 2:
            continue
        read_id = parts[0]
        chunk_start = int(parts[1])
        source = _subject_to_chr_end(row['sseqid'])
        key = (read_id, chunk_start)

        if key not in best_hits or row['bitscore'] > best_hits[key]['bitscore']:
            # Demote current best to second
            if key in best_hits and best_hits[key]['source'] != source:
                second_hits[key] = best_hits[key]
            best_hits[key] = {
                'source': source,
                'pident': row['pident'],
                'bitscore': row['bitscore'],
            }
        elif key not in second_hits or row['bitscore'] > second_hits[key]['bitscore']:
            if source != best_hits[key]['source']:
                second_hits[key] = {
                    'source': source,
                    'pident': row['pident'],
                    'bitscore': row['bitscore'],
                }

    # Merge best + second into final structure
    results = {}
    for (read_id, chunk_start), info in best_hits.items():
        sec = second_hits.get((read_id, chunk_start), {})
        info['second_source'] = sec.get('source', '')
        info['second_pident'] = sec.get('pident', 0.0)
        info['second_bitscore'] = sec.get('bitscore', 0)
        results.setdefault(read_id, []).append((chunk_start, info))

    # Sort by chunk_start
    for read_id in results:
        results[read_id].sort(key=lambda x: x[0])

    return results


def _subject_to_chr_end(sseqid):
    """Extract chr end from BLAST subject ID."""
    m = re.search(r'(chr\d+[LR])', sseqid)
    return m.group(1) if m else sseqid

# ---------------------------------------------------------------------------
# Spacer quick check (exact substring match)
# ---------------------------------------------------------------------------

def load_reference_spacer(bed_features, day0_ref_fasta):
    """Extract the reference spacer sequence for this chr end from the day0 reference.
    Returns the spacer sequence string, or '' if not found."""
    spacer_feats = [f for f in bed_features if f['feature_type'] == 'spacer']
    if not spacer_feats:
        return ''
    # Use the first (should be only) spacer feature
    spacer = spacer_feats[0]
    # Read the reference FASTA to extract the spacer region
    # We need the contig sequence from spacer['start'] to spacer['end']
    try:
        from recombination_utils import read_fasta
        ref_seqs = read_fasta(day0_ref_fasta)
        contig = spacer['contig']
        for header, seq in ref_seqs.items():
            if contig in header or header.startswith(contig):
                return seq[spacer['start']:spacer['end']].upper()
    except Exception:
        pass
    return ''


def spacer_quick_check(read_seq, reference_spacer):
    """Check if the reference spacer exists as an exact substring in the read.
    Returns True if found (100% identity, 100% length match)."""
    if not reference_spacer or len(reference_spacer) < 100:
        return False
    return reference_spacer in read_seq.upper()


# ---------------------------------------------------------------------------
# Spacer / X element analysis (chunk walk)
# ---------------------------------------------------------------------------

def analyze_chunks(read_id, chunk_hits, expected_chr_end, feature_name):
    """Walk chunks for one read, detect source switches.

    Returns dict with analysis results.
    """
    if not chunk_hits:
        return {
            f'{feature_name}_chunks_analyzed': 0,
            f'{feature_name}_source': '',
            f'{feature_name}_switch_chunk': -1,
            f'{feature_name}_best_identity': 0.0,
            f'{feature_name}_second_best_identity': 0.0,
            f'{feature_name}_confidence': 0.0,
            f'{feature_name}_recombination': 'no_data',
        }

    # Count chunks matching each source
    source_counts = {}
    source_identities = {}
    for pos, info in chunk_hits:
        src = info['source']
        source_counts[src] = source_counts.get(src, 0) + 1
        source_identities.setdefault(src, []).append(info['pident'])

    n_chunks = len(chunk_hits)
    best_source = max(source_counts, key=source_counts.get)
    best_avg_identity = sum(source_identities[best_source]) / len(source_identities[best_source])

    # Second-best source
    second_source = ''
    second_avg_identity = 0.0
    remaining = {k: v for k, v in source_counts.items() if k != best_source}
    if remaining:
        second_source = max(remaining, key=remaining.get)
        second_avg_identity = sum(source_identities[second_source]) / len(source_identities[second_source])

    # Detect switch point: walk chunks and find where source changes
    switch_chunk = -1
    sources_in_order = [(pos, info['source']) for pos, info in chunk_hits]
    if len(set(s for _, s in sources_in_order)) > 1:
        first_source = sources_in_order[0][1]
        for i, (pos, src) in enumerate(sources_in_order):
            if src != first_source:
                switch_chunk = i
                break

    # Determine recombination status
    all_match_expected = all(info['source'] == expected_chr_end for _, info in chunk_hits)
    if all_match_expected:
        recomb_status = 'no_change'
    elif switch_chunk >= 0:
        recomb_status = 'switch_detected'
    elif best_source != expected_chr_end:
        recomb_status = 'full_switch'
    else:
        recomb_status = 'partial_match'

    # Confidence: Factor 2 (separation between best and second-best)
    if second_avg_identity > 0:
        separation = (best_avg_identity - second_avg_identity) / max(best_avg_identity, 1)
    else:
        separation = 1.0

    distinctiveness = DISTINCTIVENESS.get(feature_name, 0.5)
    confidence = distinctiveness * separation

    if recomb_status == 'no_change':
        confidence = 0.95

    return {
        f'{feature_name}_chunks_analyzed': n_chunks,
        f'{feature_name}_source': best_source if recomb_status != 'no_change' else expected_chr_end,
        f'{feature_name}_switch_chunk': switch_chunk,
        f'{feature_name}_best_identity': round(best_avg_identity, 2),
        f'{feature_name}_second_best_identity': round(second_avg_identity, 2),
        f'{feature_name}_confidence': round(confidence, 4),
        f'{feature_name}_recombination': recomb_status,
    }

# ---------------------------------------------------------------------------
# Y prime analysis
# ---------------------------------------------------------------------------

def parse_y_prime_header(header):
    """Parse Y prime FASTA header like >Y_Prime_chr2L1#Short/Solo/ID4_Green-Light."""
    header = header.lstrip('>')
    name_part, class_part = (header.split('#', 1) + [''])[:2]
    y_id, color_group = '', ''
    if class_part:
        parts = class_part.split('/')
        if len(parts) >= 3:
            color_group = parts[2]
            y_id = color_group.split('_', 1)[0]
    return {'id': y_id, 'color_group': color_group, 'origin': name_part, 'full_name': header.split()[0]}


def build_y_prime_info(y_prime_lib_fasta):
    """Parse Y prime library. Returns {seq_name: info_dict}."""
    name_to_info = {}
    with open(y_prime_lib_fasta) as fh:
        for line in fh:
            if line.strip().startswith('>'):
                info = parse_y_prime_header(line.strip())
                name_to_info[info['full_name']] = info
    return name_to_info


def get_reference_y_prime_order(features, name_to_info):
    """Get day0 Y prime order for this chr end from BED + library info."""
    yp_features = [f for f in features if f['feature_type'] == 'y_prime']
    yp_features.sort(key=lambda f: f['start'])

    result = []
    for f in yp_features:
        y_id = ''
        for seq_name, info in name_to_info.items():
            feat_name = f['name']
            origin = info.get('origin', '')
            if feat_name.replace('_Y_Prime_', '') in origin or origin in feat_name:
                y_id = info.get('id', '')
                break
        if not y_id:
            id_match = re.search(r'(ID\d+)', f['name'])
            y_id = id_match.group(1) if id_match else f['name']
        result.append({'feature_name': f['name'], 'id': y_id, 'start': f['start'], 'end': f['end']})
    return result


def blast_y_primes(read_seqs, y_prime_lib, tmp_dir):
    """BLAST all reads against Y prime library. Returns per-read Y prime hits."""
    reads_fasta = os.path.join(tmp_dir, 'all_reads.fasta')
    with open(reads_fasta, 'w') as fh:
        for read_id, seq in read_seqs.items():
            fh.write(f'>{read_id}\n{seq}\n')

    blast_df = run_blast(reads_fasta, y_prime_lib, tmp_dir, label='y_prime',
                         min_pident=80.0, evalue=1e-10)

    if blast_df.empty:
        return {}

    # Convert to internal format
    rows = []
    for _, h in blast_df.iterrows():
        rows.append({
            'read_id': h['qseqid'],
            'y_prime_name': h['sseqid'],
            'match_start': min(int(h['qstart']), int(h['qend'])),
            'match_end': max(int(h['qstart']), int(h['qend'])),
            'bitscore': int(h['bitscore']),
            'pident': h['pident'],
        })

    df = pd.DataFrame(rows)

    # Merge adjacent hits from same Y prime
    merged = []
    for (rid, yname), grp in df.groupby(['read_id', 'y_prime_name']):
        grp = grp.sort_values('match_start')
        current = grp.iloc[0].to_dict()
        for i in range(1, len(grp)):
            row = grp.iloc[i]
            if row['match_start'] - current['match_end'] <= MAX_YPRIME_MERGE_GAP:
                current['match_end'] = row['match_end']
                current['bitscore'] += row['bitscore']
            else:
                merged.append(current)
                current = row.to_dict()
        merged.append(current)

    # Filter short hits
    merged = [m for m in merged if m['match_end'] - m['match_start'] >= MIN_YPRIME_HIT_LEN]

    # Deduplicate overlapping hits (keep best score)
    per_read = {}
    for hit in merged:
        per_read.setdefault(hit['read_id'], []).append(hit)

    deduped = {}
    for rid, hits in per_read.items():
        hits.sort(key=lambda h: h['bitscore'], reverse=True)
        kept = []
        for hit in hits:
            overlaps = False
            for k in kept:
                ov_start = max(hit['match_start'], k['match_start'])
                ov_end = min(hit['match_end'], k['match_end'])
                if ov_end > ov_start and (ov_end - ov_start) > 0.5 * (hit['match_end'] - hit['match_start']):
                    overlaps = True
                    break
            if not overlaps:
                kept.append(hit)
        deduped[rid] = kept

    return deduped


def analyze_y_primes(read_id, y_prime_hits, telo_side, ref_y_primes, name_to_info):
    """Analyze Y primes for one read."""
    hits = y_prime_hits or []
    hits.sort(key=lambda h: h['match_start'])
    if telo_side == 'beginning':
        hits = hits[::-1]

    # Observed Y prime array (anchor-to-telomere)
    observed_array = []
    for hit in hits:
        info = name_to_info.get(hit['y_prime_name'], {})
        observed_array.append(info.get('id', hit['y_prime_name']))

    # Position-by-position comparison with reference
    divergence_idx = -1
    expected_at_div = ''
    found_at_div = ''
    downstream_consistent = True
    status = 'No Change'

    if len(hits) < len(ref_y_primes):
        divergence_idx = len(hits)
        if divergence_idx < len(ref_y_primes):
            expected_at_div = ref_y_primes[divergence_idx]['id']
        status = "Y' Loss"
    else:
        for i in range(len(ref_y_primes)):
            if i >= len(observed_array):
                break
            if observed_array[i] != ref_y_primes[i]['id']:
                divergence_idx = i
                expected_at_div = ref_y_primes[i]['id']
                found_at_div = observed_array[i]
                if i + 1 < len(observed_array):
                    first_new = observed_array[i]
                    for j in range(i + 1, len(observed_array)):
                        if observed_array[j] != first_new:
                            downstream_consistent = False
                            break
                status = "1st Y' Change" if i == 0 else "Y' Recombination"
                break
        if status == 'No Change' and len(observed_array) > len(ref_y_primes):
            status = "Y' Gain"

    # Find compatible reference ends
    compatible_ends = find_compatible_ends(observed_array, name_to_info)

    return {
        'y_prime_count_on_read': len(hits),
        'y_prime_observed_array': ','.join(observed_array) if observed_array else '',
        'y_prime_recombination_status': status,
        'y_prime_divergence_idx': divergence_idx,
        'y_prime_expected_at_divergence': expected_at_div,
        'y_prime_found_at_divergence': found_at_div,
        'y_prime_downstream_consistent': downstream_consistent,
        'y_prime_compatible_ends': ','.join(compatible_ends),
    }


def find_compatible_ends(observed_array, name_to_info):
    """Find reference chr ends whose first Y prime ID matches the observed first Y prime."""
    if not observed_array:
        return []

    first_observed = observed_array[0]

    # Build per-chr-end first Y prime from library
    chr_end_first_yp = {}
    for seq_name, info in name_to_info.items():
        origin = info.get('origin', '')
        m = re.search(r'(chr\d+[LR])', origin)
        if not m:
            continue
        ce = m.group(1)
        pos_match = re.search(r'chr\d+[LR](\d+)', origin)
        pos = int(pos_match.group(1)) if pos_match else 0
        if ce not in chr_end_first_yp or pos < chr_end_first_yp[ce][0]:
            chr_end_first_yp[ce] = (pos, info.get('id', ''))

    compatible = [ce for ce, (_, yid) in chr_end_first_yp.items() if yid == first_observed]
    return sorted(compatible)

# ---------------------------------------------------------------------------
# Cross-feature reconciliation + confidence scoring
# ---------------------------------------------------------------------------

def reconcile_features(spacer_result, x_element_result, y_prime_result, supp_contigs, chr_end):
    """Compare results across features, compute overall confidence."""
    spacer_source = spacer_result.get('spacer_source', '')
    x_source = x_element_result.get('x_element_source', '')
    spacer_recomb = spacer_result.get('spacer_recombination', 'no_data')
    x_recomb = x_element_result.get('x_element_recombination', 'no_data')
    y_status = y_prime_result.get('y_prime_recombination_status', '')
    y_compatible = [e for e in y_prime_result.get('y_prime_compatible_ends', '').split(',') if e]

    # Parse supplementary contigs
    supp_chr_ends = []
    for contig in supp_contigs:
        m = re.search(r'chr(\d+)', contig)
        if m:
            supp_chr_ends.append(f'chr{m.group(1)}')

    spacer_has_recomb = spacer_recomb in ('switch_detected', 'full_switch')
    x_has_recomb = x_recomb in ('switch_detected', 'full_switch')
    y_has_recomb = y_status not in ('No Change', '', 'no_data')

    if not spacer_has_recomb and not x_has_recomb and not y_has_recomb:
        return {
            'recombination_source': '',
            'recombination_detected': False,
            'overall_confidence': 0.95,
            'cross_feature_consistent': True,
            'cross_feature_detail': 'all_features_match_reference',
            'is_complex_event': False,
            'confidence_factors': 'no_recombination',
        }

    # Collect source candidates
    sources = {}
    if spacer_has_recomb and spacer_source and spacer_source != chr_end:
        sources['spacer'] = spacer_source
    if x_has_recomb and x_source and x_source != chr_end:
        sources['x_element'] = x_source
    if supp_chr_ends:
        sources['supplementary'] = supp_chr_ends[0]

    unique_sources = set(sources.values())
    all_agree = len(unique_sources) <= 1

    # Y prime support check
    y_supports = False
    if unique_sources and y_compatible:
        for src in unique_sources:
            if any(src in ce for ce in y_compatible):
                y_supports = True
                break

    # Detail string
    detail_parts = []
    for feat, src in sources.items():
        detail_parts.append(f'{feat}={src}')
    if y_compatible:
        detail_parts.append(f'y_prime_compatible={",".join(y_compatible[:5])}')
    if y_has_recomb:
        detail_parts.append(f'y_prime_status={y_status}')

    # Confidence scoring
    feature_confidences = []
    if spacer_has_recomb:
        feature_confidences.append(spacer_result.get('spacer_confidence', 0))
    if x_has_recomb:
        feature_confidences.append(x_element_result.get('x_element_confidence', 0))

    base_confidence = max(feature_confidences) if feature_confidences else 0.3

    if all_agree:
        complexity_factor = 1.1 if y_supports else 1.0
    elif len(unique_sources) == 2:
        complexity_factor = PARTIAL_AGREEMENT_FACTOR
    else:
        complexity_factor = COMPLEXITY_PENALTY

    overall_confidence = min(1.0, base_confidence * complexity_factor)

    # Best source
    if unique_sources:
        from collections import Counter
        source_votes = Counter(sources.values())
        best_source = source_votes.most_common(1)[0][0]
    else:
        best_source = 'ambiguous'

    is_complex = not all_agree and len(unique_sources) > 1

    return {
        'recombination_source': best_source,
        'recombination_detected': True,
        'overall_confidence': round(overall_confidence, 4),
        'cross_feature_consistent': all_agree,
        'cross_feature_detail': '; '.join(detail_parts),
        'is_complex_event': is_complex,
        'confidence_factors': f'base={base_confidence:.2f}, complexity={complexity_factor:.2f}',
    }

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    print(f'analyze_features.py -- chr_end={args.chr_end}')

    for path in [args.reads_fasta, args.day0_bed, args.y_prime_lib]:
        if not os.path.exists(path):
            print(f'ERROR: missing input: {path}', file=sys.stderr)
            sys.exit(1)

    # Load reads
    print('  Loading reads...')
    read_seqs, read_headers = load_reads(args.reads_fasta)
    n_reads = len(read_seqs)
    print(f'  {n_reads} reads loaded')

    if n_reads == 0:
        print('  No reads -- writing empty output')
        write_results_tsv([], args.output_tsv)
        return

    # Load BED features and Y prime info
    features = load_bed_features(args.day0_bed, args.chr_end)
    name_to_info = build_y_prime_info(args.y_prime_lib)
    ref_y_primes = get_reference_y_prime_order(features, name_to_info)
    print(f'  BED: {len(features)} features, {len(ref_y_primes)} Y primes for {args.chr_end}')

    # Load supplementary alignment info from step 10
    supp_info = {}
    if os.path.exists(args.alignment_tsv) and os.path.getsize(args.alignment_tsv) > 0:
        df_aln = pd.read_csv(args.alignment_tsv, sep='\t')
        for _, row in df_aln.iterrows():
            rid = row.get('read_id', '')
            contigs = str(row.get('supplementary_contigs', ''))
            supp_info[rid] = contigs.split(';') if contigs and contigs != 'nan' else []

    with tempfile.TemporaryDirectory() as tmp_dir:
        # === Batch preprocessing ===
        print('  Chunking reads...')
        all_chunks = chunk_reads(read_seqs)
        print(f'  {len(all_chunks)} chunks from {n_reads} reads')

        # Build combined libraries from directories
        spacer_lib = build_combined_library(args.spacer_lib_dir, tmp_dir, 'spacer')
        x_element_lib = build_combined_library(args.x_element_lib_dir, tmp_dir, 'x_element')

        # 1. Spacer chunk BLAST
        spacer_hits = {}
        if spacer_lib:
            print('  BLASTing chunks against spacer library...')
            spacer_blast = batch_blast_chunks(all_chunks, spacer_lib, tmp_dir, 'spacer')
            spacer_hits = parse_chunk_results(spacer_blast)
            print(f'  Spacer: {len(spacer_blast)} hits, {len(spacer_hits)} reads with matches')
        else:
            print('  No spacer library found -- skipping spacer analysis')

        # 2. X element chunk BLAST
        x_hits = {}
        if x_element_lib:
            print('  BLASTing chunks against X element library...')
            x_blast = batch_blast_chunks(all_chunks, x_element_lib, tmp_dir, 'x_element')
            x_hits = parse_chunk_results(x_blast)
            print(f'  X element: {len(x_blast)} hits, {len(x_hits)} reads with matches')
        else:
            print('  No X element library found -- skipping X element analysis')

        # 3. Y prime BLAST
        print('  BLASTing reads against Y prime library...')
        y_prime_hits = blast_y_primes(read_seqs, args.y_prime_lib, tmp_dir)
        n_with_yp = sum(1 for hits in y_prime_hits.values() if hits)
        print(f'  Y prime: {n_with_yp} reads with hits')

    # Load reference spacer for quick check
    ref_spacer_seq = ''
    if args.day0_ref and os.path.exists(args.day0_ref):
        ref_spacer_seq = load_reference_spacer(features, args.day0_ref)
        if ref_spacer_seq:
            print(f'  Reference spacer loaded: {len(ref_spacer_seq)}bp')
        else:
            print('  No spacer feature in BED -- quick check disabled')
    else:
        print('  No day0 reference provided -- spacer quick check disabled')

    # === Per-read analysis ===
    print('  Analyzing features per read...')
    rows = []
    n_quick_check_passed = 0
    for read_id in read_seqs:
        telo_side = telo_side_from_header(read_headers.get(read_id, ''), args.chr_end)

        # Spacer quick check: exact substring match
        quick_check_skipped = False
        if ref_spacer_seq and spacer_quick_check(read_seqs[read_id], ref_spacer_seq):
            quick_check_skipped = True
            n_quick_check_passed += 1
            spacer_result = {
                'spacer_chunks_analyzed': 0,
                'spacer_source': args.chr_end,
                'spacer_switch_chunk': -1,
                'spacer_best_identity': 100.0,
                'spacer_second_best_identity': 0.0,
                'spacer_confidence': 0.95,
                'spacer_recombination': 'no_change',
                'spacer_quick_check_skipped': True,
            }
        else:
            spacer_result = analyze_chunks(read_id, spacer_hits.get(read_id, []), args.chr_end, 'spacer')
            spacer_result['spacer_quick_check_skipped'] = False
        x_result = analyze_chunks(read_id, x_hits.get(read_id, []), args.chr_end, 'x_element')
        y_result = analyze_y_primes(read_id, y_prime_hits.get(read_id, []), telo_side, ref_y_primes, name_to_info)

        reconciliation = reconcile_features(
            spacer_result, x_result, y_result, supp_info.get(read_id, []), args.chr_end)

        row = {
            'read_id': read_id,
            'chr_end': args.chr_end,
            'read_length': len(read_seqs[read_id]),
            'telo_side': telo_side,
        }
        row.update(spacer_result)
        row.update(x_result)
        row.update(y_result)
        row.update(reconciliation)
        row['qc_flags'] = ''
        rows.append(row)

    # Write output
    os.makedirs(os.path.dirname(args.output_tsv) or '.', exist_ok=True)
    write_results_tsv(rows, args.output_tsv)

    from collections import Counter
    recomb_counts = Counter(r['recombination_detected'] for r in rows)
    print(f'  Output: {args.output_tsv} ({len(rows)} reads)')
    print(f'  Recombination: {recomb_counts.get(True, 0)} detected, '
          f'{recomb_counts.get(False, 0)} no recombination')
    pct_quick = 100 * n_quick_check_passed / max(len(rows), 1)
    print(f'  Spacer quick check: {n_quick_check_passed}/{len(rows)} reads passed ({pct_quick:.1f}%)')


if __name__ == '__main__':
    main()
