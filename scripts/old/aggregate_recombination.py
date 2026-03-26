"""
Step 12/13: Aggregate per-chr-end recombination results.

Two modes:

  Per-chr-end mode (Step 12):
    Joins the three element-specific TSVs (Y prime, X prime, spacer) with the
    minimap2 breakpoints TSV into a unified per-read TSV.  Applies evidence
    hierarchy: structural (0.5) > Y prime (0.35) > spacer (0.15).

  Summary mode (Step 13):
    Reads all per-chr-end TSVs in the recombination directory and writes a
    single per-chr-end summary TSV.

Usage — per-chr-end:
  python aggregate_recombination.py \
      --breakpoints-tsv  results/{base}/recombination/{base}_{chr_end}_breakpoints.tsv \
      --y-prime-tsv      .../{chr_end}_y_prime_recomb.tsv \
      --x-prime-tsv      .../{chr_end}_x_prime_recomb.tsv \
      --spacer-tsv       .../{chr_end}_spacer_recomb.tsv \
      --chr-end  chr4R \
      --output-reads     results/{base}/recombination/{base}_{chr_end}_recombination_reads.tsv

Usage — summary:
  python aggregate_recombination.py --summarize \
      --recombination-dir  results/{base}/recombination/ \
      --base-name          {base_name} \
      --output-summary     results/{base}/recombination/{base_name}_recombination_summary.tsv
"""

import argparse
import glob as glob_module
import os
import sys

import pandas as pd

from recombination_utils import (
    Hypothesis,
    hypotheses_to_row_dict,
    is_recombination_candidate,
    has_structural_evidence,
    get_first_breakpoint_element,
    get_first_breakpoint_feature_type,
    get_post_breakpoint_contig,
    normalize_hypotheses,
    write_results_tsv,
    MIN_CONFIDENCE_TO_REPORT,
    MIN_MAPQ,
)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description='Aggregate recombination analysis results')
    p.add_argument('--summarize', action='store_true',
                   help='Run in summary mode (across all chr ends)')

    # Per-chr-end mode
    p.add_argument('--breakpoints-tsv')
    p.add_argument('--y-prime-tsv')
    p.add_argument('--x-prime-tsv')
    p.add_argument('--spacer-tsv')
    p.add_argument('--chr-end')
    p.add_argument('--output-reads')

    # Summary mode
    p.add_argument('--recombination-dir')
    p.add_argument('--base-name')
    p.add_argument('--output-summary')

    return p.parse_args()

# ---------------------------------------------------------------------------
# Safe TSV loader
# ---------------------------------------------------------------------------

def _load_tsv(path):
    """Load TSV; return empty DataFrame if file missing or empty."""
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep='\t')
        return df if not df.empty else pd.DataFrame()
    except Exception:
        return pd.DataFrame()

# ---------------------------------------------------------------------------
# Evidence hierarchy — build unified hypothesis for one read
# ---------------------------------------------------------------------------

TIER_STRUCTURAL = 0.5
TIER_YPRIME = 0.35
TIER_SPACER = 0.15

def _get_h1(df, read_id, col_prefix='h1_'):
    """Extract h1 fields from a TSV row for a given read_id."""
    if df.empty or 'read_id' not in df.columns:
        return {}
    row = df[df['read_id'] == read_id]
    if row.empty:
        return {}
    row = row.iloc[0]
    result = {}
    for suffix in ['type', 'description', 'confidence', 'source_chr_ends',
                   'breakpoint_element', 'ambiguous', 'is_compound']:
        key = col_prefix + suffix
        result[suffix] = row.get(key, '')
    return result


def _safe_conf(h1_dict):
    try:
        return float(h1_dict.get('confidence', 0.0))
    except (ValueError, TypeError):
        return 0.0


def _source_from(h1_dict):
    src = h1_dict.get('source_chr_ends', '')
    return str(src).split(';')[0] if src else ''


def build_unified_hypotheses(read_id, bp_row, yp_h1, xp_h1, sp_h1):
    """
    Combine evidence from all tiers into unified hypothesis list.

    Priority: structural > Y prime > spacer/X prime.
    """
    hypotheses = []

    alignment_class = str(bp_row.get('alignment_class', ''))
    recomb_flagged = is_recombination_candidate(bp_row)
    struct_evidence = has_structural_evidence(bp_row)
    post_breakpoint_contig = get_post_breakpoint_contig(bp_row)
    identity = float(bp_row.get('segment_1_identity', 0.0) or 0.0)
    coverage = float(bp_row.get('segment_1_coverage', 0.0) or 0.0)
    bp_element = get_first_breakpoint_element(bp_row)
    bp_feature_type = get_first_breakpoint_feature_type(bp_row)

    # --- No recombination case ---
    if not recomb_flagged:
        c_no_recomb = identity * (coverage ** 0.5)
        hypotheses.append(Hypothesis(
            h_type='no_recombination',
            description='Read matches day0 reference; no recombination detected',
            confidence=c_no_recomb,
        ))
        return normalize_hypotheses(hypotheses), 'none', True, 'none'

    # Collect evidence tiers
    evidence_tiers = []

    # --- Tier 1: Structural ---
    c_structural = 0.0
    structural_source = ''
    if struct_evidence and post_breakpoint_contig:
        c_structural = TIER_STRUCTURAL
        structural_source = post_breakpoint_contig
        evidence_tiers.append('structural')

    # --- Tier 2: Y prime ---
    yp_conf = _safe_conf(yp_h1)
    yp_type = yp_h1.get('type', '')
    yp_source = _source_from(yp_h1)
    c_yprime = 0.0

    if yp_conf > 0 and yp_type not in ('no_recombination', 'no_y_prime_change', ''):
        c_yprime = yp_conf * TIER_YPRIME
        evidence_tiers.append('y_prime')

    # --- Tier 3: Spacer / X prime ---
    sp_conf = _safe_conf(sp_h1)
    sp_source = _source_from(sp_h1)
    xp_conf = _safe_conf(xp_h1)
    xp_source = _source_from(xp_h1)
    c_spacer = 0.0
    spacer_source = ''

    if sp_conf > 0:
        c_spacer = sp_conf * TIER_SPACER
        spacer_source = sp_source
        evidence_tiers.append('spacer')
    elif xp_conf > 0 and not c_yprime and not c_structural:
        c_spacer = xp_conf * TIER_SPACER
        spacer_source = xp_source
        evidence_tiers.append('x_prime')

    # --- Determine agreed source chr end ---
    # Sources from each tier that has a result
    tier_sources = [s for s in [structural_source, yp_source, spacer_source] if s]

    if not tier_sources:
        hypotheses.append(Hypothesis(
            h_type='ambiguous',
            description=f'Recombination flagged in {bp_feature_type} but no source identified',
            confidence=0.3,
            breakpoint_element=bp_element,
            ambiguous=True,
        ))
        return normalize_hypotheses(hypotheses), ';'.join(evidence_tiers), False, ''

    # Check cross-element consistency
    unique_sources = set(tier_sources)
    cross_consistent = len(unique_sources) == 1
    agreed_source = tier_sources[0]  # use first (highest-priority) source

    if cross_consistent:
        # All tiers agree — sum confidences and apply consistency boost
        c_final = (c_structural + c_yprime + c_spacer) * 1.3
        c_final = min(1.0, c_final)

        # Determine description from best evidence tier
        if yp_conf > 0 and yp_type not in ('no_recombination', 'no_y_prime_change', ''):
            description = yp_h1.get('description', f"Recombined to {agreed_source}")
        elif sp_conf > 0:
            description = sp_h1.get('description', f"Spacer switch to {agreed_source}")
        elif struct_evidence:
            description = f"Recombined to {agreed_source} (structural evidence)"
        else:
            description = f"Recombined to {agreed_source}"

        hypotheses.append(Hypothesis(
            h_type=bp_feature_type if bp_feature_type else 'y_prime',
            description=description,
            confidence=c_final,
            breakpoint_element=bp_element,
            source_chr_ends=[agreed_source],
            ambiguous=False,
        ))

    else:
        # Tiers disagree — generate separate hypotheses per source
        seen_sources = set()
        for src, c, desc_template in [
            (structural_source, c_structural, f"Recombined to {{src}} (structural)"),
            (yp_source, c_yprime, yp_h1.get('description', "Y prime: recombined to {src}")),
            (spacer_source, c_spacer, sp_h1.get('description', "Spacer: switch to {src}")),
        ]:
            if not src or src in seen_sources or c < MIN_CONFIDENCE_TO_REPORT:
                continue
            seen_sources.add(src)
            desc = desc_template.replace('{src}', src) if '{src}' in desc_template else desc_template
            hypotheses.append(Hypothesis(
                h_type='ambiguous',
                description=desc + ' — INCONSISTENT across elements',
                confidence=c,
                breakpoint_element=bp_element,
                source_chr_ends=[src],
                ambiguous=True,
            ))

    # Y prime isolation rule: if spacer/X prime is ambiguous, promote Y prime to full weight
    if yp_conf > 0 and sp_conf < 0.4 and xp_conf < 0.4 and not struct_evidence:
        for h in hypotheses:
            if h.h_type in ('y_prime', 'ambiguous') and yp_source in h.source_chr_ends:
                h.confidence = min(1.0, h.confidence / TIER_YPRIME)
                h.description = (
                    h.description.rstrip('.')
                    + ' (Y prime evidence only; spacer ambiguous)'
                )
                break

    return (
        normalize_hypotheses(hypotheses),
        ';'.join(evidence_tiers),
        cross_consistent,
        agreed_source if cross_consistent else ';'.join(unique_sources),
    )

# ---------------------------------------------------------------------------
# Per-chr-end aggregation
# ---------------------------------------------------------------------------

def aggregate_chr_end(args):
    df_bp = _load_tsv(args.breakpoints_tsv)
    df_yp = _load_tsv(args.y_prime_tsv)
    df_xp = _load_tsv(args.x_prime_tsv)
    df_sp = _load_tsv(args.spacer_tsv)

    if df_bp.empty:
        print(f'  No reads in breakpoints TSV — writing empty output')
        write_results_tsv([], args.output_reads)
        return

    rows = []
    for _, bp_row in df_bp.iterrows():
        read_id = bp_row['read_id']

        yp_h1 = _get_h1(df_yp, read_id)
        xp_h1 = _get_h1(df_xp, read_id)
        sp_h1 = _get_h1(df_sp, read_id)

        hypotheses, evidence_tiers, cross_consistent, final_source = build_unified_hypotheses(
            read_id, bp_row, yp_h1, xp_h1, sp_h1
        )

        row = {
            'read_id': read_id,
            'chr_end': args.chr_end,
            'read_length': bp_row.get('read_length', ''),
            'telo_side': bp_row.get('telo_side', ''),
            'alignment_class': bp_row.get('alignment_class', ''),
            'segment_1_contig': bp_row.get('segment_1_contig', ''),
            'segment_1_identity': bp_row.get('segment_1_identity', ''),
            'segment_1_coverage': bp_row.get('segment_1_coverage', ''),
            'segment_2_contig': bp_row.get('segment_2_contig', ''),
            'recombination_detected': is_recombination_candidate(bp_row),
            'structural_evidence': has_structural_evidence(bp_row),
            'post_breakpoint_contig': get_post_breakpoint_contig(bp_row),
            'breakpoint_element': get_first_breakpoint_element(bp_row),
            'breakpoint_feature_type': get_first_breakpoint_feature_type(bp_row),
            'y_prime_source_chr_end': _source_from(yp_h1),
            'y_prime_validation_result': df_yp[df_yp['read_id'] == read_id].iloc[0].get(
                'y_prime_validation_result', '') if not df_yp.empty and
                'read_id' in df_yp.columns and not df_yp[df_yp['read_id'] == read_id].empty
                else '',
            'cross_element_consistent': cross_consistent,
            'evidence_tiers_used': evidence_tiers,
            'n_hypotheses': len(hypotheses),
            'qc_flags': bp_row.get('qc_flags', ''),
        }
        row.update(hypotheses_to_row_dict(hypotheses))
        rows.append(row)

    os.makedirs(os.path.dirname(args.output_reads) or '.', exist_ok=True)
    write_results_tsv(rows, args.output_reads)
    print(f'  Per-read TSV: {args.output_reads} ({len(rows)} reads)')

# ---------------------------------------------------------------------------
# Summary mode
# ---------------------------------------------------------------------------

def summarize(args):
    # Look for both new (_features.tsv) and old (_recombination_reads.tsv) formats
    pattern_new = os.path.join(args.recombination_dir, '*_features.tsv')
    pattern_old = os.path.join(args.recombination_dir, '*_recombination_reads.tsv')
    files = sorted(glob_module.glob(pattern_new))
    suffix = '_features.tsv'
    if not files:
        files = sorted(glob_module.glob(pattern_old))
        suffix = '_recombination_reads.tsv'
    print(f'  Found {len(files)} per-chr-end TSVs')

    summary_rows = []
    for fpath in files:
        df = _load_tsv(fpath)
        if df.empty:
            continue
        fname = os.path.basename(fpath)
        # Extract chr_end from filename
        parts = fname.replace(suffix, '').split('_')
        chr_end = parts[-1] if parts else fname

        total = len(df)
        if total == 0:
            continue

        def pct(mask):
            return round(mask.sum() / total * 100, 2) if total > 0 else 0.0

        # Classify by h1_type
        h1_type = df.get('h1_type', pd.Series([''] * total))
        recomb_detected = df.get('recombination_detected', pd.Series([False] * total))

        row = {
            'chr_end': chr_end,
            'total_reads': total,
            'pct_no_recombination': pct(h1_type == 'no_recombination'),
            'pct_no_y_prime_change': pct(h1_type == 'no_y_prime_change'),
            'pct_spacer_recomb': pct(h1_type == 'spacer'),
            'pct_x_prime_recomb': pct(h1_type == 'x_prime'),
            'pct_y_prime_recomb': pct(h1_type == 'y_prime'),
            'pct_compound_recomb': pct(h1_type == 'compound'),
            'pct_ambiguous': pct(h1_type == 'ambiguous'),
            'pct_recombination_flagged': pct(recomb_detected == True),
        }

        # Mean confidence of top hypothesis
        if 'h1_confidence' in df.columns:
            h1_conf = pd.to_numeric(df['h1_confidence'], errors='coerce')
            row['mean_h1_confidence'] = round(h1_conf.mean(), 4)
        else:
            row['mean_h1_confidence'] = ''

        # Most common h1 description
        if 'h1_description' in df.columns:
            row['most_common_h1_description'] = (
                df['h1_description'].value_counts().index[0]
                if not df['h1_description'].value_counts().empty else ''
            )
        else:
            row['most_common_h1_description'] = ''

        # Y prime transition counts
        if 'y_prime_recombination_status' in df.columns:
            vc = df['y_prime_recombination_status'].value_counts().to_dict()
            row['y_prime_transition_counts'] = ';'.join(
                f'{k}={v}' for k, v in vc.items()
            )
        else:
            row['y_prime_transition_counts'] = ''

        summary_rows.append(row)

    os.makedirs(os.path.dirname(args.output_summary) or '.', exist_ok=True)
    write_results_tsv(summary_rows, args.output_summary)
    print(f'  Summary TSV: {args.output_summary} ({len(summary_rows)} chr ends)')

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if args.summarize:
        print('aggregate_recombination.py — summary mode')
        if not args.recombination_dir or not args.output_summary:
            print('ERROR: --summarize requires --recombination-dir and --output-summary',
                  file=sys.stderr)
            sys.exit(1)
        summarize(args)
    else:
        print(f'aggregate_recombination.py — chr_end={args.chr_end}')
        if not args.breakpoints_tsv or not args.output_reads:
            print('ERROR: per-chr-end mode requires --breakpoints-tsv and --output-reads',
                  file=sys.stderr)
            sys.exit(1)
        aggregate_chr_end(args)


if __name__ == '__main__':
    main()
