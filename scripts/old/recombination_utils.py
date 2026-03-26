"""
Shared utilities for the telomere recombination analysis pipeline.

Used by:
  run_alignment.py
  analyze_y_prime_recombination.py
  analyze_x_prime_recombination.py
  analyze_spacer_recombination.py
  aggregate_recombination.py
"""

import os
import subprocess
import tempfile
from dataclasses import dataclass, field

import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MIN_MAPQ = 10
MIN_WHOLE_READ_IDENTITY = 0.95   # both must exceed this for "no_recombination"
MIN_WHOLE_READ_COVERAGE = 0.95
MIN_SEGMENT_IDENTITY = 0.80     # minimum identity for a segment to be considered "clean"
MAX_INTER_SEGMENT_GAP = 100     # bp — max gap between chained segments
MIN_DIVERGED_LEN = 200          # bp — minimum unaligned tail worth analyzing
BOUNDARY_BUFFER = 50            # bp — within this of a feature boundary = clean switch
BLAST_MIN_PIDENT = 75.0
BLAST_EVALUE = 1e-5
MIN_CONFIDENCE_TO_REPORT = 0.05
N_MAX_HYPOTHESES = 5
SPACER_CHUNK_SIZE = 250          # bp
SPACER_CHUNK_STEP = 125          # bp (50 % overlap)

# ---------------------------------------------------------------------------
# Shared data structures
# ---------------------------------------------------------------------------

@dataclass
class Hypothesis:
    h_type: str          # "no_recombination"|"no_y_prime_change"|"spacer"|"x_prime"|"y_prime"|"compound"|"ambiguous"
    description: str
    confidence: float    # 0–1
    breakpoint_element: str = ""
    breakpoint_pos_on_read: int = -1
    breakpoint_pos_on_ref: int = -1
    source_chr_ends: list = field(default_factory=list)
    ambiguous: bool = False
    is_compound: bool = False
    child_hypotheses: list = field(default_factory=list)


@dataclass
class BreakpointInfo:
    read_id: str
    chr_end: str
    read_length: int
    telo_side: str                 # "beginning" | "end"
    alignment_identity: float
    alignment_coverage: float
    aligned_to_contig: str
    breakpoint_element: str        # BED feature name (empty if none)
    breakpoint_feature_type: str   # "spacer"|"x_core"|"x_variable"|"y_prime"|"its"|"anchor"|""
    is_mid_element: bool
    breakpoint_pos_on_read: int    # -1 if no breakpoint
    breakpoint_pos_on_ref: int     # -1 if no breakpoint
    diverged_seq_length: int       # 0 if none
    qc_flags: list
    recombination_flagged: bool
    post_breakpoint_contig: str = ""   # from supplementary alignment (or "" for Case B)
    structural_evidence: bool = False  # True = supplementary alignment confirmed source


# ---------------------------------------------------------------------------
# BED parsing
# ---------------------------------------------------------------------------

FEATURE_TYPE_MAP = {
    '_Telomere_Repeat': 'telomere',
    '_x_variable_element': 'x_variable',
    '_x_core_element': 'x_core',
    '_space_between_anchor': 'spacer',
    '_Y_Prime_': 'y_prime',
    'ITS_': 'its',
}


def _classify_feature_type(name):
    for key, ftype in FEATURE_TYPE_MAP.items():
        if key in name:
            return ftype
    if '_anchor' in name and '_x_' not in name:
        return 'anchor'
    return 'unknown'


def load_bed_features(bed_file, chr_end):
    """Parse a simplified BED file and return features for chr_end sorted by start."""
    features = []
    with open(bed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            contig, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            # Match chr_end: feature name should start with chr_end (e.g., "chr9R_anchor")
            if not name.startswith(chr_end):
                continue
            ftype = _classify_feature_type(name)
            features.append({
                'contig': contig,
                'start': start,
                'end': end,
                'name': name,
                'feature_type': ftype,
            })
    features.sort(key=lambda f: f['start'])
    return features


def classify_breakpoint(ref_pos, features):
    """
    Given a reference position, determine which BED feature it falls in.

    Returns (feature_name, feature_type, is_mid_element).
    is_mid_element=True means the breakpoint is within a feature body (not near a boundary).
    """
    for feat in features:
        if feat['start'] <= ref_pos <= feat['end']:
            near_start = abs(ref_pos - feat['start']) <= BOUNDARY_BUFFER
            near_end = abs(ref_pos - feat['end']) <= BOUNDARY_BUFFER
            is_mid = not (near_start or near_end)
            return feat['name'], feat['feature_type'], is_mid

    # Not in any annotated feature — find flanking features
    prev_name, next_name = '', ''
    for i, feat in enumerate(features):
        if feat['end'] < ref_pos:
            prev_name = feat['name']
        elif feat['start'] > ref_pos:
            next_name = feat['name']
            break
    gap_name = f"gap_{prev_name}_{next_name}" if prev_name or next_name else "unannotated"
    return gap_name, 'gap', False


# ---------------------------------------------------------------------------
# Read orientation
# ---------------------------------------------------------------------------

def telo_side_from_strand(chr_end, strand):
    """
    Determine which end of the read the telomere is on.

    Mirrors make_y_prime_repeatmasker_tsv.py:61:
      ('L' in chr_end and strand=='AC') or ('R' in chr_end and strand=='TG') → 'beginning'
    """
    if ('L' in chr_end and strand == 'AC') or ('R' in chr_end and strand == 'TG'):
        return 'beginning'
    return 'end'


def telo_side_from_header(fasta_header, chr_end):
    """Infer telo_side from a FASTA read header containing 'AC' or 'TG'."""
    if 'AC' in fasta_header:
        strand = 'AC'
    elif 'TG' in fasta_header:
        strand = 'TG'
    else:
        return 'end'
    return telo_side_from_strand(chr_end, strand)


# ---------------------------------------------------------------------------
# BLAST helpers
# ---------------------------------------------------------------------------

BLAST_COLUMNS = [
    'qseqid', 'qlen', 'length', 'qstart', 'qend',
    'sseqid', 'slen', 'sstart', 'send', 'pident', 'bitscore', 'evalue',
]


def _build_blast_db(library_fasta, tmp_dir):
    db_path = os.path.join(tmp_dir, os.path.basename(library_fasta) + '.blastdb')
    nhr = db_path + '.nhr'
    if not os.path.exists(nhr):
        subprocess.run(
            ['makeblastdb', '-in', library_fasta, '-dbtype', 'nucl', '-out', db_path],
            check=True, capture_output=True,
        )
    return db_path


def run_blast(query_fasta, library_fasta, tmp_dir, label='blast',
              min_pident=BLAST_MIN_PIDENT, evalue=BLAST_EVALUE):
    """
    Run blastn of query_fasta against library_fasta.
    Returns a DataFrame with standard BLAST columns (empty DataFrame if no hits).
    """
    if not os.path.exists(query_fasta) or os.path.getsize(query_fasta) == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)

    db_path = _build_blast_db(library_fasta, tmp_dir)
    out_file = os.path.join(tmp_dir, f'{label}_blast.tsv')

    result = subprocess.run(
        [
            'blastn', '-query', query_fasta, '-db', db_path,
            '-perc_identity', str(min_pident), '-evalue', str(evalue),
            '-outfmt', '6 qseqid qlen length qstart qend sseqid slen sstart send pident bitscore evalue',
            '-out', out_file,
        ],
        check=True, capture_output=True,
    )

    if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)

    df = pd.read_csv(out_file, sep='\t', header=None, names=BLAST_COLUMNS)
    return df


# ---------------------------------------------------------------------------
# Hypothesis normalization
# ---------------------------------------------------------------------------

def normalize_hypotheses(hypotheses):
    """
    Normalize hypothesis confidences so they sum to ≤ 1.
    Sort descending, filter below MIN_CONFIDENCE_TO_REPORT, cap at N_MAX_HYPOTHESES.
    """
    total = sum(h.confidence for h in hypotheses)
    if total > 1.0:
        for h in hypotheses:
            h.confidence /= total
    hypotheses.sort(key=lambda h: h.confidence, reverse=True)
    return [h for h in hypotheses if h.confidence >= MIN_CONFIDENCE_TO_REPORT][:N_MAX_HYPOTHESES]


# ---------------------------------------------------------------------------
# TSV output helpers
# ---------------------------------------------------------------------------

def hypotheses_to_row_dict(hypotheses, n_slots=N_MAX_HYPOTHESES):
    """Flatten a list of Hypothesis objects into h1_* … h5_* dict entries."""
    row = {}
    for i in range(n_slots):
        prefix = f'h{i+1}_'
        if i < len(hypotheses):
            h = hypotheses[i]
            row[prefix + 'type'] = h.h_type
            row[prefix + 'description'] = h.description
            row[prefix + 'confidence'] = round(h.confidence, 4)
            row[prefix + 'breakpoint_element'] = h.breakpoint_element
            row[prefix + 'breakpoint_pos_on_read'] = h.breakpoint_pos_on_read
            row[prefix + 'breakpoint_pos_on_ref'] = h.breakpoint_pos_on_ref
            row[prefix + 'source_chr_ends'] = ';'.join(h.source_chr_ends)
            row[prefix + 'ambiguous'] = h.ambiguous
            row[prefix + 'is_compound'] = h.is_compound
        else:
            for suffix in ['type', 'description', 'confidence', 'breakpoint_element',
                           'breakpoint_pos_on_read', 'breakpoint_pos_on_ref',
                           'source_chr_ends', 'ambiguous', 'is_compound']:
                row[prefix + suffix] = ''
    return row


def write_results_tsv(rows, output_path):
    """Write a list of row dicts to a TSV file."""
    if not rows:
        # Write empty file with header derived from a dummy row
        pd.DataFrame().to_csv(output_path, sep='\t', index=False)
        return
    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep='\t', index=False)


# ---------------------------------------------------------------------------
# Step 10 column helpers (new schema with semicolon-delimited breakpoint fields)
# ---------------------------------------------------------------------------

def is_recombination_candidate(bp_row):
    """True if alignment_class is two_clean_halves or one_clean_half."""
    return str(bp_row.get('alignment_class', '')) in ('two_clean_halves', 'one_clean_half')


def has_structural_evidence(bp_row):
    """True if alignment_class is two_clean_halves (source contig known)."""
    return str(bp_row.get('alignment_class', '')) == 'two_clean_halves'


def get_first_breakpoint_feature_type(bp_row):
    val = str(bp_row.get('breakpoint_feature_types', ''))
    return val.split(';')[0] if val else ''


def get_first_breakpoint_element(bp_row):
    val = str(bp_row.get('breakpoint_elements', ''))
    return val.split(';')[0] if val else ''


def get_first_breakpoint_pos_on_read(bp_row):
    val = str(bp_row.get('breakpoint_positions_on_read', ''))
    try:
        return int(val.split(';')[0])
    except (ValueError, IndexError):
        return -1


def get_first_breakpoint_is_mid_element(bp_row):
    val = str(bp_row.get('breakpoint_is_mid_element', ''))
    return val.split(';')[0] == 'True' if val else False


def get_post_breakpoint_contig(bp_row):
    """Get the source contig from breakpoint_post_contigs or segment_2_contig."""
    val = str(bp_row.get('breakpoint_post_contigs', ''))
    first = val.split(';')[0] if val else ''
    if first:
        return first
    return str(bp_row.get('segment_2_contig', ''))


# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------

def read_fasta(fasta_path):
    """Read a FASTA file, return dict {header_without_gt: sequence}."""
    sequences = {}
    current_header = None
    current_seq = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
    if current_header is not None:
        sequences[current_header] = ''.join(current_seq)
    return sequences


def write_fasta(sequences, fasta_path):
    """Write dict {header: sequence} to a FASTA file."""
    with open(fasta_path, 'w') as fh:
        for header, seq in sequences.items():
            fh.write(f'>{header}\n{seq}\n')


def write_fasta_list(entries, fasta_path):
    """Write list of (header, sequence) tuples to a FASTA file."""
    with open(fasta_path, 'w') as fh:
        for header, seq in entries:
            fh.write(f'>{header}\n{seq}\n')
