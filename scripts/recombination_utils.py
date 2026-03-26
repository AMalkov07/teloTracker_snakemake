"""
Shared utilities for the telomere recombination analysis pipeline.

Used by:
  run_alignment.py      (step 10 — supplementary evidence)
  analyze_features.py   (step 11 — chunk-based feature analysis)
  aggregate_recombination.py  (step 12 — summary)
"""

import os
import subprocess
from dataclasses import dataclass, field

import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BLAST_MIN_PIDENT = 75.0
BLAST_EVALUE = 1e-5
SPACER_CHUNK_SIZE = 250         # bp — non-overlapping chunks

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


# ---------------------------------------------------------------------------
# Read orientation
# ---------------------------------------------------------------------------

def telo_side_from_strand(chr_end, strand):
    """
    Determine which end of the read the telomere is on.
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

    subprocess.run(
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
# TSV output helpers
# ---------------------------------------------------------------------------

def write_results_tsv(rows, output_path):
    """Write a list of row dicts to a TSV file."""
    if not rows:
        pd.DataFrame().to_csv(output_path, sep='\t', index=False)
        return
    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep='\t', index=False)


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
