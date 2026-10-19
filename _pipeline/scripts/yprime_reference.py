"""Expected Y' copies per chromosome end: the zero line of the delta-Y' plots.

The counts come from a label BED (pretelomeric_regions_<name>_simp.bed: each end's Y' rows) or a
counts table (chr_end<TAB>expected_y_primes, all 32 ends). resolve() turns the yprime_reference
config value into one of those files:

  auto    the shipped reference named by the strain; else, in recombination, the day-0 label BED;
          else none
  none    no reference: plots show absolute Y' counts per end
  <name>  a shipped reference: references/day0/<name>/ (or <name>_curated/) label BED, or
          references/yprime_expected_counts/<name>.tsv
  <path>  a label BED (.bed) or counts table (.tsv)
"""
import os
import re

CHR_ENDS = [f"chr{n}{arm}" for n in range(1, 17) for arm in ("L", "R")]

# chr4R_Y_Prime_1 is a Y' copy of chr4R; ITS rows (ITS_chr4R_Y_Prime_1-2) are not.
Y_PRIME_ROW = re.compile(r'^(chr\d+[LR])_Y_Prime_\d+$')
END_ROW = re.compile(r'^(chr\d+[LR])_')


def load_expected_counts(path):
    """{chr_end: expected Y' copies} from a label BED or a counts table. A BED leaves out ends it
    has no rows for at all (not labelled): their expected count is unknown, not zero."""
    if path.endswith('.bed'):
        return _counts_from_bed(path)
    return _counts_from_table(path)


def _counts_from_bed(path):
    labelled, counts = set(), {}
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#') or len(fields) < 4:
                continue
            end = END_ROW.match(fields[3])
            if end:
                labelled.add(end.group(1))
            yp = Y_PRIME_ROW.match(fields[3])
            if yp:
                counts[yp.group(1)] = counts.get(yp.group(1), 0) + 1
    return {end: counts.get(end, 0) for end in CHR_ENDS if end in labelled}


def _counts_from_table(path):
    counts = {}
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        if header[:2] != ['chr_end', 'expected_y_primes']:
            raise ValueError(f"{path}: the header must be 'chr_end<TAB>expected_y_primes'")
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            end, n = line.rstrip('\n').split('\t')[:2]
            if end not in CHR_ENDS:
                raise ValueError(f"{path}: unknown chromosome end {end!r} (expected chr1L ... chr16R)")
            counts[end] = int(n)
    missing = [end for end in CHR_ENDS if end not in counts]
    if missing:
        raise ValueError(f"{path}: no expected count for {', '.join(missing)}")
    return counts


def shipped_reference(name, references_dir):
    """The shipped label BED or counts table called name, or None."""
    for folder in (name, f"{name}_curated"):
        bed = os.path.join(references_dir, 'day0', folder, f'pretelomeric_regions_{folder}_simp.bed')
        if os.path.isfile(bed):
            return bed
    table = os.path.join(references_dir, 'yprime_expected_counts', f'{name}.tsv')
    return table if os.path.isfile(table) else None


def shipped_names(references_dir):
    names = []
    day0 = os.path.join(references_dir, 'day0')
    if os.path.isdir(day0):
        names += sorted(d for d in os.listdir(day0) if shipped_reference(d, references_dir))
    tables = os.path.join(references_dir, 'yprime_expected_counts')
    if os.path.isdir(tables):
        names += sorted(f[:-4] for f in os.listdir(tables) if f.endswith('.tsv'))
    return names


def resolve(value, strain, day0_bed, references_dir):
    """(file or None, label) for a yprime_reference config value. day0_bed is the day-0 label BED
    in recombination runs (None otherwise). Raises ValueError for a name or file that is not there."""
    value = 'auto' if value is None or str(value).strip() == '' else str(value).strip()
    if value.lower() == 'none':
        return None, 'none (absolute counts)'
    if value.lower() == 'auto':
        if strain:
            path = shipped_reference(str(strain), references_dir)
            if path:
                return path, f'{strain} (shipped; matched by strain)'
        if day0_bed and os.path.isfile(day0_bed):
            return day0_bed, 'day-0 label BED'
        return None, 'none (absolute counts; no reference found for this strain)'
    if os.sep in value or value.endswith(('.bed', '.tsv')):
        if not os.path.isfile(value):
            raise ValueError(f"yprime_reference file not found: {value}")
        return value, value
    path = shipped_reference(value, references_dir)
    if path is None:
        raise ValueError(f"yprime_reference '{value}' is not a shipped reference "
                         f"({', '.join(shipped_names(references_dir))}), 'auto', 'none' or a file")
    return path, f'{value} (shipped)'
