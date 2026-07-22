"""
Utility functions for telomere extension and polishing pipeline.

Import this with: from subtelomere_reference_pipeline_utils import *
"""

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools
import subprocess
import os
import tempfile
import pandas as pd


def get_75th_percentile_reads(input_tsv, output_dir, output_file):
    """
    Read TSV file and get the 75th percentile read for each chr_end based on repeat_length.

    Parameters:
    input_tsv (str): Path to the input TSV file
    output_dir (str): Path to output directory
    output_file (str): Path to output file for read IDs

    Returns:
    tuple: (dict mapping chr_end to read_id, list of all selected read_ids)
    """
    print("Reading input TSV file and selecting 75th percentile reads...")

    # Read the TSV
    df = pd.read_csv(input_tsv, sep='\t')

    # Filter for telomere reads
    df_telomere_reads = df[(df['repeat_length'] >= 30) & (df['Adapter_After_Telomere'] == True)]

    # Define all chromosome ends
    chr_ends = [f'{chr_num}{side}' for chr_num in range(1, 17) for side in ['L', 'R']]

    chr_end_to_read = {}
    all_selected_read_ids = []

    for chr_end in chr_ends:
        if chr_end not in df_telomere_reads['chr_end'].values:
            print(f"Warning: {chr_end} not found in the input file.")
            continue

        df_chr_telomere_reads = df_telomere_reads[df_telomere_reads['chr_end'] == chr_end]

        if df_chr_telomere_reads.empty:
            print(f"Warning: No telomere reads found for {chr_end}.")
            continue

        # Get mode y_prime_probe_count
        mode_y_prime = df_chr_telomere_reads['y_prime_probe_count'].mode()[0]

        # Filter for mode y_prime
        df_mode_y_chr_telomere_reads = df_chr_telomere_reads[
            df_chr_telomere_reads['y_prime_probe_count'] == mode_y_prime]

        if df_mode_y_chr_telomere_reads.empty:
            print(f"Warning: No reads found with mode y_prime_probe_count for {chr_end}.")
            continue

        # Get the read ID for the read with the 75% quantile repeat_length.
        # read_id is a deterministic tiebreak: repeat_length ties are common
        # (7172 chr4L has 2 reads at the percentile value) and pandas' default
        # quicksort is unstable, so without it the scaffold -- and therefore the
        # reference -- depends on input row order. Measured: 7172 chr4L flips
        # between SRR33298432.121462 and SRR33298432.129899 across runs.
        sorted_df = df_mode_y_chr_telomere_reads.sort_values(
            by=['repeat_length', 'read_id'], kind='mergesort')
        size_of_df = len(sorted_df)
        quantile_index = int(size_of_df * 0.75)
        if quantile_index >= size_of_df:
            quantile_index = size_of_df - 1

        quantile_read_id = sorted_df.iloc[quantile_index]['read_id']
        chr_end_to_read[chr_end] = quantile_read_id
        all_selected_read_ids.append(quantile_read_id)

        print(f"{chr_end}: Selected read {quantile_read_id} (repeat_length: {sorted_df.iloc[quantile_index]['repeat_length']})")

    print(f"\nSelected {len(chr_end_to_read)} reads for {len(chr_end_to_read)} chromosome ends")

    # Save read IDs to file
    with open(output_file, 'w') as f:
        for chr_end, read_id in chr_end_to_read.items():
            f.write(f"{chr_end}\t{read_id}\n")

    # Also save just the read IDs for easy extraction
    read_ids_only_file = output_file.replace('.txt', '_only.txt')
    with open(read_ids_only_file, 'w') as f:
        for read_id in all_selected_read_ids:
            f.write(f"{read_id}\n")

    return chr_end_to_read, all_selected_read_ids


# ---------------------------------------------------------------------------
# Consensus scaffold selection
# ---------------------------------------------------------------------------

def _load_fai(fai_path):
    idx = {}
    with open(fai_path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 5:
                idx[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
    return idx


def _fetch_seq(fh, rec):
    """Pull one sequence from an open binary FASTA handle using its .fai record."""
    length, offset, line_bases, line_width = rec
    fh.seek(offset)
    n_lines = (length + line_bases - 1) // line_bases
    return ''.join(fh.read(n_lines * line_width).decode().split())[:length]


def _extension_of(seq, anchor_info):
    """The telomere-side portion of a read, past its anchor.

    This is the part that becomes the soft clip and is grafted onto the
    reference, so it is the only part worth comparing between candidates.
    Mirrors the pipeline's own `wanted_section_of_read` convention.
    """
    if anchor_info['wanted_section_of_read'] == 'before_match_start_on_read':
        return seq[:int(anchor_info['match_start_on_read'])]
    return seq[int(anchor_info['match_end_on_read']):]


def _pairwise_agreement(records, threads, min_identity, min_coverage):
    """All-vs-all align candidate extensions; return {(a,b): (pident, coverage)}.

    Coverage is deliberately computed from the SINGLE BEST HSP rather than the
    sum of all HSPs. A read that is colinear with its peers but carries a large
    insertion then scores low coverage and is rejected -- which is what we want:
    on 7172 chr4L the selected scaffold was colinear at 99.5% yet carried a
    1,661bp insertion (52bp poly-A + 1,609bp of (ATGGA)n tandem repeat, absent
    from the S. cerevisiae genome, present in 1.16% of anchored reads). Summed
    coverage scores that read 1.00 and would have kept it; best-HSP coverage
    scores it 0.60 and drops it. Both indels and rearrangements are disqualifying
    in a scaffold, so one metric covering both is the right trade.
    """
    lens = {r[0]: len(r[1]) for r in records}
    best = {}
    with tempfile.TemporaryDirectory() as tmp:
        qf = os.path.join(tmp, 'cand.fa')
        with open(qf, 'w') as out:
            for r in records:
                out.write(f'>{r[0]}\n{r[1]}\n')
        db = os.path.join(tmp, 'db')
        subprocess.run(['makeblastdb', '-in', qf, '-dbtype', 'nucl', '-out', db],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        res = subprocess.run(
            ['blastn', '-query', qf, '-db', db, '-task', 'dc-megablast',
             '-outfmt', '6 qseqid sseqid pident length', '-num_threads', str(threads),
             '-max_target_seqs', '50'],
            capture_output=True, text=True, check=True)

    for line in res.stdout.splitlines():
        f = line.split('\t')
        if len(f) < 4 or f[0] == f[1]:
            continue
        key = tuple(sorted((f[0], f[1])))
        pid, ln = float(f[2]), int(f[3])
        cov = min(1.0, ln / max(1, min(lens[f[0]], lens[f[1]])))
        score = pid * cov
        if key not in best or score > best[key][0]:
            best[key] = (score, pid, cov)
    return {k: (v[1], v[2]) for k, v in best.items()}


def _largest_agreeing_group(ids, agreement, min_identity, min_coverage):
    """Largest set of candidates that ALL mutually agree (exact, brute force).

    Mutual agreement matters: a star-shaped 'agrees with the middle one' set
    would let a read that disagrees with the others ride in on one good pair.
    Candidate counts here are tiny (n<=11), so exhaustive search is free.
    """
    def edge(a, b):
        pid, cov = agreement.get(tuple(sorted((a, b))), (0.0, 0.0))
        return pid >= min_identity and cov >= min_coverage

    for size in range(len(ids), 1, -1):
        for combo in itertools.combinations(ids, size):
            if all(edge(a, b) for a, b in itertools.combinations(combo, 2)):
                return list(combo)
    return [ids[0]] if ids else []


def select_consensus_scaffold_reads(input_tsv, all_matches_tsv, reads_fasta,
                                    output_dir, output_file,
                                    n_candidates=5, min_agree=3,
                                    min_identity=90.0, min_coverage=0.90,
                                    max_widenings=2, threads=4,
                                    report_tsv=None):
    """Pick one scaffold read per chr_end by consensus among near-percentile reads.

    Replaces the single 75th-percentile pick. That pick trusted one read with no
    check that it was typical, and on measured day-0 data it landed on a
    structurally aberrant read at 5 of 96 chromosome ends across three strains
    (7172 4L, 7858 1R, 7858 14L, 7871 12R, 7871 14L) -- reads that then had
    20-50kb grafted onto the reference.

    Method
    ------
    Take the `n_candidates` reads bracketing the 75th-percentile position, align
    their extensions all-vs-all, and require the scaffold to mutually agree with
    at least `min_agree - 1` peers. The read with the highest anchor identity in
    that agreeing group wins.

    Note this moves the scaffold at most ends (29 of 32 on 7172), since the
    highest-identity candidate is rarely the percentile read itself. That is
    expected, not drift: every candidate in the group has been shown to be
    structurally equivalent, so the choice between them is a quality question.

    WHY ANCHOR IDENTITY ONLY WORKS *AFTER* CONSENSUS
    ------------------------------------------------
    It is a fine way to rank corroborated reads and a terrible way to filter
    uncorroborated ones -- as a filter it is inverted. Measured over 468
    candidates in three strains, median anchor identity was 99.48% for
    structurally aberrant reads vs 99.63% for normal ones. At 4 of the 6 ends
    that had an outlier, the outlier had the HIGHEST anchor identity of all five
    candidates (7172 4L 99.94%, 7858 14L 99.90%, 7871 5R 99.96%, 7871 12R
    99.97%), so "prefer >=99%" would have actively selected the bad scaffold.
    The reason is structural: the anchor is a ~5kb reference-matched window, and
    a recombined or artifact-bearing read has a perfectly normal anchor because
    the damage lies downstream of it. Identity there measures local basecall
    quality, which is orthogonal to structural normality, and is further
    confounded by anchor position -- ONT accuracy decays ~2 points per 100kb, so
    a read whose anchor sits 35kb in scores lower for reasons unrelated to its
    quality.

    So the order matters: exclude aberrant reads by consensus FIRST, then rank
    what survives by anchor identity. Never the reverse.

    Guard rails
    -----------
    Two distinct failure modes, deliberately handled differently:
      * Too few reads to vote (pool < min_agree). Measured on day-0: 7858 12R
        has 1 read, 7871 12R has 2, 7858 14L has 4. This is coverage starvation,
        not disagreement -- retrying a different window cannot help. Falls back
        to the plain percentile pick with a loud warning.
      * Reads disagree (no mutually-agreeing group of min_agree). Widen the
        window and retry up to `max_widenings` times; raise if still unresolved,
        since silently scaffolding off a read no peer corroborates is exactly
        the failure this function exists to prevent.
    """
    print('Selecting scaffold reads by consensus among near-percentile candidates...')
    print(f'  n_candidates={n_candidates}  min_agree={min_agree}  '
          f'thresholds: pident>={min_identity} coverage>={min_coverage}')

    fai = reads_fasta + '.fai'
    if not os.path.exists(fai):
        raise FileNotFoundError(
            f'{fai} not found -- consensus selection needs random access to the '
            f'reads (run: samtools faidx {reads_fasta})')
    idx = _load_fai(fai)

    df = pd.read_csv(input_tsv, sep='\t')
    am = pd.read_csv(all_matches_tsv, sep='\t')
    am = am.sort_values('bitscore', ascending=False).drop_duplicates('read_id')
    anchor = am.set_index('read_id')[
        ['match_start_on_read', 'match_end_on_read', 'wanted_section_of_read', 'pident']
    ].to_dict('index')

    df_telomere_reads = df[(df['repeat_length'] >= 30) & (df['Adapter_After_Telomere'] == True)]
    chr_ends = [f'{n}{s}' for n in range(1, 17) for s in ['L', 'R']]

    chr_end_to_read = {}
    all_selected_read_ids = []
    report = []

    with open(reads_fasta, 'rb') as fh:
        for chr_end in chr_ends:
            sub = df_telomere_reads[df_telomere_reads['chr_end'] == chr_end]
            if sub.empty:
                print(f'Warning: no telomere reads found for {chr_end}.')
                continue
            mode_y = sub['y_prime_probe_count'].mode()
            if len(mode_y) == 0:
                print(f'Warning: no mode y_prime_probe_count for {chr_end}.')
                continue
            sub = sub[sub['y_prime_probe_count'] == mode_y[0]]
            if sub.empty:
                print(f'Warning: no reads with mode y_prime_probe_count for {chr_end}.')
                continue

            # read_id tiebreak + stable sort: see get_75th_percentile_reads.
            # Without it the candidate window itself shifts between runs.
            sub = sub.sort_values(by=['repeat_length', 'read_id'], kind='mergesort')
            pool = len(sub)
            q_idx = min(int(pool * 0.75), pool - 1)
            percentile_pick = sub.iloc[q_idx]['read_id']

            # --- guard rail 1: not enough reads to hold a vote at all
            if pool < min_agree:
                print(f'{chr_end}: WARNING only {pool} candidate read(s) in the pool '
                      f'(need {min_agree} to vote) -- falling back to the plain '
                      f'percentile pick {percentile_pick}. This end is scaffolded '
                      f'off an UNCORROBORATED read.')
                chr_end_to_read[chr_end] = percentile_pick
                all_selected_read_ids.append(percentile_pick)
                report.append({'chr_end': chr_end, 'pool': pool, 'method': 'fallback_low_coverage',
                               'n_candidates': pool, 'group_size': 1,
                               'selected_read': percentile_pick,
                               'percentile_pick': percentile_pick, 'changed': False})
                continue

            # --- widen the window until a mutually-agreeing majority appears
            chosen = group = cands = None
            max_records_seen = 0
            window_n = n_candidates
            for widening in range(max_widenings + 1):
                n = n_candidates + 2 * widening
                lo = max(0, min(q_idx - n // 2, pool - n))
                picks = sub.iloc[lo:lo + n]

                records = []
                for rid in picks['read_id']:
                    if rid not in idx or rid not in anchor:
                        continue
                    ext = _extension_of(_fetch_seq(fh, idx[rid]), anchor[rid])
                    if len(ext) >= 1000:
                        records.append((rid, ext, anchor[rid]['pident']))
                max_records_seen = max(max_records_seen, len(records))
                if len(records) < min_agree:
                    continue

                ids = [r[0] for r in records]
                agreement = _pairwise_agreement(records, threads, min_identity, min_coverage)
                grp = _largest_agreeing_group(ids, agreement, min_identity, min_coverage)
                if len(grp) >= min_agree:
                    pid_of = {r[0]: r[2] for r in records}
                    # Highest anchor identity within the agreeing group. Safe here
                    # precisely because consensus has already run: the trap that
                    # makes anchor identity useless as a FILTER -- an aberrant read
                    # scoring highest, seen at 4 of 6 outlier-bearing ends -- cannot
                    # bite once those reads are excluded. Among structurally
                    # equivalent candidates it is a reasonable proxy for basecall
                    # quality. (It is only a proxy: the anchor is not itself grafted,
                    # and identity is depressed for reads whose anchor sits late,
                    # since ONT accuracy decays ~2 points per 100kb.)
                    chosen = max(grp, key=lambda r: pid_of.get(r, 0.0))
                    group, cands = grp, ids
                    window_n = n
                    if widening:
                        print(f'{chr_end}: WARNING the initial window of {n_candidates} '
                              f'reads held no {min_agree} mutually-agreeing candidates; '
                              f'widened to {n} to find a group of {len(grp)}. This should '
                              f'be rare -- inspect this end, as it means several '
                              f'near-percentile reads disagree with each other.')
                    break

            # --- guard rail 2: the window could never be filled, so a failure to
            # agree is under-sampling rather than evidence against any read.
            #
            # The agreement test (best-HSP coverage) is deliberately strict: over
            # 937 measured pairs it rejected 44/44 independently-classified
            # aberrant pairs while keeping 98.4% of normal ones. Its cost is that
            # highly repetitive extensions fragment into many HSPs and can fail it
            # even when colinear -- 7858 14L's four reads are 97-99% identical over
            # 65-79kb but split into 174-196 HSPs across an 11-element Y' array.
            # With a full window that costs at most a pair or two and a group still
            # forms; only where the pool cannot fill the window does it block. Both
            # weaker metrics tried (summed-HSP coverage, longest-colinear-chain
            # coverage) recovered those pairs but readmitted 8-29 aberrant ones, so
            # the strict test is kept and the starved case handled here instead.
            if chosen is None and max_records_seen < n_candidates:
                print(f'{chr_end}: WARNING only {max_records_seen} comparable candidate(s) '
                      f'(wanted {n_candidates}); cannot tell an aberrant read from thin '
                      f'sampling -- falling back to the plain percentile pick '
                      f'{percentile_pick}. This end is scaffolded off an UNCORROBORATED read.')
                chr_end_to_read[chr_end] = percentile_pick
                all_selected_read_ids.append(percentile_pick)
                report.append({'chr_end': chr_end, 'pool': pool, 'method': 'fallback_window_unfilled',
                               'n_candidates': max_records_seen, 'group_size': 1,
                               'selected_read': percentile_pick,
                               'percentile_pick': percentile_pick, 'changed': False})
                continue

            # --- guard rail 3: a FULL window of candidates that still disagree.
            # Nothing here is explainable by sampling, so stop rather than guess.
            if chosen is None:
                raise RuntimeError(
                    f'{chr_end}: no group of {min_agree} mutually-agreeing reads after '
                    f'{max_widenings + 1} window size(s) (pool={pool}, thresholds '
                    f'pident>={min_identity} coverage>={min_coverage}). The candidates '
                    f'around the 75th percentile do not corroborate each other, so no '
                    f'read here can be trusted as a scaffold. Inspect this end before '
                    f'rerunning.')

            chr_end_to_read[chr_end] = chosen
            all_selected_read_ids.append(chosen)
            changed = chosen != percentile_pick
            rejected = [c for c in cands if c not in group]
            note = f'  [rejected: {", ".join(rejected)}]' if rejected else ''
            flag = '  <== differs from percentile pick' if changed else ''
            print(f'{chr_end}: {chosen} (agreeing group {len(group)}/{len(cands)}, '
                  f'anchor {anchor[chosen]["pident"]:.2f}%){flag}{note}')
            report.append({'chr_end': chr_end, 'pool': pool, 'method': 'consensus',
                           'n_candidates': len(cands), 'group_size': len(group),
                           'window_n': window_n, 'widened': window_n != n_candidates,
                           'selected_read': chosen, 'percentile_pick': percentile_pick,
                           'changed': changed})

    n_changed = sum(1 for r in report if r['changed'])
    n_fallback = sum(1 for r in report if r['method'].startswith('fallback_'))
    print(f'\nSelected {len(chr_end_to_read)} reads for {len(chr_end_to_read)} chromosome ends')
    print(f'  consensus differed from the percentile pick at {n_changed} end(s)')
    print(f'  fell back to an uncorroborated read at {n_fallback} end(s)')

    with open(output_file, 'w') as f:
        for chr_end, read_id in chr_end_to_read.items():
            f.write(f'{chr_end}\t{read_id}\n')
    read_ids_only_file = output_file.replace('.txt', '_only.txt')
    with open(read_ids_only_file, 'w') as f:
        for read_id in all_selected_read_ids:
            f.write(f'{read_id}\n')
    if report_tsv:
        pd.DataFrame(report).to_csv(report_tsv, sep='\t', index=False)
        print(f'  per-end decisions written to {report_tsv}')

    return chr_end_to_read, all_selected_read_ids


def extract_selected_reads(reads_fastq, read_ids_file, output_fastq):
    """
    Extract selected reads from FASTQ file.

    Parameters:
    reads_fastq (str): Path to input FASTQ file
    read_ids_file (str): Path to file with read IDs (one per line)
    output_fastq (str): Path to output FASTQ file
    """
    print(f"Extracting selected reads from {reads_fastq}")

    # Read the read IDs
    with open(read_ids_file, 'r') as f:
        selected_read_ids = set(line.strip() for line in f if line.strip())

    # Index the FASTQ reads by ID
    print("Indexing FASTQ file...")
    fastq_index = SeqIO.to_dict(SeqIO.parse(reads_fastq, "fastq"))

    # Extract selected reads
    matched_reads = [fastq_index[read_id] for read_id in selected_read_ids if read_id in fastq_index]

    if not matched_reads:
        raise ValueError("No FASTQ reads matched the selected read IDs!")

    print(f"Writing {len(matched_reads)} selected reads to {output_fastq}")
    with open(output_fastq, "w") as out_f:
        SeqIO.write(matched_reads, out_f, "fastq")


def get_softclip_from_cigar(read, side="start"):
    """Extract soft-clipped sequence from read based on CIGAR string"""
    if not read.cigartuples:
        return ""

    if side == "start":
        if read.cigartuples[0][0] == 4:  # 4 = soft clip
            clip_length = read.cigartuples[0][1]
            return read.query_sequence[:clip_length]
    elif side == "end":
        if read.cigartuples[-1][0] == 4:  # 4 = soft clip
            clip_length = read.cigartuples[-1][1]
            return read.query_sequence[-clip_length:]

    return ""


def extend_reference_multi(bamfile, reference, read_ids_file, output_fasta, trim,
                           chr_arm_pairs_file=None):
    """
    Extend reference genome using soft-clipped bases from multiple reads

    Args:
        bamfile: Path to BAM file
        reference: Path to reference FASTA file
        read_ids_file: Path to file with read IDs (one per line)
        output_fasta: Path to output extended FASTA file
        trim: Number of bp to trim from the end of extensions to remove adapters
        chr_arm_pairs_file: "<chr_end>\\t<read_id>" pairs from select_reads. When
            supplied, each scaffold contributes ONLY to its own arm's side:
            an L-arm read's start-clip may become a prefix, an R-arm read's
            end-clip may become a suffix, and the opposite-side clip of each is
            ignored entirely.

    WHY ARM BINDING MATTERS
    ----------------------
    Each contig carries both arms, so extensions[contig] holds one 'prefix' and
    one 'suffix' list. Without arm awareness a clip is filed purely by which end
    of the CIGAR it sits on, and the winner is whichever is LONGEST. Arm and
    clip-side normally coincide geometrically (verified: all 32 arms map L->
    contig start, R-> contig end), but nothing enforces it, so a read clipped on
    BOTH sides lands in both lists and can win the other arm's slot.

    That is what corrupted 7172 chr16R: the chr16L scaffold is a fold-back whose
    spurious 82,350bp end-clip beat chr16R's legitimate 6,814bp clip, grafting
    chr16L's inverted subtelomere onto chr16R.

    Removing bad reads alone does not close this. Perfectly healthy R-arm reads
    carry a telomeric start-clip (1,848bp on 1R, 2,049bp on 9R) because the base
    reference is truncated before the telomeres -- larger than the legitimate
    prefix clips of chr7L (1,935bp) and chr15L (1,778bp). Binding by arm makes
    the wrong-side clip ineligible rather than merely out-competed.
    """
    # Read the read IDs
    with open(read_ids_file, 'r') as f:
        read_ids = [line.strip() for line in f if line.strip()]

    # read_id -> chr_end ("16L" / "16R"), used to bind clips to their own arm
    read_to_arm = {}
    if chr_arm_pairs_file:
        with open(chr_arm_pairs_file, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    read_to_arm[parts[1]] = parts[0]
        print(f"Arm binding ENABLED: {len(read_to_arm)} scaffold/arm pairs loaded")
    else:
        print("WARNING: no chr_arm pairs file supplied -- falling back to "
              "side-of-CIGAR filing (a read clipped on both sides can win the "
              "wrong arm's slot)")

    samfile = pysam.AlignmentFile(bamfile, "rb")
    ref_seqs = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))

    read_ids_set = set(read_ids)
    extensions = {}  # Store extensions per reference

    print(f"\nSearching for soft-clipped bases from {len(read_ids)} reads...")
    if trim != 0:
        print(f"Will trim {trim}bp from extension ends to remove adapters")

    # Collect all extensions from all specified reads
    for read in samfile.fetch(until_eof=True):
        if read.query_name not in read_ids_set or read.is_secondary or read.is_supplementary:
            continue

        ref_name = read.reference_name

        if ref_name not in extensions:
            extensions[ref_name] = {'prefix': [], 'suffix': []}

        # Which side is this read allowed to contribute to?
        # L arm -> contig start -> prefix only.  R arm -> contig end -> suffix only.
        # The opposite-side clip is never examined, so it cannot compete.
        arm = read_to_arm.get(read.query_name)
        if read_to_arm and arm is None:
            print(f"  WARNING: {read.query_name} has no arm assignment, skipping")
            continue
        allow_prefix = (arm is None) or arm.upper().endswith('L')
        allow_suffix = (arm is None) or arm.upper().endswith('R')

        # Soft-clipped bases at start (5' extension) -- L arms only
        if allow_prefix:
            prefix = get_softclip_from_cigar(read, side="start")
            if prefix and len(prefix) > trim:
                trimmed_prefix = prefix[:-trim] if trim > 0 else prefix
                extensions[ref_name]['prefix'].append((read.query_name, trimmed_prefix, len(prefix)))

        # Soft-clipped bases at end (3' extension) -- R arms only
        if allow_suffix:
            suffix = get_softclip_from_cigar(read, side="end")
            if suffix and len(suffix) > trim:
                trimmed_suffix = suffix[trim:] if trim > 0 else suffix
                extensions[ref_name]['suffix'].append((read.query_name, trimmed_suffix, len(suffix)))

    samfile.close()

    # Create extended reference(s)
    extended_records = []
    total_extensions = 0

    for ref_name, ext_data in extensions.items():
        orig_ref = str(ref_seqs[ref_name].seq)

        # With arm binding there is exactly one scaffold per arm, so each side
        # must have at most one candidate. More than one means something
        # upstream is broken (duplicate rows in the pairs file, or two reads
        # assigned the same arm) -- fail loudly rather than silently picking.
        if read_to_arm:
            for side in ('prefix', 'suffix'):
                if len(ext_data[side]) > 1:
                    who = ', '.join(f'{r}({len(s)}bp)' for r, s, _ in ext_data[side])
                    raise RuntimeError(
                        f'{ref_name}: {len(ext_data[side])} {side} candidates but arm '
                        f'binding permits exactly 1 -- {who}. Check '
                        f'selected_read_chr_arm_pairs.txt for duplicate arms.')

        # Use the longest prefix and suffix (after trimming) if multiple exist
        prefix = ""
        prefix_read = ""
        prefix_orig_len = 0
        if ext_data['prefix']:
            prefix_read, prefix, prefix_orig_len = max(ext_data['prefix'], key=lambda x: len(x[1]))
            print(f"Using 5' extension for {ref_name} from {prefix_read} "
                  f"[{read_to_arm.get(prefix_read, 'unbound')}]: added {len(prefix)}bp")

        suffix = ""
        suffix_read = ""
        suffix_orig_len = 0
        if ext_data['suffix']:
            suffix_read, suffix, suffix_orig_len = max(ext_data['suffix'], key=lambda x: len(x[1]))
            print(f"Using 3' extension for {ref_name} from {suffix_read} "
                  f"[{read_to_arm.get(suffix_read, 'unbound')}]: added {len(suffix)}bp")

        if prefix or suffix:
            new_ref = prefix + orig_ref + suffix

            # Create proper header
            new_id = f"{ref_name}_extended"
            description = ""
            if prefix:
                description += f" (added {prefix_orig_len}bp)"
            if suffix:
                description += f" (added {suffix_orig_len}bp)"

            new_record = SeqRecord(
                Seq(new_ref),
                id=new_id,
                description=description
            )

            extended_records.append(new_record)
            total_extensions += 1

    if extended_records:
        with open(output_fasta, "w") as out_fasta:
            SeqIO.write(extended_records, out_fasta, "fasta")

        print(f"\n{total_extensions} reference(s) extended and written to {output_fasta}")
        return True
    else:
        print(f"\nNo soft-clipped bases found near chromosome ends for any of the specified reads")
        # Copy original reference to output so pipeline can continue
        subprocess.run(["cp", reference, output_fasta], check=True)
        print(f"Copied original reference to {output_fasta}")
        return False
