"""
yprime_path.py -- explain a read's gained Y' array as an ordered path of donor
segments, using Y' IDs first and ITS lengths as the tie-breaker.

Tokens. A Y' array is a sequence of (id, its) tokens: the Y' variant ID and the
length of the interstitial telomeric sequence (ITS) that follows it, or None for
the last copy. Reference tokens come from the day-0 BED (ITS_<end>_Y_Prime_n-m
features) + the Y' library (ID per position); read tokens from the RepeatMasker
hit coordinates (gap between consecutive hits), which reproduce the reference
ITS to +-2 bp.

Templates. Every contiguous piece R[j..k] of an end's array can donate
  linear   : R[j], R[j+1], ..., R[k]           (BIR / gene conversion)
  circular : the piece excised as a circle and copied in tandem,
             R[j..k] R[j..k] ... starting at any phase; the ITS at the
             wrap-around junction is unconstrained (it is new sequence)
Two tokens match when the IDs are equal and, for every ITS that is internal to
the template, |its_read - its_ref| <= ITS_TOL.

Parse. Dynamic programme over the read tokens: cover them with the fewest
segments (lexicographic: n_segments, n_ambiguous, -n_ITS_verified), each
segment being one template match. Per segment the donor is the end that gives
that match; several ends giving the same match are reported together as an
ambiguous segment. The read's own end is a donor like any other, and a segment
explained by the own end is labelled 'self'.
"""
import re

ITS_TOL = 8          # bp; read ITS gaps reproduce the reference to +-2 bp (occasional +-10 outliers)
MAX_TOKENS = 40


# ---------------------------------------------------------------------------
# tokens
# ---------------------------------------------------------------------------

def read_tokens_from_positions(pos_str, telo_side):
    """'ID2:100-5600;ID1:5763-12400' -> [('ID2', 163), ('ID1', None)] in anchor-to-telomere order."""
    hits = []
    for p in str(pos_str or '').split(';'):
        if ':' not in p:
            continue
        yid, rng = p.rsplit(':', 1)
        m = re.match(r'(-?\d+)-(-?\d+)$', rng)
        if not m:
            continue
        hits.append((int(m.group(1)), int(m.group(2)), yid))
    return read_tokens_from_hits(hits, telo_side)


def read_tokens_from_hits(hits, telo_side):
    """hits: iterable of (start, end, id) on the read (any order)."""
    hits = sorted(hits)
    if not hits:
        return []
    gaps = [hits[i + 1][0] - hits[i][1] for i in range(len(hits) - 1)]
    ids = [h[2] for h in hits]
    if telo_side == 'beginning':          # read runs telomere -> anchor; flip to anchor -> telomere
        ids, gaps = ids[::-1], gaps[::-1]
    return [(ids[i], gaps[i] if i < len(gaps) else None) for i in range(len(ids))]


def build_reference_tokens(bed_path, location_to_id):
    """{chr_end: [(id, its), ...]} for every end with Y' in the BED.
    location_to_id: {(chr_end, pos): id} from the Y' library headers."""
    yp, its = {}, {}
    for line in open(bed_path):
        p = line.rstrip('\n').split('\t')
        if len(p) < 4:
            continue
        name = p[3]
        m = re.match(r'(chr\d+[LR])_Y_Prime_(\d+)$', name)
        if m:
            yp.setdefault(m.group(1), {})[int(m.group(2))] = None
            continue
        m = re.match(r'ITS_(chr\d+[LR]?)_Y_Prime_(\d+)-(\d+)$', name)
        if m:
            length = int(p[5]) if len(p) > 5 and p[5].isdigit() else int(p[2]) - int(p[1]) + 1
            its.setdefault(m.group(1), {})[(int(m.group(2)), int(m.group(3)))] = length
    out = {}
    for ce, positions in yp.items():
        toks = []
        n = max(positions)
        for pos in range(1, n + 1):
            yid = location_to_id.get((ce, pos), f'{ce}_Y_Prime_{pos}')
            gap = its.get(ce, {}).get((pos, pos + 1)) if pos < n else None
            toks.append((yid, gap))
        out[ce] = toks
    return out


# ---------------------------------------------------------------------------
# matching
# ---------------------------------------------------------------------------

def _its_ok(read_its, ref_its):
    if read_its is None or ref_its is None:
        return True, 0
    return abs(read_its - ref_its) <= ITS_TOL, 1


def match_template(tokens, i, ref, j, k, circular, phase=0):
    """Match tokens[i:] against piece ref[j..k] (linear) or its tandem repeat
    (circular). Returns (n_tokens_consumed, n_its_verified).

    Circular junction ITS (the ITS between one copy of the circle and the
    next). A rolling-circle product carries the SAME junction ITS at every
    repeat, so the junctions inside the read must agree with each other
    (+-ITS_TOL). A circle excised from a tandem array by recombination between
    two of its ITS keeps one of the donor's ITS as its junction, so a junction
    that equals the donor ITS flanking the piece (ref[k].its or ref[j-1].its)
    counts as verified, which lets the donor win over ends that merely carry
    the same IDs."""
    n, L = len(tokens), k - j + 1
    consumed, verified = 0, 0
    t = 0
    junctions = []
    while i + t < n:
        r = j + ((phase + t) % L) if circular else j + t
        if not circular and r > k:
            break
        rid, rits = ref[r]
        tid, tits = tokens[i + t]
        if tid != rid:
            break
        at_piece_end = (r == k)
        if i + t + 1 < n:
            if circular and at_piece_end:
                if tits is not None:
                    if junctions and abs(tits - junctions[0]) > ITS_TOL:
                        consumed = t + 1                      # inconsistent junction: the repeat stops here
                        return consumed, verified
                    junctions.append(tits)
            elif not circular and at_piece_end:
                pass                                          # last copy of a linear piece: no constraint
            else:
                ok, v = _its_ok(tits, rits)
                if not ok:
                    consumed = t + 1
                    return consumed, verified
                verified += v
        t += 1
        consumed = t
        if not circular and r == k:
            break
    if circular and junctions:
        flank = [x for x in (ref[k][1], ref[j - 1][1] if j > 0 else None) if x is not None]
        if any(abs(junctions[0] - f) <= ITS_TOL for f in flank):
            verified += len(junctions)
    return consumed, verified


# ---------------------------------------------------------------------------
# parse
# ---------------------------------------------------------------------------

def _candidates(tokens, ref_tokens):
    """Per start position: list of (L, verified, end, mode, unit_ids, piece(j,k) 1-based)."""
    n = len(tokens)
    cands = [[] for _ in range(n)]
    for i in range(n):
        for ce, ref in ref_tokens.items():
            m = len(ref)
            for j in range(m):
                L, v = match_template(tokens, i, ref, j, m - 1, False)
                if L >= 1:
                    cands[i].append((L, v, ce, 'lin', tuple(x[0] for x in ref[j:j + L]), (j + 1, j + L)))
                for k in range(j, m):
                    piece = k - j + 1
                    for phase in range(piece):            # a circle can insert starting at any of its copies
                        Lc, vc = match_template(tokens, i, ref, j, k, True, phase)
                        if Lc > piece and (piece >= 2 or Lc >= 3):   # must wrap; a 1-copy circle needs >= 3 copies
                            unit = tuple(x[0] for x in ref[j:k + 1])
                            unit = unit[phase:] + unit[:phase]
                            cands[i].append((Lc, vc, ce, 'circ', unit, (j + 1, k + 1)))
    return cands


def _best_by_length(cands_i):
    by_L = {}
    for L, v, ce, mode, unit, piece in cands_i:
        by_L.setdefault(L, []).append((v, ce, mode, unit, piece))
    # parsimony: when a linear (verbatim) match explains a segment as well as a
    # circle does, the circle explanations are dropped
    for L, lst in by_L.items():
        vmax = max(v for v, *_ in lst)
        if any(mode == 'lin' and v == vmax for v, ce, mode, unit, piece in lst):
            by_L[L] = [x for x in lst if not (x[2] == 'circ' and x[0] == vmax)]
    return by_L


def parse_path(tokens, ref_tokens, self_end):
    """Return {'segments': [...], 'n_segments', 'primary_donor', 'primary_len',
    'circles': [...], 'its_verified', 'its_checked'}.

    Segment fields: start, len, ids, donor, tag (self | donor | ambiguous |
    tentative | unk), mode (lin | circ), unit, piece (1-based copy range in the
    donor), repeats, its_verified, circle_support (strong | moderate | weak),
    alt (an equally plausible linear reading of the same tokens, if any)."""
    n = len(tokens)
    empty = {'segments': [], 'n_segments': 0, 'primary_donor': '', 'primary_len': 0, 'circles': [],
             'its_verified': 0, 'its_checked': 0}
    if n == 0 or n > MAX_TOKENS:
        return empty
    best_at = [_best_by_length(c) for c in _candidates(tokens, ref_tokens)]

    INF = (10 ** 6, 10 ** 6, 0)
    best = [INF] * (n + 1)
    choice = [None] * (n + 1)
    best[n] = (0, 0, 0)
    for i in range(n - 1, -1, -1):
        options = list(best_at[i].items()) or [(1, [(0, '?', 'unk', (tokens[i][0],), (0, 0))])]
        for L, lst in options:
            vmax = max(v for v, *_ in lst)
            donors = sorted({ce for v, ce, mode, unit, piece in lst if v == vmax})
            amb = 0 if (len(donors) == 1 or self_end in donors) else 1
            nxt = best[i + L]
            cost = (1 + nxt[0], amb + nxt[1], -vmax + nxt[2])
            if cost < best[i]:
                best[i] = cost
                choice[i] = (L, vmax, lst)

    segments, i = [], 0
    while i < n:
        L, vmax, lst = choice[i]
        top = [x for x in lst if x[0] == vmax]
        donors = sorted({ce for v, ce, mode, unit, piece in top})
        seg_ids = [t[0] for t in tokens[i:i + L]]
        if self_end in donors:
            donor, tag = self_end, 'self'
        elif len(donors) == 1:
            donor, tag = donors[0], 'donor'
        else:
            donor, tag = '|'.join(donors), 'ambiguous'
        pick = donor if tag != 'ambiguous' else donors[0]
        v, ce, mode, unit, piece = next(x for x in top if x[1] == pick)
        if tag == 'unk':
            mode, unit, piece = 'unk', (seg_ids[0],), (0, 0)
        seg = {'start': i, 'len': L, 'ids': seg_ids, 'donor': donor, 'tag': tag, 'mode': mode,
               'unit': list(unit), 'piece': piece, 'its_verified': vmax,
               'repeats': (L / len(unit)) if (mode == 'circ' and unit) else 1.0,
               'circle_support': '', 'alt': ''}
        # a piece made of ONE Y' copy (linear or circle) is only ever suggestive:
        # the same ID with a matching ITS could come from any end that carries it
        if tag == 'donor' and (len(unit) == 1 if mode == 'circ' else L == 1):
            the_id = unit[0] if mode == 'circ' else seg_ids[0]
            n_ends_with_id = sum(1 for ce2, r2 in ref_tokens.items() if the_id in [x[0] for x in r2])
            if n_ends_with_id > 1:
                seg['tag'] = 'tentative'
        if mode == 'circ' and tag in ('donor', 'self') and donor in ref_tokens:
            # how much better is the circle reading than a linear one?
            #   weak     : ONE verbatim piece of the donor explains the same copies as well
            #   strong   : more copies than the donor array holds (no single copy can do that)
            #   moderate : otherwise; a two-piece linear reading, when one exists, is
            #              reported as 'alt' (two events, or a whole-array circle)
            arr_len = len(ref_tokens[donor])
            lin = parse_path_linear_only(tokens[i:i + L], {donor: ref_tokens[donor]}, donor)
            alt = ' + '.join(f"{donor}[{a}-{b}]" if a != b else f"{donor}[{a}]" for a, b in lin['pieces'])
            if lin['n_segments'] == 1 and lin['unexplained'] == 0 and lin['its_verified'] >= vmax - 1:
                seg['circle_support'] = 'weak'
                seg['alt'] = alt
            elif L > arr_len:
                seg['circle_support'] = 'strong'
            else:
                seg['circle_support'] = 'moderate'
                if lin['n_segments'] == 2 and lin['unexplained'] == 0:
                    seg['alt'] = alt
        segments.append(seg)
        i += L
    circles = [s for s in segments if s['mode'] == 'circ']
    primary = max(segments, key=lambda s: (s['len'], s['its_verified']))
    primary_donor = (primary['donor'] if primary['tag'] == 'donor' else
                     'self' if primary['tag'] == 'self' else '')
    checked = sum(1 for t in tokens[:-1] if t[1] is not None)
    return {'segments': segments, 'n_segments': len(segments), 'primary_donor': primary_donor,
            'primary_len': primary['len'], 'circles': circles,
            'its_verified': sum(s['its_verified'] for s in segments), 'its_checked': checked}


def parse_path_linear_only(tokens, ref_tokens, self_end):
    """Fewest linear pieces (no circles) explaining tokens; used to judge circle support."""
    n = len(tokens)
    cands = [[c for c in cs if c[3] == 'lin'] for cs in _candidates(tokens, ref_tokens)]
    INF = (10 ** 6, 10 ** 6, 0)
    best = [INF] * (n + 1); choice = [None] * (n + 1); best[n] = (0, 0, 0)
    for i in range(n - 1, -1, -1):
        opts = {}
        for L, v, ce, mode, unit, piece in cands[i]:
            if v > opts.get(L, (-1, None))[0]:
                opts[L] = (v, piece)
        if not opts:
            opts = {1: (0, None)}
        for L, (v, piece) in opts.items():
            nxt = best[i + L]
            cost = (1 + nxt[0], (piece is None) + nxt[1], -v + nxt[2])
            if cost < best[i]:
                best[i] = cost; choice[i] = (L, v, piece)
    pieces, i, unexplained, verified = [], 0, 0, 0
    while i < n:
        L, v, piece = choice[i]
        if piece is None:
            unexplained += 1
        else:
            pieces.append(piece); verified += v
        i += L
    return {'n_segments': len(pieces) + unexplained, 'pieces': pieces, 'unexplained': unexplained, 'its_verified': verified}


def _donor_label(seg):
    d = 'self' if seg['tag'] == 'self' else seg['donor']
    if seg['tag'] == 'tentative':
        d += '?'
    if seg['piece'] != (0, 0) and '|' not in seg['donor']:
        a, b = seg['piece']
        d += f'[{a}]' if a == b else f'[{a}-{b}]'
    return d


def format_path(parsed):
    parts = []
    for s in parsed.get('segments', []):
        core = ','.join(s['ids'])
        if s['tag'] == 'unk':
            parts.append(f"?:{core}")
        elif s['mode'] == 'circ':
            txt = f"{_donor_label(s)}:{core}(circ x{s['repeats']:.1f}"
            if s['circle_support']:
                txt += f" {s['circle_support']}"
            if s['alt']:
                txt += f" | alt {s['alt']}"
            parts.append(txt + ')')
        else:
            parts.append(f"{_donor_label(s)}:{core}")
    return ' > '.join(parts)


def format_circles(parsed):
    return ';'.join(f"{_donor_label(s)}:{','.join(s['unit'])}x{s['repeats']:.1f}:{s['circle_support'] or 'na'}"
                    for s in parsed.get('circles', []) if s['tag'] != 'tentative')
