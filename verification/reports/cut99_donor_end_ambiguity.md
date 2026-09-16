# The donor column names an element, not a confirmed source — this matters for the counts

Good catch: the donor tables so far (`cut99_end_mismatch_counts.md`, and the earlier
recipient/donor pattern analyses) counted `donor_element`'s chromosome end as if it were the
source. It is not. `donor_element` is simply the single best-scoring member of whatever cut99
group the read's Y' matched — and **126 of 160 donor calls (79%) point to a group with more
than one member**, so a specific end was named in cases where the real answer is "one of
several possible ends."

## How often is this a problem, concretely

| | mismatched copies |
|---|---|
| donor group has exactly 1 end (name is safe) | 34 |
| donor group spans >1 end (name is arbitrary) | **126** |

Example: a `6991_day0` read's Y' was named `chr7R-1` as the donor. But at that sample's cut99
grouping, `chr7R-1` sits in a 17-end group with `chr12R-2..7`, `chr14L-1/2`, `chr15R-1`,
`chr16L-1` and `chr4R-1..7`. The evidence supports "the donor is something in this group of 18
elements," not "the donor is chr7R specifically."

## Corrected donor-end counts

Two more defensible countings, both excluding the 388 chr14L-1 reference-defect rows:

* **fractional credit** — split each mismatch's donor credit evenly across every end sharing
  the matched group (so the total still sums to 160, and an end that is *only ever* a
  co-member gets partial rather than full credit)
* **possible-donor count** — every end that could have been the source gets full credit (this
  double- and triple-counts a single event across all its candidate ends, so it does not sum to
  160, but it answers "how often is this end a plausible source at all")

| end | named donor (old, wrong) | fractional credit | possible-donor count |
|---|---|---|---|
| chr14L | 29 | 23.3 | 64 |
| chr2L | 25 | 20.0 | 40 |
| chr6L | 15 | 20.0 | 40 |
| chr13L | 15 | 18.9 | 40 |
| chr8R | 14 | 14.0 | 14 |
| chr9L | 13 | 11.0 | 22 |
| chr10L | 9 | 11.0 | 22 |
| chr7R | 9 | 4.8 | 27 |
| chr12R | 5 | 9.1 | 30 |
| chr14R | 5 | 5.0 | 5 |
| chr16L | 7 | 4.8 | 27 |
| chr5R | 4 | 4.0 | 4 |
| chr15R | 3 | 4.1 | 25 |
| chr4R | 1 | 4.1 | 25 |
| chr16R | 2 | 2.0 | 2 |
| chr5L | 2 | 2.0 | 2 |
| chr8L | 2 | 2.0 | 2 |

**Read the fractional-credit column as the corrected version of the earlier donor table.**
chr14L, chr2L, chr6L and chr13L remain the leading donors under this correction too, though
their margins narrow (chr14L 29→23.3, chr2L 25→20.0) and the ranking among the top few
tightens.

**Read the possible-donor count as a different question**: how often is this end *even in
consideration*. chr14L jumps to 64 here because it is a 5-copy end that co-occurs with the
large chr4R/chr12R/chr7R/chr16L group in several samples — its high count partly reflects
being a member of a big shared group, not necessarily being picked often.

## Which counts, unambiguously, either way

`chr14R-1`, `chr5R-1`, `chr16R-1`, `chr5L-1`, `chr8L-1` are all singletons in their groups in
every sample examined — their donor counts (5, 4, 2, 2, 2) are genuine, not group artefacts.

## Recommendation

The recipient-end counts (`as_recipient`) are unaffected by this issue — a read's own anchor
end is directly observed, not inferred through a group. Only the donor side needs the
correction above. Anywhere a specific donor chromosome end is quoted from these tables going
forward, it should carry the group's other members alongside it, or use the fractional-credit
number rather than the raw named-donor count.
