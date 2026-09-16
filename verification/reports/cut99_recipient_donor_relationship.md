# Is there a fixed recipient -> donor relationship, or does it vary?

**There is a clear, fixed relationship for the high-volume recipients, and it varies for the
low-volume ones.** All 10 samples, excluding the 388 chr14L-1 reference-defect rows; "donor
group" uses the chr-end signature (stable across samples whose group labels differ).

## Concentration by recipient (sorted by sample size)

| recipient | n | dominant donor group | concentration |
|---|---|---|---|
| **chr16R** | 13 | {chr13L, chr14L} | **100%** -- every single instance |
| **chr6L** | 16 | {chr13L, chr14L} | 81% |
| **chr2L** | 12 | {chr8R} | 83% |
| **chr13L** | 41 | {chr2L, chr6L} | 71% |
| chr10L | 12 | {chr12R,chr14L,chr15R,chr16L,chr4R,chr7R} | 75% |
| chr5R | 7 | {chr10L, chr9L} | 71% |
| chr8R | 6 | {chr2L, chr6L} | 67% |
| chr14R | 24 | {chr10L, chr9L} | only 38% -- no dominant donor |
| chr12L, chr5L, chr16L, chr8L, chr7R, chr15R, chr14L | <=7 | scattered | 33-60%, too few events to read a pattern |

The four recipients with the most data (chr16R, chr6L, chr2L, chr13L; n=13-41) all have one
donor group taking 71-100% of their events. That is not what random mis-assignment across a
12-group partition would look like -- with 12 groups available, a recipient hitting the same
one 71-100% of the time is a real preference, not noise.

## A reciprocal pair: chr13L <-> chr6L/chr2L

chr13L's dominant donor is {chr2L, chr6L} (29 of 41). chr6L's dominant donor is {chr13L,
chr14L} (13 of 16). At the specific-element level this is even tighter than the group view
shows:

* chr13L's donor is named `chr2L-1` 15 times and `chr6L-1` 14 times -- almost exactly even
  between the two members of that group, not favouring one.
* chr6L's donor is named `chr14L-5` 9 of 16 times, specifically -- not spread evenly across
  {chr13L-1, chr14L-3/4/5}.

So the relationship is not simply "these two groups are near-identical and get confused
symmetrically" -- chr6L's donor calls concentrate on one specific element (chr14L-5) within its
target group, which a purely symmetric identity-confusion model would not produce on its own.

## chr16R: the cleanest case

All 13 chr16R events point to the {chr13L, chr14L} group, split as chr13L-1 (5), chr14L-5 (4),
chr14L-3 (3), chr14L-4 (1) -- spread across the group's members but never landing outside it,
across multiple independent samples. This is the strongest single-recipient signal in the
dataset.

## chr2L -> chr8R-1 specifically

10 of 12 chr2L events name `chr8R-1` exactly (not spread across a larger group -- chr8R-1 sits
alone in a singleton group in most samples, so this is an unambiguous, specific relationship,
not a group-level approximation).

## What this does not establish

* **Directionality is not resolved by this table.** A recipient's dominant donor group and that
  group's own dominant "donor" (when it appears as a recipient) need not be the same relationship
  read in reverse -- these counts are recipient-centric, not symmetric confusion matrices.
* **Sample size varies widely** (2 to 41), so the low-n recipients' "no clear pattern" may
  simply reflect too little data, not genuine randomness.
* **This does not distinguish clonal recurrence from a true hotspot** -- the earlier
  heterogeneity analysis showed at least one of these relationships (chr13L<->chr2L) recurs
  across independent sequencing runs, consistent with a real, possibly standing, recombination
  relationship rather than symmetric noise between similar sequences.
