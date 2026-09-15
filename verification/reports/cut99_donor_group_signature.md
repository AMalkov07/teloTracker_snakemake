# Donor group counts, pooled by chromosome-end signature (all samples)

Cut99 group LABELS (G1, G2 ...) are local to each sample and not comparable. But the SET
of chromosome ends in a group is a stable signature -- if two samples both group
{chr7R, chr14L, chr16L} together, that is the same real grouping even if one calls it G8
and the other G9. This pools donor counts by that signature across all ten samples.

Excludes the 388 6991_day0_with_selection chr14L-1 reference-defect rows.

| chr-end signature (donor group) | times used as donor | seen in N samples |
|---|---|---|
| {chr2L, chr6L} | 40 | 8 |
| {chr13L, chr14L} | 37 | 10 |
| {chr12R, chr14L, chr15R, chr16L, chr4R, chr7R} | 22 | 6 |
| {chr10L, chr9L} | 22 | 6 |
| {chr8R} | 14 | 7 |
| {chr12R} | 5 | 2 |
| {chr14R} | 5 | 2 |
| {chr5R} | 4 | 3 |
| {chr12R, chr13L, chr14L, chr15R, chr16L, chr4R, chr7R} | 3 | 1 |
| {chr14L, chr16L, chr7R} | 2 | 1 |
| {chr5L} | 2 | 1 |
| {chr16R} | 2 | 2 |
| {chr8L} | 2 | 2 |