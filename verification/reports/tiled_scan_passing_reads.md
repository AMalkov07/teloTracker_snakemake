# Directional tiled-window scan: detail on the 12 reads that pass

The 12 of 31 mid-Y' candidates (10 of 24 in 6991, 1 of 2 in 7172, 1 of 5 in 7302) that clear
the tiled-window scan's final confirmation bar (both halves >=90% global identity), with the
exact split point and the two half-alignments that drove the call. Method recap: non-overlapping
300bp windows from the Y' anchor boundary, each compared at the same absolute offset to the
recipient reference; once 3 consecutive windows score <99% against the recipient, an open
search across the library (excluding same-cut99-group elements) checks whether one single
element beats the recipient on all 3 -- if so, that is the split point (end of the 3rd window)
and that element is the donor tested below. First-half identity is the read's own arm globally
aligned (`edlib`, `NW`) against a same-length slice from the **start** of the recipient
reference; second-half identity is the read's remaining arm globally aligned against a
same-length slice from the **end** of the donor reference.

| strain | chr_end | read_id | own (recipient) | split point | read length | % native | first half vs recipient | donor found | second half vs donor | whole-read donor (for comparison) |
|---|---|---|---|---|---|---|---|---|---|---|
| 6991 | chr10L | SRR33298384.498901 | chr10L-1 | 6000 bp | 6594 bp | 91% | 90.00% | chr12R-2 | 98.48% | chr7R-1 |
| 6991 | chr10L | SRR33298461.57209 | chr10L-1 | 900 bp | 6582 bp | 14% | 92.22% | chr14R-1 | 93.80% | chr7R-1 |
| 6991 | chr10L | SRR33298373.545529 | chr10L-1 | 6000 bp | 6582 bp | 91% | 90.00% | chr12R-2 | 91.07% | chr16L-1 |
| 6991 | chr10L | SRR33298373.491903 | chr10L-1 | 900 bp | 6593 bp | 14% | 92.22% | chr14R-1 | 94.26% | chr16L-1 |
| 6991 | chr10L | SRR33298373.122658 | chr10L-1 | 900 bp | 6627 bp | 14% | 90.11% | chr14R-1 | 91.50% | chr12R-2 |
| 6991 | chr16R | SRR33298377.534605 | chr16R-1 | 2100 bp | 5176 bp | 41% | 98.05% | chr12L-1 | 97.01% | chr13L-1 |
| 6991 | chr16R | SRR33298384.248800 | chr16R-1 | 3300 bp | 5189 bp | 64% | 99.76% | chr12L-1 | 98.62% | chr14L-3 |
| 6991 | chr16R | SRR33298373.72733 | chr16R-1 | 4200 bp | 5183 bp | 81% | 94.43% | chr12L-1 | 91.96% | chr14L-3 |
| 6991 | chr5L | SRR33298434.226995 | chr5L-1 | 5700 bp | 6307 bp | 90% | 92.91% | chr14R-1 | 90.94% | chr12R-1 |
| 6991 | chr7R | SRR33298434.164214 | chr7R-1 | 5700 bp | 6587 bp | 87% | 93.05% | chr14R-1 | 98.08% | chr14R-1 (agrees) |
| 7172 | chr16R | SRR33298432.177252 | chr16R-1 | 4200 bp | 5194 bp | 81% | 94.55% | chr12L-1 | 98.49% | chr14L-3 |
| 7302 | chr14R | SRR33298452.246202 | chr14R-1 | 5400 bp | 6623 bp | 82% | 95.72% | chr12R-2 | 98.20% | chr14L-1 |

**Two distinct shapes among these 12**, as noted when this method was first run:

* **8 "late" breakpoints** (41-91% native) -- most of the Y' is native, a short telomere-distal
  tail switches to a specific donor. The more biologically plausible shape: a short localized
  gene-conversion tract.
* **4 "early" breakpoints, all at exactly 900 bp** (14% native) -- `SRR33298461.57209`,
  `SRR33298373.491903`, `SRR33298373.122658` (all chr10L), and none from 7172/7302. Given the
  99%-per-window trigger already fires this early on pure ONT sequencing noise in a third of
  native controls (see the main 6991 report), this shape should be read with more caution even
  though it does clear the final 90% bar on both sides.

**Only 1 of 12 (`SRR33298434.164214`, chr7R -> chr14R-1) agrees with the donor the whole-read
global-alignment test flagged earlier.** The other 11 land on a different element -- most often
`chr12L-1` (4x) or `chr14R-1` (3x) or `chr12R-2` (2x) -- because this test is asking what best
explains the short trailing segment specifically, not the whole molecule on average; a real
chimeric read can have its whole-molecule average pulled toward one group while a short embedded
tract points somewhere else entirely.
