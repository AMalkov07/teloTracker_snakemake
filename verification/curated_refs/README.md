# Curated references (copied from TeloTracker/references/<strain>_features/ on Argon)

One local correction:

* `7172_features/repeatmasker_7172_all_y_primes.fasta`: header `>Y_Prime_chr5R11` -> `>Y_Prime_chr5R1`.
  chr5R carries exactly one Y' in both the curated 7172 BED and our 7172 day-0 reference, and the
  same entry is spelled `chr5R1` in the 6991 and 7302 curated libraries. Left as `chr5R11` the
  library declares a Y' at chr5R position 11 and none at position 1, which makes the strict-library
  check abort (a BED Y' that cannot be resolved silently turns every read at that end into a
  "1st Y' Change"). The file on Argon is unchanged.
