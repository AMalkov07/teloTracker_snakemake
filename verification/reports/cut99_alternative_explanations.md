# What else explains a Y' that misses its group but shows weak or no recombination evidence?

Across 6991, 7172 and 7302 (excluding the 388 reference-defect rows): 61 strong, 27 weak,
32 FAILS, 40 no-junction. Below, each candidate explanation with whether the data supports it.

## Supported by the data

**1. RepeatMasker / BLAST disagreement — the best-supported explanation for FAILS.**
The group call comes from RepeatMasker; BLAST is an independent check. Re-aligning the whole
Y' region with BLAST:

| category | whole-Y' BLAST favours the donor | favours the expected element |
|---|---|---|
| FAILS | 18 | **13** |
| no junction | 32 | 6 |

For **13 of 31 FAILS, BLAST says the read is simply normal** — and by clear margins
(97.66 % vs 94.76 %, 96.42 % vs 94.37 %, 98.79 % vs 96.45 %). Those reads are not recombinant
and not ambiguous; RepeatMasker mis-assigned the copy and BLAST disagrees. A methodological
artefact, not biology.

**2. Whole-element conversion, or a junction that cannot be localised.**
32 of 38 no-junction reads have BLAST agreeing the donor wins across the whole element. So the
donor call is real; what is missing is a junction. That is either genuine conversion of the
entire Y', or a junction hidden because donor and recipient are too alike to place it.

**3. Reference defects — the largest single cause overall.**
388 of 6991's 530 mismatches. A mis-assembled or mis-bounded reference element splits into its
own group and then mis-sorts every read at that end. Confirmed twice: the chr16L boundary
over-extension (77 bp of X-element sequence) and chr14L-1 in `6991_day0_with_selection`
(5,720 bp vs 6,654).

**4. Read quality — a modest contributor, not a main one.**
Best identity to any reference: strong median 99.24 %, FAILS median 98.74 %. Below 97 %:
6 of 32 FAILS against 3 of 61 strong. Real but small.

## Tested and NOT supported

**5. "The junction is near an element end, so one flank is too short to measure."**
Not supported. Median smaller half: strong 2,368 bp, weak 2,400 bp — essentially identical.
Only 2 weak reads have a flank under 500 bp. Weak calls are not short-flank calls.

**6. "Donor and recipient are too near-identical to localise a junction."**
Not supported as a general explanation. Median pair-homology identity: strong 97.12 %,
no-junction 97.95 %, weak 97.88 % — barely separated. It was true for one specific case
(7302 chr2L-1/chr8R-1 at 99.14 %, behind both of that strain's no-junction reads) but it does
not generalise.

## Not tested — remaining possibilities

* **Population heterogeneity.** Day-0 populations are not clonal. A subpopulation carrying a
  Y' variant absent from the reference would leave neither candidate correct, which is what a
  FAILS looks like. This is the most plausible untested explanation for the 18 FAILS where
  BLAST does favour the donor but no junction resolves.
* **Library chimerism.** A prep chimera joins at homology and would mimic recombination.
* **Mis-anchoring.** If a read were assigned to the wrong chromosome end, its "expected"
  element would be wrong by construction.
* **Recombination in the spacer or X element.** Never tested — reads were confirmed to have a
  native anchor, but the region between anchor and Y' was not examined.
* **The internal tandem repeat.** Every Y' element carries a repeat ~1,610 bp from its
  telomere-distal end, which fragments alignments locally and can let a wrong reference win.

## Practical reading

Of the 99 non-strong, non-defect rows, roughly **13 are demonstrable method artefacts**,
**32 are real donor calls lacking a localisable junction**, and the rest are unresolved. The
61 "strong" calls are the only set I would treat as established mid-Y' recombination.
