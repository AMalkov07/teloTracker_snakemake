#!/usr/bin/env python3
"""
filter_reads.py — Filter FASTQ reads by minimum length and minimum qscore.

Reads are kept only if:
  - Length >= MIN_LENGTH
  - qs:f:<value> tag in the header is >= MIN_QSCORE (if the tag is present)

Reads whose headers don't carry a qs:f: tag (e.g. older Dorado basecalls or
FASTQs that lost the tag during a BAM conversion) are NOT filtered by qscore
— they pass the qscore check and are still subject to the length filter.
A summary of how many reads lacked the tag is printed at the end so the
missing metadata isn't silent.

Usage:
    python filter_reads.py <input.fastq> <output.fastq> [min_length] [min_qscore]

Defaults:
    min_length = 2000
    min_qscore = 10
"""

import sys

input_file  = sys.argv[1]
output_file = sys.argv[2]
min_length  = int(sys.argv[3])   if len(sys.argv) > 3 else 2000
min_qscore  = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0

kept = 0
skipped_length = 0
skipped_qscore = 0
missing_qs_tag = 0

with open(input_file) as infile, open(output_file, "w") as outfile:
    while True:
        header = infile.readline()
        if not header:
            break
        seq  = infile.readline()
        plus = infile.readline()
        qual = infile.readline()

        # Parse qs:f: tag from header
        qs_value = None
        for tag in header.split():
            tag = tag.strip()
            if tag.startswith("qs:f:"):
                try:
                    qs_value = float(tag[5:])
                except ValueError:
                    pass
                break

        # Filter by qscore (only when a tag was present)
        if qs_value is None:
            missing_qs_tag += 1
        elif qs_value < min_qscore:
            skipped_qscore += 1
            continue

        # Filter by length
        read_len = len(seq.strip())
        if read_len < min_length:
            skipped_length += 1
            continue

        outfile.write(header)
        outfile.write(seq)
        outfile.write(plus)
        outfile.write(qual)
        kept += 1

total = kept + skipped_length + skipped_qscore
print(f"  Total reads:        {total}")
print(f"  Kept:               {kept}")
print(f"  Skipped (qscore):   {skipped_qscore}  [qs < {min_qscore}]")
print(f"  Skipped (length):   {skipped_length}  [len < {min_length}bp]")
if missing_qs_tag:
    print(f"  WARNING: {missing_qs_tag} reads had no qs:f: tag in their header — "
          f"they bypassed the qscore filter. Re-basecall with a Dorado version that "
          f"emits qs tags if strict quality filtering is required.")
