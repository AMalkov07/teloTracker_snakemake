"""Input files the pipeline needs before it starts, checked up front.

Shared by run_pipeline.py (before it calls Snakemake) and the Snakefile (for direct `snakemake`
runs), so both find a sample's input file the same way and report the same missing day-0 files.
"""
import glob
import os

# A sample's input file is <input dir>/<sample><extension>; the first match wins, in this order.
READ_INPUTS = [(".bam", "bam"), (".fastq", "fastq"), (".fq", "fastq"),
               (".fastq.gz", "fastqgz"), (".fq.gz", "fastqgz")]


def find_sample_input(input_dir, sample, candidates=READ_INPUTS):
    """(path, type) of the sample's input file, or None when there is none."""
    for ext, kind in candidates:
        path = os.path.join(input_dir, sample + ext)
        if os.path.exists(path):
            return path, kind
    return None


def tried_inputs(input_dir, sample, candidates=READ_INPUTS):
    return ", ".join(os.path.join(input_dir, sample + ext) for ext, _ in candidates)


LABEL_HINT = "made by label_regions"


def missing_day0_inputs(day0_ref, day0_bed, y_prime_lib, spacer_lib_dir, x_element_lib):
    """[(what, path, how it is made)] for every day-0 file recombination cannot run without.

    Without this check a missing reference, label BED or Y' library surfaces only as a Snakemake
    MissingInputException once a sample's run has started, and a missing spacer or X library not
    at all: every end is then marked not_in_library and the run still succeeds."""
    problems = []
    if not os.path.isfile(day0_ref):
        problems.append(("day-0 reference", day0_ref,
                         "made by create_ref, or set reference_fasta to a curated genome"))
    for what, path, hint in (("day-0 label BED", day0_bed, LABEL_HINT),
                             ("day-0 Y' library", y_prime_lib, LABEL_HINT + ", or y_prime_lib_override"),
                             ("day-0 X element library", x_element_lib, LABEL_HINT)):
        if not os.path.isfile(path):
            problems.append((what, path, hint))
    if not glob.glob(os.path.join(spacer_lib_dir, "*.fasta")):
        problems.append(("day-0 spacer library (a folder of .fasta files)", spacer_lib_dir, LABEL_HINT))
    return problems
