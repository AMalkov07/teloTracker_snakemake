#!/usr/bin/env python3

import os
import subprocess
import logging
from multiprocessing import Pool
from functools import partial
import argparse


def setup_logger(log_file="repeatmasker.log"):
    logger = logging.getLogger("RepeatMaskerRunner")
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger


def run_repeatmasker_single(chr_num, side, base_name_chr_anchor_outputs_dir, base_name,
                             anchor_set, repeatmasker_y_prime_dir, strain_number, logger):
    try:
        tag = f"chr{chr_num}{side}"
        logger.info(f"[{tag}] Starting RepeatMasker run")

        query = os.path.join(
            base_name_chr_anchor_outputs_dir,
            f"{base_name}_blasted_{anchor_set}_{tag}_anchor_reads"
        )
        database = f"repeatmasker_{strain_number}_all_y_primes"
        repeatmasker_file_name = os.path.join(
            repeatmasker_y_prime_dir,
            f"{base_name}_{tag}_repeatmasker_results.ssv"
        )
        fasta_file = f"{query}.fasta"
        out_file = f"{fasta_file}.out"

        subprocess.run([
            "RepeatMasker", fasta_file,
            "-lib", f"{database}.fasta",
            "-s", "-pa", "8",
            "--cutoff", "1000",
            "-no_is", "-norna", "-gff",
            "-dir", repeatmasker_y_prime_dir
        ], check=True)

        header = (
            "SW_score divergence_percent deletion_percent insertion_percent read_id "
            "match_start_on_read match_end_on_read leftover_on_read strand y_prime_id "
            "y_prime_group match_start_on_y_prime match_end_on_y_prime leftover_on_y_prime "
            "match_id sub_match"
        )
        with open(repeatmasker_file_name, "w") as out_f:
            out_f.write(header + "\n")

        filtered_file = os.path.join(repeatmasker_y_prime_dir, f"{tag}_filtered.out")
        with open(out_file, "r") as infile, open(filtered_file, "w") as filt:
            for line in infile:
                if "Y_Prime" in line:
                    filt.write(line)

        with open(filtered_file, "r") as filt, open(repeatmasker_file_name, "a") as out_f:
            out_f.writelines(filt)

        logger.info(f"[{tag}] Completed successfully")
        return (chr_num, side, "OK")

    except subprocess.CalledProcessError as e:
        logger.error(f"[{tag}] RepeatMasker failed: {e}")
        return (chr_num, side, "FAIL")
    except Exception as e:
        logger.exception(f"[{tag}] Unexpected error")
        return (chr_num, side, "FAIL")


def run_repeatmasker_parallel(base_name_chr_anchor_outputs_dir, base_name,
                              anchor_set, repeatmasker_y_prime_dir,
                              strain_number, max_parallel_jobs=4,
                              log_file="repeatmasker.log"):
    logger = setup_logger(log_file)

    jobs = [(chr_num, side) for chr_num in range(1, 18) for side in ['L', 'R']]

    fn = partial(
        run_repeatmasker_single,
        base_name_chr_anchor_outputs_dir=base_name_chr_anchor_outputs_dir,
        base_name=base_name,
        anchor_set=anchor_set,
        repeatmasker_y_prime_dir=repeatmasker_y_prime_dir,
        strain_number=strain_number,
        logger=logger
    )

    logger.info(f"Launching RepeatMasker jobs in parallel using {max_parallel_jobs} workers...")

    with Pool(processes=max_parallel_jobs) as pool:
        results = pool.starmap(fn, jobs)

    failed = [f"{chr}{side}" for chr, side, status in results if status != "OK"]
    if failed:
        logger.warning(f"RepeatMasker failed on: {', '.join(failed)}")
    else:
        logger.info("All RepeatMasker jobs completed successfully.")


def main():
    parser = argparse.ArgumentParser(description="Run RepeatMasker in parallel across chr1–chr17 (L/R) combinations.")
    parser.add_argument("--base-dir", required=True,
                        help="Path to base_name_chr_anchor_outputs_dir")
    parser.add_argument("--base-name", required=True,
                        help="Prefix for input/output files")
    parser.add_argument("--anchor-set", required=True,
                        help="Anchor set name")
    parser.add_argument("--repeatmasker-dir", required=True,
                        help="Output directory for RepeatMasker results")
    parser.add_argument("--strain-number", required=True,
                        help="Strain number used to build database")
    parser.add_argument("--jobs", type=int, default=10,
                        help="Number of parallel jobs to run")
    parser.add_argument("--log", default="repeatmasker.log",
                        help="Log file name")

    args = parser.parse_args()

    run_repeatmasker_parallel(
        base_name_chr_anchor_outputs_dir=args.base_dir,
        base_name=args.base_name,
        anchor_set=args.anchor_set,
        repeatmasker_y_prime_dir=args.repeatmasker_dir,
        strain_number=args.strain_number,
        max_parallel_jobs=args.jobs,
        log_file=args.log
    )


if __name__ == "__main__":
    main()

