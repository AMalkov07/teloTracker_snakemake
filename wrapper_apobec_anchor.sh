#!/bin/bash
#$ -q UI,TELOMERE2
#$ -pe smp 16
#$ -j y
#$ -cwd
#$ -N apobec_anchor
# APOBEC motif detection in ANCHOR regions, A3A vs EV (EV = control).
#
# Design (per arm):
#   1. map that arm's DAY-0 anchored reads to the reference anchors -> build a
#      per-anchor DAY-0 CONSENSUS (majority vote) = strain-adjusted "ancestral"
#      baseline; report substitution differences vs the reference anchors.
#   2. map the arm's SURVIVOR (4c) anchored reads to that day-0 baseline.
#   3. tally survivor mismatches vs baseline and classify the APOBEC motif two
#      ways (RAW per-read + RECURRENCE-supported), reporting both.
# The A3A-4c vs EV-4c comparison is the readout (shared ONT error cancels).
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate consensus
export LANG=C LC_ALL=C

ANCHORS="_pipeline/references/telomerase_shutoff_anchors.fasta"
OUT="apobec_anchor_analysis"
mkdir -p "$OUT"
samtools faidx "$ANCHORS"

MPILEUP() { samtools mpileup -B -q 20 -Q 10 -d 100000 \
    --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY -f "$1" "$2"; }

run_arm() {
    local tag="$1" d0="$2" surv="$3"
    local d0reads="results/${d0}/_pipeline/${d0}_all_chr_anchored_reads.fasta"
    local survreads="results/${surv}/_pipeline/${surv}_all_chr_anchored_reads.fasta"
    echo "===== ARM $tag : baseline=$d0  survivor=$surv ====="

    # 1. day-0 -> reference anchors -> consensus baseline
    minimap2 -ax map-ont --secondary=no -t 16 "$ANCHORS" "$d0reads" 2>/dev/null \
        | samtools sort -@ 8 -o "$OUT/${tag}_d0.bam" -
    samtools index "$OUT/${tag}_d0.bam"
    MPILEUP "$ANCHORS" "$OUT/${tag}_d0.bam" > "$OUT/${tag}_d0.pileup"
    ndiff=$(python3 apobec_build_consensus.py "$ANCHORS" "$OUT/${tag}_d0.pileup" \
                "$OUT/${tag}_baseline.fa" "$OUT/${tag}_day0_diffs.tsv")
    echo "[$tag] day-0 anchor differences vs reference anchors: $ndiff  (see ${tag}_day0_diffs.tsv)"
    samtools faidx "$OUT/${tag}_baseline.fa"

    # 2. survivor -> day-0 baseline
    minimap2 -ax map-ont --secondary=no -t 16 "$OUT/${tag}_baseline.fa" "$survreads" 2>/dev/null \
        | samtools sort -@ 8 -o "$OUT/${tag}_surv.bam" -
    samtools index "$OUT/${tag}_surv.bam"
    MPILEUP "$OUT/${tag}_baseline.fa" "$OUT/${tag}_surv.bam" > "$OUT/${tag}_surv.pileup"
}

run_arm A3A dorado_AM7575_A3A-3-D0_072326 dorado_AM7575_A3A-3-4c
run_arm EV  dorado_AM7575_EV-3-D0_072226  dorado_AM7575_EV-3-4c

echo
echo "############## APOBEC MOTIF RESULTS (anchors) ##############"
python3 apobec_classify.py "$OUT/A3A_baseline.fa" "$OUT/A3A_surv.pileup" "A3A-4c (APOBEC3A)"  | tee "$OUT/A3A_result.txt"
python3 apobec_classify.py "$OUT/EV_baseline.fa"  "$OUT/EV_surv.pileup"  "EV-4c (control)"     | tee "$OUT/EV_result.txt"
echo "DONE"
