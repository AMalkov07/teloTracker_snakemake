#!/bin/bash
#$ -q UI,TELOMERE2
#$ -pe smp 16
#$ -j y
#$ -cwd
#$ -N apobec_xelem
# APOBEC motif detection in X-ELEMENTS (core+variable as one feature), A3A vs EV.
#
# X-elements are a conserved family across chr ends, so reads could cross-map to
# the wrong end's X. We AVOID that structurally: every read is already anchored
# to a specific end, so we map each end's anchored reads ONLY to that end's own X
# sequence (single-reference minimap2) -> a read can never land on another end's
# X. Per-arm day-0 subtraction is kept as verification/strain-correction.
#
# Design per arm: build day-0 pileup (ancestral) + survivor pileup, both per-end;
# call survivor mismatches vs the X reference and keep as ACQUIRED only those
# absent in day-0 (removes strain diffs + any standing/paralog variation).
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate consensus
export LANG=C LC_ALL=C

# 32 per-end X sequences (core+variable), named chr<N><L/R>_x_ends
XREF_SRC="results/dorado_AM7575_A3A-3-D0_072326/_pipeline/pretelomeric_labels/x_element_clustering_7575/x_elements_working.fasta"
OUT="apobec_xelement_analysis"
mkdir -p "$OUT" "$OUT/xref"
cp "$XREF_SRC" "$OUT/x_reference.fasta"
samtools faidx "$OUT/x_reference.fasta"

# pre-extract each end's X sequence as its own single-sequence reference
while read -r xname _; do
    samtools faidx "$OUT/x_reference.fasta" "$xname" > "$OUT/xref/${xname}.fa"
    samtools faidx "$OUT/xref/${xname}.fa"
done < "$OUT/x_reference.fasta.fai"

MPILEUP() { samtools mpileup -B -q 20 -Q 10 -d 100000 \
    --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY -f "$OUT/x_reference.fasta" "$1"; }

# map a dataset's per-end anchored reads, each end -> its own X, merge
map_perend() {
    local base="$1" outbam="$2" tmp="$OUT/tmp_$(basename $outbam .bam)"
    rm -rf "$tmp"; mkdir -p "$tmp"; local list="$tmp/list.txt"; : > "$list"
    while read -r xname _; do
        local end="${xname%_x_ends}"
        local rf="results/${base}/_pipeline/blast/chr_anchor_reads/${base}_blasted_telomerase_shutoff_anchors_${end}_anchor_reads.fasta"
        [ -s "$rf" ] || continue
        minimap2 -ax map-ont --secondary=no -t 16 "$OUT/xref/${xname}.fa" "$rf" 2>/dev/null \
            | samtools sort -@ 4 -o "$tmp/${xname}.bam" -
        echo "$tmp/${xname}.bam" >> "$list"
    done < "$OUT/x_reference.fasta.fai"
    samtools merge -f -b "$list" "$outbam"
    samtools index "$outbam"
    echo "  [$base] $(samtools view -c -F 0x904 $outbam) primary alignments across $(wc -l < $list) ends"
    rm -rf "$tmp"
}

run_arm() {
    local tag="$1" d0="$2" surv="$3"
    echo "===== ARM $tag : baseline(day0)=$d0  survivor=$surv ====="
    map_perend "$d0"   "$OUT/${tag}_d0.bam"
    map_perend "$surv" "$OUT/${tag}_surv.bam"
    MPILEUP "$OUT/${tag}_d0.bam"   > "$OUT/${tag}_d0.pileup"
    MPILEUP "$OUT/${tag}_surv.bam" > "$OUT/${tag}_surv.pileup"
    ndiff=$(python3 apobec_build_consensus.py "$OUT/x_reference.fasta" "$OUT/${tag}_d0.pileup" \
                "$OUT/${tag}_day0_consensus.fa" "$OUT/${tag}_day0_diffs.tsv")
    echo "[$tag] day-0 X differences vs X reference: $ndiff  (see ${tag}_day0_diffs.tsv)"
}

run_arm A3A dorado_AM7575_A3A-3-D0_072326 dorado_AM7575_A3A-3-4c
run_arm EV  dorado_AM7575_EV-3-D0_072226  dorado_AM7575_EV-3-4c

echo
echo "############## APOBEC MOTIF RESULTS (X-elements) ##############"
python3 apobec_classify_subtract.py "$OUT/x_reference.fasta" "$OUT/A3A_surv.pileup" "$OUT/A3A_d0.pileup" "A3A-4c (APOBEC3A)" | tee "$OUT/A3A_result.txt"
python3 apobec_classify_subtract.py "$OUT/x_reference.fasta" "$OUT/EV_surv.pileup"  "$OUT/EV_d0.pileup"  "EV-4c (control)"    | tee "$OUT/EV_result.txt"
echo "DONE"
