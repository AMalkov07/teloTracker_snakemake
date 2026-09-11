#!/bin/bash
#$ -q UI,TELOMERE2
#$ -pe smp 16
#$ -j y
#$ -cwd
#$ -N apobec_spacer
# APOBEC motif detection in SPACERS (space_between_anchor: the variable region
# between the anchor and the X-element), A3A vs EV. Same design as the X-element
# run (per-end read mapping + day-0 subtraction), but the spacer reference is
# built by ALIGNMENT (step 0) because the label bed is stale vs the rebuilt ref.
#
# Step 0: place each end's anchor and X-element on the A3A-D0 reference with
# minimap2, take the gap between them as that end's spacer, extract it. Spacers
# can be long (up to ~30 kb); depth (not length) is what matters for consensus
# mutation calling, and per-end single-ref mapping still prevents cross-mapping.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate consensus
export LANG=C LC_ALL=C

REF="results/dorado_AM7575_A3A-3-D0_072326/_pipeline/assembly_7575/assembly_7575_dorado_reference.fasta"
ANCHORS="_pipeline/references/telomerase_shutoff_anchors.fasta"
XELEM="results/dorado_AM7575_A3A-3-D0_072326/_pipeline/pretelomeric_labels/x_element_clustering_7575/x_elements_working.fasta"
OUT="apobec_spacer_analysis"
mkdir -p "$OUT" "$OUT/xref"
samtools faidx "$REF"

# ---- Step 0: build per-end spacer reference by alignment ----
minimap2 -x map-ont --secondary=no -t 16 "$REF" "$ANCHORS" 2>/dev/null > "$OUT/anchors.paf"
minimap2 -x map-ont --secondary=no -t 16 "$REF" "$XELEM"   2>/dev/null > "$OUT/xelem.paf"
python3 compute_spacer_intervals.py "$OUT/anchors.paf" "$OUT/xelem.paf" > "$OUT/spacer_regions.txt"
: > "$OUT/spacer_reference.fasta"
while read -r end region; do
    samtools faidx "$REF" "$region" | awk -v e="$end" 'NR==1{print ">"e"_spacer"; next}{print}' >> "$OUT/spacer_reference.fasta"
done < "$OUT/spacer_regions.txt"
samtools faidx "$OUT/spacer_reference.fasta"
echo "built $(grep -c '^>' $OUT/spacer_reference.fasta) per-end spacer sequences"

# pre-extract each end's spacer as its own single-sequence reference
while read -r name _; do
    samtools faidx "$OUT/spacer_reference.fasta" "$name" > "$OUT/xref/${name}.fa"
    samtools faidx "$OUT/xref/${name}.fa"
done < "$OUT/spacer_reference.fasta.fai"

MPILEUP() { samtools mpileup -B -q 20 -Q 10 -d 100000 \
    --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY -f "$OUT/spacer_reference.fasta" "$1"; }

# map a dataset's per-end anchored reads, each end -> its own spacer, merge
map_perend() {
    local base="$1" outbam="$2" tmp="$OUT/tmp_$(basename $outbam .bam)"
    rm -rf "$tmp"; mkdir -p "$tmp"; local list="$tmp/list.txt"; : > "$list"
    while read -r name _; do
        local end="${name%_spacer}"
        local rf="results/${base}/_pipeline/blast/chr_anchor_reads/${base}_blasted_telomerase_shutoff_anchors_${end}_anchor_reads.fasta"
        [ -s "$rf" ] || continue
        minimap2 -ax map-ont --secondary=no -t 16 "$OUT/xref/${name}.fa" "$rf" 2>/dev/null \
            | samtools sort -@ 4 -o "$tmp/${name}.bam" -
        echo "$tmp/${name}.bam" >> "$list"
    done < "$OUT/spacer_reference.fasta.fai"
    samtools merge -f -b "$list" "$outbam"; samtools index "$outbam"
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
    ndiff=$(python3 apobec_build_consensus.py "$OUT/spacer_reference.fasta" "$OUT/${tag}_d0.pileup" \
                "$OUT/${tag}_day0_consensus.fa" "$OUT/${tag}_day0_diffs.tsv")
    echo "[$tag] day-0 spacer differences vs reference: $ndiff  (see ${tag}_day0_diffs.tsv)"
}

run_arm A3A dorado_AM7575_A3A-3-D0_072326 dorado_AM7575_A3A-3-4c
run_arm EV  dorado_AM7575_EV-3-D0_072226  dorado_AM7575_EV-3-4c

echo
echo "############## APOBEC MOTIF RESULTS (spacers) ##############"
python3 apobec_classify_subtract.py "$OUT/spacer_reference.fasta" "$OUT/A3A_surv.pileup" "$OUT/A3A_d0.pileup" "A3A-4c (APOBEC3A)" | tee "$OUT/A3A_result.txt"
python3 apobec_classify_subtract.py "$OUT/spacer_reference.fasta" "$OUT/EV_surv.pileup"  "$OUT/EV_d0.pileup"  "EV-4c (control)"    | tee "$OUT/EV_result.txt"
echo "DONE"
