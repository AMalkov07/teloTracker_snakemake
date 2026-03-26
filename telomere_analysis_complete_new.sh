#$ -q TELOMERE2,UI
#$ -pe smp 108
#$ -j y
# -o /dev/null
#$ -cwd

# Base name of the file/run wanting to analyze
#Example to put after submission = dorado_2191_repeat1_old-detect-only
base_name=$1

# Anchor set you want to blast for chromosome anchors
# Default is new_all_best_anchors
### The new default is test_anchors
anchor_set=test_anchors


#######################################################################################
#######################################################################################

# Other paths that will be used - Don't touch
telomere_analysis_dir=./
base_name_dir=results/$base_name/
strain_number=${base_name%%_day*}
strain_number=${strain_number#dorado_}
base_name_chr_anchor_outputs_dir=results/$base_name/chr_anchor_included_individual_files/
blast_y_primes_output_dir=results/$base_name/y_prime_blast_results_from_chr_anchored_reads/
base_name_repeat_split_dir=results/$base_name/split_to_repeat_outputs/
base_name_porechop_dir=results/$base_name/porechop_trim/

adapter_seq_file=references/offical_nanopore_adapter_seq+trunc.txt

porechop_abi_log_file=$base_name_dir/$base_name\_all_chrs_split_to_telomere_porechopped.log

bam_dir=../samples_dorado_basecalled/

#######################################################################################
#######################################################################################

echo "Starting $base_name"

conda activate consensus

# Step 0
# Prep files
echo "Step 0. Prep files by moving/creating files from nanopore_analysis to telomere_analysis"

mkdir -p $base_name_dir

# Original bam->fastq conversion (commented out - using pre-existing fastq instead)
#samtools fastq -@ 108 -T '*' $bam_dir/$base_name.bam > $base_name_dir/$base_name.fastq

# Copy fastq from samples_dorado_basecalled2 directory
cp samples_dorado_basecalled2/$base_name.fastq $base_name_dir/$base_name.fastq

seqtk seq -A $base_name_dir/$base_name.fastq > $base_name_dir/$base_name.fasta

# Step 1
# Blast for the Reads with Anchors

echo "Step 1. Blast for the Reads with Anchors"

query=$base_name_dir/$base_name
database=references/$anchor_set
blast_output=$base_name_dir/$base_name\_blasted_$anchor_set.tsv

echo -e \
"read_id\ttotal_read_length\tread_bp_used_for_match\tmatch_start_on_read\tmatch_end_on_read\tanchor_name\ttotal_anchor_length\tmatch_start_on_anchor\tmatch_end_on_anchor\tpident\tbitscore\tevalue" \
> $blast_output

blastn -query $query.fasta -db $database.fasta -task dc-megablast -perc_identity 85 -min_raw_gapped_score 3000 -num_threads 108 \
-outfmt "6 qseqid qlen length qstart qend sseqid slen sstart send pident bitscore evalue" >> $blast_output


#######################################################################################

# Step 2
# Filter the Blast Results for the Reads with Anchors

echo "Step 2. Filter the Blast Results for the Reads with chromosome Anchors"

python scripts/filter_for_reads_with_anchors.py $base_name $anchor_set


#######################################################################################

# Step 3
# Once Finished Run the Command Below

echo "Step 3. Grabs/makes a new file with the chromosome anchored reads for each chromosome with labels of AC/TG and the chromosome end name"

#python scripts/split_and_label_all_reads.py $base_name $anchor_set

python scripts/split_and_label_all_reads_include_anchor.py $base_name $anchor_set

#######################################################################################

mkdir -p $base_name_porechop_dir


############# cat $base_name_repeat_split_dir/chr*_split_to_telomere_all_reads_best_pre_trim.fasta > $all_chrs_split_to_telomere_base\_pre_trim.fasta
cat $base_name_chr_anchor_outputs_dir/$base_name\_blasted_$anchor_set\_chr*_anchor_reads.fasta > $base_name_dir/$base_name\_all_chromosome_anchored_reads_pre_trim.fasta

# Check for adapters and tags

echo "Starting porechop_abi"

porechop_abi -t 108 --format fastq -ddb -v 3 --no_split -cap $adapter_seq_file -i $base_name_dir/$base_name\_all_chromosome_anchored_reads_pre_trim.fasta \
             -o $base_name_porechop_dir/all_chromosome_anchored_read_trimmed.fastq > $porechop_abi_log_file

sed -n '1~4s/^@/>/p;2~4p' $base_name_porechop_dir/all_chromosome_anchored_read_trimmed.fastq > $base_name_porechop_dir/all_chromosome_anchored_read_trimmed.fasta

python scripts/check_for_adapters.py $base_name $porechop_abi_log_file

# Compare tags vs adapters

#dorado summary $bam_dir/$base_name.bam > $base_name_dir/$base_name\_sequencing_summary.tsv

python scripts/make_summary_file_for_comparison.py $base_name

python scripts/compare_adapter_callers_dorado.py $base_name $anchor_set

# Fine trim telomere for adapter and tag sequences

python scripts/fine_telomere_trimming.py $base_name $anchor_set


#######################################################################################
#######################################################################################
# Step 7
# Blast Anchored Reads for Y' sequences/probes

echo "Step 7. Blast Anchored Reads for pre-telomere repeats"


# Blasting the Y' probe sequence

echo "Blasting Y' probe sequence"

blast_y_primes_output_dir_probe=$blast_y_primes_output_dir/probes/

mkdir -p $blast_y_primes_output_dir_probe

for chr_num in {1..16}
do
    for side in {"L","R"}
    do
        echo "Blasting $chr_num$side"

        query=$base_name_chr_anchor_outputs_dir/$base_name\_blasted_$anchor_set\_chr$chr_num$side\_anchor_reads
        database="references/probe"
        blasting_file_name=$blast_y_primes_output_dir_probe/$base_name\_$chr_num$side\_anchor_reads_blasted_probe.tsv

        echo -e \
        "read_id\ttotal_read_length\tread_bp_used_for_match\tmatch_start_on_read\tmatch_end_on_read\tanchor_name\ttotal_anchor_length\tmatch_start_on_anchor\tmatch_end_on_anchor\tpident\tbitscore\tevalue" \
        > $blasting_file_name

        blastn -query $query.fasta -db $database.fasta -perc_identity 90 -num_threads 108 \
        -outfmt "6 qseqid qlen length qstart qend sseqid slen sstart send pident bitscore evalue" >> $blasting_file_name

    done
done


echo "Done with Step 7."

#######################################################################################

# Step 8.
# Running Y Prime Analysis

echo "Step 8. Run Y Prime Analysis"

#python scripts/filter_for_y_primes_in_reads.py $base_name

# Read through the y prime probe blast results, make figures for Y primes, and the outputs directory final tsv
python scripts/y_prime_analysis.py $base_name $anchor_set

echo "Done with Step 8."

#######################################################################################


# Step 9

echo "Step 9. Making Telomere Graphs"

mkdir -p results/outputs/

python scripts/combined_plot_telomeres.py $base_name results/outputs/$base_name\_post_y_prime_probe.tsv results/outputs/

echo "Done with Step 9."

#######################################################################################
#######################################################################################

# ============================================================================
# NEW STEPS 10–13: Recombination analysis against day0 reference
# (replaces the old RepeatMasker + pairing approach)
# ============================================================================

STRAIN_ID="${strain_number}"
THREADS=108
THREADS_PER_JOB=4

# --- Paths derived from the day0 run for this strain ---
DAY0_BASE_NAME="dorado_${STRAIN_ID}_day0_PromethION_no_tag_yes_rejection"
DAY0_REF="results/${DAY0_BASE_NAME}/assembly_${STRAIN_ID}/assembly_${STRAIN_ID}_dorado_reference.fasta"
DAY0_BED="results/${DAY0_BASE_NAME}/pretelomeric_labels/pretelomeric_regions_${STRAIN_ID}_simp.bed"
Y_PRIME_LIB="references/extracted_yprimes_${STRAIN_ID}.fasta"
SPACER_LIB_DIR="references/pairings_for_spacers/${STRAIN_ID}_pairings/"
X_ELEMENT_LIB_DIR="references/pairings_for_x_element_ends/${STRAIN_ID}_pairings/"
CHR_ANCHOR_DIR="results/${base_name}/chr_anchor_included_individual_files"
RECOMB_DIR="results/${base_name}/recombination"

# --- Validate day0 reference files ---
for f in "${DAY0_REF}" "${DAY0_BED}" "${Y_PRIME_LIB}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: Required day0 file missing: ${f}"
        echo "Please run create_ref.sh and label_regions.sh for strain ${STRAIN_ID} (day0) first."
        exit 1
    fi
done

mkdir -p "${RECOMB_DIR}"

MAX_PARALLEL=$(( THREADS / THREADS_PER_JOB ))

# ============================================================================
# Step 10: minimap2 alignment (one job per chr end, parallelized)
# ============================================================================

echo "Step 10. Running minimap2 alignment for all 32 chromosome ends"

for chr_num in {1..16}; do
    for side in L R; do
        chr_end="chr${chr_num}${side}"
        reads="${CHR_ANCHOR_DIR}/${base_name}_blasted_${anchor_set}_${chr_end}_anchor_reads.fasta"

        [ ! -f "${reads}" ] && continue

        # Simple semaphore: wait if already at max parallel jobs
        while [ "$(jobs -r | wc -l)" -ge "${MAX_PARALLEL}" ]; do
            sleep 1
        done

        python scripts/run_alignment.py \
            --reads-fasta "${reads}" \
            --day0-ref    "${DAY0_REF}" \
            --day0-bed    "${DAY0_BED}" \
            --chr-end     "${chr_end}" \
            --output-tsv  "${RECOMB_DIR}/${base_name}_${chr_end}_breakpoints.tsv" \
            --output-bam  "${RECOMB_DIR}/${base_name}_${chr_end}.bam" \
            --threads     "${THREADS_PER_JOB}" &
    done
done
wait

echo "Done with Step 10."

#######################################################################################

# ============================================================================
# Step 11: Element-specific recombination analysis (parallelized per chr end)
# All three element scripts run sequentially within each chr-end job,
# but chr ends are processed in parallel.
# ============================================================================

echo "Step 11. Running element-specific recombination analysis"

for chr_num in {1..16}; do
    for side in L R; do
        chr_end="chr${chr_num}${side}"
        bp_tsv="${RECOMB_DIR}/${base_name}_${chr_end}_breakpoints.tsv"
        div_fasta="${RECOMB_DIR}/${base_name}_${chr_end}_diverged_tails.fasta"
        reads="${CHR_ANCHOR_DIR}/${base_name}_blasted_${anchor_set}_${chr_end}_anchor_reads.fasta"

        [ ! -f "${bp_tsv}" ] && continue

        while [ "$(jobs -r | wc -l)" -ge "${MAX_PARALLEL}" ]; do
            sleep 1
        done

        (
            # 11a: Y prime analysis (runs on full reads)
            python scripts/analyze_y_prime_recombination.py \
                --breakpoints-tsv "${bp_tsv}" \
                --reads-fasta     "${reads}" \
                --diverged-fasta  "${div_fasta}" \
                --day0-bed        "${DAY0_BED}" \
                --y-prime-lib     "${Y_PRIME_LIB}" \
                --chr-end         "${chr_end}" \
                --strain          "${STRAIN_ID}" \
                --output-tsv      "${RECOMB_DIR}/${base_name}_${chr_end}_y_prime_recomb.tsv" \
                --threads         "${THREADS_PER_JOB}"

            # 11b: X prime analysis (runs on diverged tails)
            python scripts/analyze_x_prime_recombination.py \
                --breakpoints-tsv   "${bp_tsv}" \
                --diverged-fasta    "${div_fasta}" \
                --x-element-lib-dir "${X_ELEMENT_LIB_DIR}" \
                --chr-end           "${chr_end}" \
                --strain            "${STRAIN_ID}" \
                --output-tsv        "${RECOMB_DIR}/${base_name}_${chr_end}_x_prime_recomb.tsv" \
                --threads           "${THREADS_PER_JOB}"

            # 11c: Spacer analysis (chunk + tail hybrid)
            python scripts/analyze_spacer_recombination.py \
                --breakpoints-tsv "${bp_tsv}" \
                --diverged-fasta  "${div_fasta}" \
                --reads-fasta     "${reads}" \
                --day0-ref        "${DAY0_REF}" \
                --spacer-lib-dir  "${SPACER_LIB_DIR}" \
                --chr-end         "${chr_end}" \
                --strain          "${STRAIN_ID}" \
                --output-tsv      "${RECOMB_DIR}/${base_name}_${chr_end}_spacer_recomb.tsv" \
                --threads         "${THREADS_PER_JOB}"
        ) &
    done
done
wait

echo "Done with Step 11."

#######################################################################################

# ============================================================================
# Step 12: Aggregate per-chr-end results into unified per-read TSVs
# ============================================================================

echo "Step 12. Aggregating per-read recombination results"

for chr_num in {1..16}; do
    for side in L R; do
        chr_end="chr${chr_num}${side}"
        bp_tsv="${RECOMB_DIR}/${base_name}_${chr_end}_breakpoints.tsv"

        [ ! -f "${bp_tsv}" ] && continue

        python scripts/aggregate_recombination.py \
            --breakpoints-tsv  "${bp_tsv}" \
            --y-prime-tsv      "${RECOMB_DIR}/${base_name}_${chr_end}_y_prime_recomb.tsv" \
            --x-prime-tsv      "${RECOMB_DIR}/${base_name}_${chr_end}_x_prime_recomb.tsv" \
            --spacer-tsv       "${RECOMB_DIR}/${base_name}_${chr_end}_spacer_recomb.tsv" \
            --chr-end          "${chr_end}" \
            --output-reads     "${RECOMB_DIR}/${base_name}_${chr_end}_recombination_reads.tsv"
    done
done

echo "Done with Step 12."

#######################################################################################

# ============================================================================
# Step 13: Write cross-chr-end summary TSV
# ============================================================================

echo "Step 13. Writing recombination summary"

python scripts/aggregate_recombination.py --summarize \
    --recombination-dir "${RECOMB_DIR}" \
    --base-name         "${base_name}" \
    --output-summary    "${RECOMB_DIR}/${base_name}_recombination_summary.tsv"

echo "Done with Step 13."

echo ""
echo "========================================================================"
echo "Recombination analysis complete for ${base_name}"
echo "Results in: ${RECOMB_DIR}"
echo "  - 32x *_breakpoints.tsv (minimap2 alignment results)"
echo "  - 32x *_diverged_tails.fasta (post-breakpoint sequences)"
echo "  - 32x *.bam (aligned reads)"
echo "  - 32x *_y_prime_recomb.tsv"
echo "  - 32x *_x_prime_recomb.tsv"
echo "  - 32x *_spacer_recomb.tsv"
echo "  - 32x *_recombination_reads.tsv (unified per-read)"
echo "  - 1x  *_recombination_summary.tsv (per-chr-end summary)"
echo "========================================================================"
