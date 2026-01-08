#$ -q TELOMERE2,UI
#$ -pe smp 16
#$ -j y
# -o /dev/null
#$ -cwd

# Base name of the file/run wanting to analyze
#Example to put after submission = dorado_2191_repeat1_old-detect-only
#base_name=$1
base_name=dorado_6991_day0_repeat2_PromethION_no_tag_yes_rejection

# Anchor set you want to blast for chromosome anchors
# Default is new_all_best_anchors
### The new default is test_anchors
#anchor_set=test_anchors
#anchor_set=tlc1_mph1_anchor_sequences_from_day1
anchor_set=test_anchors


#######################################################################################
#######################################################################################

# Other paths that will be used - Don't touch
telomere_analysis_dir=./
base_name_dir=$base_name/
strain_number=${base_name%%_re*}
strain_number=${strain_number#dorado_}
base_name_chr_anchor_outputs_dir=$base_name/chr_anchor_included_individual_files/
blast_y_primes_output_dir=$base_name/y_prime_blast_results_from_chr_anchored_reads/
base_name_repeat_split_dir=$base_name/split_to_repeat_outputs/
base_name_porechop_dir=$base_name/porechop_trim/

adapter_seq_file=offical_nanopore_adapter_seq+trunc.txt

porechop_abi_log_file=$base_name_dir/$base_name\_all_chrs_split_to_telomere_porechopped.log

bam_dir=../samples_dorado_basecalled2/

#######################################################################################
#######################################################################################

: "
Order of steps/commands/scripts run:

"

#######################################################################################
#######################################################################################

echo "Starting $base_name"

#conda activate consensus

# Step 0
# Prep files
echo "Step 0. Prep files by moving/creating files from nanopore_analysis to telomere_analysis"

#mkdir -p $base_name_dir

#samtools fastq -@ 56 -T '*' $bam_dir/$base_name.bam > $base_name_dir/$base_name.fastq
#seqtk seq -A $base_name_dir/$base_name.fastq > $base_name_dir/$base_name.fasta

# Step 1
# Blast for the Reads with Anchors

echo "Step 1. Blast for the Reads with Anchors"

query=$base_name_dir/$base_name
database=$anchor_set

echo -e \
"read_id\ttotal_read_length\tread_bp_used_for_match\tmatch_start_on_read\tmatch_end_on_read\tanchor_name\ttotal_anchor_length\tmatch_start_on_anchor\tmatch_end_on_anchor\tpident\tbitscore\tevalue" \
> $query\_blasted_$database.tsv

#blastn -query $query.fasta -db $database.fasta -perc_identity 85 -min_raw_gapped_score 3000 -num_threads 56 \
#-outfmt "6 qseqid qlen length qstart qend sseqid slen sstart send pident bitscore evalue" >> $query\_blasted_$database.tsv


# New testing blast command

#mv $base_name_dir/$base_name\_blasted_1003_anchors.tsv $base_name_dir/$base_name\_blasted_1003_anchors-old.tsv
#mv $base_name_dir/all_matches_$base_name\_blasted_1003_anchors.tsv $base_name_dir/all_matches_$base_name\_blasted_1003_anchors-old.tsv
#mv $base_name_dir/top_matches_$base_name\_blasted_1003_anchors.fasta $base_name_dir/top_matches_$base_name\_blasted_1003_anchors-old.fasta

blastn -query $query.fasta -db $database.fasta -task megablast -perc_identity 85 -min_raw_gapped_score 3000 -num_threads 56 \
-outfmt "6 qseqid qlen length qstart qend sseqid slen sstart send pident bitscore evalue" >> $query\_blasted_$database.tsv


#######################################################################################

# Step 2
# Filter the Blast Results for the Reads with Anchors

echo "Step 2. Filter the Blast Results for the Reads with chromosome Anchors"

python filter_for_reads_with_anchors.py $base_name $anchor_set


#######################################################################################

# Step 3
# Once Finished Run the Command Below

echo "Step 3. Grabs/makes a new file with the chromosome anchored reads for each chromosome with labels of AC/TG and the chromosome end name"

#python split_and_label_all_reads.py $base_name $anchor_set

python split_and_label_all_reads_include_anchor.py $base_name $anchor_set

#######################################################################################

mkdir -p $base_name_porechop_dir


############# cat $base_name_repeat_split_dir/chr*_split_to_telomere_all_reads_best_pre_trim.fasta > $all_chrs_split_to_telomere_base\_pre_trim.fasta
cat $base_name_chr_anchor_outputs_dir/$base_name\_blasted_$anchor_set\_chr*_anchor_reads.fasta > $base_name/$base_name\_all_chromosome_anchored_reads_pre_trim.fasta

# Check for adapters and tags

echo "Starting porechop_abi"

porechop_abi -t 56 --format fastq -ddb -v 3 --no_split -cap $adapter_seq_file -i $base_name/$base_name\_all_chromosome_anchored_reads_pre_trim.fasta \
             -o $base_name_porechop_dir/all_chromosome_anchored_read_trimmed.fastq > $porechop_abi_log_file

sed -n '1~4s/^@/>/p;2~4p' $base_name_porechop_dir/all_chromosome_anchored_read_trimmed.fastq > $base_name_porechop_dir/all_chromosome_anchored_read_trimmed.fasta

python check_for_adapters.py $base_name $porechop_abi_log_file

# Compare tags vs adapters

#dorado summary $bam_dir/$base_name.bam > $base_name_dir/$base_name\_sequencing_summary.tsv

python make_summary_file_for_comparison.py $base_name

python compare_adapter_callers_dorado.py $base_name $anchor_set

# Fine trim telomere for adapter and tag sequences

python fine_telomere_trimming.py $base_name $anchor_set


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
        database="probe"
        blasting_file_name=$blast_y_primes_output_dir_probe/$base_name\_$chr_num$side\_anchor_reads_blasted_$database.tsv

        echo -e \
        "read_id\ttotal_read_length\tread_bp_used_for_match\tmatch_start_on_read\tmatch_end_on_read\tanchor_name\ttotal_anchor_length\tmatch_start_on_anchor\tmatch_end_on_anchor\tpident\tbitscore\tevalue" \
        > $blasting_file_name

        blastn -query $query.fasta -db $database.fasta -perc_identity 90 -num_threads 56 \
        -outfmt "6 qseqid qlen length qstart qend sseqid slen sstart send pident bitscore evalue" >> $blasting_file_name

    done
done


echo "Done with Step 7."

#######################################################################################

# Step 8.
# Running Y Prime Analysis

echo "Step 8. Run Y Prime Analysis"

#python filter_for_y_primes_in_reads.py $base_name

# Read through the y prime probe blast results, make figures for Y primes, and the outputs directory final tsv
python y_prime_analysis.py $base_name $anchor_set

echo "Done with Step 8."

#######################################################################################


# Step 9

echo "Step 9. Making Telomere Graphs"

cd outputs/

python combined_plot_telomeres.py $base_name\_post_y_prime_probe.tsv

cd ../

echo "Done with Step 9."

#######################################################################################

