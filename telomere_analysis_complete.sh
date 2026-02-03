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

bam_dir=samples_dorado_basecalled/

#######################################################################################
#######################################################################################

: "
Order of steps/commands/scripts run:

"

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
cp $bam_dir/$base_name.fastq $base_name_dir/$base_name.fastq

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

: "


# Step 10

# Blasting the Y' sequences - without redundant Y' sequences

echo "Step 10. Finding Y' sequences by RepeatMasker"

# mamba activate repeatmasker  # Using conda activate consensus at top

repeatmasker_y_prime_dir=$base_name_dir/read_repeatmasker_results/

mkdir -p $repeatmasker_y_prime_dir

for chr_num in {1..17}
do
    for side in {"L","R"}
    do

        echo "Searching $chr_num$side reads"

        query=$base_name_chr_anchor_outputs_dir/$base_name\_blasted_$anchor_set\_chr$chr_num$side\_anchor_reads
        database="references/extracted_yprimes_${strain_number}"
        repeatmasker_file_name=$repeatmasker_y_prime_dir/$base_name\_chr$chr_num$side\_repeatmasker_results.ssv

        RepeatMasker $query.fasta -lib $database.fasta -s -pa 14 --cutoff 1000 -no_is -norna -gff -dir $repeatmasker_y_prime_dir

        echo -e \
        "SW_score divergence_percent deletion_percent insertion_percent read_id match_start_on_read match_end_on_read leftover_on_read strand y_prime_id y_prime_group match_start_on_y_prime match_end_on_y_prime leftover_on_y_prime match_id sub_match" \
        > $repeatmasker_file_name

        filtered_repeatmasker_data_file=$repeatmasker_y_prime_dir/chr$chr_num$side\_filtered.out

        grep "Y_Prime" $repeatmasker_y_prime_dir/$base_name\_blasted_$anchor_set\_chr$chr_num$side\_anchor_reads.fasta.out \
        > $filtered_repeatmasker_data_file

        cat $filtered_repeatmasker_data_file >> $repeatmasker_file_name

    done
done

# mamba activate base  # Using conda activate consensus at top

python scripts/make_y_prime_repeatmasker_tsv.py $base_name

python scripts/get_summary_stats_for_y_prime_repeatmasker.py $base_name $anchor_set

python scripts/make_pairings_from_y_primes.py $base_name $anchor_set


echo "Done with Step 10."

#######################################################################################


# Step 11. Find any x element switches

echo "Step 11. Find any x element end switches"


# RepeatMasker for x element for y prime pairings

echo "Finding x element end sequences by RepeatMasker for y prime pairings"


# mamba activate repeatmasker  # Using conda activate consensus at top

input_reads_dir=results/$base_name/paired_by_y_prime_reads/

# Get all file names in the directory with the ".fasta" suffix
fasta_files=("$input_reads_dir"/*.fasta)

# List to store file names without the ".fasta" suffix
parings_to_run=()

# Iterate over the file names
for file in "${fasta_files[@]}"; do
    # Remove the ".fasta" suffix and append to the list
    parings_to_run+=("$(basename "$file" .fasta)")
done

repeatmasker_dir=$base_name_dir/paired_x_element_ends_repeatmasker_results/

mkdir -p $repeatmasker_dir

for pairing in "${parings_to_run[@]}"; do

    echo "Starting $pairing reads"

    query=$input_reads_dir/$pairing
    database=references/pairings_for_x_element_ends/$strain_number\_pairings/$strain_number\_paired_$pairing

    echo "Running RepeatMasker on reads $query.fasta with reference $database.fasta"

    RepeatMasker $query.fasta -lib $database.fasta -s -pa 14 --cutoff 500 -no_is -norna -gff -dir $repeatmasker_dir


    filtered_repeatmasker_data_file=$repeatmasker_dir/$pairing\_x_element_ends_filtered.out
    grep "x_ends" $repeatmasker_dir/$pairing.fasta.out \
    > $filtered_repeatmasker_data_file

    repeatmasker_file_name=$repeatmasker_dir/$base_name\_$pairing\_x_element_ends_repeatmasker_results.ssv
    echo -e \
    "SW_score divergence_percent deletion_percent insertion_percent read_id match_start_on_read match_end_on_read leftover_on_read strand x_element_ends section_number match_start_on_chr_end_section match_end_on_chr_end_section leftover_on_chr_end_section match_id sub_match" \
    > $repeatmasker_file_name
    cat $filtered_repeatmasker_data_file >> $repeatmasker_file_name

done

# mamba activate base  # Using conda activate consensus at top

python scripts/make_x_element_ends_pairs_repeatmasker_tsv.py $base_name


echo "Done with Step 11."

#######################################################################################

# Step 12
# Determine when an end switches from one chromosome arm to another
# Finding 250bp tracts of spacer sequences by RepeatMasker for spacer sequences for y prime pairings

echo "Finding 250bp tracts of spacer sequences by RepeatMasker for y prime pairings"

# mamba activate repeatmasker  # Using conda activate consensus at top

input_reads_dir=results/$base_name/paired_by_y_prime_reads/

# Get all file names in the directory with the ".fasta" suffix
fasta_files=("$input_reads_dir"/*.fasta)

# List to store file names without the ".fasta" suffix
parings_to_run=()

# Iterate over the file names
for file in "${fasta_files[@]}"; do
    # Remove the ".fasta" suffix and append to the list
    parings_to_run+=("$(basename "$file" .fasta)")
done

repeatmasker_dir=$base_name_dir/paired_spacer_repeatmasker_results/

mkdir -p $repeatmasker_dir

for pairing in "${parings_to_run[@]}"; do

    echo "Starting $pairing reads"

    query=$input_reads_dir/$pairing
    database=references/pairings_for_spacers/$strain_number\_pairings/$strain_number\_paired_$pairing

    echo "Running RepeatMasker on reads $query.fasta with reference $database.fasta"

    RepeatMasker $query.fasta -lib $database.fasta -s -pa 14 --cutoff 500 -no_is -norna -gff -dir $repeatmasker_dir


    filtered_repeatmasker_data_file=$repeatmasker_dir/$pairing\_spacer_filtered.out
    grep "_from_repeat_to_plus_50kb" $repeatmasker_dir/$pairing.fasta.out \
    > $filtered_repeatmasker_data_file

    repeatmasker_file_name=$repeatmasker_dir/$base_name\_$pairing\_spacer_repeatmasker_results.ssv
    echo -e \
    "SW_score divergence_percent deletion_percent insertion_percent read_id match_start_on_read match_end_on_read leftover_on_read strand chr_end_tract section_number match_start_on_chr_end_section match_end_on_chr_end_section leftover_on_chr_end_section match_id sub_match" \
    > $repeatmasker_file_name
    cat $filtered_repeatmasker_data_file >> $repeatmasker_file_name

done


# mamba activate base  # Using conda activate consensus at top

python scripts/make_spacer_pairs_repeatmasker_tsv.py $base_name


echo "Done with Step 12."

#######################################################################################

# Step 13

"



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

: "
Order of steps/commands/scripts run:




python scripts/repeatmasker_runner.py \
  --base-dir $base_name_chr_anchor_outputs_dir \
  --base-name $base_name \
  --anchor-set $anchor_set \
  --repeatmasker-dir $repeatmasker_y_prime_dir \
  --strain-number $strain_number \
  --jobs 10





"

#######################################################################################
#######################################################################################

echo "Starting $base_name"

# mamba activate telotracker  # Using conda activate consensus at top

# Rerun this part

# Read through the y prime probe blast results, make figures for Y primes, and the outputs directory final tsv
python scripts/y_prime_analysis.py $base_name $anchor_set

python scripts/fine_telomere_trimming.py $base_name $anchor_set



# Step 10

# Blasting the Y' sequences - without redundant Y' sequences

echo "Step 10. Finding Y' sequences by RepeatMasker"

repeatmasker_y_prime_dir=$base_name_dir/read_repeatmasker_results/

mkdir -p $repeatmasker_y_prime_dir

for chr_num in {1..16}
do
    for side in {"L","R"}
    do

        echo "Searching $chr_num$side reads"

        query=$base_name_chr_anchor_outputs_dir/$base_name\_blasted_$anchor_set\_chr$chr_num$side\_anchor_reads
        database="references/extracted_yprimes_${strain_number}"
        repeatmasker_file_name=$repeatmasker_y_prime_dir/$base_name\_chr$chr_num$side\_repeatmasker_results.ssv

        RepeatMasker $query.fasta -lib $database.fasta -s -pa 12 --cutoff 1000 -no_is -norna -gff -dir $repeatmasker_y_prime_dir

        echo -e \
        "SW_score divergence_percent deletion_percent insertion_percent read_id match_start_on_read match_end_on_read leftover_on_read strand y_prime_id y_prime_group match_start_on_y_prime match_end_on_y_prime leftover_on_y_prime match_id sub_match" \
        > $repeatmasker_file_name

        filtered_repeatmasker_data_file=$repeatmasker_y_prime_dir/chr$chr_num$side\_filtered.out

        grep "Y_Prime" $repeatmasker_y_prime_dir/$base_name\_blasted_$anchor_set\_chr$chr_num$side\_anchor_reads.fasta.out \
        > $filtered_repeatmasker_data_file

        cat $filtered_repeatmasker_data_file >> $repeatmasker_file_name

    done
done

python scripts/make_y_prime_repeatmasker_tsv.py $base_name

python scripts/get_summary_stats_for_y_prime_repeatmasker.py $base_name $anchor_set $strain_number

python scripts/get_stats_of_recombination.py $base_name

##################### 
# Only for 0 --> 1 Y's
#python scripts/make_pairings_from_y_primes.py $base_name $anchor_set $strain_number

# For all change in first Y's
python scripts/make_pairings_from_y_primes_all_ends.py $base_name 
#####################

echo "Done with Step 10."

#######################################################################################


# Step 11. Find any x element switches

echo "Step 11. Find any x element end switches"


# RepeatMasker for x element for y prime pairings

echo "Finding x element end sequences by RepeatMasker for y prime pairings"


input_reads_dir=results/$base_name/paired_by_y_prime_reads/

# Get all file names in the directory with the ".fasta" suffix
fasta_files=("$input_reads_dir"/*.fasta)

# List to store file names without the ".fasta" suffix
parings_to_run=()

# Iterate over the file names
for file in "${fasta_files[@]}"; do
    # Remove the ".fasta" suffix and append to the list
    parings_to_run+=("$(basename "$file" .fasta)")
done

repeatmasker_dir=$base_name_dir/paired_x_element_ends_repeatmasker_results/

mkdir -p $repeatmasker_dir

for pairing in "${parings_to_run[@]}"; do

    echo "Starting $pairing reads"

    query=$input_reads_dir/$pairing
    database=references/pairings_for_x_element_ends/$strain_number\_pairings/$strain_number\_paired_$pairing

    echo "Running RepeatMasker on reads $query.fasta with reference $database.fasta"

    RepeatMasker $query.fasta -lib $database.fasta -s -pa 12 --cutoff 500 -no_is -norna -gff -dir $repeatmasker_dir


    filtered_repeatmasker_data_file=$repeatmasker_dir/$pairing\_x_element_ends_filtered.out
    grep "x_ends" $repeatmasker_dir/$pairing.fasta.out \
    > $filtered_repeatmasker_data_file

    repeatmasker_file_name=$repeatmasker_dir/$base_name\_$pairing\_x_element_ends_repeatmasker_results.ssv
    echo -e \
    "SW_score divergence_percent deletion_percent insertion_percent read_id match_start_on_read match_end_on_read leftover_on_read strand x_element_ends section_number match_start_on_chr_end_section match_end_on_chr_end_section leftover_on_chr_end_section match_id sub_match" \
    > $repeatmasker_file_name
    cat $filtered_repeatmasker_data_file >> $repeatmasker_file_name

done


python scripts/make_x_element_ends_pairs_repeatmasker_tsv.py $base_name $strain_number


echo "Done with Step 11."

#######################################################################################

# Step 12
# Determine when an end switches from one chromosome arm to another
# Finding 250bp tracts of spacer sequences by RepeatMasker for spacer sequences for y prime pairings

echo "Finding 250bp tracts of spacer sequences by RepeatMasker for y prime pairings"


input_reads_dir=results/$base_name/paired_by_y_prime_reads/

# Get all file names in the directory with the ".fasta" suffix
fasta_files=("$input_reads_dir"/*.fasta)

# List to store file names without the ".fasta" suffix
parings_to_run=()

# Iterate over the file names
for file in "${fasta_files[@]}"; do
    # Remove the ".fasta" suffix and append to the list
    parings_to_run+=("$(basename "$file" .fasta)")
done

repeatmasker_dir=$base_name_dir/paired_spacer_repeatmasker_results/

mkdir -p $repeatmasker_dir

for pairing in "${parings_to_run[@]}"; do

    echo "Starting $pairing reads"

    query=$input_reads_dir/$pairing
    database=references/pairings_for_spacers/$strain_number\_pairings/$strain_number\_paired_$pairing

    echo "Running RepeatMasker on reads $query.fasta with reference $database.fasta"

    RepeatMasker $query.fasta -lib $database.fasta -s -pa 12 --cutoff 500 -no_is -norna -gff -dir $repeatmasker_dir


    filtered_repeatmasker_data_file=$repeatmasker_dir/$pairing\_spacer_filtered.out
    grep "_from_repeat_to_plus_50kb" $repeatmasker_dir/$pairing.fasta.out \
    > $filtered_repeatmasker_data_file

    repeatmasker_file_name=$repeatmasker_dir/$base_name\_$pairing\_spacer_repeatmasker_results.ssv
    echo -e \
    "SW_score divergence_percent deletion_percent insertion_percent read_id match_start_on_read match_end_on_read leftover_on_read strand chr_end_tract section_number match_start_on_chr_end_section match_end_on_chr_end_section leftover_on_chr_end_section match_id sub_match" \
    > $repeatmasker_file_name
    cat $filtered_repeatmasker_data_file >> $repeatmasker_file_name

done


python scripts/make_spacer_pairs_repeatmasker_tsv.py $base_name $strain_number


echo "Done with Step 12."

#######################################################################################

# Step 13
