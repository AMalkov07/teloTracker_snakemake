#$ -q TELOMERE2,UI
#$ -pe smp 12
#$ -j y
# -o /dev/null
#$ -cwd

# Base name of the file/run wanting to analyze
#Example to put after submission = dorado_2191_repeat1_old-detect-only
base_name=dorado_6991_day0_repeat2_PromethION_no_tag_yes_rejection

# Anchor set you want to blast for chromosome anchors
# Default is new_all_best_anchors
### The new default is test_anchors
anchor_set=test_anchors



#######################################################################################
#######################################################################################

# Other paths that will be used - Don't touch
telomere_analysis_dir=./
base_name_dir=$base_name/
strain_number=${base_name%%_day*}
strain_number=${strain_number#dorado_}
base_name_chr_anchor_outputs_dir=$base_name/chr_anchor_included_individual_files/
blast_y_primes_output_dir=$base_name/y_prime_blast_results_from_chr_anchored_reads/
base_name_repeat_split_dir=$base_name/split_to_repeat_outputs/
base_name_porechop_dir=$base_name/porechop_trim/

adapter_seq_file=offical_nanopore_adapter_seq+trunc.txt

porechop_abi_log_file=$base_name_dir/$base_name\_all_chrs_split_to_telomere_porechopped.log

bam_dir=../samples_dorado_basecalled/


#######################################################################################
#######################################################################################

#Order of steps/commands/scripts run:




#######################################################################################
#######################################################################################

echo "Starting $base_name"

#mamba activate telotracker

# Rerun this part

# Read through the y prime probe blast results, make figures for Y primes, and the outputs directory final tsv

#python fine_telomere_trimming.py $base_name $anchor_set

#python y_prime_analysis.py $base_name $anchor_set


# Step 10

# Blasting the Y' sequences - without redundant Y' sequences

echo "Step 10. Finding Y' sequences by RepeatMasker"

repeatmasker_y_prime_dir=$base_name_dir/read_repeatmasker_results/

mkdir -p $repeatmasker_y_prime_dir


# tldr of loop is filter out any reads that have some y primes in them
for chr_num in {1..16}
do
    for side in {"L","R"}
    do

        echo "Searching $chr_num$side reads"

        query=$base_name_chr_anchor_outputs_dir/$base_name\_blasted_$anchor_set\_chr$chr_num$side\_anchor_reads # figure out exactly what these files are and when htey are created
        database="repeatmasker_${strain_number}_all_y_primes"
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



python make_y_prime_repeatmasker_tsv.py $base_name

#python get_summary_stats_for_y_prime_repeatmasker.py $base_name #?? needs more inputs # don't think I need this line
python get_stats_of_recombination.py $base_name #script outputs *_recombination.tsv file and it seems to rely on the make_y_prime_repeatmasker_tsv.py script but not the script right above this one

#####################
# Only for 0 --> 1 Y's
#python make_pairings_from_y_primes.py $base_name $anchor_set $strain_number

# For all change in first Y's
python make_pairings_from_y_primes_all_ends.py $base_name #also relies on the make_y_prime_repeatmasker_tsv.py script, but not the other 2 python scripts above this one
#####################

echo "Done with Step 10."




##########################
#########################
######################
#########################



#echo "Starting $base_name"

##mamba activate telotracker

## Rerun this part

## Read through the y prime probe blast results, make figures for Y primes, and the outputs directory final tsv

#python fine_telomere_trimming.py $base_name $anchor_set

#python y_prime_analysis.py $base_name $anchor_set


## Step 10

## Blasting the Y' sequences - without redundant Y' sequences

#echo "Step 10. Finding Y' sequences by RepeatMasker"

#repeatmasker_y_prime_dir=$base_name_dir/read_repeatmasker_results/

#mkdir -p $repeatmasker_y_prime_dir


#for chr_num in {1..16}
#do
    #for side in {"L","R"}
    #do

        #echo "Searching $chr_num$side reads"

        #query=$base_name_chr_anchor_outputs_dir/$base_name\_blasted_$anchor_set\_chr$chr_num$side\_anchor_reads
        #database="repeatmasker_${strain_number}_all_y_primes"
        #repeatmasker_file_name=$repeatmasker_y_prime_dir/$base_name\_chr$chr_num$side\_repeatmasker_results.ssv

        #RepeatMasker $query.fasta -lib $database.fasta -s -pa 12 --cutoff 1000 -no_is -norna -gff -dir $repeatmasker_y_prime_dir

        #echo -e \
        #"SW_score divergence_percent deletion_percent insertion_percent read_id match_start_on_read match_end_on_read leftover_on_read strand y_prime_id y_prime_group match_start_on_y_prime match_end_on_y_prime leftover_on_y_prime match_id sub_match" \
        #> $repeatmasker_file_name

        #filtered_repeatmasker_data_file=$repeatmasker_y_prime_dir/chr$chr_num$side\_filtered.out

        #grep "Y_Prime" $repeatmasker_y_prime_dir/$base_name\_blasted_$anchor_set\_chr$chr_num$side\_anchor_reads.fasta.out \
        #> $filtered_repeatmasker_data_file

        #cat $filtered_repeatmasker_data_file >> $repeatmasker_file_name

    #done
#done


#python make_y_prime_repeatmasker_tsv.py $base_name

#python get_summary_stats_for_y_prime_repeatmasker.py $base_name

#python get_stats_of_recombination.py $base_name

######################
## Only for 0 --> 1 Y's
##python make_pairings_from_y_primes.py $base_name $anchor_set $strain_number

## For all change in first Y's
#python make_pairings_from_y_primes_all_ends.py $base_name
######################

#echo "Done with Step 10."

#######################################################################################


# Step 11. Find any x element switches

echo "Step 11. Find any x element end switches"


# RepeatMasker for x element for y prime pairings

echo "Finding x element end sequences by RepeatMasker for y prime pairings"


input_reads_dir=$base_name/paired_by_y_prime_reads/

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
    database=pairings_for_x_element_ends/$strain_number\_pairings/$strain_number\_paired_$pairing

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


python make_x_element_ends_pairs_repeatmasker_tsv.py $base_name $strain_number


echo "Done with Step 11."

#######################################################################################

# Step 12
# Determine when an end switches from one chromosome arm to another
# Finding 250bp tracts of spacer sequences by RepeatMasker for spacer sequences for y prime pairings

echo "Finding 250bp tracts of spacer sequences by RepeatMasker for y prime pairings"


input_reads_dir=$base_name/paired_by_y_prime_reads/

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
    database=pairings_for_spacers/$strain_number\_pairings/$strain_number\_paired_$pairing

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


python make_spacer_pairs_repeatmasker_tsv.py $base_name $strain_number


echo "Done with Step 12."

#######################################################################################

# Step 13



