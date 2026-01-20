import os
import glob

# Load configuration
configfile: "config.yaml"

# Variables from config
BASE = config["base_name"]
ANCHOR = config["anchor_set"]
BAM_IN = f"{config['bam_dir']}/{BASE}.bam"
FASTA_IN = f"results/{BASE}/{BASE}.fasta"
# Use strain from config if provided, otherwise extract from base name
STRAIN = config.get("strain", BASE.split('_')[1] if 'dorado_' in BASE else BASE.split('_')[0])

# Define the 16 yeast chromosomes and their two sides (L and R)
CHROMS = [str(i) for i in range(1, 17)]
SIDES = ["L", "R"]

# This creates a list like ['chr1L', 'chr1R', 'chr2L', ..., 'chr17R']
CHROM_SIDES = [f"chr{c}{s}" for c in CHROMS for s in SIDES]

# Paths to label_regions.sh outputs (created by running label_regions.sh before this pipeline)
# These are located in the assembly directory created by create_ref.sh
ASSEMBLY_DIR = f"results/{BASE}/assembly_{STRAIN}"
PRETELOMERIC_LABELS_DIR = f"{ASSEMBLY_DIR}/pretelomeric_labels"
REPEATMASKER_YPRIMES_FASTA = f"{PRETELOMERIC_LABELS_DIR}/extracted_yprimes_{STRAIN}.fasta"
FEATURES_BED = f"{PRETELOMERIC_LABELS_DIR}/pretelomeric_regions_{STRAIN}_simp.bed"

rule all:
    input:
        f"results/{BASE}/{BASE}_post_y_prime_probe.tsv",
        f"results/{BASE}/{BASE}_stats_y_prime.txt",
        f"results/{BASE}/y_prime_blast/all_{BASE}_probe_matches.tsv",
        directory(f"results/{BASE}/figures_for_y_primes"),
        f"results/{BASE}/{BASE}_y_prime_recombination.tsv",
        # X element and spacer pairing analysis - requires strain-specific reference pairings
        # Uncomment when pairings_for_x_element_ends/{STRAIN}_pairings directory is created
        # f"results/{BASE}/{BASE}_paired_x_element_ends_repeatmasker.tsv",
        # f"results/{BASE}/{BASE}_good_x_element_ends_paired_repeatmasker.tsv",
        # f"results/{BASE}/{BASE}_good_gained_y_x_element_ends_paired_repeatmasker.tsv",
        # Spacer analysis - requires strain-specific reference pairings and BED files
        # f"results/{BASE}/{BASE}_paired_spacer_repeatmasker.tsv",
        # f"results/{BASE}/{BASE}_paired_good_spacer_repeatmasker.tsv",
        # f"results/{BASE}/{BASE}_paired_good_gained_spacer_repeatmasker.tsv"

# --- STEP 0: Automatic Detection ---
# If the FASTA doesn't exist, we create this rule to generate it from BAM.
if not os.path.exists(FASTA_IN):
    rule prep_from_bam:
        input:
            bam = BAM_IN
        output:
            fastq = f"results/{BASE}/{BASE}.fastq",
            fasta = FASTA_IN
        threads: 16
        shell:
            """
            samtools fastq -@ {threads} -T '*' {input.bam} > {output.fastq}
            seqtk seq -A {output.fastq} > {output.fasta}
            """

# New rule specifically for indexing
rule index_fasta:
    input:
        fasta = FASTA_IN
    output:
        fai = f"{FASTA_IN}.fai"
    shell:
        "samtools faidx {input.fasta}"


# --- STEP 1: Blast for Reads with Anchors ---
rule blast_anchors:
    input:
        # Snakemake handles the dependency whether it comes from 'prep_from_bam' 
        # or an existing file on disk.
        query = FASTA_IN,
        db = config["references"]["anchors"]
    output:
        tsv = f"results/{BASE}/{BASE}_blasted_{ANCHOR}.tsv"
    threads: 56
    shell:
        """
        echo -e "read_id\ttotal_read_length\tread_bp_used_for_match\tmatch_start_on_read\tmatch_end_on_read\tanchor_name\ttotal_anchor_length\tmatch_start_on_anchor\tmatch_end_on_anchor\tpident\tbitscore\tevalue" > {output.tsv}
        
        blastn -query {input.query} -db {input.db} -task megablast \
            -perc_identity 85 -min_raw_gapped_score 3000 -num_threads {threads} \
            -outfmt "6 qseqid qlen length qstart qend sseqid slen sstart send pident bitscore evalue" \
            >> {output.tsv}
        """

# --- STEP 2: Filter Blast Results ---
rule filter_anchors:
    input:
        tsv = rules.blast_anchors.output.tsv
    output:
        # Adjust this filename if your python script produces something different
        all_matches = f"results/{BASE}/all_matches_{BASE}_blasted_{ANCHOR}.tsv",
        top_matches = f"results/{BASE}/top_matches_{BASE}_blasted_{ANCHOR}.tsv"
    shell:
        "python scripts/filter_for_reads_with_anchors.py {input.tsv} {output.all_matches} {output.top_matches}"

# --- STEP 3: Split and Label ---
rule split_and_label:
    input:
        tsv = rules.filter_anchors.output.top_matches,
        fasta = FASTA_IN,
        fai = f"{FASTA_IN}.fai"
    output:
        # Instead of directory(), we list all 34 expected files
        reads = expand("results/{{BASE}}/chr_anchor_included_individual_files/{{BASE}}_blasted_{{ANCHOR}}_{cs}_anchor_reads.fasta", 
                       cs=CHROM_SIDES)
    shell:
        """
        mkdir -p results/{BASE}/chr_anchor_included_individual_files/
        python scripts/split_and_label_all_reads_include_anchor.py {input.tsv} {input.fasta} results/{BASE}/chr_anchor_included_individual_files/ {BASE} {ANCHOR}
        """

rule aggregate_chromosomes:
    input:
        # This tells Snakemake to wait until all chromosome files from Step 3 are complete
        reads = expand(f"results/{BASE}/chr_anchor_included_individual_files/{BASE}_blasted_{ANCHOR}_{{cs}}_anchor_reads.fasta",
                       cs=CHROM_SIDES)
    output:
        merged = f"results/{BASE}/{BASE}_all_chromosome_anchored_reads_pre_trim.fasta"
    shell:
        # Concatenate all chromosome-specific FASTA files
        "cat {input.reads} > {output.merged}"

rule porechop_trim:
    input:
        fasta = rules.aggregate_chromosomes.output.merged,
        adapters = config["references"]["adapters"]
    output:
        fastq = f"results/{BASE}/porechop_trim/all_chromosome_anchored_read_trimmed.fastq",
        fasta = f"results/{BASE}/porechop_trim/all_chromosome_anchored_read_trimmed.fasta",
        log = f"results/{BASE}/{BASE}_all_chrs_split_to_telomere_porechopped.log"
    threads: 56
    shell:
        """
        # Run Porechop
        porechop_abi -t {threads} --format fastq -ddb -v 3 --no_split \
            -cap {input.adapters} -i {input.fasta} \
            -o {output.fastq} > {output.log}

        # Convert FASTQ to FASTA using sed as per original script
        sed -n '1~4s/^@/>/p;2~4p' {output.fastq} > {output.fasta}
        """

# --- STEP 6: Create Sequencing Summary ---
# (Optimized: replaced Python script with seqkit for 10-100x speed improvement)
rule sequence_summary:
    input:
        fasta = FASTA_IN
    output:
        summary = f"results/{BASE}/{BASE}_sequencing_summary.tsv"
    shell:
        """
        echo -e "read_id\\tsequence_length_template" > {output.summary}
        seqkit fx2tab -n -l -i {input.fasta} >> {output.summary}
        """

# --- STEP 7: Parse Porechop Logs ---
# (From check_for_adapters.py)
rule parse_porechop_log:
    input:
        log = f"results/{BASE}/{BASE}_all_chrs_split_to_telomere_porechopped.log"
    output:
        tsv = f"results/{BASE}/{BASE}_porechopped_results.tsv"
    shell:
        "python scripts/check_for_adapters.py {BASE} {input.log} {output.tsv}"

# --- STEP 8: Compare Callers (QC Stats) ---
# (From compare_adapter_callers_dorado.py)
rule compare_callers:
    input:
        summary = rules.sequence_summary.output.summary,
        porechop = rules.parse_porechop_log.output.tsv,
        blast_raw = f"results/{BASE}/{BASE}_blasted_{ANCHOR}.tsv",
        blast_top = f"results/{BASE}/top_matches_{BASE}_blasted_{ANCHOR}.tsv"
    output:
        stats = f"results/{BASE}/{BASE}_adapter_trimming_check.stats",
        table = f"results/{BASE}/{BASE}_adapter_trimming_check.tsv"
    shell:
        """
        python scripts/compare_adapter_callers_dorado.py \
            {input.summary} \
            {input.porechop} \
            {input.blast_raw} \
            {input.blast_top} \
            {output.stats} \
            {output.table} \
            {ANCHOR}
        """

# Create a dedicated rule for indexing the merged fasta file
rule index_merged_fasta:
    input:
        rules.aggregate_chromosomes.output.merged
    output:
        fai = f"{rules.aggregate_chromosomes.output.merged}.fai"
    shell:
        "samtools faidx {input}"

# --- STEP 9: Fine Telomere Trimming ---
# (From fine_telomere_trimming.py)
rule fine_trimming:
    input:
        fasta = rules.aggregate_chromosomes.output.merged,
        fai = f"{rules.aggregate_chromosomes.output.merged}.fai",
        adapter_info = rules.compare_callers.output.table,
        best_anchor = f"results/{BASE}/top_matches_{BASE}_blasted_{ANCHOR}.tsv"
    output:
        main_tsv = f"results/{BASE}/{BASE}_post_telo_trimming.tsv",
        fasta_trimmed = f"results/{BASE}/{BASE}_fine_trimmed.fasta",
        trim_dir = directory(f"results/{BASE}/repeat_trim_files/")
    shell:
        """
        mkdir -p {output.trim_dir}
        python scripts/fine_telomere_trimming.py \
            {input.fasta} \
            {input.adapter_info} \
            {input.best_anchor} \
            {output.main_tsv} \
            {output.trim_dir} \
            {BASE} \
            {output.fasta_trimmed} 
        """

# optional
#rule make_probe_db:
    #input:
        #fasta = config["references"]["probe_fasta"]
    #output:
        ## BLAST creates several files; we track the main database file
        #db_file = f"{config['references']['probe_fasta']}.nin"
    #shell:
        #"makeblastdb -in {input.fasta} -dbtype nucl"

rule blast_y_primes:
    input:
        # This relies on the individual FASTA files created in Step 3
        query = "results/{BASE}/chr_anchor_included_individual_files/{BASE}_blasted_" + ANCHOR + "_{chrom_side}_anchor_reads.fasta",
        # Assuming your probe database is in a references folder
        db = config["references"]["probe"]
    output:
        tsv = "results/{BASE}/y_prime_blast/{BASE}_{chrom_side}_blasted_probe.tsv"
    threads: 4  # Running 14 jobs at once (14 * 4 = 56 cores) is much faster than a loop!
    shell:
        """
        # Create the TSV header exactly like your original script
        echo -e "read_id\\ttotal_read_length\\tread_bp_used_for_match\\tmatch_start_on_read\\tmatch_end_on_read\\tanchor_name\\ttotal_anchor_length\\tmatch_start_on_anchor\\tmatch_end_on_anchor\\tpident\\tbitscore\\tevalue" > {output.tsv}
        
        # Run BLAST
        blastn -query {input.query} -db {input.db} -perc_identity 90 -num_threads {threads} \
            -outfmt "6 qseqid qlen length qstart qend sseqid slen sstart send pident bitscore evalue" \
            >> {output.tsv}
        """

rule aggregate_y_prime_blasts:
    input:
        # Wait for all 32 chromosome-side BLAST results
        expand("results/{BASE}/y_prime_blast/{BASE}_{chrom_side}_blasted_probe.tsv",
               BASE=BASE, chrom_side=CHROM_SIDES)
    output:
        # Create the probe directory as a marker that all blasts are done
        probe_dir = directory(f"results/{BASE}/y_prime_blast")
    shell:
        "echo 'All Y prime BLAST jobs complete'"

rule y_prime_analysis:
    input:
        best_anchor = f"results/{BASE}/top_matches_{BASE}_blasted_{ANCHOR}.tsv",
        telo_results = f"results/{BASE}/{BASE}_post_telo_trimming.tsv",
        # Directly depend on all 32 blast files
        probe_blasts = expand("results/{BASE}/y_prime_blast/{BASE}_{chrom_side}_blasted_probe.tsv",
                             BASE=BASE, chrom_side=CHROM_SIDES)
    output:
        main_tsv = f"results/{BASE}/{BASE}_post_y_prime_probe.tsv",
        stats = f"results/{BASE}/{BASE}_stats_y_prime.txt",
        combined_probe = f"results/{BASE}/y_prime_blast/all_{BASE}_probe_matches.tsv",
        figures = directory(f"results/{BASE}/figures_for_y_primes")
    params:
        base_name = BASE,
        anchor_set = ANCHOR,
        probe_dir = f"results/{BASE}/y_prime_blast"  # Pass as parameter instead
    shell:
        """
        python scripts/y_prime_analysis.py {params.base_name} {params.anchor_set} \
            --top-anchor-blast {input.best_anchor} \
            --telomere-repeat-results {input.telo_results} \
            --probe-blast-dir {params.probe_dir} \
            --output-tsv {output.main_tsv} \
            --output-stats {output.stats} \
            --figure-dir {output.figures} \
            --combined-probe-output {output.combined_probe}
        """

rule repeatmasker_y_primes:
    input:
        query = "results/{base}/chr_anchor_included_individual_files/{base}_blasted_{anchor}_{chrom_side}_anchor_reads.fasta",
        database = REPEATMASKER_YPRIMES_FASTA  # From label_regions.sh output
    output:
        out_file = "results/{base}/read_repeatmasker_results/{base}_blasted_{anchor}_{chrom_side}_anchor_reads.fasta.out",
        filtered = "results/{base}/read_repeatmasker_results/{base}_{anchor}_{chrom_side}_filtered.out",
        ssv = "results/{base}/read_repeatmasker_results/{base}_{anchor}_{chrom_side}_repeatmasker_results.ssv"
    params:
        outdir = "results/{base}/read_repeatmasker_results/"
    threads: 12
    shell:
        """
        mkdir -p {params.outdir}
        
        RepeatMasker {input.query} -lib {input.database} -s -pa {threads} \
            --cutoff 1000 -no_is -norna -gff -dir {params.outdir}
        
        # Create header
        echo -e "SW_score\\tdivergence_percent\\tdeletion_percent\\tinsertion_percent\\tread_id\\tmatch_start_on_read\\tmatch_end_on_read\\tleftover_on_read\\tstrand\\ty_prime_id\\ty_prime_group\\tmatch_start_on_y_prime\\tmatch_end_on_y_prime\\tleftover_on_y_prime\\tmatch_id\\tsub_match" > {output.ssv}
        
        # Filter for Y_Prime and append
        grep "Y_Prime" {output.out_file} > {output.filtered} || touch {output.filtered}
        cat {output.filtered} >> {output.ssv}
        """

rule aggregate_repeatmasker_y_primes:
    input:
        expand("results/{base}/read_repeatmasker_results/{base}_{anchor}_{chrom_side}_repeatmasker_results.ssv",
               base=BASE, anchor=ANCHOR, chrom_side=CHROM_SIDES)
    output:
        marker = touch(f"results/{BASE}/read_repeatmasker_results/.done")
    shell:
        "echo 'All RepeatMasker Y prime jobs complete'"

rule make_y_prime_repeatmasker_tsv:
    input:
        repeatmasker_done = f"results/{BASE}/read_repeatmasker_results/.done",
        y_prime_probe = f"results/{BASE}/{BASE}_post_y_prime_probe.tsv"
    output:
        all_repeatmasker = f"results/{BASE}/{BASE}_repeatmasker.tsv",
        good_end = f"results/{BASE}/{BASE}_good_end_y_repeatmasker.tsv",
        gained_y = f"results/{BASE}/{BASE}_gained_y_repeatmasker.tsv"
    params:
        repeatmasker_dir = f"results/{BASE}/read_repeatmasker_results/"
    shell:
        """
        python scripts/make_y_prime_repeatmasker_tsv.py \
            {params.repeatmasker_dir} \
            {input.y_prime_probe} \
            {output.all_repeatmasker} \
            {output.good_end} \
            {output.gained_y}
        """

rule get_stats_of_recombination:
    input:
        good_end_y = f"results/{BASE}/{BASE}_good_end_y_repeatmasker.tsv",
        y_prime_probe = f"results/{BASE}/{BASE}_post_y_prime_probe.tsv"
    output:
        recomb = f"results/{BASE}/{BASE}_y_prime_recombination.tsv"
    params:
        strain = STRAIN
    shell:
        """
        python scripts/get_stats_of_recombination.py \
            {input.good_end_y} \
            {input.y_prime_probe} \
            {params.strain} \
            {output.recomb}
        """

checkpoint make_pairings_from_y_primes_all_ends:
    input:
        good_end_y = f"results/{BASE}/{BASE}_good_end_y_repeatmasker.tsv",
        chr_reads = expand(f"results/{BASE}/chr_anchor_included_individual_files/{BASE}_blasted_{ANCHOR}_{{cs}}_anchor_reads.fasta",
                          cs=CHROM_SIDES)
    output:
        pairings_dir = directory(f"results/{BASE}/paired_by_y_prime_reads/")
    params:
        strain = STRAIN,
        anchor = ANCHOR,
        base_name = BASE
    shell:
        """
        mkdir -p {output.pairings_dir}
        python scripts/make_pairings_from_y_primes_all_ends.py \
            {input.good_end_y} \
            results/{params.base_name}/chr_anchor_included_individual_files/ \
            {output.pairings_dir} \
            {params.strain} \
            {params.anchor} \
            {params.base_name}
        """

def get_pairing_names(wildcards):
    """Get list of pairing file basenames (without .fasta extension) from checkpoint"""
    checkpoint_output = checkpoints.make_pairings_from_y_primes_all_ends.get(**wildcards).output.pairings_dir
    # List all .fasta files in the directory
    pairing_files = glob.glob(f"{checkpoint_output}/*.fasta")
    # Return basenames without .fasta extension
    return [os.path.basename(f).replace('.fasta', '') for f in pairing_files]

rule repeatmasker_x_elements:
    input:
        query = "results/{base}/paired_by_y_prime_reads/{pairing}.fasta",
        database = f"references/pairings_for_x_element_ends/{STRAIN}_pairings/{STRAIN}_paired_{{pairing}}.fasta"
    output:
        out_file = "results/{base}/paired_x_element_ends_repeatmasker_results/{pairing}.fasta.out",
        filtered = "results/{base}/paired_x_element_ends_repeatmasker_results/{pairing}_x_element_ends_filtered.out",
        ssv = "results/{base}/paired_x_element_ends_repeatmasker_results/{base}_{pairing}_x_element_ends_repeatmasker_results.ssv"
    params:
        outdir = "results/{base}/paired_x_element_ends_repeatmasker_results/"
    threads: 12
    shell:
        """
        mkdir -p {params.outdir}
        
        RepeatMasker {input.query} -lib {input.database} -s -pa {threads} \
            --cutoff 500 -no_is -norna -gff -dir {params.outdir}
        
        # Create header
        echo -e "SW_score\\tdivergence_percent\\tdeletion_percent\\tinsertion_percent\\tread_id\\tmatch_start_on_read\\tmatch_end_on_read\\tleftover_on_read\\tstrand\\tx_element_ends\\tsection_number\\tmatch_start_on_chr_end_section\\tmatch_end_on_chr_end_section\\tleftover_on_chr_end_section\\tmatch_id\\tsub_match" > {output.ssv}
        
        # Filter for x_ends and append
        grep "x_ends" {output.out_file} > {output.filtered} || touch {output.filtered}
        cat {output.filtered} >> {output.ssv}
        """

def get_pairing_names_with_refs(wildcards):
    """Get list of pairings that have corresponding reference databases"""
    pairing_names = get_pairing_names(wildcards)
    # Filter to only pairings that have reference files
    valid_pairings = []
    for pairing in pairing_names:
        ref_file = f"references/pairings_for_x_element_ends/{STRAIN}_pairings/{STRAIN}_paired_{pairing}.fasta"
        if os.path.exists(ref_file):
            valid_pairings.append(pairing)
        else:
            print(f"Warning: No reference file for pairing {pairing}, skipping")
    return valid_pairings

def aggregate_x_element_inputs(wildcards):
    """Aggregate all x_element RepeatMasker results based on checkpoint output"""
    pairing_names = get_pairing_names_with_refs(wildcards)
    return expand("results/{base}/paired_x_element_ends_repeatmasker_results/{base}_{pairing}_x_element_ends_repeatmasker_results.ssv",
                  base=wildcards.base,
                  pairing=pairing_names)

rule aggregate_x_elements:
    input:
        aggregate_x_element_inputs
    output:
        marker = touch("results/{base}/paired_x_element_ends_repeatmasker_results/.done")
    shell:
        "echo 'All X element RepeatMasker jobs complete'"
    
rule make_x_element_ends_pairs_repeatmasker_tsv:
    input:
        x_elements_done = "results/{base}/paired_x_element_ends_repeatmasker_results/.done",
        pairings_dir = "results/{base}/paired_by_y_prime_reads/",
        y_prime_probe = "results/{base}/{base}_post_y_prime_probe.tsv"
    output:
        all_tsv = "results/{base}/{base}_paired_x_element_ends_repeatmasker.tsv",
        good_tsv = "results/{base}/{base}_good_x_element_ends_paired_repeatmasker.tsv",
        gained_tsv = "results/{base}/{base}_good_gained_y_x_element_ends_paired_repeatmasker.tsv"
    params:
        strain = STRAIN,
        repeatmasker_dir = "results/{base}/paired_x_element_ends_repeatmasker_results/"
    shell:
        """
        python scripts/make_x_element_ends_pairs_repeatmasker_tsv.py \
            {params.repeatmasker_dir} \
            {params.strain} \
            {input.y_prime_probe} \
            {output.all_tsv} \
            {output.good_tsv} \
            {output.gained_tsv}
        """

#######################################################################################
# Step 12: Find 250bp tracts of spacer sequences
#######################################################################################

rule repeatmasker_spacers:
    input:
        query = "results/{base}/paired_by_y_prime_reads/{pairing}.fasta",
        database = f"references/pairings_for_spacers/{STRAIN}_pairings/{STRAIN}_paired_{{pairing}}.fasta"
    output:
        out_file = "results/{base}/paired_spacer_repeatmasker_results/{pairing}.fasta.out",
        filtered = "results/{base}/paired_spacer_repeatmasker_results/{pairing}_spacer_filtered.out",
        ssv = "results/{base}/paired_spacer_repeatmasker_results/{base}_{pairing}_spacer_repeatmasker_results.ssv"
    params:
        outdir = "results/{base}/paired_spacer_repeatmasker_results/"
    threads: 12
    shell:
        """
        mkdir -p {params.outdir}

        RepeatMasker {input.query} -lib {input.database} -s -pa {threads} \
            --cutoff 500 -no_is -norna -gff -dir {params.outdir}

        # Create header
        echo -e "SW_score\\tdivergence_percent\\tdeletion_percent\\tinsertion_percent\\tread_id\\tmatch_start_on_read\\tmatch_end_on_read\\tleftover_on_read\\tstrand\\tchr_end_tract\\tsection_number\\tmatch_start_on_chr_end_section\\tmatch_end_on_chr_end_section\\tleftover_on_chr_end_section\\tmatch_id\\tsub_match" > {output.ssv}

        # Filter for spacer sequences and append
        grep "_from_repeat_to_plus_50kb" {output.out_file} > {output.filtered} || touch {output.filtered}
        cat {output.filtered} >> {output.ssv}
        """

def aggregate_spacer_inputs(wildcards):
    """Aggregate all spacer RepeatMasker results based on checkpoint output"""
    pairing_names = get_pairing_names_with_refs(wildcards)
    return expand("results/{base}/paired_spacer_repeatmasker_results/{base}_{pairing}_spacer_repeatmasker_results.ssv",
                  base=wildcards.base,
                  pairing=pairing_names)

rule aggregate_spacers:
    input:
        aggregate_spacer_inputs
    output:
        marker = touch("results/{base}/paired_spacer_repeatmasker_results/.done")
    shell:
        "echo 'All spacer RepeatMasker jobs complete'"

rule make_spacer_pairs_repeatmasker_tsv:
    input:
        spacers_done = "results/{base}/paired_spacer_repeatmasker_results/.done",
        pairings_dir = "results/{base}/paired_by_y_prime_reads/",
        y_prime_probe = "results/{base}/{base}_post_y_prime_probe.tsv",
        anchors_bed = f"references/{STRAIN}_anchors_and_distances.bed",
        features_bed = FEATURES_BED  # From label_regions.sh output (simplified BED)
    output:
        all_tsv = "results/{base}/{base}_paired_spacer_repeatmasker.tsv",
        good_tsv = "results/{base}/{base}_paired_good_spacer_repeatmasker.tsv",
        gained_tsv = "results/{base}/{base}_paired_good_gained_spacer_repeatmasker.tsv"
    params:
        strain = STRAIN,
        repeatmasker_dir = "results/{base}/paired_spacer_repeatmasker_results/"
    shell:
        """
        python scripts/make_spacer_pairs_repeatmasker_tsv.py \
            {params.repeatmasker_dir} \
            {params.strain} \
            {input.anchors_bed} \
            {input.features_bed} \
            {input.y_prime_probe} \
            {output.all_tsv} \
            {output.good_tsv} \
            {output.gained_tsv}
        """