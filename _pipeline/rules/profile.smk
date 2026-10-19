# Read profiling (steps 0-6), included by _pipeline/Snakefile. Uses only the shared settings
# defined there (BASE, ANCHOR, paths, INPUT_FILE / INPUT_TYPE, the profiling options); needs no
# strain or day-0 files.

# ===========================================================================
# CORE ANALYSIS (Steps 0-6)
# ===========================================================================

if INPUT_KIND == "contigs":

    # -----------------------------------------------------------------------
    # Step 0 (contigs): each contig cut into two halves, one per chromosome end
    # -----------------------------------------------------------------------

    rule split_contig_ends:
        input: INPUT_FILE
        output: f"{RESULTS}/{BASE}.fasta"
        threads: 1
        shell: "python _pipeline/scripts/split_contig_ends.py {input} {output}"

else:

    # ---------------------------------------------------------------------------
    # Step 0a: Convert input to raw FASTQ
    # ---------------------------------------------------------------------------

    rule prep_raw_fastq:
        input: INPUT_FILE
        output: f"{RESULTS}/{BASE}_raw.fastq"
        threads: 1
        run:
            if INPUT_TYPE == "bam":
                shell("samtools fastq -T '*' {input} > {output}")
            elif INPUT_TYPE == "fastqgz":
                shell("zcat {input} > {output}")
            else:
                shell("cp {input} {output}")

    # ---------------------------------------------------------------------------
    # Step 0b: Filter reads (>= min_read_length bp, qscore >= min_qscore; defaults 2000, 10)
    # ---------------------------------------------------------------------------

    rule filter_reads:
        input: f"{RESULTS}/{BASE}_raw.fastq"
        output: f"{RESULTS}/{BASE}.fastq"
        threads: 1
        shell: "python _pipeline/scripts/filter_reads.py {input} {output} {MIN_READ_LENGTH} {MIN_QSCORE}"

    # ---------------------------------------------------------------------------
    # Step 0c: Convert filtered FASTQ to FASTA
    # ---------------------------------------------------------------------------

    rule convert_to_fasta:
        input: f"{RESULTS}/{BASE}.fastq"
        output: f"{RESULTS}/{BASE}.fasta"
        threads: 1
        shell: "sed -n '1~4s/^@/>/p;2~4p' {input} > {output}"

# ---------------------------------------------------------------------------
# Step 1: BLAST reads for chromosome anchors
# ---------------------------------------------------------------------------

rule blast_anchors:
    input:
        fasta = f"{RESULTS}/{BASE}.fasta",
        db    = ANCHOR_DB,
    output: f"{BLAST_DIR}/{BASE}_blasted_{ANCHOR}.tsv"
    threads: workflow.cores
    shell:
        """
        mkdir -p {BLAST_DIR}
        echo -e "read_id\\ttotal_read_length\\tread_bp_used_for_match\\tmatch_start_on_read\\tmatch_end_on_read\\tanchor_name\\ttotal_anchor_length\\tmatch_start_on_anchor\\tmatch_end_on_anchor\\tpident\\tbitscore\\tevalue" > {output}
        blastn -query {input.fasta} -db {input.db} \
            -task dc-megablast -perc_identity 85 -min_raw_gapped_score {MIN_GAPPED_SCORE} \
            -num_threads {threads} \
            -outfmt "6 qseqid qlen length qstart qend sseqid slen sstart send pident bitscore evalue" \
            >> {output}
        """

# ---------------------------------------------------------------------------
# Step 2: Filter BLAST results
# ---------------------------------------------------------------------------

rule detect_foldbacks:
    """Flag hairpin reads: the same anchor twice, in opposite orientations.

    Re-BLASTs each ANCHORED read against ONLY its assigned anchor, without
    -min_raw_gapped_score (that flag filters on the preliminary gapped score,
    which systematically hides the degraded second copy of a fold-back).
    Detection only -- it cannot reassign a read or anchor a new one.
    """
    input:
        blast   = f"{BLAST_DIR}/{BASE}_blasted_{ANCHOR}.tsv",
        fasta   = f"{RESULTS}/{BASE}.fasta",
        anchors = ANCHOR_DB,   # db prefix is the FASTA path itself
    output:
        ids    = f"{BLAST_DIR}/{BASE}_foldback_read_ids.txt",
        detail = f"{BLAST_DIR}/{BASE}_foldback_detail.tsv",
    threads: workflow.cores
    shell:
        """
        samtools faidx {input.fasta}
        python _pipeline/scripts/detect_foldback_reads.py \
            --blast-tsv     {input.blast} \
            --reads-fasta   {input.fasta} \
            --anchors-fasta {input.anchors} \
            --output        {output.ids} \
            --detail-tsv    {output.detail} \
            --threads       {threads}
        """


rule filter_anchors:
    input:
        blast    = f"{BLAST_DIR}/{BASE}_blasted_{ANCHOR}.tsv",
        foldback = f"{BLAST_DIR}/{BASE}_foldback_read_ids.txt",
    output: f"{BLAST_DIR}/top_matches_{BASE}_blasted_{ANCHOR}.tsv"
    threads: 1
    shell:
        "python _pipeline/scripts/filter_for_reads_with_anchors.py "
        "{BASE} {ANCHOR} {RESULTS} {input.foldback}"

# ---------------------------------------------------------------------------
# Step 3: Split and label reads per chromosome end
# ---------------------------------------------------------------------------

rule split_and_label:
    input: f"{BLAST_DIR}/top_matches_{BASE}_blasted_{ANCHOR}.tsv"
    output:
        fasta = expand(f"{CHR_DIR}/{BASE}_blasted_{ANCHOR}_{{ce}}_anchor_reads.fasta", ce=CHROM_ENDS),
        tsv   = expand(f"{CHR_DIR}/{BASE}_blasted_{ANCHOR}_{{ce}}_anchor_reads.tsv",   ce=CHROM_ENDS),
    threads: 1
    shell:
        """
        python _pipeline/scripts/split_and_label_all_reads_include_anchor.py {BASE} {ANCHOR} {RESULTS}
        for f in {output.fasta} {output.tsv}; do touch -a "$f"; done
        """

# ---------------------------------------------------------------------------
# Step 4a: Concatenate per-chr-end reads
# ---------------------------------------------------------------------------

rule aggregate_chromosomes:
    input: expand(f"{CHR_DIR}/{BASE}_blasted_{ANCHOR}_{{ce}}_anchor_reads.fasta", ce=CHROM_ENDS)
    output: f"{RESULTS}/{BASE}_all_chr_anchored_reads.fasta"
    threads: 1
    shell: "cat {input} > {output}"

if ADAPTER_CHECK:

    # ---------------------------------------------------------------------------
    # Step 4b: Porechop adapter trimming
    # ---------------------------------------------------------------------------

    rule porechop_trim:
        input: f"{RESULTS}/{BASE}_all_chr_anchored_reads.fasta"
        output:
            fastq = f"{PCHOP_DIR}/{BASE}_trimmed.fastq",
            fasta = f"{PCHOP_DIR}/{BASE}_trimmed.fasta",
            log   = f"{PCHOP_DIR}/{BASE}_porechop.log",
        threads: workflow.cores
        shell:
            """
            mkdir -p {PCHOP_DIR}
            porechop_abi -t {threads} --format fastq -ddb -v 3 --no_split \
                -cap {ADAPTER_FILE} \
                -i {input} -o {output.fastq} > {output.log}
            sed -n '1~4s/^@/>/p;2~4p' {output.fastq} > {output.fasta}
            """

    # ---------------------------------------------------------------------------
    # Step 4c: Adapter check + telomere trimming
    # ---------------------------------------------------------------------------

    rule parse_porechop_and_trim:
        input:
            log   = f"{PCHOP_DIR}/{BASE}_porechop.log",
            fasta = f"{RESULTS}/{BASE}.fasta",
        output:
            adapter_tsv = f"{RESULTS}/{BASE}_adapter_trimming_check.tsv",
            telo_tsv    = f"{RESULTS}/{BASE}_post_telo_trimming.tsv",
        threads: 1
        shell:
            """
            python _pipeline/scripts/check_for_adapters.py {BASE} {input.log} {RESULTS}
            python _pipeline/scripts/compare_adapter_callers_dorado.py {BASE} {ANCHOR} {RESULTS} {input.fasta}
            python _pipeline/scripts/fine_telomere_trimming.py {BASE} {ANCHOR} {RESULTS}
            """

else:

    # -----------------------------------------------------------------------
    # Step 4 (adapter check off): no porechop; every anchored read goes to
    # telomere trimming, and its adapter columns are left blank (not checked)
    # -----------------------------------------------------------------------

    rule skip_adapter_check_and_trim:
        input:
            anchored = f"{RESULTS}/{BASE}_all_chr_anchored_reads.fasta",
            top      = f"{BLAST_DIR}/top_matches_{BASE}_blasted_{ANCHOR}.tsv",
        output:
            adapter_tsv = f"{RESULTS}/{BASE}_adapter_trimming_check.tsv",
            telo_tsv    = f"{RESULTS}/{BASE}_post_telo_trimming.tsv",
        threads: 1
        shell:
            """
            python _pipeline/scripts/no_adapter_check.py {BASE} {ANCHOR} {RESULTS} {input.anchored}
            python _pipeline/scripts/fine_telomere_trimming.py {BASE} {ANCHOR} {RESULTS}
            """

# ---------------------------------------------------------------------------
# Step 5: BLAST anchored reads for Y' probe (per chr end)
# ---------------------------------------------------------------------------

rule blast_y_primes:
    input: f"{CHR_DIR}/{BASE}_blasted_{ANCHOR}_{{chr_end}}_anchor_reads.fasta"
    output: f"{YP_BLAST}/{BASE}_{{chr_end}}_blasted_probe.tsv"
    threads: workflow.cores
    shell:
        """
        mkdir -p {YP_BLAST}
        echo -e "read_id\\ttotal_read_length\\tread_bp_used_for_match\\tmatch_start_on_read\\tmatch_end_on_read\\tanchor_name\\ttotal_anchor_length\\tmatch_start_on_anchor\\tmatch_end_on_anchor\\tpident\\tbitscore\\tevalue" > {output}
        # An end with no anchored reads has an empty FASTA: keep the header-only table. Any
        # other blastn failure fails the rule rather than reading as zero Y' copies.
        if [ -s {input} ]; then
            blastn -query {input} -db {PROBE_DB} \
                -perc_identity 90 -num_threads {threads} \
                -outfmt "6 qseqid qlen length qstart qend sseqid slen sstart send pident bitscore evalue" \
                >> {output}
        fi
        """

# ---------------------------------------------------------------------------
# Step 6: Y prime analysis + plots
# ---------------------------------------------------------------------------

rule y_prime_analysis:
    input:
        probes = expand(f"{YP_BLAST}/{BASE}_{{ce}}_blasted_probe.tsv", ce=CHROM_ENDS),
        telo   = f"{RESULTS}/{BASE}_post_telo_trimming.tsv",
    output: f"{RESULTS}/{BASE}_post_y_prime_probe.tsv"
    params:
        # expected Y' copies per end (zero line of the delta plots); empty = none, absolute counts.
        # --opt=value form: :q renders an empty value as nothing, which would drop the argument.
        expected = config.get("yprime_reference_path", ""),
        label    = config.get("yprime_reference_label", ""),
    threads: 1
    shell:
        """
        python _pipeline/scripts/y_prime_analysis.py {BASE} {RESULTS} \
            --expected={params.expected:q} --expected-label={params.label:q} {ADAPTER_FILTER_ARG}
        """

rule single_sample_plots:
    input: f"{RESULTS}/{BASE}_post_y_prime_probe.tsv"
    output: touch(f"{GRAPH_DIR}/.plots_done")
    threads: 1
    shell:
        """
        mkdir -p {GRAPH_DIR}
        python _pipeline/scripts/single_sample_plots.py {BASE} {RESULTS} {ADAPTER_FILTER_ARG}
        """

# Per-end anchored vs qualifying counts. Qualifying uses the same filter as the
# 500bp/facet plots, so this file states in numbers what those graphs show.
rule read_summary:
    input: f"{RESULTS}/{BASE}_post_y_prime_probe.tsv"
    output: f"{RESULTS}/{BASE}_read_summary.tsv"
    threads: 1
    shell: "python _pipeline/scripts/read_summary.py {BASE} {RESULTS} {ADAPTER_FILTER_ARG}"

# The settings this sample was profiled with (shown in the sample report).
rule profile_settings:
    output: f"{RESULTS}/{BASE}_profile_settings.tsv"
    params: settings = PROFILE_SETTINGS
    threads: 1
    run:
        with open(output[0], "w") as fh:
            fh.write("setting\tvalue\n")
            for key, value in params.settings:
                fh.write(f"{key}\t{value}\n")
