# ===========================================================================
# RECOMBINATION ANALYSIS (Steps 9.5-13)
# ===========================================================================
# Included by _pipeline/Snakefile only when recombination is on. Every rule here works against the
# day-0 sample: its reference, its labels and the Y' / spacer / X libraries label_regions builds.

from prerequisites import missing_day0_inputs

if not config.get("strain"):
    print("ERROR: recombination is on but the config has no strain (it names the day-0 files)", file=sys.stderr)
    sys.exit(1)
STRAIN = config["strain"]

# Day-0 reference + labels.  For the create_ref / label_regions steps these
# default to results/<base_name>/.  For the recombination step on a later
# time point, run_pipeline.py rewrites base_name to the time-point sample but
# writes day0_ref / day0_bed pointing back at the day-0 directory.
DAY0_REF  = config["references"].get(
    "day0_ref",
    f"results/{BASE}/_pipeline/assembly_{STRAIN}/assembly_{STRAIN}_dorado_reference.fasta",
)
DAY0_BED  = config["references"].get(
    "day0_bed",
    f"results/{BASE}/_pipeline/pretelomeric_labels/pretelomeric_regions_{STRAIN}_simp.bed",
)

# Use curated Y prime library override if specified, otherwise use extracted version
Y_PRIME_LIB    = config["references"].get("y_prime_lib_override",
                    config["references"]["y_prime_lib"].replace("{strain}", STRAIN))
SPACER_LIB_DIR = config["references"]["spacer_lib_dir"].replace("{strain}", STRAIN)
X_ELEM_LIB     = config["references"]["x_element_lib"].replace("{strain}", STRAIN)

_missing = missing_day0_inputs(DAY0_REF, DAY0_BED, Y_PRIME_LIB, SPACER_LIB_DIR, X_ELEM_LIB)
if _missing:
    print("ERROR: recombination needs these day-0 files, which are missing:", file=sys.stderr)
    for _what, _path, _hint in _missing:
        print(f"  {_what}: {_path}  ({_hint})", file=sys.stderr)
    sys.exit(1)

# ---------------------------------------------------------------------------
# Step 9.5: Filter reads to telomere-confirmed only
# ---------------------------------------------------------------------------

rule filter_telomere_reads:
    input:
        telo_tsv  = f"{RESULTS}/{BASE}_post_telo_trimming.tsv",
        chr_reads = expand(f"{CHR_DIR}/{BASE}_blasted_{ANCHOR}_{{ce}}_anchor_reads.fasta", ce=CHROM_ENDS),
    output: expand(f"{TELO_DIR}/{BASE}_{{ce}}_telomere_reads.fasta", ce=CHROM_ENDS)
    threads: 1
    shell:
        """
        python _pipeline/scripts/filter_telomere_reads.py \
            --telo-tsv    {input.telo_tsv} \
            --input-dir   {CHR_DIR} \
            --output-dir  {TELO_DIR} \
            --base-name   {BASE} \
            --anchor-set  {ANCHOR}
        """

# ---------------------------------------------------------------------------
# Step 10: minimap2 alignment (supplementary evidence only)
# ---------------------------------------------------------------------------

rule recombination_alignment:
    input:
        reads    = f"{TELO_DIR}/{BASE}_{{chr_end}}_telomere_reads.fasta",
        day0_ref = DAY0_REF,
    output:
        tsv = f"{RECOMB_DIR}/{BASE}_{{chr_end}}_alignment.tsv",
        bam = f"{RECOMB_DIR}/{BASE}_{{chr_end}}.bam",
    threads: workflow.cores
    shell:
        """
        mkdir -p {RECOMB_DIR}
        python _pipeline/scripts/run_alignment.py \
            --reads-fasta {input.reads} \
            --day0-ref    {input.day0_ref} \
            --chr-end     {wildcards.chr_end} \
            --output-tsv  {output.tsv} \
            --output-bam  {output.bam} \
            --threads     {threads}
        """

# ---------------------------------------------------------------------------
# Step 11: Combined feature analysis (per chr end)
# ---------------------------------------------------------------------------

rule recombination_analyze:
    input:
        reads       = f"{TELO_DIR}/{BASE}_{{chr_end}}_telomere_reads.fasta",
        alignment   = f"{RECOMB_DIR}/{BASE}_{{chr_end}}_alignment.tsv",
        anchor_tsv  = f"{CHR_DIR}/{BASE}_blasted_{ANCHOR}_{{chr_end}}_anchor_reads.tsv",
        day0_bed    = DAY0_BED,
        day0_ref    = DAY0_REF,
        y_prime_lib = Y_PRIME_LIB,
        telo_tsv    = f"{RESULTS}/{BASE}_post_telo_trimming.tsv",
        probe_tsv   = f"{RESULTS}/{BASE}_post_y_prime_probe.tsv",
    output:
        tsv = f"{RECOMB_DIR}/{BASE}_{{chr_end}}_features.tsv",
    params:
        chr_end           = "{chr_end}",
        strain            = STRAIN,
        spacer_lib_dir    = SPACER_LIB_DIR,
        x_element_lib     = X_ELEM_LIB,
        attribution_mode  = config.get("attribution_mode", "v2"),
        y_prime_id_level  = config.get("y_prime_id_level", "family"),
    threads: workflow.cores
    shell:
        """
        python _pipeline/scripts/analyze_features.py \
            --reads-fasta       {input.reads} \
            --alignment-tsv     {input.alignment} \
            --anchor-tsv        {input.anchor_tsv} \
            --day0-bed          {input.day0_bed} \
            --day0-ref          {input.day0_ref} \
            --y-prime-lib       {input.y_prime_lib} \
            --spacer-lib-dir    {params.spacer_lib_dir} \
            --x-element-lib     {params.x_element_lib} \
            --telo-tsv          {input.telo_tsv} \
            --probe-tsv         {input.probe_tsv} \
            --attribution-mode  {params.attribution_mode} \
            --y-prime-id-level  {params.y_prime_id_level} \
            --chr-end           {params.chr_end} \
            --strain            {params.strain} \
            --output-tsv        {output.tsv} \
            --threads           {threads}
        """

# ---------------------------------------------------------------------------
# Step 12: Summary across all chr ends
# ---------------------------------------------------------------------------

rule recombination_summary:
    input: expand(f"{RECOMB_DIR}/{BASE}_{{ce}}_features.tsv", ce=CHROM_ENDS)
    output: f"{RECOMB_DIR}/{BASE}_recombination_summary.tsv"
    params:
        recomb_dir = RECOMB_DIR,
        base_name  = BASE,
    threads: 1
    shell:
        """
        python _pipeline/scripts/aggregate_recombination.py \
            --recombination-dir {params.recomb_dir} \
            --base-name         {params.base_name} \
            --output-summary    {output}
        """

# ---------------------------------------------------------------------------
# Step 13a: Filter event-only TSVs (slim + full) for downstream investigation
# ---------------------------------------------------------------------------

rule extract_recombination_events:
    input: expand(f"{RECOMB_DIR}/{BASE}_{{ce}}_features.tsv", ce=CHROM_ENDS)
    output:
        marker = f"{RESULTS}/recombination_events/.events_extracted",
    params:
        recomb_dir = RECOMB_DIR,
        base_name  = BASE,
        events_dir = f"{RESULTS}/recombination_events",
    threads: 1
    shell:
        """
        mkdir -p {params.events_dir}
        python _pipeline/scripts/extract_recombination_events.py \
            {params.recomb_dir} {params.base_name} {params.events_dir}
        touch {output.marker}
        """

# ---------------------------------------------------------------------------
# Step 13b: Per-chr_end track plots (only chr_ends with at least one event)
# ---------------------------------------------------------------------------

# "Onion skin": how each gained Y' array was assembled from donor pieces -- circles (a donor
# piece copied in tandem) and the same donor contributing several pieces. Aggregates the
# y_prime_path / y_prime_path_circles columns analyze_features.py already writes, per end,
# collapses reads into distinct gained arrays (the events table, also the input of the
# across-timepoint step run_pipeline.py runs after the last sample), and draws a
# not-to-scale schematic of its onion-skin and possible onion-skin reads, one PNG per end
# that has any.
rule onion_skin:
    input: expand(f"{RECOMB_DIR}/{BASE}_{{ce}}_features.tsv", ce=CHROM_ENDS)
    output:
        summary = f"{RESULTS}/recombination_events/{BASE}_onion_skin_summary.tsv",
        events  = f"{RESULTS}/recombination_events/{BASE}_onion_skin_events.tsv",
        marker  = f"{GRAPH_DIR}/onion_skin/.plots_done",
    params:
        recomb_dir = RECOMB_DIR,
        base_name  = BASE,
        plots_dir  = f"{GRAPH_DIR}/onion_skin",
        read_ids   = f"{GRAPH_DIR}/onion_skin/{BASE}_gain_like_read_ids.txt",
    threads: 1
    shell:
        """
        mkdir -p {params.plots_dir}
        python _pipeline/scripts/onion_skin_summary.py \
            {params.base_name} {params.recomb_dir} {output.summary} {params.read_ids} \
            --events-out {output.events}
        python _pipeline/scripts/plot_yprime_copies.py \
            {params.recomb_dir} {params.base_name} {params.plots_dir} \
            --schematic --reads {params.read_ids} --tiers {output.events} --onion-only
        touch {output.marker}
        """

# Recombinant Y' variants: copies built from pieces of two or more day-0 Y's (e.g. the 5' part
# of one element and the rest of another). RepeatMasker names such a copy after one parent; this
# step re-reads every Y' copy against all of the day-0 Y's and reports the recombinants that recur
# (consensus-confirmed). Report-only: no existing output changes.
rule yprime_variants:
    input:
        features = expand(f"{RECOMB_DIR}/{BASE}_{{ce}}_features.tsv", ce=CHROM_ENDS),
        reads    = expand(f"{TELO_DIR}/{BASE}_{{ce}}_telomere_reads.fasta", ce=CHROM_ENDS),
        day0_ref = DAY0_REF,
        day0_bed = DAY0_BED,
    output:
        tsv    = f"{RESULTS}/recombination_events/{BASE}_yprime_variants.tsv",
        copies = f"{RESULTS}/recombination_events/{BASE}_yprime_variant_copies.tsv.gz",
        fasta  = f"{RESULTS}/recombination_events/{BASE}_yprime_variants.fasta",
    params:
        recomb_dir = RECOMB_DIR,
        telo_dir   = TELO_DIR,
        base_name  = BASE,
        out_dir    = f"{RESULTS}/recombination_events",
    threads: workflow.cores
    shell:
        """
        python _pipeline/scripts/yprime_variants.py \
            {params.base_name} {params.recomb_dir} {params.telo_dir} \
            --day0-ref {input.day0_ref} --day0-bed {input.day0_bed} \
            --out-dir {params.out_dir} --threads {threads}
        """

rule recombination_track_plots:
    input: expand(f"{RECOMB_DIR}/{BASE}_{{ce}}_features.tsv", ce=CHROM_ENDS)
    output:
        marker = f"{GRAPH_DIR}/recombination_tracks/.plots_done",
    params:
        recomb_dir = RECOMB_DIR,
        base_name  = BASE,
        plots_dir  = f"{GRAPH_DIR}/recombination_tracks",
    threads: 1
    shell:
        """
        mkdir -p {params.plots_dir}
        python _pipeline/scripts/plot_recombination_tracks.py \
            {params.recomb_dir} {params.base_name} {params.plots_dir}
        """
