# Expected Y′ copies per chromosome end

Each `<name>.tsv` gives, for all 32 ends, the number of Y′ copies a strain carries
(`chr_end<TAB>expected_y_primes`). Profiling draws its delta-Y′ plots against these numbers:
a read's delta is its Y′ probe count minus the expected count of its end.

Pick one with `yprime_reference: <name>` in `pipeline_config.yaml`. With the default
`yprime_reference: auto`, a table is used when its name equals the config's `strain`.

The 6991, 7172 and 7302 counts are not here: they come from the curated day-0 label BEDs in
`../day0/` (`yprime_reference: 6991_curated`, or `auto` with `strain: 6991_curated`), which give the
same numbers as the old tables.

## Provenance

Copied unchanged from the per-strain tables that `_pipeline/scripts/y_prime_analysis.py` held up to
commit e7cb48b, where they were picked by matching text in the sample name. They have not been
re-checked against assemblies.

| File | Old match | Total Y′ | Note |
|---|---|---|---|
| `7021.tsv` | `7021` | 34 | |
| `7093.tsv` | `7093` | 29 | |
| `7154.tsv` | `7154` | 21 | |
| `7174.tsv` | `7174` | 36 | |
| `7250.tsv` | `7250` | 39 | |
| `7321.tsv` | `7321` | 43 | |
| `7323.tsv` | `7323` | 29 | |
| `7324.tsv` | `7324` | 56 | |
| `7372.tsv` | `7372` | 34 | identical to 6991 |
| `KRLT3_7637.tsv` | `KRLT3` or `7637` | 33 | the old match also caught KRLT36, KRLT37, KRLT39 |
| `tlc1.tsv` | `tlc1` | 9 | |

To add a strain, write a table in the same format (all 32 ends), or point `yprime_reference` at its
label BED.
