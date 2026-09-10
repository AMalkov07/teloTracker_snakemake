# Day-0 reference verification: 7302_day0_with_selection (strain 7302)

- curated bed: `verification/curated_refs/7302_features/7302_final_features.bed`
- curated lib: `verification/curated_refs/7302_features/repeatmasker_7302_all_y_primes.fasta`
- our bed: `verification/snapshot/7302_day0_with_selection/pretelomeric_labels/pretelomeric_regions_7302_day0_with_selection_simp.bed`
- our lib: `verification/snapshot/7302_day0_with_selection/pretelomeric_labels/extracted_yprimes_7302_day0_with_selection.fasta`

## Verdicts

| check | result | detail |
|---|---|---|
| A1 per-end Y' counts | **PASS** | 0 end(s) differ |
| A2 feature coordinates | **PASS** | 0 major, 3 minor difference(s) |
| A3 Y' grouping | **FAIL** | ARI vs curated variants=0.543, vs curated families=0.604 over 36 elements (8 our clusters vs 12 curated variants / 7 families); 4 concordant, 26 split_in_curated, 2 merged_in_curated, 4 CONFLICT |
| A4 library provenance | **PASS** | 0 failing check(s) |

## A1 -- per-end Y' counts (mismatches only)

_(none)_


Full table: `7302_day0_with_selection_counts.tsv`

## A2 -- features outside tolerance or missing

| chr_end | feature | ftype | length_curated | length_ours | dlen | offset_curated | offset_ours | doffset | verdict | severity |
|---|---|---|---|---|---|---|---|---|---|---|
| chr5L | ITS_Y_Prime_0-1 | its | 25 |  |  | 747 |  |  | missing_in_ours | minor |
| chr12R | ITS_Y_Prime_1-2 | its | 193 | 160.0 | -33.0 | 12460 | 12463.0 | 3.0 | length_diff | minor |
| chr13R | anchor | anchor | 5039 | 5045.0 | 6.0 | -5038 | -5044.0 | -6.0 | length_diff | minor |


## A3 -- Y' grouping concordance

Adjusted Rand Index vs curated **variants** (ID + colour shade) = **0.543**; vs curated **families** (ID number only) = **0.604** (1.0 = identical partition; pass >= 0.85 at the variant level).

Verdicts: `split_in_curated` = curated is finer (our cluster = union of several curated variants); `merged_in_curated` = ours is finer; `CONFLICT` = the element is grouped with different partners in the two schemes.

### Contingency (rows = our ID, cols = curated variant)

| our_id | ID1_Gray | ID2_Red-Dark | ID2_Red-Light | ID3_Orange | ID4_Green-Dark | ID4_Green-Light | ID5_Blue-Dark | ID5_Blue-Light | ID6_Purple-Dark | ID6_Purple-Light | ID6_Purple-Neutral | ID7_Yellow |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ID1 | 1 | 1 | 15 | 0 | 0 | 0 | 2 | 1 | 0 | 0 | 0 | 0 |
| ID2 | 3 | 0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 |
| ID3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |
| ID4 | 0 | 0 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
| ID5 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ID6 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 |
| ID7 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ID8 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |

### Cluster mapping ours -> curated variant

| our_id | curated_id | overlap | size | purity | curated_ids_seen |
|---|---|---|---|---|---|
| ID1 | ID2_Red-Light | 15 | 20 | 0.75 | ID2_Red-Light:15,ID5_Blue-Dark:2,ID1_Gray:1,ID2_Red-Dark:1,ID5_Blue-Light:1 |
| ID2 | ID3_Orange | 3 | 8 | 0.375 | ID3_Orange:3,ID1_Gray:3,ID7_Yellow:2 |
| ID3 | ID6_Purple-Dark | 2 | 2 | 1.0 | ID6_Purple-Dark:2 |
| ID4 | ID4_Green-Dark | 1 | 2 | 0.5 | ID4_Green-Dark:1,ID4_Green-Light:1 |
| ID5 | ID1_Gray | 1 | 1 | 1.0 | ID1_Gray:1 |
| ID6 | ID6_Purple-Neutral | 1 | 1 | 1.0 | ID6_Purple-Neutral:1 |
| ID7 | ID1_Gray | 1 | 1 | 1.0 | ID1_Gray:1 |
| ID8 | ID6_Purple-Light | 1 | 1 | 1.0 | ID6_Purple-Light:1 |

### Curated variants that our clustering merged (pairwise identity of representatives)

| our_id | curated_variant_a | curated_variant_b | pident | coverage_pct | size_a | size_b |
|---|---|---|---|---|---|---|
| ID4 | ID4_Green-Dark | ID4_Green-Light | 99.048 | 100.2 | Short | Short |
| ID1 | ID1_Gray | ID2_Red-Dark | 98.287 | 99.4 | Long | Long |
| ID1 | ID1_Gray | ID2_Red-Light | 98.377 | 99.4 | Long | Long |
| ID1 | ID1_Gray | ID5_Blue-Dark | 98.415 | 100.8 | Long | Long |
| ID1 | ID1_Gray | ID5_Blue-Light | 98.753 | 99.4 | Long | Long |
| ID1 | ID2_Red-Dark | ID2_Red-Light | 99.955 | 100.0 | Long | Long |
| ID1 | ID2_Red-Dark | ID5_Blue-Dark | 99.219 | 100.1 | Long | Long |
| ID1 | ID2_Red-Dark | ID5_Blue-Light | 99.249 | 100.1 | Long | Long |
| ID1 | ID2_Red-Light | ID5_Blue-Dark | 99.309 | 100.1 | Long | Long |
| ID1 | ID2_Red-Light | ID5_Blue-Light | 99.279 | 100.1 | Long | Long |
| ID1 | ID5_Blue-Dark | ID5_Blue-Light | 99.97 | 100.0 | Long | Long |
| ID2 | ID1_Gray | ID3_Orange | 98.48 | 99.5 | Short | Short |
| ID2 | ID1_Gray | ID7_Yellow | 98.554 | 99.6 | Short | Short |
| ID2 | ID3_Orange | ID7_Yellow | 99.599 | 99.9 | Short | Short |

### Elements not concordant

| chr_end | pos | element | our_id | our_size | our_group_n | curated_variant | curated_family | curated_size | curated_group_n | size_class_agrees | inherited_variant_preclustering | inherited_matches_curated | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chr2L | 1 | chr2L_Y_Prime_1 | ID4 | Short | 2 | ID4_Green-Light | ID4 | Short | 1 | True | ID4_Green-Light | True | split_in_curated |
| chr4R | 1 | chr4R_Y_Prime_1 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 2 | chr4R_Y_Prime_2 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 3 | chr4R_Y_Prime_3 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 4 | chr4R_Y_Prime_4 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 5 | chr4R_Y_Prime_5 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 6 | chr4R_Y_Prime_6 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 7 | chr4R_Y_Prime_7 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Light | True | split_in_curated |
| chr5R | 1 | chr5R_Y_Prime_1 | ID1 | Long | 20 | ID1_Gray | ID1 | Long | 6 | True | ID1_Gray | True | CONFLICT |
| chr6L | 1 | chr6L_Y_Prime_1 | ID4 | Short | 2 | ID4_Green-Dark | ID4 | Short | 1 | True | ID4_Green-Light | False | split_in_curated |
| chr7R | 1 | chr7R_Y_Prime_1 | ID1 | Long | 20 | ID5_Blue-Dark | ID5 | Long | 2 | True | ID5_Blue-Dark | True | split_in_curated |
| chr8L | 1 | chr8L_Y_Prime_1 | ID2 | Short | 8 | ID1_Gray | ID1 | Short | 6 | True | ID1_Gray | True | CONFLICT |
| chr8R | 1 | chr8R_Y_Prime_1 | ID2 | Short | 8 | ID1_Gray | ID1 | Short | 6 | True | ID1_Gray | True | CONFLICT |
| chr12L | 1 | chr12L_Y_Prime_1 | ID2 | Short | 8 | ID1_Gray | ID1 | Short | 6 | True | ID1_Gray | True | CONFLICT |
| chr12R | 1 | chr12R_Y_Prime_1 | ID5 | Long | 1 | ID1_Gray | ID1 | Long | 6 | True | ID1_Gray | True | merged_in_curated |
| chr12R | 2 | chr12R_Y_Prime_2 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Light | True | split_in_curated |
| chr12R | 3 | chr12R_Y_Prime_3 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Light | True | split_in_curated |
| chr12R | 4 | chr12R_Y_Prime_4 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Light | True | split_in_curated |
| chr12R | 5 | chr12R_Y_Prime_5 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Light | True | split_in_curated |
| chr12R | 6 | chr12R_Y_Prime_6 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Dark | False | split_in_curated |
| chr13L | 1 | chr13L_Y_Prime_1 | ID2 | Short | 8 | ID7_Yellow | ID7 | Short | 2 | True | ID1_Gray | False | split_in_curated |
| chr13L | 2 | chr13L_Y_Prime_2 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Dark | False | split_in_curated |
| chr13L | 3 | chr13L_Y_Prime_3 | ID2 | Short | 8 | ID7_Yellow | ID7 | Short | 2 | True | ID1_Gray | False | split_in_curated |
| chr13L | 4 | chr13L_Y_Prime_4 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Dark | False | split_in_curated |
| chr14L | 1 | chr14L_Y_Prime_1 | ID1 | Long | 20 | ID5_Blue-Light | ID5 | Long | 1 | True | ID5_Blue-Light | True | split_in_curated |
| chr14L | 2 | chr14L_Y_Prime_2 | ID1 | Long | 20 | ID2_Red-Light | ID2 | Long | 15 | True | ID2_Red-Light | True | split_in_curated |
| chr14L | 3 | chr14L_Y_Prime_3 | ID2 | Short | 8 | ID3_Orange | ID3 | Short | 3 | True | ID3_Orange | True | split_in_curated |
| chr14L | 4 | chr14L_Y_Prime_4 | ID2 | Short | 8 | ID3_Orange | ID3 | Short | 3 | True | ID3_Orange | True | split_in_curated |
| chr14L | 5 | chr14L_Y_Prime_5 | ID2 | Short | 8 | ID3_Orange | ID3 | Short | 3 | True | ID3_Orange | True | split_in_curated |
| chr15R | 1 | chr15R_Y_Prime_1 | ID1 | Long | 20 | ID2_Red-Dark | ID2 | Long | 1 | True | ID2_Red-Dark | True | split_in_curated |
| chr16L | 1 | chr16L_Y_Prime_1 | ID1 | Long | 20 | ID5_Blue-Dark | ID5 | Long | 2 | True | ID5_Blue-Dark | True | split_in_curated |
| chr16R | 1 | chr16R_Y_Prime_1 | ID7 | Short | 1 | ID1_Gray | ID1 | Short | 6 | True | ID1_Gray | True | merged_in_curated |


Inherited (pre-clustering, 6991-curated-scheme) ID equals this strain's curated ID for 83% of elements.

### Sequence confirmation (our unique variants BLASTed vs curated variants)

min pident 99.98, min coverage 99.4%; 10 variant(s) whose best hit disagrees with the cluster mapping

| our_header | our_id | our_size | best_hit | hit_curated_id | hit_size | pident | coverage_pct | bitscore | mapped_curated_id | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| Y_Prime_chr13L1,3#Short/Tandem/ID2_Red | ID2 | Short | Y_Prime_chr13L3#Short/Solo/ID7_Yellow | ID7_Yellow | Short | 100.0 | 100.0 | 10129.0 | ID3_Orange | hit_id_differs_from_mapping |
| Y_Prime_chr12L1#Short/Solo/ID2_Red | ID2 | Short | Y_Prime_chr12L1#Short/Solo/ID1_Gray | ID1_Gray | Short | 100.0 | 100.0 | 9588.0 | ID3_Orange | hit_id_differs_from_mapping |
| Y_Prime_chr14L1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr14L1#Long/Solo/ID5_Blue-Light | ID5_Blue-Light | Long | 100.0 | 100.0 | 12292.0 | ID2_Red-Light | hit_id_differs_from_mapping |
| Y_Prime_chr15R1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr15R1#Long/Solo/ID2_Red-Dark | ID2_Red-Dark | Long | 99.985 | 100.0 | 12273.0 | ID2_Red-Light | hit_id_differs_from_mapping |
| Y_Prime_chr16L1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr16L1#Long/Solo/ID5_Blue-Dark | ID5_Blue-Dark | Long | 100.0 | 100.0 | 12434.0 | ID2_Red-Light | hit_id_differs_from_mapping |
| Y_Prime_chr2L1#Short/Solo/ID4_Orange | ID4 | Short | Y_Prime_chr2L1#Short/Solo/ID4_Green-Light | ID4_Green-Light | Short | 100.0 | 100.0 | 11036.0 | ID4_Green-Dark | hit_id_differs_from_mapping |
| Y_Prime_chr5R1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr5R1#Long/Solo/ID1_Gray | ID1_Gray | Long | 100.0 | 100.0 | 12364.0 | ID2_Red-Light | hit_id_differs_from_mapping |
| Y_Prime_chr7R1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr7R1#Long/Solo/ID5_Blue-Dark | ID5_Blue-Dark | Long | 100.0 | 100.0 | 12290.0 | ID2_Red-Light | hit_id_differs_from_mapping |
| Y_Prime_chr8L1#Short/Solo/ID2_Red | ID2 | Short | Y_Prime_chr8L1#Short/Solo/ID1_Gray | ID1_Gray | Short | 100.0 | 100.0 | 9347.0 | ID3_Orange | hit_id_differs_from_mapping |
| Y_Prime_chr8R1#Short/Solo/ID2_Red | ID2 | Short | Y_Prime_chr8R1#Short/Solo/ID1_Gray | ID1_Gray | Short | 100.0 | 99.4 | 10039.0 | ID3_Orange | hit_id_differs_from_mapping |


## A4 -- library provenance

| check | value | ok |
|---|---|---|
| run_config | verification/snapshot/7302_day0_with_selection/run_config.yaml | True |
| config_strain | 7302_day0_with_selection | True |
| y_prime_lib_override_present | False | True |
| y_prime_lib_resolved | results/7302_day0_with_selection/_pipeline/pretelomeric_labels/extracted_yprimes_7302_day0_with_selection.fasta | True |
| ours_lib_md5 | 1814ac330402a2f0c37409b9048b89cb | True |
| bed_yprimes_not_in_lib | - | True |
| lib_yprimes_not_in_bed | - | True |
| n_bed_yprimes | 36 | True |
| n_lib_yprime_locations | 36 | True |
| n_lib_unique_variants | 23 | True |
