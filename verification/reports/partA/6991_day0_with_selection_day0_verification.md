# Day-0 reference verification: 6991_day0_with_selection (strain 6991)

- curated bed: `verification/curated_refs/6991_features/6991_final_features.bed`
- curated lib: `verification/curated_refs/6991_features/repeatmasker_6991_all_y_primes.fasta`
- our bed: `verification/snapshot/6991_day0_with_selection/pretelomeric_labels/pretelomeric_regions_6991_day0_with_selection_simp.bed`
- our lib: `verification/snapshot/6991_day0_with_selection/pretelomeric_labels/extracted_yprimes_6991_day0_with_selection.fasta`

## Verdicts

| check | result | detail |
|---|---|---|
| A1 per-end Y' counts | **PASS** | 0 end(s) differ |
| A2 feature coordinates | **FAIL** | 9 major, 3 minor difference(s) |
| A3 Y' grouping | **FAIL** | ARI vs curated variants=0.540, vs curated families=0.599 over 34 elements (8 our clusters vs 10 curated variants / 6 families); 6 concordant, 21 split_in_curated, 2 merged_in_curated, 5 CONFLICT |
| A4 library provenance | **PASS** | 0 failing check(s) |

## A1 -- per-end Y' counts (mismatches only)

_(none)_


Full table: `6991_day0_with_selection_counts.tsv`

## A2 -- features outside tolerance or missing

| chr_end | feature | ftype | length_curated | length_ours | dlen | offset_curated | offset_ours | doffset | verdict | severity |
|---|---|---|---|---|---|---|---|---|---|---|
| chr5L | ITS_Y_Prime_0-1 | its | 25 |  |  | 747 |  |  | missing_in_ours | minor |
| chr12R | ITS_Y_Prime_1-2 | its | 193 | 160.0 | -33.0 | 12461 | 12460.0 | -1.0 | length_diff | minor |
| chr12R | ITS_Y_Prime_3-4 | its | 205 | 172.0 | -33.0 | 26097 | 26096.0 | -1.0 | length_diff | minor |
| chr14L | Y_Prime_1 | y_prime | 6657 | 5721.0 | -936.0 | 10296 | 10297.0 | 1.0 | length_diff | major |
| chr14L | ITS_Y_Prime_1-2 | its | 169 | 1943.0 | 1774.0 | 16951 | 16018.0 | -933.0 | length+offset_diff | major |
| chr14L | Y_Prime_2 | y_prime | 6655 | 6653.0 | -2.0 | 17120 | 17961.0 | 841.0 | offset_diff | major |
| chr14L | ITS_Y_Prime_2-3 | its | 172 | 172.0 | 0.0 | 23773 | 24614.0 | 841.0 | offset_diff | major |
| chr14L | Y_Prime_3 | y_prime | 5488 | 5486.0 | -2.0 | 23945 | 24786.0 | 841.0 | offset_diff | major |
| chr14L | ITS_Y_Prime_3-4 | its | 175 | 175.0 | 0.0 | 29431 | 30272.0 | 841.0 | offset_diff | major |
| chr14L | Y_Prime_4 | y_prime | 5488 | 5486.0 | -2.0 | 29606 | 30447.0 | 841.0 | offset_diff | major |
| chr14L | ITS_Y_Prime_4-5 | its | 175 | 175.0 | 0.0 | 35092 | 35933.0 | 841.0 | offset_diff | major |
| chr14L | Y_Prime_5 | y_prime | 5488 | 5485.0 | -3.0 | 35267 | 36108.0 | 841.0 | offset_diff | major |


## A3 -- Y' grouping concordance

Adjusted Rand Index vs curated **variants** (ID + colour shade) = **0.540**; vs curated **families** (ID number only) = **0.599** (1.0 = identical partition; pass >= 0.85 at the variant level).

Verdicts: `split_in_curated` = curated is finer (our cluster = union of several curated variants); `merged_in_curated` = ours is finer; `CONFLICT` = the element is grouped with different partners in the two schemes.

### Contingency (rows = our ID, cols = curated variant)

| our_id | ID1_Gray | ID2_Red-Dark | ID2_Red-Light | ID3_Orange | ID4_Green-Light | ID5_Blue-Dark | ID5_Blue-Light | ID6_Purple-Dark | ID6_Purple-Light | ID6_Purple-Neutral |
|---|---|---|---|---|---|---|---|---|---|---|
| ID1 | 1 | 1 | 14 | 0 | 0 | 2 | 1 | 0 | 0 | 0 |
| ID2 | 4 | 0 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 |
| ID3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 |
| ID4 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 0 | 0 |
| ID5 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ID6 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| ID7 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ID8 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 0 |

### Cluster mapping ours -> curated variant

| our_id | curated_id | overlap | size | purity | curated_ids_seen |
|---|---|---|---|---|---|
| ID1 | ID2_Red-Light | 14 | 19 | 0.737 | ID2_Red-Light:14,ID5_Blue-Dark:2,ID1_Gray:1,ID2_Red-Dark:1,ID5_Blue-Light:1 |
| ID2 | ID1_Gray | 4 | 7 | 0.571 | ID1_Gray:4,ID3_Orange:3 |
| ID3 | ID6_Purple-Dark | 2 | 2 | 1.0 | ID6_Purple-Dark:2 |
| ID4 | ID4_Green-Light | 2 | 2 | 1.0 | ID4_Green-Light:2 |
| ID5 | ID1_Gray | 1 | 1 | 1.0 | ID1_Gray:1 |
| ID6 | ID6_Purple-Neutral | 1 | 1 | 1.0 | ID6_Purple-Neutral:1 |
| ID7 | ID1_Gray | 1 | 1 | 1.0 | ID1_Gray:1 |
| ID8 | ID6_Purple-Light | 1 | 1 | 1.0 | ID6_Purple-Light:1 |

### Curated variants that our clustering merged (pairwise identity of representatives)

| our_id | curated_variant_a | curated_variant_b | pident | coverage_pct | size_a | size_b |
|---|---|---|---|---|---|---|
| ID1 | ID1_Gray | ID2_Red-Dark | 98.302 | 99.4 | Long | Long |
| ID1 | ID1_Gray | ID2_Red-Light | 98.242 | 100.0 | Long | Long |
| ID1 | ID1_Gray | ID5_Blue-Dark | 98.414 | 100.8 | Long | Long |
| ID1 | ID1_Gray | ID5_Blue-Light | 98.738 | 99.4 | Long | Long |
| ID1 | ID2_Red-Dark | ID2_Red-Light | 99.88 | 100.0 | Long | Long |
| ID1 | ID2_Red-Dark | ID5_Blue-Dark | 99.219 | 100.1 | Long | Long |
| ID1 | ID2_Red-Dark | ID5_Blue-Light | 99.249 | 100.1 | Long | Long |
| ID1 | ID2_Red-Light | ID5_Blue-Dark | 99.175 | 100.1 | Long | Long |
| ID1 | ID2_Red-Light | ID5_Blue-Light | 99.145 | 100.1 | Long | Long |
| ID1 | ID5_Blue-Dark | ID5_Blue-Light | 99.97 | 100.0 | Long | Long |
| ID2 | ID1_Gray | ID3_Orange | 99.581 | 99.9 | Short | Short |

### Elements not concordant

| chr_end | pos | element | our_id | our_size | our_group_n | curated_variant | curated_family | curated_size | curated_group_n | size_class_agrees | inherited_variant_preclustering | inherited_matches_curated | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| chr4R | 1 | chr4R_Y_Prime_1 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 2 | chr4R_Y_Prime_2 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 3 | chr4R_Y_Prime_3 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 4 | chr4R_Y_Prime_4 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 5 | chr4R_Y_Prime_5 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 6 | chr4R_Y_Prime_6 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 7 | chr4R_Y_Prime_7 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr5R | 1 | chr5R_Y_Prime_1 | ID1 | Long | 19 | ID1_Gray | ID1 | Long | 7 | True | ID1_Gray | True | CONFLICT |
| chr7R | 1 | chr7R_Y_Prime_1 | ID1 | Long | 19 | ID5_Blue-Dark | ID5 | Long | 2 | True | ID5_Blue-Dark | True | split_in_curated |
| chr8L | 1 | chr8L_Y_Prime_1 | ID2 | Short | 7 | ID1_Gray | ID1 | Short | 7 | True | ID1_Gray | True | CONFLICT |
| chr8R | 1 | chr8R_Y_Prime_1 | ID2 | Short | 7 | ID1_Gray | ID1 | Short | 7 | True | ID1_Gray | True | CONFLICT |
| chr12L | 1 | chr12L_Y_Prime_1 | ID2 | Short | 7 | ID1_Gray | ID1 | Short | 7 | True | ID1_Gray | True | CONFLICT |
| chr12R | 1 | chr12R_Y_Prime_1 | ID5 | Long | 1 | ID1_Gray | ID1 | Long | 7 | True | ID1_Gray | True | merged_in_curated |
| chr12R | 2 | chr12R_Y_Prime_2 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr12R | 3 | chr12R_Y_Prime_3 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr12R | 4 | chr12R_Y_Prime_4 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr12R | 5 | chr12R_Y_Prime_5 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr12R | 6 | chr12R_Y_Prime_6 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr12R | 7 | chr12R_Y_Prime_7 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr13L | 1 | chr13L_Y_Prime_1 | ID2 | Short | 7 | ID1_Gray | ID1 | Short | 7 | True | ID1_Gray | True | CONFLICT |
| chr14L | 1 | chr14L_Y_Prime_1 | ID1 | Short | 19 | ID5_Blue-Light | ID5 | Long | 1 | False | ID5_Blue-Light | True | split_in_curated |
| chr14L | 2 | chr14L_Y_Prime_2 | ID1 | Long | 19 | ID2_Red-Light | ID2 | Long | 14 | True | ID2_Red-Light | True | split_in_curated |
| chr14L | 3 | chr14L_Y_Prime_3 | ID2 | Short | 7 | ID3_Orange | ID3 | Short | 3 | True | ID3_Orange | True | split_in_curated |
| chr14L | 4 | chr14L_Y_Prime_4 | ID2 | Short | 7 | ID3_Orange | ID3 | Short | 3 | True | ID3_Orange | True | split_in_curated |
| chr14L | 5 | chr14L_Y_Prime_5 | ID2 | Short | 7 | ID3_Orange | ID3 | Short | 3 | True | ID3_Orange | True | split_in_curated |
| chr15R | 1 | chr15R_Y_Prime_1 | ID1 | Long | 19 | ID2_Red-Dark | ID2 | Long | 1 | True | ID2_Red-Dark | True | split_in_curated |
| chr16L | 1 | chr16L_Y_Prime_1 | ID1 | Long | 19 | ID5_Blue-Dark | ID5 | Long | 2 | True | ID5_Blue-Dark | True | split_in_curated |
| chr16R | 1 | chr16R_Y_Prime_1 | ID7 | Short | 1 | ID1_Gray | ID1 | Short | 7 | True | ID1_Gray | True | merged_in_curated |


Inherited (pre-clustering, 6991-curated-scheme) ID equals this strain's curated ID for 100% of elements.

### Sequence confirmation (our unique variants BLASTed vs curated variants)

min pident 99.98, min coverage 100.0%; 7 variant(s) whose best hit disagrees with the cluster mapping

| our_header | our_id | our_size | best_hit | hit_curated_id | hit_size | pident | coverage_pct | bitscore | mapped_curated_id | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| Y_Prime_chr14L3,4#Short/Tandem/ID2_Red | ID2 | Short | Y_Prime_chr14L3,4,5#Short/Tandem/ID3_Orange | ID3_Orange | Short | 100.0 | 100.0 | 10133.0 | ID1_Gray | hit_id_differs_from_mapping |
| Y_Prime_chr14L1#Short/Solo/ID1_Gray | ID1 | Short | Y_Prime_chr14L1#Long/Solo/ID5_Blue-Light | ID5_Blue-Light | Long | 99.983 | 100.0 | 10560.0 | ID2_Red-Light | hit_id_differs_from_mapping |
| Y_Prime_chr14L5#Short/Solo/ID2_Red | ID2 | Short | Y_Prime_chr14L3,4,5#Short/Tandem/ID3_Orange | ID3_Orange | Short | 99.982 | 100.0 | 10126.0 | ID1_Gray | hit_id_differs_from_mapping |
| Y_Prime_chr15R1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr15R1#Long/Solo/ID2_Red-Dark | ID2_Red-Dark | Long | 100.0 | 100.0 | 12279.0 | ID2_Red-Light | hit_id_differs_from_mapping |
| Y_Prime_chr16L1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr16L1#Long/Solo/ID5_Blue-Dark | ID5_Blue-Dark | Long | 100.0 | 100.0 | 12432.0 | ID2_Red-Light | hit_id_differs_from_mapping |
| Y_Prime_chr5R1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr5R1#Long/Solo/ID1_Gray | ID1_Gray | Long | 100.0 | 100.0 | 12362.0 | ID2_Red-Light | hit_id_differs_from_mapping |
| Y_Prime_chr7R1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr16L1#Long/Solo/ID5_Blue-Dark | ID5_Blue-Dark | Long | 99.985 | 100.0 | 12286.0 | ID2_Red-Light | hit_id_differs_from_mapping |


## A4 -- library provenance

| check | value | ok |
|---|---|---|
| run_config | verification/snapshot/6991_day0_with_selection/run_config.yaml | True |
| config_strain | 6991_day0_with_selection | True |
| y_prime_lib_override_present | False | True |
| y_prime_lib_resolved | results/6991_day0_with_selection/_pipeline/pretelomeric_labels/extracted_yprimes_6991_day0_with_selection.fasta | True |
| ours_lib_md5 | 4d07d11601d0df2184a151bbb95bebe3 | True |
| bed_yprimes_not_in_lib | - | True |
| lib_yprimes_not_in_bed | - | True |
| n_bed_yprimes | 34 | True |
| n_lib_yprime_locations | 34 | True |
| n_lib_unique_variants | 23 | True |
