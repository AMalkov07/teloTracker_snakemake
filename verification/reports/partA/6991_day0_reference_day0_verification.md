# Day-0 reference verification: 6991_day0_reference (strain 6991)

- curated bed: `verification/curated_refs/6991_features/6991_final_features.bed`
- curated lib: `verification/curated_refs/6991_features/repeatmasker_6991_all_y_primes.fasta`
- our bed: `verification/snapshot/6991_day0_reference/pretelomeric_labels/pretelomeric_regions_6991_day0_reference_simp.bed`
- our lib: `verification/snapshot/6991_day0_reference/pretelomeric_labels/extracted_yprimes_6991_day0_reference.fasta`

## Verdicts

| check | result | detail |
|---|---|---|
| A1 per-end Y' counts | **FAIL** | 1 end(s) differ |
| A2 feature coordinates | **FAIL** | 2 major, 3 minor difference(s) |
| A3 Y' grouping | **FAIL** | ARI vs curated variants=0.521, vs curated families=0.582 over 33 elements (8 our clusters vs 10 curated variants / 6 families); 6 concordant, 20 split_in_curated, 2 merged_in_curated, 5 CONFLICT |
| A4 library provenance | **PASS** | 0 failing check(s) |

## A1 -- per-end Y' counts (mismatches only)

| chr_end | n_curated | n_ours | n_probe_ours | diff_ours_minus_curated | verdict |
|---|---|---|---|---|---|
| chr12R | 7 | 6 | 6 | -1 | MISMATCH |


Full table: `6991_day0_reference_counts.tsv`

## A2 -- features outside tolerance or missing

| chr_end | feature | ftype | length_curated | length_ours | dlen | offset_curated | offset_ours | doffset | verdict | severity |
|---|---|---|---|---|---|---|---|---|---|---|
| chr5L | ITS_Y_Prime_0-1 | its | 25 |  |  | 747 |  |  | missing_in_ours | minor |
| chr12R | ITS_Y_Prime_1-2 | its | 193 | 160.0 | -33.0 | 12461 | 12463.0 | 2.0 | length_diff | minor |
| chr12R | ITS_Y_Prime_3-4 | its | 205 | 172.0 | -33.0 | 26097 | 26099.0 | 2.0 | length_diff | minor |
| chr12R | ITS_Y_Prime_6-7 | its | 170 |  |  | 46572 |  |  | missing_in_ours | major |
| chr12R | Y_Prime_7 | y_prime | 6654 |  |  | 46742 |  |  | missing_in_ours | major |


## A3 -- Y' grouping concordance

Adjusted Rand Index vs curated **variants** (ID + colour shade) = **0.521**; vs curated **families** (ID number only) = **0.582** (1.0 = identical partition; pass >= 0.85 at the variant level).

Verdicts: `split_in_curated` = curated is finer (our cluster = union of several curated variants); `merged_in_curated` = ours is finer; `CONFLICT` = the element is grouped with different partners in the two schemes.

### Contingency (rows = our ID, cols = curated variant)

| our_id | ID1_Gray | ID2_Red-Dark | ID2_Red-Light | ID3_Orange | ID4_Green-Light | ID5_Blue-Dark | ID5_Blue-Light | ID6_Purple-Dark | ID6_Purple-Light | ID6_Purple-Neutral |
|---|---|---|---|---|---|---|---|---|---|---|
| ID1 | 1 | 1 | 13 | 0 | 0 | 2 | 1 | 0 | 0 | 0 |
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
| ID1 | ID2_Red-Light | 13 | 18 | 0.722 | ID2_Red-Light:13,ID5_Blue-Dark:2,ID2_Red-Dark:1,ID1_Gray:1,ID5_Blue-Light:1 |
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
| chr4R | 1 | chr4R_Y_Prime_1 | ID1 | Long | 18.0 | ID2_Red-Light | ID2 | Long | 13.0 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 2 | chr4R_Y_Prime_2 | ID1 | Long | 18.0 | ID2_Red-Light | ID2 | Long | 13.0 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 3 | chr4R_Y_Prime_3 | ID1 | Long | 18.0 | ID2_Red-Light | ID2 | Long | 13.0 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 4 | chr4R_Y_Prime_4 | ID1 | Long | 18.0 | ID2_Red-Light | ID2 | Long | 13.0 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 5 | chr4R_Y_Prime_5 | ID1 | Long | 18.0 | ID2_Red-Light | ID2 | Long | 13.0 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 6 | chr4R_Y_Prime_6 | ID1 | Long | 18.0 | ID2_Red-Light | ID2 | Long | 13.0 | True | ID2_Red-Light | True | split_in_curated |
| chr4R | 7 | chr4R_Y_Prime_7 | ID1 | Long | 18.0 | ID2_Red-Light | ID2 | Long | 13.0 | True | ID2_Red-Light | True | split_in_curated |
| chr5R | 1 | chr5R_Y_Prime_1 | ID1 | Long | 18.0 | ID1_Gray | ID1 | Long | 7.0 | True | ID1_Gray | True | CONFLICT |
| chr7R | 1 | chr7R_Y_Prime_1 | ID1 | Long | 18.0 | ID5_Blue-Dark | ID5 | Long | 2.0 | True | ID5_Blue-Dark | True | split_in_curated |
| chr8L | 1 | chr8L_Y_Prime_1 | ID2 | Short | 7.0 | ID1_Gray | ID1 | Short | 7.0 | True | ID1_Gray | True | CONFLICT |
| chr8R | 1 | chr8R_Y_Prime_1 | ID2 | Short | 7.0 | ID1_Gray | ID1 | Short | 7.0 | True | ID1_Gray | True | CONFLICT |
| chr12L | 1 | chr12L_Y_Prime_1 | ID2 | Short | 7.0 | ID1_Gray | ID1 | Short | 7.0 | True | ID1_Gray | True | CONFLICT |
| chr12R | 1 | chr12R_Y_Prime_1 | ID5 | Long | 1.0 | ID1_Gray | ID1 | Long | 7.0 | True | ID1_Gray | True | merged_in_curated |
| chr12R | 2 | chr12R_Y_Prime_2 | ID1 | Long | 18.0 | ID2_Red-Light | ID2 | Long | 13.0 | True | ID2_Red-Light | True | split_in_curated |
| chr12R | 3 | chr12R_Y_Prime_3 | ID1 | Long | 18.0 | ID2_Red-Light | ID2 | Long | 13.0 | True | ID2_Red-Light | True | split_in_curated |
| chr12R | 4 | chr12R_Y_Prime_4 | ID1 | Long | 18.0 | ID2_Red-Light | ID2 | Long | 13.0 | True | ID2_Red-Light | True | split_in_curated |
| chr12R | 5 | chr12R_Y_Prime_5 | ID1 | Long | 18.0 | ID2_Red-Light | ID2 | Long | 13.0 | True | ID2_Red-Light | True | split_in_curated |
| chr12R | 6 | chr12R_Y_Prime_6 | ID1 | Long | 18.0 | ID2_Red-Light | ID2 | Long | 13.0 | True | ID2_Red-Light | True | split_in_curated |
| chr13L | 1 | chr13L_Y_Prime_1 | ID2 | Short | 7.0 | ID1_Gray | ID1 | Short | 7.0 | True | ID1_Gray | True | CONFLICT |
| chr14L | 1 | chr14L_Y_Prime_1 | ID1 | Long | 18.0 | ID5_Blue-Light | ID5 | Long | 1.0 | True | ID5_Blue-Light | True | split_in_curated |
| chr14L | 2 | chr14L_Y_Prime_2 | ID1 | Long | 18.0 | ID2_Red-Light | ID2 | Long | 13.0 | True | ID2_Red-Light | True | split_in_curated |
| chr14L | 3 | chr14L_Y_Prime_3 | ID2 | Short | 7.0 | ID3_Orange | ID3 | Short | 3.0 | True | ID3_Orange | True | split_in_curated |
| chr14L | 4 | chr14L_Y_Prime_4 | ID2 | Short | 7.0 | ID3_Orange | ID3 | Short | 3.0 | True | ID3_Orange | True | split_in_curated |
| chr14L | 5 | chr14L_Y_Prime_5 | ID2 | Short | 7.0 | ID3_Orange | ID3 | Short | 3.0 | True | ID3_Orange | True | split_in_curated |
| chr15R | 1 | chr15R_Y_Prime_1 | ID1 | Long | 18.0 | ID2_Red-Dark | ID2 | Long | 1.0 | True | ID2_Red-Dark | True | split_in_curated |
| chr16L | 1 | chr16L_Y_Prime_1 | ID1 | Long | 18.0 | ID5_Blue-Dark | ID5 | Long | 2.0 | True | ID5_Blue-Dark | True | split_in_curated |
| chr16R | 1 | chr16R_Y_Prime_1 | ID7 | Short | 1.0 | ID1_Gray | ID1 | Short | 7.0 | True | ID1_Gray | True | merged_in_curated |
| chr12R | 7 | chr12R_Y_Prime_7 |  |  |  | ID2_Red-Light |  | Long |  |  |  |  | only_in_curated |


Inherited (pre-clustering, 6991-curated-scheme) ID equals this strain's curated ID for 100% of elements.

### Sequence confirmation (our unique variants BLASTed vs curated variants)

min pident 99.98, min coverage 100.0%; 6 variant(s) whose best hit disagrees with the cluster mapping

| our_header | our_id | our_size | best_hit | hit_curated_id | hit_size | pident | coverage_pct | bitscore | mapped_curated_id | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| Y_Prime_chr14L3,4,5#Short/Tandem/ID2_Red | ID2 | Short | Y_Prime_chr14L3,4,5#Short/Tandem/ID3_Orange | ID3_Orange | Short | 100.0 | 100.0 | 10133.0 | ID1_Gray | hit_id_differs_from_mapping |
| Y_Prime_chr14L1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr14L1#Long/Solo/ID5_Blue-Light | ID5_Blue-Light | Long | 100.0 | 100.0 | 12292.0 | ID2_Red-Light | hit_id_differs_from_mapping |
| Y_Prime_chr15R1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr15R1#Long/Solo/ID2_Red-Dark | ID2_Red-Dark | Long | 99.985 | 100.0 | 12273.0 | ID2_Red-Light | hit_id_differs_from_mapping |
| Y_Prime_chr16L1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr16L1#Long/Solo/ID5_Blue-Dark | ID5_Blue-Dark | Long | 100.0 | 100.0 | 12432.0 | ID2_Red-Light | hit_id_differs_from_mapping |
| Y_Prime_chr5R1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr5R1#Long/Solo/ID1_Gray | ID1_Gray | Long | 99.985 | 100.0 | 12357.0 | ID2_Red-Light | hit_id_differs_from_mapping |
| Y_Prime_chr7R1#Long/Solo/ID1_Gray | ID1 | Long | Y_Prime_chr7R1#Long/Solo/ID5_Blue-Dark | ID5_Blue-Dark | Long | 100.0 | 100.0 | 12290.0 | ID2_Red-Light | hit_id_differs_from_mapping |


## A4 -- library provenance

| check | value | ok |
|---|---|---|
| run_config | verification/snapshot/6991_day0_reference/run_config.yaml | True |
| config_strain | 6991_day0_reference | True |
| y_prime_lib_override_present | False | True |
| y_prime_lib_resolved | results/6991_day0_reference/_pipeline/pretelomeric_labels/extracted_yprimes_6991_day0_reference.fasta | True |
| ours_lib_md5 | 5c5412b5785e086661719b3e1279fb37 | True |
| bed_yprimes_not_in_lib | - | True |
| lib_yprimes_not_in_bed | - | True |
| n_bed_yprimes | 33 | True |
| n_lib_yprime_locations | 33 | True |
| n_lib_unique_variants | 22 | True |
