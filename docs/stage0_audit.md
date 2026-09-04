# Stage 0 audit: what each rung actually ships

Gene universe: **18,481** genes (protein-coding, gnomAD v2.1.1 LOEUF, GTEx brain expression).


## Rung 1 — GTEx v10 brain

| rung | cell | genes_tested | in_universe | egenes_qval | egenes_bonf | frac_of_tested | frac_bonf_of_tested | frac_of_universe |
|---|---|---|---|---|---|---|---|---|
| bulk tissue (rung 1) | Brain_Amygdala | 24991 | 16019 | 4047 | - | 0.253 | - | 0.219 |
| bulk tissue (rung 1) | Brain_Anterior_cingulate_cortex_BA24 | 25420 | 16033 | 5967 | - | 0.372 | - | 0.323 |
| bulk tissue (rung 1) | Brain_Caudate_basal_ganglia | 25768 | 16070 | 8126 | - | 0.506 | - | 0.440 |
| bulk tissue (rung 1) | Brain_Cerebellar_Hemisphere | 26170 | 15944 | 9039 | - | 0.567 | - | 0.489 |
| bulk tissue (rung 1) | Brain_Cerebellum | 26696 | 16107 | 9151 | - | 0.568 | - | 0.495 |
| bulk tissue (rung 1) | Brain_Cortex | 26134 | 16195 | 8226 | - | 0.508 | - | 0.445 |
| bulk tissue (rung 1) | Brain_Frontal_Cortex_BA9 | 25884 | 16114 | 7956 | - | 0.494 | - | 0.430 |
| bulk tissue (rung 1) | Brain_Hippocampus | 25468 | 16111 | 5877 | - | 0.365 | - | 0.318 |
| bulk tissue (rung 1) | Brain_Hypothalamus | 26212 | 16321 | 5705 | - | 0.350 | - | 0.309 |
| bulk tissue (rung 1) | Brain_Nucleus_accumbens_basal_ganglia | 25926 | 16069 | 7664 | - | 0.477 | - | 0.415 |
| bulk tissue (rung 1) | Brain_Putamen_basal_ganglia | 24915 | 15899 | 6792 | - | 0.427 | - | 0.368 |
| bulk tissue (rung 1) | Brain_Spinal_cord_cervical_c-1 | 25721 | 16227 | 5499 | - | 0.339 | - | 0.298 |
| bulk tissue (rung 1) | Brain_Substantia_nigra | 25195 | 16078 | 3920 | - | 0.244 | - | 0.212 |


## Rung 2 — PsychENCODE

- rows: 674,626
- unique genes: 8,190 (5,252 in universe)
- **supplies a denominator: no.** Significant pairs only. Every gene in this file is an eGene, so the file cannot supply the set of genes tested-but-not-significant.


## Rungs 3–4 — SingleBrain

| rung | cell | genes_tested | in_universe | egenes_qval | egenes_bonf | frac_of_tested | frac_bonf_of_tested | frac_of_universe |
|---|---|---|---|---|---|---|---|---|
| subtype (rung 4) | Ast1 | 12358 | 11310 | 8892 | 2671 | 0.786 | 0.236 | 0.481 |
| subtype (rung 4) | Ast2 | 12674 | 11601 | 8227 | 2012 | 0.709 | 0.173 | 0.445 |
| subtype (rung 4) | Ast3 | 15245 | 13461 | 4477 | 1121 | 0.333 | 0.083 | 0.242 |
| subtype (rung 4) | Ast4 | 12676 | 11619 | 10063 | 1007 | 0.866 | 0.087 | 0.545 |
| major (rung 3) | Ast | 12257 | 11257 | 9451 | 4143 | 0.840 | 0.368 | 0.511 |
| major (rung 3) | End | 15158 | 13539 | 3843 | 988 | 0.284 | 0.073 | 0.208 |
| subtype (rung 4) | Ext1 | 11893 | 10916 | 10238 | 4868 | 0.938 | 0.446 | 0.554 |
| subtype (rung 4) | Ext2 | 12048 | 11060 | 10066 | 4245 | 0.910 | 0.384 | 0.545 |
| subtype (rung 4) | Ext3 | 12567 | 11409 | 6878 | 1275 | 0.603 | 0.112 | 0.372 |
| subtype (rung 4) | Ext4 | 11988 | 11008 | 10033 | 4312 | 0.911 | 0.392 | 0.543 |
| subtype (rung 4) | Ext5 | 12099 | 11084 | 7626 | 1612 | 0.688 | 0.145 | 0.413 |
| subtype (rung 4) | Ext6 | 12309 | 11229 | 6321 | 1460 | 0.563 | 0.130 | 0.342 |
| subtype (rung 4) | Ext7 | 11887 | 10897 | 8157 | 2466 | 0.749 | 0.226 | 0.441 |
| subtype (rung 4) | Ext8 | 11867 | 10891 | 9600 | 3630 | 0.881 | 0.333 | 0.519 |
| major (rung 3) | Ext | 12000 | 11011 | 10700 | 7471 | 0.972 | 0.679 | 0.579 |
| subtype (rung 4) | IN1 | 12399 | 11335 | 6025 | 1554 | 0.532 | 0.137 | 0.326 |
| subtype (rung 4) | IN2 | 12094 | 11130 | 7644 | 1984 | 0.687 | 0.178 | 0.414 |
| subtype (rung 4) | IN3 | 13044 | 11826 | 6473 | 835 | 0.547 | 0.071 | 0.350 |
| subtype (rung 4) | IN4 | 11864 | 10929 | 9266 | 3420 | 0.848 | 0.313 | 0.501 |
| subtype (rung 4) | IN5 | 12987 | 11820 | 6081 | 1175 | 0.514 | 0.099 | 0.329 |
| subtype (rung 4) | IN6 | 12032 | 11098 | 8392 | 2385 | 0.756 | 0.215 | 0.454 |
| subtype (rung 4) | IN7 | 12224 | 11219 | 9115 | 2937 | 0.812 | 0.262 | 0.493 |
| major (rung 3) | IN | 12226 | 11227 | 10305 | 5389 | 0.918 | 0.480 | 0.558 |
| subtype (rung 4) | MG1 | 14158 | 12838 | 6138 | 1551 | 0.478 | 0.121 | 0.332 |
| subtype (rung 4) | MG2 | 14127 | 12823 | 5965 | 1308 | 0.465 | 0.102 | 0.323 |
| subtype (rung 4) | MG3 | 14651 | 13098 | 6832 | 926 | 0.522 | 0.071 | 0.370 |
| subtype (rung 4) | MG4 | 14858 | 13238 | 5966 | 791 | 0.451 | 0.060 | 0.323 |
| major (rung 3) | MG | 12320 | 11405 | 7796 | 2583 | 0.684 | 0.226 | 0.422 |
| excluded | MiGA3 | 21059 | 14287 | 13934 | 3026 | 0.975 | 0.212 | 0.754 |
| subtype (rung 4) | OD1 | 12025 | 10978 | 6705 | 1807 | 0.611 | 0.165 | 0.363 |
| subtype (rung 4) | OD2 | 11350 | 10472 | 8135 | 3222 | 0.777 | 0.308 | 0.440 |
| subtype (rung 4) | OD3 | 11939 | 10939 | 8230 | 2884 | 0.752 | 0.264 | 0.445 |
| major (rung 3) | OD | 11505 | 10593 | 9126 | 4554 | 0.862 | 0.430 | 0.494 |
| subtype (rung 4) | OPC1 | 12055 | 11044 | 8448 | 2777 | 0.765 | 0.251 | 0.457 |
| subtype (rung 4) | OPC2 | 12681 | 11539 | 9674 | 910 | 0.838 | 0.079 | 0.523 |
| major (rung 3) | OPC | 12029 | 11031 | 8569 | 3187 | 0.777 | 0.289 | 0.464 |


## Power contrast — Bryois

| file | arm | rows_chr1 | genes_chr1 | in_universe_chr1 | id_parsed_ok | min_p | has_per_gene_correction |
|---|---|---|---|---|---|---|---|
| Astrocytes.1.gz | cell type | 3676698 | 1112 | 954 | 1.000 | 0.000 | False |
| pb.1.gz | pseudobulk | 5422234 | 1624 | 1416 | 1.000 | 0.000 | False |
