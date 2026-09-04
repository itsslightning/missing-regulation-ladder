# Stage 0 report — scope lock, data audit, and what the report's assumptions got wrong

Date: 2026-09-04. Generated tables: [`stage0_audit.md`](stage0_audit.md).
Decisions: [`../DECISIONS.md`](../DECISIONS.md).

---

## 1. Headline

All four rungs are reachable, but **not on the terms the project report assumed**.
Three of its data assumptions do not hold, one of them seriously enough to change
what the recovery curve should measure. Nothing here blocks Stage 1; two things
need a decision from you before it starts.

The single most important finding is that **the outcome measure saturates at the
single-nucleus rungs unless the denominator is fixed**, and that fixing it is
what makes the whole comparison interpretable.

---

## 2. What is usable, rung by rung

| Rung | Source | Usable | Note |
|---|---|---|---|
| 1 — bulk tissue | GTEx v10, 13 brain tissues | **Yes** | Full eGenes tables: every gene tested, with permutation q-values. Best case. |
| 2 — bulk brain | PsychENCODE DER-08b | **Partly** | Significant pairs only — 8,190 genes, all of them eGenes. No denominator in the file. |
| 2 — bulk brain | MetaBrain | **No** | No public file index; access gated behind a Google Form. Cannot be automated. |
| 3 — sn major | SingleBrain, 7 cell types | **Yes** | One row per gene, q-values, CC-BY-4.0. |
| 4 — sn subtype | SingleBrain, 28 subtypes | **Yes** | Same format. |
| Power contrast | Bryois, 8 cell types + pseudobulk | **Yes** | Nominal p-values only; needs its own correction. |

Downloaded and hash-stamped so far: 47 files, 0.29 GB, all recorded in
`logs/provenance.json` with URL, UTC date and sha256.

---

## 3. Where the report's data assumptions do not hold

### 3.1 gnomAD's `brain_expression` column is empty

The report proposes matching control genes on expression level, and gnomAD
v2.1.1 has a `brain_expression` field that looks like exactly the right
covariate. **It is `NA` for all 19,704 rows.**

Resolved: expression now comes from the GTEx v10 median-TPM matrix
(`GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz`), which also yields a
breadth-of-expression covariate (number of the 13 brain tissues with TPM ≥ 1).
No timeline impact.

### 3.2 MetaBrain cannot be scripted

`download.metabrain.nl` serves no public index and returns 404 for every
guessable path; the site routes all access through a Google Form. This is a
human-in-the-loop step, not an engineering problem.

Consequence: rung 2 currently rests on PsychENCODE alone. See **D-011**.

### 3.3 PsychENCODE supplies a numerator but not a denominator

`DER-08b_hg38_eQTL.bonferroni.txt` contains 674,626 significant variant–gene
pairs across **8,190 unique genes — every one of which is an eGene by
construction**. The set of genes tested-but-not-significant is not in the file.

You cannot compute "fraction of genes with a detectable eQTL" from a file that
only contains genes with detectable eQTLs. Two ways out:

- score PsychENCODE against the fixed gene universe (what the pipeline now
  does), accepting that its denominator is assumed rather than observed; or
- pull `Full_hg19_cis-eQTL.txt.gz` (3.3 GB, hg19 only — needs a liftOver that
  the other rungs do not).

The first is implemented. The second is available if rung 2 turns out to be
load-bearing.

---

## 4. Format mismatches between datasets

| | Gene ID | Build / GENCODE | Significance shipped |
|---|---|---|---|
| gnomAD v2.1.1 | `ENSG`, unversioned | GRCh37 / v19 era | n/a |
| SCHEMA | `ENSG`, unversioned | GRCh38 | **p-values only, no q** |
| GTEx v10 | `ENSG.version` | GRCh38 / v39 | Storey q on beta-approximated permutation p |
| PsychENCODE | `ENSG.version` | hg38 (lifted) | FDR, significant rows only |
| SingleBrain | `ENSG.version` | GRCh38 / **v38** | Storey q on within-gene corrected p |
| Bryois | `SYMBOL_ENSG` compound | GRCh37 | **nominal p only** |

Bare unversioned ENSG is the only join key reaching all six. Cost of that join,
measured rather than assumed:

- constraint → expression: 19,183 → **18,481** genes (702 lost to GENCODE drift)
- universe → SCHEMA: **17,630 of 18,481** carry SCHEMA statistics (851 do not)

Two mismatches need real handling rather than a join:

- **Bryois ships no per-gene correction.** Nominal p-values across all
  SNP–gene pairs, one file per chromosome. Scoring it the way GTEx and
  SingleBrain are scored means recomputing a per-gene correction locally. This
  is part of **D-001**.
- **SCHEMA ships no q-value.** Defining "a SCHEMA gene" therefore requires
  choosing a p-value column *and* a correction. This is part of **D-009**.

One mismatch turned out **not** to be real. I initially flagged SingleBrain's
`qval` as under-corrected relative to GTEx's, which would have invalidated the
cross-rung comparison. Checking directly: SingleBrain `qval` is a rank-preserving
transform of `Random_FDR` (Spearman ρ = 1.000000), which is a within-gene
variant correction — the same two-stage structure as GTEx's Storey-on-permutation
q. **The two are comparable.** Recorded here because the opposite conclusion
would have been a serious error.

---

## 5. The finding that matters: the outcome measure saturates

Fraction of tested genes called eGenes at q ≤ 0.05, SingleBrain major cell types:

| Cell type | of genes tested | of fixed universe |
|---|---|---|
| Excitatory neurons (Ext) | **0.972** | 0.579 |
| Inhibitory neurons (IN) | 0.918 | 0.558 |
| Astrocytes (Ast) | 0.840 | 0.511 |
| Oligodendrocytes (OD) | 0.860 | 0.511 |
| Microglia (MG) | 0.684 | 0.422 |
| Endothelial (End) | 0.284 | 0.208 |
| GTEx Brain_Cortex (rung 1) | 0.508 | **0.445** |

Read the middle column and rung 3 has no headroom: 97% of the genes excitatory
neurons tested are eGenes. A constrained-gene recovery fraction measured that way
would rise to ≈1 at rungs 3–4 **whatever the biology**, and the curve would
"support" the power account by construction.

Read the right-hand column — everything scored against one fixed 18,481-gene
universe — and the ladder becomes interpretable: bulk cortex 0.445 → best
single-nucleus cell type 0.579. A real but modest 13-percentage-point step.

The pipeline therefore fixes the denominator in `s01_gene_universe.py`. This was
not in the report's design and it is the difference between a measurable
question and a foregone conclusion.

**Two cautions that follow immediately, both for Stage 1:**

1. **The `qval` vs within-gene-Bonferroni choice moves rungs 3–4 by a factor of
   2–4.** Ext is 0.972 under `qval` and 0.679 under `Fixed_bonf`; MG4 is 0.451
   versus 0.060. This is **D-001**, and it is worth more to the final number
   than almost anything else in the project.
2. **eGene yield tracks cell abundance, not resolution.** Ext > IN > Ast > OD >
   MG > End is very close to the ordering of how many nuclei each cell type
   contributes. And GTEx *cerebellum* (0.568 of tested) essentially matches
   SingleBrain excitatory neurons (0.579 of universe) — a bulk tissue matching
   the best single-nucleus cell type. This is the power confound, visible in the
   data before any constraint stratification. It is exactly why the Bryois
   pseudobulk arm (**D-007**) matters.

Nothing above is stratified by constraint yet. These are whole-transcriptome
rates and say nothing directly about the hypotheses — but they set the ceiling
and the confound structure Stage 1 has to work inside.

---

## 6. Licensing position (early, per the Stage 3 requirement)

| Source | Licence | Redistributable | Handling |
|---|---|---|---|
| SingleBrain (Zenodo 14908182) | CC-BY-4.0 | Yes, with attribution | `data/raw/` |
| Bryois (Zenodo 7276971) | CC-BY-4.0 | Yes, with attribution | `data/raw/` |
| gnomAD v2.1.1 / v4.1 | No restrictions, citation requested | Yes | `data/raw/` |
| GTEx v10 | Open-access summary statistics | Yes | `data/raw/` |
| SCHEMA browser | Summary results freely available | Yes | `data/raw/` |
| PsychENCODE | Not yet verified | **Held as No** | `data/restricted/` |
| MetaBrain | Registration form, no redistribution grant | **No** | not obtained |
| PGC3 SCZ | Data-use agreement | **No** | `data/restricted/` when it arrives |

`data/restricted/` is excluded by `.gitignore` as a whole directory, so a new
file cannot be committed merely because no pattern matched it. The two headline
single-nucleus resources are both CC-BY-4.0, which means the derived recovery
tables can be published without difficulty.

---

## 7. Revised timeline

The report's Stage 1 estimate of weeks 2–4 still holds. Three changes:

- **Gained:** SingleBrain `top_assoc` files are ~1.6 MB each, not gigabytes.
  All 36 downloaded in minutes. Rungs 3 and 4 are essentially free (D-006).
- **Lost:** SingleBrain `full_assoc` files are 4–10 GB each (~150 GB for the
  set). Stage 2 colocalization at rungs 3–4 needs targeted per-gene extraction,
  not a bulk download. Worth planning for now rather than discovering in week 5.
- **Added:** the Bryois pseudobulk arm, if accepted (D-007), costs ~4 GB and one
  extra analysis arm, and is the only clean resolution-vs-power test available.

PGC3 remains unblocking until Stage 2, as intended.

---

## 8. Gene sets built

### 8.1 Constrained genes and matched controls (D-002)

Matching on log GTEx cortex TPM and coding-exon count, caliper 0.25 SD, **with
replacement**, seed 20260903.

| | |
|---|---|
| Constrained genes (LOEUF < 0.35) | 2,937 |
| Matched | **2,767 (94%)** |
| Unique control genes | 1,294 (serving 2,767 case-slots) |
| Effective n, control side | 513 |

Balance after matching:

| Covariate | Constrained | Control | SMD |
|---|---|---|---|
| log TPM (cortex) | 1.020 | 1.014 | 0.009 |
| coding exons | 14.617 | 14.586 | 0.003 |
| brain tissues expressed | 11.939 | 11.842 | 0.030 |
| log gene length | 4.793 | 4.443 | **0.670** (deliberately unmatched) |

Case representativeness: retained vs all constrained genes, LOEUF SMD 0.027.

Matching **with** replacement was not the default and the reason is recorded in
full under D-002: 1:1 without replacement dropped 40% of constrained genes, and
the dropped 40% were systematically the most constrained, longest and
most-expressed — the genes the hypothesis is most about.

### 8.2 SCHEMA sets (D-009)

| Set | Genes in universe |
|---|---|
| **Primary** — Singh et al. 2022, FDR < 0.05 | **32** |
| Singh et al. 2022, exome-wide p < 2.14×10⁻⁶ | 10 |
| Sensitivity — browser 2026-08-21, local BH FDR < 0.05 | 50 |

The pipeline reproduces the published result exactly: SETD1A, CUL1, XPO7, TRIO,
CACNA1G, SP4, GRIA3, GRIN2A, HERC1, RB1CC1 as the exome-wide ten.

**Two flags for Stage 1.** First, the two releases disagree more than expected —
only **12 of the published 32** are also significant in the browser release (20
published-only, 38 browser-only). They are not interchangeable and the
discordance will need explaining. Second, **32 genes is a small set**;
SCHEMA-specific recovery curves will carry wide confidence intervals. The
LOEUF-constrained set (2,937 genes) is the statistical workhorse and SCHEMA is
the sharper, smaller overlay.

---

## 9. Decisions taken, and what is still open

**Settled in Stage 0** (full rationale in `DECISIONS.md`): D-002 control
matching, D-005 constraint source, D-006 SingleBrain file scope, D-007 Bryois
pseudobulk arm accepted, D-009 SCHEMA sets, D-010 rung 1 = `Brain_Cortex`,
D-012 fixed gene universe as denominator. D-008 (PsychENCODE licence) is
provisional pending a terms read before Stage 3.

**Open, needed at Stage 1 kickoff:**

- **D-001** — what counts as a detectable eQTL. Section 5 is the evidence: this
  choice is worth a factor of 2–4 at rungs 3–4 and is the largest single lever
  in the project.
- **D-003** — multiple-testing correction across genes and rungs.

**Open, lower stakes:** **D-011**, whether rung 2 is PsychENCODE alone, waits
for MetaBrain, or is dropped for a three-rung ladder.

**Open, Stage 2:** **D-004** colocalization priors. Also blocked on PGC3 access,
which remains un-blocking until then as intended.
