# The constrained-gene eQTL deficit is a property of bulk tissue RNA-seq, not of sample size or cell-type resolution

**Toby Butler**
Code and data: https://github.com/itsslightning/missing-regulation-ladder

---

## Headline

> **Moving from bulk cortex to single-nucleus data closes 59.9% [95% CI
> 49.5–70.5%] of the cis-eQTL deficit at constrained genes. Neither donor
> count nor cell-type resolution accounts for that closure: a 6.8× increase in
> donors within bulk tissue closes −4.1% [−17.4, +7.7] of it, and varying
> cell-type granularity in two independent studies moves it by −0.007 [−0.021,
> +0.005]. What closes it is the assay. The residual ~40% survives every
> assay, resolution and sample size tested here.**

## Background

Mostafavi *et al.* (2023) showed that cis-eQTLs are depleted near constrained
genes while GWAS signal is enriched there. They proposed a **selection**
account: common regulatory variants of large effect are removed at
dosage-sensitive genes, so the missing regulation is real biology and should
persist at any resolution. Rosen *et al.* (2026, AJHG; bioRxiv
2025.08.05.668745) counter with a **power** account. Larger eQTL studies find
weaker eQTLs, and those look progressively more GWAS-like, so the gap is
substantially a measurement problem that should shrink as power and resolution
increase.

The question had been examined in blood and GTEx, but not in brain, and never
against the SCHEMA schizophrenia gene set. This project builds a resolution
ladder in brain and measures which account the data supports.

*Provisional-source note.* Rosen *et al.* and the November-2025 brain/blood
snRNA preprint have not completed peer review. **No number from either is
load-bearing here.** They supply the hypothesis being tested, not any input or
threshold, so the analysis would stand unchanged if their reported effect sizes
were revised.

## Data

Public summary statistics only; no genotype-level data.

| Source | Role | Donors | Licence |
|---|---|---|---|
| GTEx v10 cortex | bulk tissue, rung 1 | 205 | open access |
| PsychENCODE DER-08 | bulk brain, rung 2 | 1,387 | citation required |
| SingleBrain (Jang 2026) | snRNA-seq, 7 types / 28 subtypes | 983 | CC-BY-4.0 |
| Bryois 2022 | snRNA-seq, 8 types + pseudobulk | 192 | CC-BY-4.0 |
| gnomAD v2.1.1 | LOEUF / pLI constraint | — | open |
| SCHEMA (Singh 2022) | schizophrenia exome gene set | — | publisher supp. |
| PGC3 SCZ wave 3 (European) | GWAS outcome for SMR | 52,017 / 75,889 | CC-BY-4.0 |

Every file is stamped with its URL, UTC download date, byte count and sha256 in
`logs/downloads.json`. MetaBrain was dropped because access is gated behind a
form and can't be scripted.

## Methods

**Gene universe.** One fixed denominator of **18,481** protein-coding genes with
a gnomAD LOEUF value and GTEx brain expression, used for every rung. This choice
is load-bearing rather than bookkeeping. Scored against genes-tested instead,
97% of genes in SingleBrain excitatory neurons are eGenes, the measure
saturates, and constrained-gene recovery would rise to ≈1 at high resolution
*regardless of biology*.

**Case and control sets.** Cases are LOEUF < 0.35 (n = 2,937). Controls are
LOEUF ≥ 1.0, matched on log cortex TPM and coding-exon count within a 0.25 SD
caliper, **with replacement** (seeded). Without replacement, 40% of cases went
unmatched, and the unmatched ones were systematically the most constrained,
longest and most-expressed: the genes the hypothesis is most about. With
replacement retains 2,767 of 2,937 cases at SMD 0.009/0.003. The 1,294 unique
controls carry frequency weights, and control intervals are cluster-bootstrapped
over control genes.

Gene length is deliberately **not** matched on. It is a partial proxy for
regulatory-landscape complexity, which is part of the mechanism under test
rather than a nuisance. The residual imbalance (SMD 0.670) is reported.

**Detection.** Each rung ships a different significance column, and they are not
on a common evidential scale. A Storey q depends on π₀ estimated within each
study, so q ≤ 0.05 is a *more permissive* bar in a better-powered study. That
loosens the threshold precisely where power is highest. Detection is instead
recomputed identically everywhere: per-gene Bonferroni over the cis variants
tested, then Benjamini–Hochberg across all gene × cell-type tests within a rung,
so a rung with 28 cell types pays for its 28 opportunities. Absolute eGene
counts fall below every published figure as a result. That is the price of
comparability, and these numbers should not be compared to a paper's headline
eGene count.

**Causal inference.** SMR (Zhu *et al.* 2016) against PGC3 European. coloc was
excluded as primary on a data fact: it needs full regional statistics, which
exist for only two of five arms. Variants are harmonised on rsID and on the
allele pair, because an undetected orientation flip would silently invert the
sign of every causal estimate.

## Results

**The recovery curve.** The constrained-vs-control gap, primary detection rule:

| Rung | Assay | Donors | Gap [95% CI] |
|---|---|---|---|
| GTEx cortex | bulk | 205 | 0.357 [0.311, 0.402] |
| PsychENCODE | bulk | 1,387 | 0.371 [0.328, 0.413] |
| SingleBrain, 7 types | snRNA | 983 | **0.143** [0.105, 0.180] |
| SingleBrain, 28 subtypes | snRNA | 983 | 0.166 [0.124, 0.204] |

Closure from bulk cortex to single-nucleus is **59.9% [49.5, 70.5]**. Across
every structural variant tested (gnomAD v4.1 instead of v2.1.1, pLI instead of
LOEUF, Benjamini–Yekutieli instead of BH) closure spans 47.4–59.9%, and the
residual gap excludes zero in all of them.

**Donor count is excluded.** PsychENCODE has the most donors on the ladder and
the largest gap. GTEx → PsychENCODE is ×6.77 donors for **−4.1% [−17.4, +7.7]**
closure. PsychENCODE → SingleBrain is ×0.71 donors for **+61.4% [52.3, 70.5]**.

**Cell-type resolution is excluded.** Two tests vary resolution in opposite
directions with donors held fixed. Splitting SingleBrain classes into their own
subtypes gives **−0.007 [−0.021, +0.005]**, pooled over 6 classes, with every
per-class interval spanning zero. Pooling Bryois nuclei into pseudobulk against
its 8 cell types gives **+0.005 [−0.035, +0.048]**. A matched-gene-set control
confirms the two studies don't disagree. Read as a bound, resolution accounts
for at most ~10% of the observed closure.

**Assay survives.** All five arms on one gene set (2,489 constrained, 1,174
control):

| Arm | Assay | Donors | Gap |
|---|---|---|---|
| GTEx cortex | bulk tissue | 205 | 0.367 |
| PsychENCODE | bulk tissue | 1,387 | 0.367 |
| Bryois pseudobulk | snRNA-seq | 192 | 0.144 |
| Bryois, 8 cell types | snRNA-seq | 192 | 0.150 |
| SingleBrain, 7 cell types | snRNA-seq | 983 | 0.114 |

Bulk spread is **0.000** across 6.8× donors, single-nucleus spread 0.035 across
5.1×, and the between-class separation is **0.217**. The sharpest comparison is
Bryois pseudobulk, single-nucleus with all nuclei pooled and the smallest study
here, at 0.144 against PsychENCODE bulk at 0.367 with seven times the donors.

**Causal explanation of GWAS loci.** Of 184 genome-wide significant SCZ loci
(MHC excluded), bulk explains 19 (GTEx) and 27 (PsychENCODE); single-nucleus
explains 59. **41 loci gain a causal eQTL explanation only in single-nucleus
data.** The SMR significance *rate* is near-constant across rungs (6.9–8.2%),
which locates the gain precisely. Single-nucleus data doesn't make eQTLs more
likely to be causal; it makes more genes testable at all.

As a positive control, **C4A** is the top-ranked non-MHC gene with b_xy = +0.189.
Higher expression raising risk is the best-established causal direction in
schizophrenia genetics, and the sign balance across significant results is 0.501,
indicating no systematic orientation error.

**Targets.** 516 genes are causally implicated only in single-nucleus data. Of
those, 175 are small-molecule tractable and 24 are at Phase 1 or beyond,
including the calcium-channel family (CACNA1C, CACNA1D, CACNA1I, CACNB2).
Tractability rates are near-identical across visibility groups (~32–34%), so
single-nucleus data is finding more genes at the same druggable rate rather than
enriching for druggable ones.

## Interpretation

Neither hypothesis the project set out to adjudicate survives intact. Selection
predicts the gap persists at any resolution, and it doesn't. Power predicts it
closes with sample size, and it doesn't do that either. The gap is set by **what
is measured**: nuclei or whole tissue.

Two candidate mechanisms, stated as hypotheses rather than findings. snRNA-seq
captures largely nascent nuclear transcript while bulk captures mature,
stability-buffered cytoplasmic mRNA, so post-transcriptional buffering of
dosage-sensitive genes would damp genotype effects on steady-state mRNA while
leaving them visible in nascent transcription. Alternatively, cell-composition
variance between bulk samples may specifically mask eQTLs at broadly expressed
genes. Fractionation data could distinguish the two.

The residual ~40% closes under no condition tested and remains consistent with
selection. CHRM4 illustrates its practical cost. It sits in this study's own
constrained set (LOEUF 0.265, pLI 0.974), is well expressed in cortex (11.14
TPM), has been an approved drug target since 2024, and has no detectable
cis-eQTL at any assay or resolution here. An eQTL-based target-discovery
pipeline would never surface it through regulatory evidence.

## Limitations

Bulk and single-nucleus studies differ in ancestry, brain region, pipeline and
GENCODE vintage, not only in assay. The evidence that this is assay rather than
a study-level accident is that two bulk studies agree exactly despite 6.8×
different N and different consortia, and three single-nucleus arms agree across
5.1× different N and two consortia. A confound would have to track assay class
across four independent datasets.

Resolution and per-context power are intrinsically coupled in single-cell data,
so no test here isolates resolution alone. Two tests vary it in opposite
directions and bracket it.

SMR cannot separate a shared causal variant from linkage, so the Stage 2 counts
inflate. The HEIDI and coloc sensitivity arms have not been run. Locus
definition is distance-based rather than LD-based.

This is an observational re-analysis. No claim of causality follows from SMR
alone.

## Data and code availability

All analysis code, the decision log and the derived tables are at
https://github.com/itsslightning/missing-regulation-ladder, where
`uv sync && uv run python scripts/run_all.py` reproduces every table and figure.
Non-redistributable sources are excluded from the repository, and only derived
summaries that cannot reconstruct them are published.
