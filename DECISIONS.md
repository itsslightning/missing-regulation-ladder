# Decisions log

Every non-trivial methodological choice, the alternatives that were actually on
the table, and why. Written as the choice is made, not reconstructed afterwards.
This is the raw material for the thesis methods section.

Machine-readable twin: `logs/decisions.jsonl` (appended by
`pipeline/decisions.py` whenever a decision is recorded in code, and replayed at
import so a choice survives across runs).

**Status key** — `OPEN`: needs a human call before dependent code can run.
`SET`: decided, with rationale. `PROVISIONAL`: decided but expected to be
revisited when a blocker clears.

---

## Load-bearing decisions — the four sentinels

These can move the headline conclusion on their own, so `pipeline/decisions.py`
holds them as sentinels that raise on use rather than defaulting. Analysis code
that depends on one cannot run until it is recorded.

### D-001 — What counts as a "detectable" cis-eQTL? — **SET** (2026-09-04)

**Delegated to Claude by the project owner** after being held open through
Stage 0. Because it was not chosen by the person whose thesis this becomes, the
rejected options run as sensitivity arms and appear in every figure.

**Chosen:** `uniform_recomputed_fdr`, implemented as **per-gene Bonferroni over
the cis variants tested, applied identically at every rung, then one common
across-gene FDR** (D-003).

**Why not `study_native_threshold`.** GTEx and SingleBrain both ship Storey
q-values, and a Storey q depends on π₀ — the estimated fraction of true nulls —
computed *within each study*. A better-powered study has a lower π₀, which makes
q ≤ 0.05 a **more permissive** bar there than in a weaker study. The threshold
therefore loosens exactly where power is highest, which is precisely the
direction that manufactures the recovery the power account predicts. Measured:
this arm reports **88.2% gap closure vs 59.9%** under the uniform rule. Using it
as primary would have let the significance calling produce the answer.

**Why not `fixed_nominal_p`.** The number of cis variants tested per gene differs
across rungs, so a fixed nominal threshold rewards rungs with denser coverage.

**Why not `effect_size_floor` as primary.** It targets the power confound most
directly, but the rungs do not ship a common effect-size scale: GTEx slopes are
on inverse-normal-transformed expression, SingleBrain betas on
scaled/quantile-normalised expression, and Bryois ships betas without standard
errors. Retained as a sensitivity arm where the scales are close enough to be
informative (it gives 53.1% closure).

**How the uniform rule is computable everywhere:**

| Rung | Per-gene Bonferroni from |
|---|---|
| GTEx | `pval_nominal` × `num_var` |
| SingleBrain | `Fixed_bonf` (already exactly this construction) |
| PsychENCODE | `nominal_pval` × `number_of_SNPs_tested` |
| Bryois | min nominal p × variants tested per gene |

**Known cost, stated not hidden:** Bonferroni over cis variants ignores LD, so it
is stricter than the permutation p-value GTEx would use, and absolute eGene
counts fall below every published figure. That is accepted because it is
conservative *by the same construction at every rung*, which is what a
cross-rung comparison requires.

**Sensitivity arms shipped:** `native`, `effect`. See `docs/figures/fig2`.

### D-002 — Control-gene matching — **SET** (2026-09-04)

**Chosen:** `caliper_expression_exons`, caliper **0.25 SD**, **with replacement**.
Covariates: log GTEx brain-cortex median TPM and gnomAD coding-exon count.
Controls drawn from LOEUF ≥ 1.0; cases are LOEUF < 0.35. Recorded in
`logs/decisions.jsonl`; parameters in `pipeline/s03_control_matching.py`.

**Alternatives rejected:** adding gene length to the covariates; a propensity
score over all covariates; no matching with LOEUF-decile stratification instead.

**Why these covariates:** matches the report's stated caveat exactly and carries
the fewest assumptions. Gene length is deliberately *not* matched on — it is a
partial proxy for regulatory-landscape complexity, which is part of Mostafavi's
proposed mechanism rather than a nuisance covariate, so matching on it would
risk conditioning on a mediator and regressing away the effect being measured.

**Why with replacement — this was not the obvious default.** 1:1 matching
without replacement at 0.25 SD left 40% of constrained genes unmatched, and the
unmatched 40% were not random: they had more coding exons (SMD −1.37), higher
expression (−0.61) and *lower* LOEUF (+0.40) than the matched ones. The genes
the hypothesis is most about were the ones being dropped, because
highly-expressed many-exon genes are rare among unconstrained genes — which is
exactly the confound the matching exists to handle. Measured comparison:

| Strategy | Cases kept | SMD expr | SMD exons | Case bias (LOEUF) |
|---|---|---|---|---|
| 0.25 SD, no replacement | 1,753 / 2,937 | 0.026 | 0.026 | 0.160 |
| 1.0 SD, no replacement | 2,274 / 2,937 | 0.202 | 0.220 | 0.102 |
| **0.25 SD, with replacement** | **2,767 / 2,937** | **0.009** | **0.003** | **0.027** |

**Cost, and how it is handled:** 1,294 unique control genes serve 2,767 cases,
so control observations are reused and correlated. Reuse is mild (median 1×,
90th percentile 4×, max 21×; 828 controls used once). Every control-side
statistic therefore carries a frequency weight, confidence intervals on the
control curve must be cluster-robust by control gene (1,294 clusters), and
effective sample size (513 for a weighted mean) is reported next to raw counts.

**Known residual imbalance:** log gene length, SMD 0.670. Expected and accepted
per the mediator argument above; reported, not hidden.

### D-003 — Multiple-testing correction across genes and rungs — **SET** (2026-09-04)

**Delegated to Claude by the project owner**, same caveat as D-001.

**Chosen:** `bh_within_rung` — Benjamini-Hochberg across **all gene × cell-type
tests within a rung**, then a gene counts as detected at that rung if any of its
cell types is significant.

**The decisive argument is coherence, not conservatism.** Correcting across the
full gene × rung grid would make a gene's eQTL status at rung 1 depend on how
many rungs the analysis happens to include — adding a 29th SingleBrain subtype
would change whether GTEx cortex is called as having an eQTL for gene X. The
number of rungs is *my design choice, not a property of the data*, so it must
not enter the per-rung detection call. That rules out both `bh_across_all` and
`by_across_all` as the calling rule.

**Why BH is applied across gene × cell-type pairs, not per cell type.** Rung 4
has 28 cell types and rung 3 has 7. Calling per cell type and then taking the
union would hand rung 4 four times as many chances purely for having more
columns. Correcting over all pairs within the rung makes each rung pay for its
own opportunities. This is load-bearing: under it, **rung 4 detects fewer genes
than rung 3** (9,066 vs 9,906), whereas under study-native calling the order
reverses (13,961 vs 13,211).

**Why not `by_across_all`.** Beyond the coherence problem, the BY penalty over a
grid of ~18,481 genes × ~40 cell-type columns is a factor of roughly 12. It
would suppress detection at every rung severely enough that the curve would
flatten toward zero and read as support for selection through sheer
conservatism.

**Why not `permutation_null`.** Not rejected on merit — **deferred**. It is the
most defensible option and is the natural upgrade if the headline gap turns out
marginal. It needs permutation of gene labels within matched strata at every
rung: substantial compute and substantially more code to get right. The gap is
currently far from marginal (0.143, CI [0.105, 0.180]), so it is not needed yet.

**Cross-rung multiplicity** is handled where it belongs: the constrained-vs-control
gap is reported per rung with bootstrap CIs, and the small family of rung-level
gap contrasts carries its own BH correction (`gap_q` in
`recovery_by_rung.parquet`).

**Sensitivity arm run 2026-09-04 (`s09_robustness.py`).** `by_across_all` — the
most conservative correction available and the one rejected above — still gives
**47.9% closure [34.4, 60.3]**, with the residual gap excluding zero. So the
headline is not an artefact of a permissive multiple-testing rule. It is the
variant that moves the number most, which is worth stating.

### D-004 — Colocalization method and thresholds — **SET** (2026-09-06)

**Delegated to Claude by the project owner**, same caveat as D-001 and D-003:
the rejected options run as sensitivity arms wherever the data allows.

**Chosen:** `smr_heidi` — SMR at every rung from top-SNP statistics, with BH
within rung at 0.05 (matching D-003). HEIDI as a **supplementary** filter where
regional data exists. coloc as a **sensitivity** arm on the same two arms.

**Why coloc is excluded as primary — a data fact, not a preference.** coloc
needs full regional summary statistics for both traits at every gene. Those
exist for only two of five arms:

| Arm | Full regional data? |
|---|---|
| PsychENCODE | **Yes** — the complete 3.3 GB association file |
| Bryois | **Yes** — all nominal SNP–gene pairs, 4.5 GB |
| GTEx v10 | No — the archive publishes eGenes only, one row per gene |
| SingleBrain | No — `full_assoc` is 4–10 GB per cell type, ~150 GB for the set |
| Bryois pseudobulk | Yes |

At the ~650 KB/s this project has measured against Zenodo, SingleBrain's full
associations are about **64 hours** of download, and they are not indexed for
range extraction, so even a targeted per-gene pull would require fetching each
file whole. Running coloc only where the data happens to be complete would mean
comparing rungs scored by **different methods** — precisely the error D-001 was
chosen to avoid.

**Why SMR works everywhere.** It needs only the top eQTL SNP's effect and
standard error plus the GWAS effect at that variant, and every rung ships it:

| Arm | Effect | Standard error |
|---|---|---|
| GTEx | `slope` | `slope_se` |
| SingleBrain | `fixed_beta` | `fixed_sd` |
| PsychENCODE | `regression_slope` | derived from beta and nominal p |
| Bryois | beta | derived from beta and nominal p |

**Known cost, stated rather than hidden.** SMR alone cannot separate a shared
causal variant from linkage between two distinct causal variants, so absolute
counts are inflated relative to a true colocalization. Accepted on the same
logic as the Bonferroni conservatism in D-001: the bias is of similar
construction at every rung, so the cross-rung comparison survives it even
though the absolute numbers should not be quoted against published coloc counts.

**Supplementary and sensitivity arms:**

- **HEIDI** (p > 0.05, the conventional non-rejection of the single-variant
  model) on PsychENCODE and Bryois, to estimate what fraction of the SMR count
  is linkage. Explicitly not part of the uniform primary, because it needs
  regional data.
- **coloc** at p12 = 1e-5 and 1e-6, PP4 ≥ 0.8, on the same two arms. If SMR and
  coloc disagree sharply where both are computable, the cross-rung comparison
  is reported with that caveat attached.

**Consequence for Stage 2 engineering:** the join key becomes the **variant**,
not the gene. rsID is the common key for three of four rungs (SingleBrain's
`variant_id` *is* an rsID; GTEx ships `rs_id_dbSNP155_GRCh38p13`; Bryois is
rsID-native). PsychENCODE gives `chr:pos` on hg19 and needs a lookup — Bryois's
`snp_pos.txt.gz` supplies exactly that mapping. This is the variant-level
analogue of the unversioned-ENSG problem from Stage 0 and will need the same
loss accounting.

---

## Decisions arising from Stage 0

### D-005 — gnomAD v2.1.1 as the primary constraint source — **SET** (2026-09-04)

**Chosen:** v2.1.1 `oe_lof_upper` (LOEUF) and `pLI` as primary; v4.1 as a
robustness check.

**Alternatives:** v4.1 as primary (newer, ~730k exomes); v2.1.1 only.

**Why:** v2.1.1 is the vintage Mostafavi et al. 2023 and the SCHEMA papers used,
so the constrained-gene definition matches the literature being adjudicated. It
ships `oe_lof_upper_bin` (deciles) and the matching covariates in one table.

**Caveat found:** its `brain_expression` column — the obvious candidate for the
expression covariate — is `NA` for all 19,704 rows. Expression is taken from the
GTEx v10 median-TPM matrix instead.

**Robustness check run 2026-09-04 (`s09_robustness.py`): passes.** Rebuilding
the case set on v4.1 LOEUF < 0.35 more than halves it (2,767 → 1,292 genes,
because v4.1 shifts LOEUF values) and reorders genes, yet the headline closure
is **59.0% [45.5, 72.3]** against the primary's **59.9% [49.1, 70.0]**. The
conclusion does not depend on the gnomAD vintage, so v2.1.1 stays primary on
the original grounds.

**Also checked:** pLI ≥ 0.9 instead of LOEUF as the constraint metric gives
52.2% [38.3, 64.4] on a 3,000-gene case set.

### D-006 — SingleBrain `top_assoc` only for Stages 0–1 — **SET** (2026-09-04)

**Chosen:** the 36 `*_top_assoc.tsv.gz` files (~1.6 MB each). Defer
`*_full_assoc.tsv.gz` to Stage 2, fetched per gene of interest.

**Why:** `top_assoc` is one row per gene with the best variant and a q-value —
exactly and completely what the recovery curve needs. Full associations are
4–10 GB each (~150 GB for the set) and are needed only for colocalization.

**Consequence:** Stage 2 colocalization at rungs 3–4 needs a targeted
extraction strategy, not a bulk download. Noted now, not in week 5.

### D-007 — Bryois pseudobulk arm as a within-study power control — **SET** (2026-09-04)

**Chosen:** accepted. Bryois pseudobulk (`pb.[1-22].gz`) versus Bryois cell
types enters as a fifth, parallel comparison — not a rung on the main ladder.

**Why:** the four rungs differ in donor count, ancestry, normalisation and eQTL
pipeline as well as resolution, so the main curve partly measures power. The
Bryois February 2023 update added a "tissue-like" pseudobulk analysis
aggregating reads across all nuclei per individual — **same donors, same
pipeline, same normalisation** as the eight cell-type files, differing only in
resolution. It is the only comparison in this design where resolution is
isolated from sample size.

**Alternative:** rely only on reporting recovery against donor N, and accept
"resolution-plus-power" framing with no clean within-study test.

**Cost:** ~4 GB additional download, one extra analysis arm.

### D-008 — PsychENCODE held as non-redistributable — **PROVISIONAL** (2026-09-04)

**Chosen:** written to `data/restricted/`.

**Why:** the files download without a click-through, which is not the same as a
licence to rehost. Held restricted until the resource.psychencode.org terms are
read and quoted here. Derived counts are publishable either way.

**Revisit:** before Stage 3.

### D-009 — SCHEMA gene set: both releases, published primary — **SET** (2026-09-04)

**Chosen:** Singh et al. 2022 Supplementary Table 5 as the **primary** set;
the browser release 2026-08-21 as a **sensitivity** arm.

**Alternatives:** published only (fewer genes, wider CIs); browser only (more
cases but not citable, and requires inventing an FDR).

**Why:** the published table ships `Q meta`, the FDR the paper's own thresholds
are defined on, so the primary set needs no invented significance rule. The
browser release ships p-values only, so any gene set drawn from it requires a
locally computed FDR — a judgment layered on top of the gene-set definition.
Running both shows the conclusion is robust to the choice.

**Validation:** the pipeline reproduces the published result exactly — **32
genes at FDR < 0.05** and **10 at exome-wide significance** (SETD1A, CUL1,
XPO7, TRIO, CACNA1G, SP4, GRIA3, GRIN2A, HERC1, RB1CC1).

**Flag for Stage 1 — the two releases disagree substantially.** Browser release
under local BH FDR < 0.05 gives 50 genes, of which only **12 overlap** the
published 32 (20 published-only, 38 browser-only). Do not present them as
interchangeable, and expect to explain the discordance.

**Second flag:** 32 genes is a small set. SCHEMA-specific recovery curves will
have wide confidence intervals. The LOEUF-constrained set (2,937 genes) is the
statistical workhorse; SCHEMA is the sharper, smaller overlay.

**Cost:** the Singh supplementary table is keyed by gene symbol, not Ensembl ID,
so it needs a symbol → ENSG hop through gnomAD. 17,740 of 18,320 mapped; 580
lost to symbol drift.

### D-010 — Rung 1 is GTEx `Brain_Cortex` alone — **SET** (2026-09-04)

**Chosen:** `Brain_Cortex` only.

**Alternatives:** `Brain_Frontal_Cortex_BA9`; the union of all 13 brain tissues;
cortex primary with union as sensitivity.

**Why:** closest to the report's "GTEx v10 cortex (bulk)", and the cleanest
anatomical comparison to PsychENCODE prefrontal cortex and SingleBrain
neocortical nuclei. A union of 13 tissues would have far more power and would
inflate rung 1, flattening the very curve being measured.

**Supporting observation:** GTEx cerebellum reaches 0.568 of tested genes,
essentially matching SingleBrain excitatory neurons (0.579 of universe) — a bulk
tissue matching the best single-nucleus cell type. Tissue choice at rung 1 is
worth as much as a rung of the ladder, which is why it is recorded rather than
assumed.

### D-012 — Recovery is scored against a fixed gene universe — **SET** (2026-09-04)

**Chosen:** one fixed 18,481-gene denominator (protein-coding, has gnomAD v2.1.1
LOEUF, has GTEx brain expression), used for every rung.

**Alternatives:** score each rung against the genes that rung tested; score
against all protein-coding genes regardless of annotation coverage.

**Why this is load-bearing, not bookkeeping.** Measured against genes tested,
97% of the genes SingleBrain excitatory neurons tested are eGenes — the measure
has no headroom, and a constrained-gene recovery fraction would rise to ≈1 at
rungs 3–4 *whatever the biology*, "supporting" the power account by
construction. A per-rung denominator also rewards a cell type for testing fewer
genes, and constrained genes are more broadly expressed, so they would be
systematically advantaged.

Against the fixed universe the ladder becomes interpretable: bulk cortex 0.445 →
best single-nucleus cell type 0.579.

**Cost:** genes not tested at a rung count as "no detectable eQTL" there, which
conflates "tested and null" with "not expressed in this cell type". That is the
right default for this question — a gene with no eQTL detectable in a cell type
is missing regulation in that cell type either way — but it must be stated, and
the per-rung tested counts are kept in `docs/stage0_audit.md` so the alternative
can be computed.

### D-013 — Derived outputs are split on licence, not convenience — **SET** (2026-09-04)

**Chosen:** `s04` writes two tables. `data/processed/schema_gene_sets.parquet`
(published) carries the boolean set memberships and the browser-derived
columns. `data/restricted/schema_gene_sets_full.parquet` (local, gitignored)
adds `schema_p_published` and `schema_q_published`.

**Why:** those two columns are Singh et al.'s Supplementary Table 5 reproduced
for 17,740 genes. Crossref reports the paper under **Springer Nature
text-and-data-mining terms, not CC-BY**, so shipping them is republishing a
substantial part of a restricted table — which the project brief explicitly
rules out — rather than publishing a derived result.

The boolean memberships are a different matter: the 32-gene and 10-gene sets
are the paper's own headline result, freely citable, and reproducible from the
flags alone. Nothing analytically necessary is lost.

**Alternative:** publish the full table and rely on the data being "public
anyway". Rejected — freely downloadable is not the same as licensed to rehost,
which is the distinction `sources.py` exists to enforce.

**Generalises to:** any future derived table touching PsychENCODE (D-008) or
PGC3. The rule is that a derived output may ship if it cannot be used to
reconstruct a restricted source table.

---

## Still open

### D-011 — Rung 2 source, given MetaBrain is inaccessible — **OPEN**

MetaBrain has no public file index; access is gated behind a Google Form, so it
cannot be automated.

| Option | Consequence |
|---|---|
| PsychENCODE alone as rung 2 | Proceeds now. Rung 2 is one bulk-brain study rather than the meta-analysed resource the report assumed, and its file supplies no denominator of its own (see the Stage 0 report). |
| Wait for MetaBrain access | Matches the report; unknown delay, human action required. |
| Drop rung 2, use a 3-rung ladder | Simplest; loses the bulk-tissue → bulk-brain step where the report expected the first increment. |
