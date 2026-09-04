# Decisions log

Every non-trivial methodological choice, the alternatives that were actually on
the table, and why. Written as the choice is made, not reconstructed afterwards.
This is the raw material for the thesis methods section.

Machine-readable twin: `logs/decisions.jsonl` (appended by
`pipeline/decisions.py` whenever a decision is recorded in code).

**Status key** — `OPEN`: needs a human call before dependent code can run.
`PROPOSED`: recommended, awaiting sign-off. `SET`: decided, with rationale.
`PROVISIONAL`: decided but expected to be revisited when a blocker clears.

---

## Load-bearing decisions — held OPEN by design

These four can move the headline conclusion on their own, so `pipeline/decisions.py`
holds them as sentinels that raise on use rather than defaulting. Analysis code
that depends on one cannot run until it is recorded.

### D-001 — What counts as a "detectable" cis-eQTL? — **OPEN**

*Needed by: Stage 1.*

This is the y-axis of the recovery curve. A permissive definition makes every
rung look successful and flatters H2 (power); a strict one suppresses the small
effects single-nucleus data exists to reveal and flatters H1 (selection).

Complication found in Stage 0: **the rungs do not ship the same statistics.**
GTEx v10 and SingleBrain both provide per-gene permutation q-values. Bryois
provides nominal p-values only, with no per-gene correction. Any definition must
say how that is bridged.

| Option | Meaning |
|---|---|
| `study_native_threshold` | Each study's own call (q ≤ 0.05); reproduce an equivalent per-gene correction for Bryois. Comparable to published eGene counts; correction differs slightly per rung. |
| `uniform_recomputed_fdr` | One identical per-gene FDR recomputed across all rungs. Maximum internal comparability; numbers no longer match any published count. |
| `fixed_nominal_p` | One fixed nominal p threshold on the top variant per gene. Transparent, but favours rungs with denser variant coverage. |
| `effect_size_floor` | Significance plus a minimum effect size, so "detected" means the same magnitude everywhere. Targets the power confound directly; discards the real small-effect eQTLs that are the point of higher resolution. |

### D-002 — Control-gene matching tolerance and covariates — **OPEN**

*Needed by: Stage 0 (final step).*

The control curve is the entire comparison. Constrained genes are longer, more
highly expressed and have more exons, and each of those independently predicts
eQTL discovery power. Loose matching leaves a gap that is really a length or
expression artefact. Over-tight matching on covariates that are themselves
*consequences* of constraint regresses away the effect being measured —
regulatory-landscape complexity is part of Mostafavi's proposed mechanism, not a
nuisance variable.

| Option | Meaning |
|---|---|
| `caliper_expression_exons` | Nearest-neighbour within a caliper on expression decile and coding-exon count only, per the report's caveat. Fewest assumptions; leaves gene length and TSS density unbalanced. |
| `caliper_plus_length` | Adds gene length and number of cis variants tested. Closes the obvious power confounds; risks conditioning on a mediator. |
| `propensity_score` | One propensity model over all covariates. Uses everything at once; harder to defend in a viva, hides which covariate does the work. |
| `stratified_no_matching` | No matched set; compare across LOEUF deciles with covariate adjustment in regression. Discards no genes, adjustment is explicit; loses the simple two-curve plot. |

### D-003 — Multiple-testing correction across genes and rungs — **OPEN**

*Needed by: Stage 1.*

Tests are dependent in two directions at once: the same gene is retested at
every rung, and neighbouring cell subtypes share donors and nuclei. BH assumes a
dependence structure this design does not have.

| Option | Meaning |
|---|---|
| `bh_within_rung` | BH across genes within each rung, nothing across rungs. Each rung its own experiment; cross-rung comparison becomes descriptive. |
| `by_across_all` | BY across the full gene × rung grid. Valid under arbitrary dependence; substantially more conservative, will lower apparent recovery everywhere. |
| `bh_across_all` | BH across the full grid. More powerful; independence/PRDS assumption not satisfied here. |
| `permutation_null` | Permute gene labels within matched strata and calibrate empirically. No parametric dependence assumption, most defensible; more compute and more code to get right. |

### D-004 — Colocalization priors and posterior threshold — **OPEN**

*Needed by: Stage 2 (also blocked on PGC3 access).*

coloc's `p12` prior sets how readily a shared causal variant is declared, and the
Stage 2 headline — how many loci gain an eQTL explanation per rung — moves with
it. Subtype rungs test more features, compounding the effect.

| Option | Meaning |
|---|---|
| `coloc_default_priors` | p1=p2=1e-4, p12=1e-5, PP4 ≥ 0.8. Comparable to most papers; widely argued to over-declare. |
| `coloc_conservative_p12` | p12=1e-6, PP4 ≥ 0.8. Guards against over-calling; lowers absolute counts at every rung, and the cross-rung contrast is what matters. |
| `coloc_sensitivity_band` | Report across p12 ∈ [1e-6, 1e-5] as a band. Most honest; complicates the single headline statement. |
| `smr_heidi` | SMR + HEIDI instead. Needs only top-SNP statistics, so it works where full associations are too large; HEIDI rejects linkage rather than confirming a shared variant — a subtly different question. |

---

## Decisions arising from Stage 0

### D-005 — gnomAD v2.1.1 as the primary constraint source — **SET** (2026-09-04)

**Chosen:** v2.1.1 `oe_lof_upper` (LOEUF) and `pLI` as primary; v4.1 downloaded
as a robustness check.

**Alternatives:** v4.1 as primary (newer, ~730k exomes, better-powered
constraint estimates); v2.1.1 only (no robustness check).

**Why:** v2.1.1 is the vintage Mostafavi et al. 2023 and the SCHEMA papers used,
so the constrained-gene definition matches the literature the project is
adjudicating. It also ships `oe_lof_upper_bin` (the decile field) and every
covariate the control matching needs — `num_coding_exons`, `cds_length`,
`gene_length`, `brain_expression` — in a single table. v4.1 recomputes LOEUF and
shifts gene rankings; if the recovery curve changes qualitatively under v4.1
that is a finding to report, not a bug to hide.

**Revisit if:** the v4.1 robustness check materially changes the conclusion.

### D-006 — SingleBrain `top_assoc` only for Stages 0–1 — **SET** (2026-09-04)

**Chosen:** download the 36 `*_top_assoc.tsv.gz` files (~1.6 MB each) for the
recovery curve. Defer `*_full_assoc.tsv.gz` to Stage 2, fetched per gene of
interest rather than in bulk.

**Alternatives:** download all full associations up front (~150 GB across 36
files, 4–10 GB each).

**Why:** `top_assoc` is one row per gene with the best variant and a q-value,
which is exactly and completely what the recovery curve needs. Full associations
are needed only for colocalization. Pulling 150 GB to compute a statistic that
uses one row per gene would be wasteful and would put the project's storage
footprint beyond what a laptop and a public repo can carry.

**Consequence to be aware of:** Stage 2 colocalization at rungs 3–4 will need a
targeted extraction strategy, not a bulk download. Noted now so it is not a
surprise later.

### D-007 — Add a Bryois pseudobulk arm as a within-study power control — **PROPOSED**

*Not in the original report. Recommended addition; needs sign-off because it
changes the design.*

**Proposed:** treat Bryois pseudobulk (`pb[1-22].gz`) versus Bryois cell types
as a fifth, parallel comparison — not a rung on the main ladder.

**Why:** the report names cross-study comparability as the project's main
technical risk, and it is a real one: the four rungs differ in donor count,
ancestry, normalisation and eQTL pipeline as well as in resolution, so the main
curve partly measures power. The February 2023 update to the Bryois Zenodo
record added a "tissue-like" pseudobulk analysis aggregating reads across all
nuclei per individual — **same donors, same pipeline, same normalisation** as
the eight cell-type files, differing only in resolution.

That is the single comparison available anywhere in this design where resolution
is isolated from sample size. It converts the project's biggest weakness into a
controlled contrast, and it is the natural answer to the interview question
"how do you know this isn't just a power difference?".

**Cost:** ~4 GB of additional download and one extra analysis arm.

**Alternative:** rely only on reporting recovery against donor N (mitigation
already planned), and accept "resolution-plus-power" as the framing with no
clean within-study test.

### D-008 — PsychENCODE held as non-redistributable pending a licence read — **PROVISIONAL** (2026-09-04)

**Chosen:** write PsychENCODE files to `data/restricted/`.

**Why:** the files download without a click-through, which is not the same as a
licence to rehost. Held restricted until the resource.psychencode.org terms are
read and quoted here. Derived counts are publishable either way, so this costs
nothing now.

**Revisit:** before Stage 3, per the project's licensing check.

---

## Open questions escalated to the project owner

These emerged in Stage 0, are conclusion-shaping, and are **not** being defaulted.

### D-009 — Which SCHEMA release defines the gene set? — **OPEN**

The SCHEMA browser currently serves a release dated **2026-08-21** carrying
**87,959 cases / 150,587 controls** — roughly 3.6× the 24,248 cases of the
published Singh et al. 2022 analysis. It also ships **no q-value column**, only
per-gene p-values, so an FDR would have to be computed locally.

This is not a minor version bump. It changes which genes are "SCHEMA genes",
how many there are, and whether the gene set is citable to a peer-reviewed paper.

| Option | Consequence |
|---|---|
| Published Singh et al. 2022 set | Citable, peer-reviewed, matches the literature being adjudicated; fewer genes (10 exome-wide significant, 32 at FDR < 0.05), so less power in the recovery curve. |
| Browser release 2026-08-21 | ~3.6× the cases, more genes, better-powered curve; not a citable peer-reviewed set, needs a locally computed FDR, and a thesis examiner will ask why the published set was not used. |
| Both | Primary analysis on one, sensitivity analysis on the other. Most defensible; roughly doubles Stage 1 work. |

### D-010 — What defines rung 1 (GTEx brain)? — **OPEN**

GTEx v10 has 13 brain tissues. "Cortex" alone, a single frontal-cortex region,
or the union across brain regions are materially different rung-1 baselines —
a union has far more power, which would flatten the very curve being measured.
Because rung 1 is the baseline the whole ladder is measured against, this choice
propagates into every downstream number.

| Option | Consequence |
|---|---|
| `Brain_Cortex` only | Closest to the report's "GTEx v10 cortex (bulk)"; single tissue, modest N, cleanest comparison to PsychENCODE prefrontal cortex. |
| `Brain_Frontal_Cortex_BA9` only | Best anatomical match to PsychENCODE and SingleBrain (neocortex); slightly smaller N. |
| Union of all 13 brain tissues | Most power, most eGenes; inflates rung 1 and makes the ladder's first step look artificially high. |
| Cortex primary + union as sensitivity | Defensible; extra work. |

### D-011 — Rung 2 source, given MetaBrain is inaccessible — **OPEN**

MetaBrain has no public file index; access is gated behind a Google Form the
maintainers distribute links through, so it cannot be automated. See the Stage 0
report.

| Option | Consequence |
|---|---|
| PsychENCODE alone as rung 2 | Proceeds now, no blocker. Rung 2 is then one bulk-brain study rather than the meta-analysed resource the report assumed. |
| Wait for MetaBrain access | Matches the report; unknown delay, human action required. |
| Drop rung 2, use a 3-rung ladder | Simplest; loses the bulk-tissue → bulk-brain step, which is where the report expected the first increment. |
