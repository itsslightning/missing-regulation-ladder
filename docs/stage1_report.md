# Stage 1 report — the recovery curve

Date: 2026-09-04. Figures: [`figures/`](figures/). Decisions:
[`../DECISIONS.md`](../DECISIONS.md). Tables: `data/processed/recovery_by_rung.parquet`.

---

## 1. Headline

Going from bulk cortex to single-nucleus major cell types, the
constrained-gene eQTL gap falls from **0.357 to 0.143**, a closure of

> **59.9% of the gap, 95% CI [49.5%, 70.5%].**
>
> The residual **40.1% [29.5%, 50.5%]** does not close, and its interval
> excludes zero at every rung under every detection rule tested.

**And the closure is caused by neither of the things this project set out to
adjudicate.** §4 tests all three candidate explanations directly and excludes
two of them:

| Explanation | Test | Verdict |
|---|---|---|
| Donor count (Rosen's power account) | ×6.8 donors within bulk tissue | **excluded** — gap unchanged (−4.1% closure) |
| Cell-type resolution | splitting *and* pooling, two studies | **excluded** — bounded at ≤10% |
| **Assay: nuclei vs whole tissue** | bulk vs snRNA-seq, matched genes | **survives** — 0.217 separation |

The single sharpest number: **Bryois pseudobulk**, single-nucleus data with all
nuclei pooled and 192 donors, the smallest study here, has a gap of **0.144**,
against PsychENCODE bulk tissue at **1,387 donors** with a gap of **0.367**.
Seven times the donors, and more than twice the gap.

So the honest decomposition is: **~60% of the constrained-gene deficit is a
property of bulk tissue RNA-seq that single-nucleus RNA-seq does not share; the
residual ~40% survives every assay, resolution and sample size tested here and
is the part consistent with selection.**

## 2. The curve

Primary arm: uniform per-gene Bonferroni + BH within rung (D-001, D-003).
Constrained = LOEUF < 0.35 (2,767 matched genes); controls = 1,294 unique
matched unconstrained genes, frequency-weighted (D-002).

| Rung | Donors | Constrained | Control | Gap [95% CI] |
|---|---|---|---|---|
| GTEx cortex (bulk tissue) | 205 | 0.201 | 0.557 | **0.357** [0.311, 0.402] |
| PsychENCODE (bulk brain) | 1,387 | 0.324 | 0.695 | 0.371 [0.328, 0.413] |
| SingleBrain, 7 major types | 983 | 0.624 | 0.767 | **0.143** [0.105, 0.180] |
| SingleBrain, 28 subtypes | 983 | 0.542 | 0.708 | 0.166 [0.124, 0.204] |

Every gap is significant at the bootstrap resolution limit (p = 0.0005 = 1/2000
replicates, BH-adjusted across the rung contrasts).

**The SCHEMA genes move further than the constrained set as a whole.** Of the 32
published SCHEMA genes (FDR < 0.05), **4 have a detectable cis-eQTL in bulk
cortex; 19 do at single-nucleus major-cell-type resolution**, 12.5% → 59.4%.

Splitting the 32 three ways, now that rung 2 has a real denominator (full
table: [`schema_switch.md`](schema_switch.md)):

- **13 switch on**: invisible in bulk, detectable in single-nucleus data. This
  is the Stage 2 colocalization target list, and the cell types matter for
  which SingleBrain full-association files to pull (4–10 GB each, so they are
  fetched per gene, per D-006):

  | Gene | Switches on in | | Gene | Switches on in |
  |---|---|---|---|---|
  | TRIO | End, Ext, OD, OPC | | SRRM2 | Ast |
  | DNM3 | Ast, IN, MG, OPC | | CUL1 | MG |
  | FAM178A | Ext, IN, OD, OPC | | SP4 | Ext |
  | XPO7 | Ast, Ext, OD | | ZMYM2 | IN |
  | STAG1 | Ext, MG | | NR3C2 | OD |
  | KDM6B | Ext, IN | | ZNF136 | Ext |
  | SV2A | IN, OD | | | |

- **6 were already visible** in bulk: MAGI2, AKAP11, GRIN2A, ANKRD12, PREP,
  CACNA1G.
- **10 remain undetected** in any assay at any resolution: ASH1L, RB1CC1,
  FAM120A, GRIA3, HCN4, HERC1, SLC22A11, OR4P4, HIST1H1E, MAGEC1. This is the
  residual missing regulation, and the most interesting group for the selection
  hypothesis. It includes RB1CC1 and HERC1, two of the ten exome-wide
  significant genes.

**A correction from the provisional rung 2.** With the significant-only
PsychENCODE file this table read 15 switching on and GRIN2A among them, and the
report drew attention to GRIN2A appearing in oligodendrocyte lineage rather
than neurons. With a real bulk denominator, GRIN2A and CACNA1G are detectable
in bulk after all. The oligodendrocyte observation was an artefact of the broken
rung and has been withdrawn.

What survives is that several genes switch on in a **single** cell type: SRRM2
in astrocytes, CUL1 in microglia, SP4 in excitatory neurons, NR3C2 in
oligodendrocytes. That is genuine cell-type-specific regulatory signal, and it
is worth holding alongside §4's finding that granularity does not move the
*aggregate* gap: knowing *where* a gene is regulated is useful for a drug
target even when it is not what makes the gene detectable in the first place.

## 3. Robustness to the decision I was asked to make

D-001 and D-003 were delegated to me rather than chosen by the project owner,
so the rejected options run alongside the chosen one in every figure:

| D-001 arm | Gap, bulk → sn_major | Closure [95% CI] |
|---|---|---|
| **uniform Bonferroni + BH** (primary) | 0.357 → 0.143 | **59.9%** [49.5, 70.5] |
| each study's own q ≤ 0.05 | 0.306 → 0.036 | 88.2% [80.8, 95.4] |
| uniform + \|β\| ≥ 0.1 | 0.359 → 0.168 | 53.1% [42.7, 63.9] |

**The direction is robust; the magnitude is not.** Every arm shows substantial
partial closure with a residual gap excluding zero. But the study-native arm
puts closure at 88%, nearly double the effect-size arm's 53%.

This is exactly why the native arm was rejected as primary. Both GTEx and
SingleBrain ship Storey q-values, and a Storey q depends on π₀ estimated
*within each study*, so q ≤ 0.05 is a **more permissive** bar in a
better-powered study. It loosens precisely where power is highest, which
inflates apparent recovery. Anyone quoting "single-nucleus resolution closes
~90% of the gap" from published eGene counts is, in part, quoting that artefact.

### The other three choices the headline rests on

D-001 is the largest lever, but not the only one. The remaining structural
choices were re-run end to end, each rebuilding its own matched control set
under the same D-002 rules and seed, since a different case set needs its own
controls:

| Family | Variant | n cases | Gap bulk → sn | Closure [95% CI] |
|---|---|---|---|---|
| Constraint source | **gnomAD v2.1.1 LOEUF < 0.35** (primary) | 2,767 | 0.357 → 0.143 | **59.9%** [49.1, 70.0] |
| | gnomAD v4.1 LOEUF < 0.35 | 1,292 | 0.306 → 0.125 | 59.0% [45.5, 72.3] |
| Constraint metric | gnomAD v2.1.1 pLI ≥ 0.9 | 3,000 | 0.220 → 0.105 | 52.2% [38.3, 64.4] |
| Multiple testing | **BH within rung** (primary) | 2,767 | 0.357 → 0.143 | **59.9%** [49.1, 70.0] |
| | BY across the whole grid | 2,767 | 0.323 → 0.170 | 47.4% [34.3, 59.9] |
| SCHEMA release | published 32 genes (primary) | 32 | recovery 0.125 → 0.594 | — |
| | browser 50 genes | 50 | recovery 0.140 → 0.620 | — |

The primary row appears in both tables with slightly different intervals,
[49.5, 70.5] above, [49.1, 70.0] here. That is Monte Carlo noise: the two
tables draw independent 2,000-replicate bootstraps of the same quantity. A
difference of ~0.4 percentage points is the resolution of the bootstrap, not a
discrepancy, and it is worth knowing that is the precision on offer.

**Closure spans 47.4%–59.9% across all of these, and the residual gap excludes
zero in every one.** Two results are worth calling out:

- **The v4.1 check, promised in D-005, passes.** Recomputing constraint on
  ~730k exomes reorders genes and more than halves the case set (2,767 → 1,292,
  because v4.1 LOEUF values shift), yet closure is 59.0% against 59.9%. The
  conclusion does not depend on the gnomAD vintage.
- **BY across the whole grid, the most conservative correction available and
  the one I rejected on coherence grounds, still gives 47.4% closure.** So the
  headline is not an artefact of a permissive multiple-testing rule. It is,
  however, the variant that moves the number most, which is worth knowing.

## 4. The power confound, which qualifies everything above

**The rungs differ in donor count as much as in resolution.** GTEx cortex has
205 donors; SingleBrain has 983. That is a 4.8× increase in N accompanying the
resolution increase, and the report flagged this as the project's main technical
risk. It is present, and it is not small.

Four pieces of internal evidence bear on it, and they point the same way:

**(a) Splitting a cell class into its subtypes does not measurably change the
gap** (figure 4). This is the tightest test currently available. SingleBrain
reports each major class both pooled and split (Ast against Ast1–Ast4, Ext
against Ext1–Ext8, and so on) from the **same donors, the same nuclei, the
same pipeline**. The only thing that changes is the grouping.

| Class | Pooled gap | Split gap | Δ [95% CI] |
|---|---|---|---|
| Ast (4 sub) | 0.142 | 0.113 | −0.030 [−0.067, +0.006] |
| Ext (8 sub) | 0.154 | 0.183 | +0.029 [−0.005, +0.059] |
| IN (7 sub) | 0.206 | 0.170 | −0.036 [−0.077, +0.003] |
| MG (4 sub) | 0.061 | 0.034 | −0.027 [−0.060, +0.004] |
| OD (3 sub) | 0.104 | 0.114 | +0.010 [−0.016, +0.035] |
| OPC (2 sub) | 0.086 | 0.096 | +0.010 [−0.009, +0.028] |
| **Pooled** | | | **−0.007 [−0.021, +0.005]** |

Δ intervals come from a *paired* bootstrap: both arms are evaluated on the same
replicate and the difference taken within it, so the strong correlation between
them cancels rather than inflating the interval.

**Every per-class interval spans zero**, so no individual class shows an
effect. An earlier version of this report said the gap "narrows in 3 of 6
classes and widens in 3", which was reporting the sign of noise and has been
corrected. The pooled estimate is what carries the claim: **−0.007 [−0.021,
+0.005]**.

Read that as a bound, not a null. The data **rule out any change larger than
about 0.021 in either direction**, against a bulk-to-single-nucleus closure of
0.214. So cell-type resolution accounts for **at most ~10% of the observed
closure**, with a point estimate near 3%. It does not prove resolution
contributes nothing; it puts a ceiling on how much it could contribute.

**(b) The same holds at ladder scale.** Rungs 3 and 4 differ only in grouping
(7 classes vs 28 subtypes) at identical donor count, and the gap goes 0.143 →
0.166, slightly *wider*.

**(c) A bulk tissue already matches the best cell type.** From the Stage 0
audit, GTEx *cerebellum* calls eGenes for 0.568 of tested genes, essentially
matching SingleBrain excitatory neurons at 0.579. A bulk tissue with good N
reaches single-nucleus territory.

**(d) eGene yield across SingleBrain cell types tracks cell abundance.**
Ext > IN > Ast > OD > MG > End is close to the ordering of how many nuclei each
type contributes. That is a power ordering, not a biological one.

### The caveat that limits (a) and (b), and probably limits any such test

Splitting a class gives each subtype fewer nuclei, so per-context power falls
even though donor count is unchanged. Both arms of both tests therefore trade
resolution against reads-per-context rather than isolating resolution.

That is not a flaw in this design; it is a property of the data. **In
single-cell data, resolution and per-context power are intrinsically coupled**:
at fixed sequencing depth you cannot resolve more contexts without putting
fewer reads in each. The Bryois pseudobulk arm has the identical confound
running the other way: pooling all nuclei raises reads-per-context while
lowering resolution.

So Bryois will not be a clean isolation of resolution either. What the two
directions together *can* do is bracket it: if pooling (Bryois) and splitting
(this test) both leave the gap unmoved, resolution is not the active ingredient
in either direction, and the closure on the main ladder must be coming from
donor count. That is a bounded, defensible conclusion, and it is the one this
project can actually reach.

### The Bryois arm, complete, and it agrees

All 198 files (4.5 GB, 22 chromosomes) are downloaded and size-verified, so
D-007 now runs at full coverage on **2,489 constrained genes**.

| Arm | Constrained | Control | Gap [95% CI] |
|---|---|---|---|
| Bryois pseudobulk (all nuclei pooled) | 0.063 | 0.207 | 0.144 [0.103, 0.185] |
| Bryois 8 cell types | 0.141 | 0.291 | 0.150 [0.102, 0.193] |
| **Paired difference** | | | **+0.005 [−0.035, +0.048]** |

**The chr1–2 negative was noise, exactly as flagged.** At 2 of 22 chromosomes
this read −0.037 [−0.075, +0.004] and looked like it might contradict the
SingleBrain result. At full coverage it is **+0.005**, as close to zero as the
data can put it. This is why that number was not quoted at the time.

**Matched-gene-set control:** SingleBrain restricted to the same 2,489 genes
gives **+0.019 [−0.015, +0.052]**, overlapping Bryois. No evidence the studies
disagree.

### The three tests converge

| Test | Direction of change | Result | Bound, as % of the 0.214 closure |
|---|---|---|---|
| SingleBrain, split 6 classes into subtypes | more resolution | **−0.007** [−0.021, +0.005] | **≤10%** |
| Bryois, pseudobulk → 8 cell types | more resolution | **+0.005** [−0.035, +0.048] | ≤22% |
| SingleBrain on the Bryois gene set | more resolution | +0.019 [−0.015, +0.052] | ≤24% |

Two independent studies, with different donors, pipelines and normalisations,
varying resolution in **opposite directions** (Bryois by pooling nuclei
together, SingleBrain by splitting them apart) and all three intervals centre
on zero.

That is the bracket the design was built to produce. If cell-type resolution
were what closes the gap, splitting should have narrowed it and pooling should
have widened it, and the two tests should have disagreed in sign. They do not.

**Conclusion: cell-type resolution does not detectably close the
constrained-gene gap, bounded at ~10% of the observed closure by the tighter
test.** The ~60% closure on the main ladder is therefore attributable to donor
count, 205 to 983 and a 4.8× increase in power, not to seeing cell types
separately.

### CORRECTION: it is not donor count either

An earlier version of this report concluded that the closure was "attributable
to the 4.8× increase in donor count". **That was wrong**, and it was wrong
because rung 2 was still the provisional significant-only file, the single
data point that could separate donor count from assay was the broken one.

With rung 2 repaired, PsychENCODE has the **highest donor count on the ladder**
(1,387, more than SingleBrain's 983) and the **largest gap**:

| From → to | Donor change | Gap | Closure [95% CI] |
|---|---|---|---|
| GTEx cortex → PsychENCODE | **×6.77** | 0.357 → 0.371 | **−4.1%** [−17.4%, +7.7%] |
| PsychENCODE → SingleBrain | **×0.71** | 0.371 → 0.143 | **+61.4%** [52.3%, 70.5%] |

Nearly a sevenfold increase in donors closes **nothing**. A 29% *decrease* in
donors, accompanied by a change of assay, closes **61%**.

### What actually closes it: the assay

Scoring all five arms on one common gene set, the 2,489 constrained and 1,174
control genes Bryois tested, so no arm is advantaged by which genes it covers
(figure 5):

| Arm | Assay | Donors | Gap [95% CI] |
|---|---|---|---|
| GTEx cortex | bulk tissue RNA-seq | 205 | 0.367 [0.324, 0.414] |
| PsychENCODE | bulk tissue RNA-seq | 1,387 | 0.367 [0.322, 0.410] |
| Bryois pseudobulk | snRNA-seq, nuclei pooled | 192 | 0.144 [0.105, 0.187] |
| Bryois, 8 cell types | snRNA-seq | 192 | 0.150 [0.104, 0.197] |
| SingleBrain, 7 cell types | snRNA-seq | 983 | 0.114 [0.077, 0.149] |

- **Bulk arms: spread 0.000** across a 6.8× range of donors. The two bulk
  studies agree to three decimal places.
- **Single-nucleus arms: spread 0.035** across a 5.1× range of donors.
- **Between classes: 0.217**, six times the larger within-class spread.

The decisive single comparison is **Bryois pseudobulk**: single-nucleus data
with all nuclei pooled per donor, the *smallest* study on the ladder at 192
donors, and it shows the closed gap (0.144), less than half PsychENCODE's
0.367 at seven times the donors.

**So all three candidate explanations can now be tested, and two are excluded:**

| Explanation | Test | Verdict |
|---|---|---|
| Donor count | ×6.8 donors within bulk | **excluded** — gap unchanged |
| Cell-type granularity | splitting and pooling, both directions | **excluded** — bounded ≤10% |
| Assay (nuclei vs tissue) | bulk vs snRNA-seq at matched genes | **survives** — 0.217 separation |

### My read

The constrained-gene eQTL deficit is **a property of bulk tissue RNA-seq that
largely disappears in single-nucleus RNA-seq**, and neither sample size nor
cell-type resolution accounts for it.

That is a third account, distinct from both hypotheses this project set out to
adjudicate. Mostafavi's selection account predicts the gap persists at any
resolution; it does not. Rosen's power account predicts it closes with sample
size; it does not do that either. What closes it is changing what you measure.

**Candidate mechanisms, stated as hypotheses rather than findings.** I cannot distinguish these
and am not claiming to:

1. **Nuclear vs cytoplasmic RNA.** snRNA-seq captures largely nascent,
   unspliced transcript; bulk captures mature, stability-buffered cytoplasmic
   mRNA. Constrained genes are dosage-sensitive, and post-transcriptional
   buffering would damp genotype effects on steady-state mRNA while leaving
   them visible in nascent transcription. This is the mechanistically most
   interesting possibility and it is testable with nuclear/cytoplasmic
   fractionation data.
2. **Cell-composition noise.** Variation in cellular composition between bulk
   samples adds variance that could specifically mask eQTLs at broadly
   expressed genes, which constrained genes are.
3. **Normalisation and quantification differences** between the two assay
   families.

**What I cannot rule out.** The bulk and single-nucleus studies differ in more
than assay: ancestry, brain region, pipeline, GENCODE vintage. The evidence
that this is assay rather than a study-level accident is that the two bulk
studies agree exactly despite 6.8× different N and different consortia, and the
three single-nucleus arms agree despite 5.1× different N and two different
consortia. A confound would have to track assay class across four independent
datasets.

**What would change this conclusion**, stated so it can be checked:

- a bulk brain eQTL study at ~1,000 donors closing the gap much *less* than
  SingleBrain does at the same N. That would restore a resolution effect;
- a resolution axis this design cannot reach, e.g. spatial or activity-state
  contexts rather than cell-type labels;
- the residual gap disappearing under a detection rule that is uniform across
  rungs but more sensitive than per-gene Bonferroni.

## 5. What was provisional, and how it resolved

**Rung 2 (PsychENCODE) is complete.** It was scored for most of this stage from
the Bonferroni-filtered release, which contains only significant pairs: 8,190
genes, every one an eGene by construction. Genes absent from that file cannot be
told apart from tested-and-null, so its rates were deflated by an unknown amount
and its position on the curve meant nothing. The full association file (3.3 GB,
hg19, no liftOver needed because gene-level detection uses no coordinates) has
since been downloaded and verified against its byte count, and every number in
this report is scored from it. The provisional shading has been removed from the
figures.

**The curiosity did not survive, and something better replaced it.** On the
filtered file PsychENCODE had the highest donor count on the ladder (1,387) and
the *lowest* constrained recovery (0.154), which looked like evidence against a
pure-N account. With the real denominator its constrained recovery is 0.324,
above GTEx cortex at 0.201, so the ordering was an artefact of the filtered file
exactly as suspected.

What the repair produced instead is stronger. The constrained-gene *gap* is
0.357 at 205 donors and 0.371 at 1,387: unchanged across a 6.8-fold difference
in sample size, in two studies from different consortia. That is the observation
§4 and the assay contrast are built on, and it is a cleaner argument against the
power account than the artefact it replaced, because it rests on the quantity
the matched controls were designed to protect.

## 6. Method notes worth carrying into the thesis

**The unadjusted LOEUF-decile view is misleading and must not be shown alone**
(figure 3). Across LOEUF deciles 0 → 9, median cortex TPM falls from 12.71 to
0.18 (a 70-fold drop), coding exons from 17.9 to 3.4, and brain tissues
expressed from 12.1 to 4.5. The least-constrained decile is barely expressed, so
it has few detectable eQTLs for reasons unrelated to selection. The confound is
strong enough that **at single-nucleus rungs the unadjusted view inverts, making
the most-constrained decile look best recovered**. Holding expression fixed
restores the expected direction at bulk rungs and shows a nearly flat curve at
single-nucleus rungs. This is the clearest possible demonstration of why the
matched control set (D-002) exists.

**Constrained genes are detected in fewer cell types when detected at all**:
mean 1.56 vs 1.96 major cell types, and 2.41 vs 4.01 subtypes, against matched
controls. Consistent with more restricted, more context-specific regulation,
which is the direction both hypotheses predict, so it discriminates neither, but
it is a real feature of the data.

**Absolute eGene counts are below every published figure**, deliberately.
Bonferroni over cis variants ignores LD and is stricter than the permutation
p-values GTEx reports. It is applied identically at every rung, which is what a
cross-rung comparison needs, and it is the reason these numbers should not be
compared to a paper's headline eGene count.

## 7. Status

- All four rungs: complete, rung 2 on the full PsychENCODE association file.
- Within-SingleBrain resolution test (fig 4): complete.
- Bryois power-control arm (D-007): complete, and it carries more weight than
  planned. It is the arm that separates assay from donor count.
- D-004: decided (SMR with HEIDI, recorded 2026-09-06). Stage 2 is complete;
  the HEIDI and coloc sensitivity arms are not yet run.
- Robustness battery (constraint source, constraint metric, multiple-testing
  rule, SCHEMA release): complete. Closure spans 47.4%-59.9%.
