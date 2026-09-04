# Stage 1 report — the recovery curve

Date: 2026-09-04. Figures: [`figures/`](figures/). Decisions:
[`../DECISIONS.md`](../DECISIONS.md). Tables: `data/processed/recovery_by_rung.parquet`.

---

## 1. Headline

Going from bulk cortex to single-nucleus major cell types, the
constrained-gene eQTL gap falls from **0.357 to 0.143** — a closure of

> **59.9% of the gap, 95% CI [49.5%, 70.5%].**
>
> The residual **40.1% [29.5%, 50.5%]** does not close, and its interval
> excludes zero at every rung under every detection rule tested.

That is a decomposition, not a verdict, and it is the honest shape of the
answer. But there is a serious qualification in §4 that has to travel with it.

## 2. The curve

Primary arm: uniform per-gene Bonferroni + BH within rung (D-001, D-003).
Constrained = LOEUF < 0.35 (2,767 matched genes); controls = 1,294 unique
matched unconstrained genes, frequency-weighted (D-002).

| Rung | Donors | Constrained | Control | Gap [95% CI] |
|---|---|---|---|---|
| GTEx cortex (bulk tissue) | 205 | 0.201 | 0.557 | **0.357** [0.311, 0.402] |
| PsychENCODE (bulk brain) † | 1,387 | 0.154 | 0.502 | 0.347 [0.301, 0.391] |
| SingleBrain, 7 major types | 983 | 0.624 | 0.767 | **0.143** [0.105, 0.180] |
| SingleBrain, 28 subtypes | 983 | 0.542 | 0.708 | 0.166 [0.124, 0.204] |

† provisional — see §5.

Every gap is significant at the bootstrap resolution limit (p = 0.0005 = 1/2000
replicates, BH-adjusted across the rung contrasts).

**The SCHEMA genes move further than the constrained set as a whole.** Of the 32
published SCHEMA genes (FDR < 0.05), **4 have a detectable cis-eQTL in bulk
cortex; 19 do at single-nucleus major-cell-type resolution** — 12.5% → 59.4%.
Fifteen genes that carry exome-wide schizophrenia evidence acquire regulatory
evidence that bulk tissue could not see. That is the Stage 2 target list.

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
*within each study* — so q ≤ 0.05 is a **more permissive** bar in a
better-powered study. It loosens precisely where power is highest, which
inflates apparent recovery. Anyone quoting "single-nucleus resolution closes
~90% of the gap" from published eGene counts is, in part, quoting that artefact.

## 4. The power confound — the thing that qualifies everything above

**The rungs differ in donor count as much as in resolution.** GTEx cortex has
205 donors; SingleBrain has 983. That is a 4.8× increase in N accompanying the
resolution increase, and the report flagged this as the project's main technical
risk. It is present, and it is not small.

Four pieces of internal evidence bear on it, and they point the same way:

**(a) Splitting a cell class into its subtypes does nothing to the gap**
(figure 4). This is the tightest test currently available. SingleBrain reports
each major class both pooled and split — Ast against Ast1–Ast4, Ext against
Ext1–Ext8, and so on — from the **same donors, the same nuclei, the same
pipeline**. The only thing that changes is the grouping. Result: the gap moves
by a mean of **−0.007**, narrowing in 3 of 6 classes and widening in 3.

| Class | Pooled gap | Split gap | Change |
|---|---|---|---|
| Ast (4 sub) | 0.142 | 0.113 | −0.030 |
| Ext (8 sub) | 0.154 | 0.183 | +0.029 |
| IN (7 sub) | 0.206 | 0.170 | −0.036 |
| MG (4 sub) | 0.061 | 0.034 | −0.027 |
| OD (3 sub) | 0.104 | 0.114 | +0.010 |
| OPC (2 sub) | 0.086 | 0.096 | +0.010 |

If cell-type resolution were the mechanism, splitting should have narrowed the
gap systematically. It did not.

**(b) The same holds at ladder scale.** Rungs 3 and 4 differ only in grouping
(7 classes vs 28 subtypes) at identical donor count, and the gap goes 0.143 →
0.166 — slightly *wider*.

**(c) A bulk tissue already matches the best cell type.** From the Stage 0
audit, GTEx *cerebellum* calls eGenes for 0.568 of tested genes, essentially
matching SingleBrain excitatory neurons at 0.579. A bulk tissue with good N
reaches single-nucleus territory.

**(d) eGene yield across SingleBrain cell types tracks cell abundance.**
Ext > IN > Ast > OD > MG > End is close to the ordering of how many nuclei each
type contributes — a power ordering, not a biological one.

### The caveat that limits (a) and (b), and probably limits any such test

Splitting a class gives each subtype fewer nuclei, so per-context power falls
even though donor count is unchanged. Both arms of both tests therefore trade
resolution against reads-per-context rather than isolating resolution.

That is not a flaw in this design; it is a property of the data. **In
single-cell data, resolution and per-context power are intrinsically coupled** —
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

**My read.** The data are consistent with a large power component, a real
residual, and — on four independent internal checks — **no evidence that
cell-type resolution itself is what closes the gap**.

That is a sharper claim than "Rosen is right", and it is the most interesting
thing Stage 1 has produced. It splits Rosen's account in two. The *power* half
looks well supported here. The *resolution* half has no support in this data at
all: every time resolution is varied with donor count held fixed, the gap does
not move. The natural implication — testable by someone with the data, not by
me — is that **bulk cortex at 983 donors would close about as much of the gap as
SingleBrain does**. Nobody has run that.

**I cannot settle this from the rungs alone** and am not going to claim
otherwise. The Bryois pseudobulk arm (D-007) is still the best remaining
evidence and is downloading (198 files, ~4.8 GB). Per the caveat above it will
not isolate resolution cleanly either, but it varies it in the opposite
direction, and agreement between the two would bracket the conclusion.

Until it lands, the defensible framing is **"resolution-plus-power closes ~60%
of the gap, and what evidence there is points at the power half"**.

## 5. What is provisional

**Rung 2 (PsychENCODE) should not be read yet.** The available file is the
Bonferroni-filtered release, which contains only significant pairs — 8,190
genes, every one an eGene by construction. Genes absent from it cannot be
distinguished from tested-and-null, so its rates are deflated by an unknown
amount and its position on the curve is not meaningful. The full association
file (3.3 GB, hg19) is downloading; it needs no liftOver because gene-level
detection uses no coordinates. It is shaded in every figure.

**A curiosity to revisit once rung 2 is fixed:** PsychENCODE has the highest
donor count on the ladder (1,387) and currently shows the *lowest* constrained
recovery. If that survives the full file, it is evidence against a pure-N
account and worth a hard look. If it does not survive, it was an artefact of the
filtered file. Either way it should not be interpreted now.

## 6. Method notes worth carrying into the thesis

**The unadjusted LOEUF-decile view is misleading and must not be shown alone**
(figure 3). Across LOEUF deciles 0 → 9, median cortex TPM falls from 12.71 to
0.18 — a 70-fold drop — coding exons from 17.9 to 3.4, and brain tissues
expressed from 12.1 to 4.5. The least-constrained decile is barely expressed, so
it has few detectable eQTLs for reasons unrelated to selection. The confound is
strong enough that **at single-nucleus rungs the unadjusted view inverts, making
the most-constrained decile look best recovered**. Holding expression fixed
restores the expected direction at bulk rungs and shows a nearly flat curve at
single-nucleus rungs. This is the clearest possible demonstration of why the
matched control set (D-002) exists.

**Constrained genes are detected in fewer cell types when detected at all**:
mean 1.56 vs 1.96 major cell types, and 2.41 vs 4.01 subtypes, against matched
controls. Consistent with more restricted, more context-specific regulation —
which is the direction both hypotheses predict, so it discriminates neither, but
it is a real feature of the data.

**Absolute eGene counts are below every published figure**, deliberately.
Bonferroni over cis variants ignores LD and is stricter than the permutation
p-values GTEx reports. It is applied identically at every rung, which is what a
cross-rung comparison needs, and it is the reason these numbers should not be
compared to a paper's headline eGene count.

## 7. Status

- Rungs 1, 3, 4: complete.
- Rung 2: provisional pending the full PsychENCODE file (downloading).
- Within-SingleBrain resolution test (fig 4): complete.
- Bryois power-control arm (D-007): downloading.
- D-004 (colocalization priors) remains open; Stage 2, also awaiting PGC3.
- `by_across_all` sensitivity arm for D-003: available, not yet run.
