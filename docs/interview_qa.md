# Interview Q&A — The Missing Regulation Ladder

Answers grounded in what the pipeline found, with the numbers. Reasoning traces
to [`../DECISIONS.md`](../DECISIONS.md); every figure is reproducible with
`uv run python scripts/run_all.py`.

---

### "Give me the project in thirty seconds."

Disease genes are the ones eQTL studies are worst at finding regulatory
variants for. Two published accounts disagree about why: Mostafavi says
selection removed the variants, Rosen says we lack power to see them. Nobody
had tested it in brain. I built a resolution ladder — bulk cortex, bulk brain,
single-nucleus cell types, subtypes — and measured the gap between constrained
genes and matched controls at each rung.

The gap closes 59.9% [49.5–70.5%]. But the cause is neither hypothesis. Donor
count doesn't do it and cell-type resolution doesn't do it. **The assay does** —
bulk tissue versus single-nucleus. That's a third account, and the residual 40%
still looks like selection.

---

### "How do you know this isn't just a power/sample-size difference in disguise?"

*This was the project's main flagged risk, so I tested it three ways rather than
caveating it.*

1. **Within bulk tissue, donors do nothing to the gap.** PsychENCODE has 1,387
   donors against GTEx's 205 — 6.8× — and the gap is 0.367 in both, agreeing to
   three decimal places. Closure across that step is **−4.1% [−17.4, +7.7]**.
2. **The winning arm has *fewer* donors.** PsychENCODE → SingleBrain is a 29%
   *decrease* in donors and a **+61.4%** closure.
3. **The smallest study shows the closed gap.** Bryois pseudobulk, 192 donors —
   the smallest arm here — has a gap of 0.144, against PsychENCODE's 0.367 at
   seven times the donors.

Bulk spread is 0.000 across 6.8× donors; single-nucleus spread is 0.035 across
5.1×; the separation between assay classes is 0.217. If sample size drove this,
the gap would track N within each class. It doesn't track it at all.

**Where power *does* matter:** absolute counts. Stage 2 loci go 19 → 27 within
bulk for 6.8× donors. More donors mean more genes clear the instrument
threshold. They just don't close the constrained-gene *deficit*.

---

### "Why brain and not blood?"

The hypotheses are about schizophrenia genes, and SCHEMA genes are expressed and
act in brain. Blood eQTLs would test whether the phenomenon exists, not whether
it explains the disease. It also isn't new: FastGxC and OneK1K/CIGMA have done
the constrained-gene question in blood — the brain, and SCHEMA specifically,
was the open part.

The honest limitation: brain is where the sample sizes are worst, which is why
the power confound had to be handled properly rather than waved at.

---

### "Why these four rungs and not more or fewer?"

The four were pre-registered from the report, and they span the axis being
tested: bulk single tissue → bulk meta-analysed brain → 7 major cell types → 28
subtypes.

What I'd change with hindsight: the ladder is the wrong *shape* for the answer
it produced. Rungs 3 and 4 differ only in granularity, which turned out not to
matter, while rungs 1–2 versus 3–4 differ in assay, which turned out to be
everything. If I rebuilt it, the primary axis would be assay with donor count
as a covariate, and granularity as a secondary arm.

I also added a rung the report didn't specify — Bryois pseudobulk — because it
varies resolution in the *opposite* direction (pooling rather than splitting)
with donors fixed. That arm did more work than any of the original four.

---

### "What would change your conclusion?"

Four concrete things, in order of how much they'd move me:

1. **A bulk brain eQTL study at ~1,000 donors closing the gap as much as
   SingleBrain does.** That would restore the power account and kill the assay
   claim. Nobody has run it; it's the obvious next experiment.
2. **A confound tracking assay class across all four datasets.** The assay
   classes also differ in ancestry, brain region and pipeline. My defence is
   that two bulk studies from different consortia agree exactly, and three
   single-nucleus arms from two consortia agree — a confound would have to
   track assay across four independent datasets. Not impossible.
3. **The residual disappearing under a more sensitive uniform rule.** My
   per-gene Bonferroni is conservative by design. If a permutation-calibrated
   uniform rule closed the residual 40%, selection would be in trouble.
4. **HEIDI or coloc disagreeing sharply with SMR** where both are computable.
   That wouldn't touch Stage 1, but it would undercut the Stage 2 locus counts.

---

### "What counts as a 'detectable' eQTL, and why does it matter?"

It's the single largest lever in the project, which is why it's decision D-001
and why the dashboard exposes it as a toggle.

Each study ships a different significance column and they aren't on a common
scale. GTEx and SingleBrain both ship **Storey q-values**, and a Storey q
depends on π₀ — the estimated fraction of true nulls — computed *within each
study*. A better-powered study has lower π₀, so q ≤ 0.05 is a **more permissive
bar** there. It loosens exactly where power is highest, which is the direction
that manufactures the recovery the power account predicts.

Concretely: study-native thresholds report **88%** closure. My uniform rule
reports **60%**. Anyone quoting "single-nucleus closes ~90% of the gap" from
published eGene counts is partly quoting that artefact.

So I recompute identically everywhere: per-gene Bonferroni over cis variants,
then BH within rung. It's conservative — absolute counts fall below every
published figure — but conservative *by the same construction at every rung*,
which is what a cross-rung comparison needs.

---

### "Your control genes — how do you know the comparison is fair?"

Constrained genes are longer, more expressed and have more exons, and every one
of those independently predicts eQTL discovery. So I match on log cortex TPM and
coding-exon count within a 0.25 SD caliper.

Two things I'd want to be asked about:

**Why matching with replacement.** 1:1 without replacement left 40% of cases
unmatched — and not a random 40%: they had more exons (SMD −1.37), higher
expression (−0.61) and *lower* LOEUF (+0.40). It was discarding the genes the
hypothesis is most about. With replacement retains 94% at better balance (SMD
0.009/0.003) and a representative case set. The cost is reused controls, so
they carry frequency weights and intervals are cluster-bootstrapped.

**Why I deliberately *don't* match on gene length.** Length is a partial proxy
for regulatory-landscape complexity, which is part of Mostafavi's proposed
mechanism — matching on it would condition on a mediator and regress away the
effect being measured. The residual imbalance (SMD 0.670) is reported rather
than hidden.

---

### "Show me you didn't just p-hack the denominator."

Fair, because the denominator nearly did decide the answer.

Scored against *genes tested*, 97% of genes in SingleBrain excitatory neurons
are eGenes. The measure saturates, and constrained-gene recovery would rise to
≈1 at high resolution **regardless of biology**. Scored against one fixed
18,481-gene universe, it becomes measurable: 0.445 → 0.579.

That's decision D-012, and it was made before the recovery curve was computed.

The related trap: the raw LOEUF-decile view is dominated by expression — median
cortex TPM falls 12.71 → 0.18 across deciles — and at single-nucleus rungs it
*inverts*, making the most constrained genes look best recovered. That figure
is in the repo with all three panels, because it's the clearest demonstration
of why matched controls exist.

---

### "You were given two hypotheses and returned a third. Isn't that convenient?"

It's the opposite of convenient — it's the result I'd least like to defend,
because it means neither of the papers I set out to referee is right.

I'd point at the falsifiability. The assay claim makes a sharp prediction: bulk
brain at ~1,000 donors should *not* close the gap. That's testable and would
sink the claim. And I found it by repairing a broken rung, not by looking for
it: PsychENCODE was provisional for most of the project because its available
file contained only significant genes. When I got the full 3.3 GB file, the
highest-donor bulk study turned out to have the largest gap — which overturned
a conclusion I'd already committed and written up.

The correction is in the git history and in the report.

---

### "What did you get wrong?"

Three, all in the repo:

1. **I concluded the closure was donor count** and committed it. Rung 2 was
   broken at the time — the one data point that could separate donors from
   assay. Repairing it inverted the conclusion.
2. **I reported the sign of noise.** The resolution test compared point
   estimates against an arbitrary cutoff and I wrote up "narrows in 3 of 6
   classes." With paired-bootstrap intervals attached, all six span zero. The
   corrected version is a *stronger* claim — a bound of ≤10% — but it was wrong
   as first written.
3. **A figure with a false title.** I labelled the LOEUF-decile figure
   "monotone" when it plainly wasn't, and the underlying pattern turned out to
   be an expression confound worth a three-panel figure of its own.

I'd rather show these than a clean narrative — the audit trail is the point.

---

### "The CHRM4 connection — did it work out?"

No, and that's the better answer.

No muscarinic receptor reaches SMR significance. CHRM4 has **zero detectable
cis-eQTLs at any assay or resolution**. But look at what CHRM4 is: LOEUF 0.265
(decile 0 — in this study's own constrained case set), pLI 0.974, 11.14 TPM in
cortex so not a low-expression artefact, and an approved drug target since
Cobenfy in 2024.

CHRM4 is the thesis in one gene. It's a validated target with no findable
regulatory variation, so a genetics-led pipeline built on eQTLs would never have
surfaced it — my other repo found it by a route that didn't depend on eQTLs.
It sits in the residual 40% that no assay closes.

That's the practical cost of missing regulation: eQTL-based target discovery has
a systematic blind spot exactly at the dosage-sensitive genes most likely to
matter.

---

### "Why SMR rather than coloc?"

A data constraint, not a preference. coloc needs full regional summary
statistics at every gene; those exist for only two of five arms. GTEx v10
publishes eGenes only, and SingleBrain's full associations are ~150 GB — about
64 hours at the throughput I measured, and not range-indexed. Running coloc only
where the data happens to be complete would compare rungs scored by different
methods, which is the exact error D-001 exists to avoid.

SMR needs only the top SNP's effect and standard error, which every rung ships.
The cost: it can't separate a shared causal variant from linkage, so counts
inflate — but identically at every rung, and the question is a cross-rung
comparison. **These counts shouldn't be quoted against published coloc
figures**, and HEIDI/coloc sensitivity arms on the two complete arms are the
next thing I'd run.

---

### "How do I know your pipeline is right?"

A positive control chosen before I looked: **C4A** is the top non-MHC gene with
b_xy = +0.189 — higher expression raises risk, the best-established causal
direction in schizophrenia. CACNA1C, FURIN, KCTD13, BTN3A2 are all in the top
15. The b_xy sign balance is 0.501, so no systematic orientation error.

That last one matters more than it sounds. Harmonisation reported a 99.7% allele
flip, which looked like a serious bug. It isn't: PGC3 defines A1 as the
*reference* allele while GTEx and SingleBrain report effects against the
*alternate*. Because I match on actual allele letters rather than an assumed
convention, the flip rate is confirmation. Had I got it wrong, every causal
estimate would have flipped sign — a risk-increasing gene reading as protective,
silently and plausibly.

---

### "What would you do next with more time?"

1. **HEIDI and coloc** on PsychENCODE and Bryois, to bound how much of the SMR
   count is linkage.
2. **Nuclear vs cytoplasmic fractionation eQTLs**, to test the mechanism the
   assay result implies — that constrained genes are post-transcriptionally
   buffered, so genotype effects survive in nascent transcript and are damped in
   steady-state mRNA.
3. **A permutation-calibrated detection rule**, the one option in D-003 I
   deferred rather than rejected, to check the residual 40% isn't an artefact of
   Bonferroni conservatism.
4. **Push someone to run bulk brain at 1,000 donors**, which is the experiment
   that would settle the assay claim either way.
