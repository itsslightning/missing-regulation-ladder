"""Record the Stage 1 decisions, D-001 and D-003.

Both were delegated to me by the project owner on 2026-09-04, having been held
open through Stage 0. Because they are conclusion-shaping and were not chosen
by the person whose thesis this becomes, the analysis runs the rejected options
as sensitivity arms rather than only the chosen one. If the headline conclusion
flips under a rejected option, that is a finding to report, not a detail to
bury.

Run once. `decide()` refuses to overwrite without `force=True`.
"""

from __future__ import annotations

from pipeline import decisions


D001_RATIONALE = (
    "Chosen 2026-09-04, delegated by the project owner. "
    "study_native_threshold is rejected as PRIMARY for a specific reason found "
    "in Stage 0: both GTEx and SingleBrain ship Storey q-values, and a Storey "
    "q depends on pi0, the estimated fraction of true nulls, computed WITHIN "
    "each study. A better-powered study has a lower pi0, which makes q<=0.05 a "
    "MORE permissive evidential bar there than in a weaker study. The threshold "
    "is therefore not constant across rungs -- it loosens exactly where power "
    "is higher, which is precisely the direction that would manufacture the "
    "recovery the power account predicts. Using it as primary would let the "
    "significance calling produce the answer. "
    "fixed_nominal_p is rejected because the number of cis variants tested per "
    "gene differs across rungs, so a fixed nominal threshold rewards rungs with "
    "denser variant coverage. "
    "effect_size_floor is attractive in principle -- it targets the power "
    "confound most directly -- but is rejected as primary because the rungs do "
    "not ship a common effect-size scale that survives scrutiny: GTEx slopes "
    "are on inverse-normal-transformed expression, SingleBrain betas on "
    "scaled/centred/quantile-normalised expression, and Bryois ships betas with "
    "no accompanying standard errors in the summary files. It is retained as a "
    "sensitivity arm where the scales are close enough to be informative. "
    "uniform_recomputed_fdr is chosen and implemented as: per-gene Bonferroni "
    "over the cis variants tested for that gene, applied identically at every "
    "rung, followed by one common across-gene FDR (D-003). This is computable "
    "at every rung from what each ships -- GTEx pval_nominal x num_var, "
    "SingleBrain's Fixed_bonf column which is exactly this construction, "
    "PsychENCODE nominal_pval x number_of_SNPs_tested, and Bryois min nominal p "
    "x variants per gene. It is conservative relative to permutation because it "
    "ignores LD between cis variants, but it is conservative BY THE SAME "
    "CONSTRUCTION at every rung, which is what a cross-rung comparison "
    "requires. Absolute eGene counts will be lower than any published figure; "
    "that is expected and is the price of comparability. "
    "study_native_threshold and effect_size_floor both run as sensitivity arms."
)

D003_RATIONALE = (
    "Chosen 2026-09-04, delegated by the project owner. "
    "The decisive argument is coherence, not conservatism. Correcting across "
    "the full gene x rung grid makes a gene's eQTL status at rung 1 depend on "
    "how many rungs the analysis happens to include -- adding a 29th SingleBrain "
    "subtype would change whether GTEx cortex is called as having an eQTL for "
    "gene X. The number of rungs is a design choice of mine, not a property of "
    "the data, so it must not enter the per-rung detection call. That rules out "
    "bh_across_all and by_across_all as the primary calling rule. "
    "by_across_all is additionally rejected because the BY penalty over a grid "
    "this size (~18,481 genes x ~40 rung-columns) is a factor of roughly 12, "
    "which would suppress detection at every rung so severely that the curve "
    "would flatten toward zero and read as support for selection through sheer "
    "conservatism. "
    "permutation_null is the most defensible option and is not rejected on "
    "merit -- it is deferred. It requires permuting gene labels within matched "
    "strata across every rung, which is substantial compute and substantially "
    "more code to get right, and it is the natural upgrade if the headline gap "
    "turns out to be marginal. "
    "bh_within_rung is chosen: Benjamini-Hochberg across genes within each "
    "rung, which is the family every eQTL study corrects within, and which is "
    "identical in construction at every rung given D-001. Cross-rung "
    "multiplicity is then handled where it actually belongs -- in the "
    "inferential statements about the gap, not in the detection call. The "
    "constrained-vs-control gap is reported per rung with confidence intervals, "
    "and the small family of rung-level gap contrasts carries its own BH "
    "correction. by_across_all runs as a sensitivity arm."
)


def main() -> None:
    for decision, option, rationale in (
        (decisions.DETECTABLE_EQTL, "uniform_recomputed_fdr", D001_RATIONALE),
        (decisions.MULTIPLE_TESTING, "bh_within_rung", D003_RATIONALE),
    ):
        if decision.decided:
            print(f"{decision.key} already recorded as {decision.value!r}")
            continue
        decision.decide(option, rationale=rationale)
        print(f"{decision.key} recorded: {option}")

    print("\nDecision status:")
    print(decisions.summary())


if __name__ == "__main__":
    main()
