"""Record D-004, the colocalization method and thresholds.

Delegated to Claude by the project owner on 2026-09-06, having been held open
since Stage 0. As with D-001 and D-003, the rejected options run as sensitivity
arms wherever the data allows, rather than being argued away.

Run once. `decide()` refuses to overwrite without `force=True`.
"""

from __future__ import annotations

from pipeline import decisions


D004_RATIONALE = (
    "Chosen 2026-09-06, delegated by the project owner. "
    "The governing principle is the one that decided D-001, D-003 and D-012: "
    "the quantity of interest is a COMPARISON ACROSS RUNGS, how many "
    "additional GWAS loci gain an eQTL explanation as the assay changes, not "
    "an absolute count that matches any published figure. A method must "
    "therefore be computable IDENTICALLY at every rung, and that requirement "
    "eliminates coloc on this data. "
    "coloc needs full regional summary statistics for both traits at every "
    "gene tested. Those exist here for only two of the five arms. PsychENCODE "
    "ships a complete 3.3 GB association file and Bryois ships all nominal "
    "SNP-gene pairs, but GTEx v10 publishes only eGenes (one row per gene) in "
    "the archive used here, and SingleBrain's full_assoc files are 4-10 GB "
    "each, roughly 150 GB for the set, at the ~650 KB/s this project has "
    "measured against Zenodo that is about 64 hours of download, and they are "
    "not indexed for range extraction, so a targeted per-gene pull would still "
    "require fetching each file whole. Running coloc only where the data "
    "happens to be complete would mean comparing rungs scored by different "
    "methods, which is exactly the error D-001 was chosen to avoid. "
    "SMR needs only the top eQTL SNP's effect and standard error, plus the "
    "GWAS effect at that same variant. Every rung ships this: GTEx has slope "
    "and slope_se, SingleBrain has fixed_beta and fixed_sd, PsychENCODE and "
    "Bryois give beta and a nominal p from which the standard error follows. "
    "So SMR is computable everywhere, identically. "
    "The known cost is stated rather than hidden: SMR alone cannot separate a "
    "shared causal variant from linkage between two distinct causal variants, "
    "so absolute counts will be inflated relative to a true colocalization. "
    "That is accepted on the same logic as the Bonferroni conservatism in "
    "D-001: the bias is in the same direction and of similar construction at "
    "every rung, so the cross-rung comparison survives it. The multiplicity "
    "asymmetry between rungs with many cell types and rungs with one is "
    "handled as in D-003, by BH across all gene x cell-type tests within a "
    "rung. "
    "HEIDI is the standard remedy for the linkage ambiguity, but it needs "
    "regional data and so is available only for PsychENCODE and Bryois. It "
    "therefore runs as a SUPPLEMENTARY filter on those two arms, explicitly "
    "not as part of the uniform primary, and is used to estimate how much of "
    "the SMR count is linkage rather than to change the headline. "
    "coloc runs as a sensitivity arm on those same two arms, at p12=1e-5 and "
    "p12=1e-6, to check that SMR is not over-calling relative to a proper "
    "colocalization where both can be computed. If SMR and coloc disagree "
    "sharply there, the SMR-based cross-rung comparison is reported with that "
    "caveat attached. "
    "Thresholds: SMR significance by Benjamini-Hochberg within rung at 0.05, "
    "matching D-003; HEIDI p > 0.05 as the conventional non-rejection of the "
    "single-variant model; coloc PP4 >= 0.8."
)


def main() -> None:
    d = decisions.COLOC_PRIORS
    if d.decided:
        print(f"{d.key} already recorded as {d.value!r}\n  {d.rationale}")
        return
    d.decide("smr_heidi", rationale=D004_RATIONALE)
    print(f"{d.key} recorded: {d.value}")
    print("\nDecision status:")
    print(decisions.summary())


if __name__ == "__main__":
    main()
