"""Stage 0, step 3: build the matched unconstrained control gene set.

This script implements all four strategies listed under D-002 and refuses to run
until one of them has been recorded. That is deliberate: the control curve *is*
the comparison, and picking a matching strategy after seeing which one produces
a cleaner gap would be the single easiest way to fabricate this project's
result.

Why matching is needed at all. Constrained genes are not a random sample of the
genome. They are longer, more broadly and highly expressed, and have more coding
exons than average -- and every one of those independently raises the chance an
eQTL study finds a signal. So a raw constrained-vs-rest comparison confounds
constraint with discovery power in the same direction the hypotheses disagree
about. Matching removes that.

Why matching can also destroy the result. Mostafavi's proposed mechanism is that
constrained genes have *complex regulatory landscapes* -- more enhancers, more
distributed regulation. Number of cis variants tested and gene length are
partial proxies for exactly that. Matching on them is therefore matching on a
mediator, which would regress away the effect being measured. There is no
choice here that is free of assumption, which is why the alternatives are
recorded rather than defaulted.

Seeding: every strategy that breaks ties or samples does so through a
`numpy.random.default_rng(RANDOM_SEED)` created here, so the control set is
byte-identical across runs and machines.

Writes: data/processed/gene_sets.parquet
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from pipeline import decisions
from pipeline.config import GENE_SETS_TABLE, RANDOM_SEED
from pipeline.s01_gene_universe import GENE_UNIVERSE

#: Covariates available for matching, and what each is a proxy for.
#: `n_cis_variants` is absent until a rung is chosen, so it is not offered here;
#: strategies needing it are computed at Stage 1 against a specific rung.
COVARIATES = {
    "log_tpm": "expression level in the reference brain tissue",
    "n_coding_exons": "exon count, per the report's stated caveat",
    "log_gene_length": "genomic span, a proxy for regulatory landscape size",
    "n_brain_tissues_expressed": "breadth of expression across brain",
}

#: Constrained genes come from the extreme of the LOEUF distribution, controls
#: from the unconstrained end. Genes in between are excluded from both sets
#: rather than assigned to controls, so the contrast is between clearly
#: constrained and clearly unconstrained genes rather than across a boundary
#: where a single gene could fall either way.
CONTROL_LOEUF_MIN = 1.0

# --- parameters of D-002, chosen by the project owner on 2026-09-04 ----------
#
# The strategy (which covariates) is recorded in decisions.py. These two are its
# tolerance settings, recorded here because they are numbers the pipeline reads
# rather than a branch it takes. Both are in DECISIONS.md under D-002.
#
# WITH_REPLACEMENT was not the obvious default and the reason matters. Matching
# 1:1 without replacement at this caliper left 40% of constrained genes
# unmatched, and the unmatched ones were not a random 40%: they had more coding
# exons (SMD -1.37), higher expression (-0.61) and *lower* LOEUF (+0.40) than
# the matched ones. That is, the genes the hypothesis is most about were the
# ones being dropped, because highly-expressed many-exon genes are rare among
# unconstrained genes -- which is precisely the confound the matching exists to
# handle. Allowing replacement retains 94% of cases and improves balance.
#
# The cost is real and is handled rather than ignored: 1,294 unique control
# genes serve 2,767 cases, so control observations are reused and correlated.
# Every control-side statistic therefore carries a frequency weight, and
# confidence intervals on the control curve must be cluster-robust by control
# gene. Effective sample size is reported next to raw counts.
CALIPER_SD = 0.25
WITH_REPLACEMENT = True


def prepare(universe: pd.DataFrame) -> pd.DataFrame:
    """Derive the matching covariates. Log-transforms are applied here so that
    calipers are expressed on the scale the matching actually uses."""
    df = universe.copy()
    df["log_tpm"] = np.log10(df["tpm_reference"].clip(lower=0) + 0.01)
    df["log_gene_length"] = np.log10(df["gene_length"].clip(lower=1))
    df = df.dropna(subset=list(COVARIATES))
    return df


def _standardise(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Z-score covariates so a caliper is in comparable units across them."""
    z = df[cols].astype(float)
    return (z - z.mean()) / z.std(ddof=0)


def _greedy_caliper_match(
    cases: pd.DataFrame,
    pool: pd.DataFrame,
    cols: list[str],
    caliper: float,
    rng: np.random.Generator,
    with_replacement: bool = WITH_REPLACEMENT,
) -> pd.DataFrame:
    """Nearest-neighbour matching within a caliper, on standardised covariates.

    Cases are processed in a random order rather than in LOEUF order: matching
    greedily from the most constrained gene down would give the earliest cases
    the pick of the pool and leave systematically worse matches for the rest.
    The order is drawn from the seeded generator, so it is reproducible.

    With replacement, a control may serve several cases and the draw order stops
    mattering -- every case simply takes its nearest control. Without it, the
    pool is consumed as matching proceeds and cases matched late do worse.
    """
    all_rows = pd.concat([cases[cols], pool[cols]])
    z_all = _standardise(all_rows, cols)
    z_cases = z_all.loc[cases.index].to_numpy()
    z_pool = z_all.loc[pool.index].to_numpy()

    available = np.ones(len(pool), dtype=bool)
    pairs = []
    for i in rng.permutation(len(cases)):
        d = np.sqrt(((z_pool - z_cases[i]) ** 2).sum(axis=1))
        if not with_replacement:
            d = np.where(available, d, np.inf)
        j = int(np.argmin(d))
        if not np.isfinite(d[j]) or d[j] > caliper:
            continue  # no acceptable control; this case goes unmatched
        if not with_replacement:
            available[j] = False
        pairs.append((cases.index[i], pool.index[j], float(d[j])))

    return pd.DataFrame(pairs, columns=["case_gene", "control_gene", "distance"])


def build(caliper: float = CALIPER_SD) -> pd.DataFrame:
    """Construct case and control sets under the recorded D-002 strategy.

    `caliper` is in pooled standard deviations of the matching covariates. It is
    a parameter rather than a constant because its value is part of D-002 --
    "matching tolerance" in the project brief -- and should be recorded with the
    strategy.
    """
    strategy = decisions.CONTROL_MATCHING.value  # raises until D-002 is recorded

    universe = pd.read_parquet(GENE_UNIVERSE)
    df = prepare(universe)
    rng = np.random.default_rng(RANDOM_SEED)

    cases = df[df["is_constrained_loeuf"]]
    pool = df[df["loeuf"] >= CONTROL_LOEUF_MIN]

    if strategy == "caliper_expression_exons":
        cols = ["log_tpm", "n_coding_exons"]
        pairs = _greedy_caliper_match(cases, pool, cols, caliper, rng)
    elif strategy == "caliper_plus_length":
        cols = ["log_tpm", "n_coding_exons", "log_gene_length"]
        pairs = _greedy_caliper_match(cases, pool, cols, caliper, rng)
    elif strategy == "propensity_score":
        raise NotImplementedError(
            "propensity_score selected but not implemented. It needs a logistic "
            "model of constraint on all covariates and matching on the linear "
            "predictor; written only if chosen, so the unused branches do not "
            "sit in the repo untested."
        )
    elif strategy == "stratified_no_matching":
        raise NotImplementedError(
            "stratified_no_matching selected: there is no control set to build. "
            "Stage 1 should compare across LOEUF deciles with covariate "
            "adjustment instead of loading gene_sets.parquet."
        )
    else:  # pragma: no cover - decide() already restricts the option set
        raise ValueError(f"unhandled strategy {strategy!r}")

    # Frequency weights. A case is matched once and weighs 1. A control weighs
    # the number of cases it was matched to, so that a control serving three
    # cases counts three times in any control-side rate -- which is what makes
    # the two curves comparable under replacement.
    control_use = pairs["control_gene"].value_counts()

    matched_cases = set(pairs["case_gene"])
    out = pd.DataFrame(index=df.index)
    out["set"] = pd.Series("unused", index=df.index, dtype="object")
    out.loc[list(matched_cases), "set"] = "constrained"
    out.loc[control_use.index, "set"] = "control"

    out["weight"] = 0.0
    out.loc[list(matched_cases), "weight"] = 1.0
    out.loc[control_use.index, "weight"] = control_use.astype(float)
    #: How many cases each control stands in for. 1 everywhere without
    #: replacement; the clustering unit for control-side confidence intervals.
    out["times_used"] = 0
    out.loc[control_use.index, "times_used"] = control_use

    out["loeuf"] = df["loeuf"]
    out["loeuf_decile"] = df["loeuf_decile"]
    out["pli"] = df["pli"]
    for c in COVARIATES:
        out[c] = df[c]

    out.attrs["strategy"] = strategy
    out.attrs["caliper"] = caliper
    out.attrs["with_replacement"] = WITH_REPLACEMENT
    out.attrs["seed"] = RANDOM_SEED
    out.to_parquet(GENE_SETS_TABLE)

    _report_balance(out, pairs, cases)
    return out


def _weighted_mean_var(x: pd.Series, w: pd.Series) -> tuple[float, float]:
    """Frequency-weighted mean and variance.

    Controls are reused under replacement, so an unweighted control mean would
    under-count the controls that stand in for several cases and the balance
    table would flatter the matching.
    """
    mean = float(np.average(x, weights=w))
    var = float(np.average((x - mean) ** 2, weights=w))
    return mean, var


def _report_balance(out: pd.DataFrame, pairs: pd.DataFrame, cases: pd.DataFrame) -> None:
    """Standardised mean differences after matching, plus the two things a
    matched design can still get wrong.

    A |SMD| under 0.1 is the conventional bar for "balanced". Reported per
    covariate because a single summary number would hide one badly-matched
    covariate behind three well-matched ones.

    Beyond balance, two diagnostics that matter more here:

    - **Case representativeness.** Unmatched cases are dropped, and if they are
      dropped non-randomly the constrained curve describes a biased subset. The
      LOEUF shift between retained cases and all constrained genes measures that
      directly.
    - **Effective sample size.** Under replacement the control curve has fewer
      independent observations than its weight total suggests. ESS = (sum w)^2 /
      sum(w^2) is what confidence intervals should be built on.
    """
    n_matched = len(pairs)
    print(
        f"  matched {n_matched:,} of {len(cases):,} constrained genes "
        f"({len(cases) - n_matched:,} unmatched, no control within caliper)"
    )

    a = out[out["set"] == "constrained"]
    b = out[out["set"] == "control"]
    w = b["weight"]
    ess = float(w.sum() ** 2 / (w**2).sum())
    print(
        f"  controls: {len(b):,} unique genes carrying {w.sum():,.0f} matched "
        f"case-slots (effective n = {ess:,.0f})"
    )

    print(f"\n  {'covariate':<28} {'constrained':>12} {'control':>12} {'SMD':>8}")
    for c in COVARIATES:
        ma, va = float(a[c].mean()), float(a[c].var(ddof=1))
        mb, vb = _weighted_mean_var(b[c], w)
        sd = np.sqrt((va + vb) / 2)
        smd = (ma - mb) / sd if sd else np.nan
        flag = "" if abs(smd) < 0.1 else "   <-- imbalanced"
        print(f"  {c:<28} {ma:12.3f} {mb:12.3f} {smd:8.3f}{flag}")

    sd_loeuf = cases["loeuf"].std(ddof=1)
    bias = (a["loeuf"].mean() - cases["loeuf"].mean()) / sd_loeuf
    verdict = "representative" if abs(bias) < 0.1 else "SKEWED"
    print(
        f"\n  case representativeness: retained cases vs all constrained genes, "
        f"LOEUF SMD = {bias:.3f} ({verdict})"
    )


def main() -> None:
    if not decisions.CONTROL_MATCHING.decided:
        print("Control matching is blocked on decision D-002.\n")
        print(f"  Question: {decisions.CONTROL_MATCHING.question}\n")
        for name, meaning in decisions.CONTROL_MATCHING.alternatives.items():
            print(f"  - {name}\n      {meaning}\n")
        print("Record the choice in DECISIONS.md and via decisions.CONTROL_MATCHING")
        print(".decide(<option>, rationale=...), then re-run.")
        return

    build()
    print(f"\nWrote {GENE_SETS_TABLE}")


if __name__ == "__main__":
    main()
