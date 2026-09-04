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
) -> pd.DataFrame:
    """One-to-one nearest-neighbour matching without replacement.

    Cases are processed in a random order rather than in LOEUF order: matching
    greedily from the most constrained gene down would give the earliest cases
    the pick of the pool and leave systematically worse matches for the rest.
    The order is drawn from the seeded generator, so it is reproducible.
    """
    all_rows = pd.concat([cases[cols], pool[cols]])
    z_all = _standardise(all_rows, cols)
    z_cases = z_all.loc[cases.index].to_numpy()
    z_pool = z_all.loc[pool.index].to_numpy()

    available = np.ones(len(pool), dtype=bool)
    pairs = []
    for i in rng.permutation(len(cases)):
        d = np.sqrt(((z_pool - z_cases[i]) ** 2).sum(axis=1))
        d[~available] = np.inf
        j = int(np.argmin(d))
        if not np.isfinite(d[j]) or d[j] > caliper:
            continue  # no acceptable control; this case goes unmatched
        available[j] = False
        pairs.append((cases.index[i], pool.index[j], float(d[j])))

    return pd.DataFrame(pairs, columns=["case_gene", "control_gene", "distance"])


def build(caliper: float = 0.25) -> pd.DataFrame:
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

    matched_cases = set(pairs["case_gene"])
    out = pd.DataFrame(index=df.index)
    out["set"] = pd.Series("unused", index=df.index, dtype="object")
    out.loc[list(matched_cases), "set"] = "constrained"
    out.loc[list(pairs["control_gene"]), "set"] = "control"
    out["loeuf"] = df["loeuf"]
    out["loeuf_decile"] = df["loeuf_decile"]
    out["pli"] = df["pli"]
    for c in COVARIATES:
        out[c] = df[c]

    out.attrs["strategy"] = strategy
    out.attrs["caliper"] = caliper
    out.attrs["seed"] = RANDOM_SEED
    out.to_parquet(GENE_SETS_TABLE)

    _report_balance(out, pairs, cases)
    return out


def _report_balance(out: pd.DataFrame, pairs: pd.DataFrame, cases: pd.DataFrame) -> None:
    """Standardised mean differences after matching.

    A |SMD| under 0.1 is the conventional bar for "balanced". Reported per
    covariate because a single summary number would hide one badly-matched
    covariate behind three well-matched ones.
    """
    print(
        f"  matched {len(pairs):,} of {len(cases):,} constrained genes "
        f"({len(cases) - len(pairs):,} unmatched, no control within caliper)"
    )
    a = out[out["set"] == "constrained"]
    b = out[out["set"] == "control"]
    print(f"\n  {'covariate':<28} {'constrained':>12} {'control':>12} {'SMD':>8}")
    for c in COVARIATES:
        ma, mb = a[c].mean(), b[c].mean()
        sd = np.sqrt((a[c].var(ddof=1) + b[c].var(ddof=1)) / 2)
        smd = (ma - mb) / sd if sd else np.nan
        flag = "" if abs(smd) < 0.1 else "   <-- imbalanced"
        print(f"  {c:<28} {ma:12.3f} {mb:12.3f} {smd:8.3f}{flag}")


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
