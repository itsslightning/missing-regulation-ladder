"""Stage 1, step 4: does splitting a cell type help, holding everything else fixed?

The recovery curve cannot separate resolution from donor count, because
SingleBrain arrives with 4.8x GTEx's donors. The Bryois pseudobulk arm (D-007)
is the decisive test and is still downloading. This is the tightest test
available in the meantime, and it is a genuinely clean one.

SingleBrain reports each major cell class BOTH as one class and as its
constituent subtypes:

    Ast  <->  Ast1 Ast2 Ast3 Ast4
    Ext  <->  Ext1 ... Ext8
    IN   <->  IN1 ... IN7
    MG   <->  MG1 ... MG4
    OD   <->  OD1 OD2 OD3
    OPC  <->  OPC1 OPC2

Same donors, same study, same pipeline, same normalisation, same nuclei -- the
ONLY difference is whether those nuclei were pooled into one class or split.
So this isolates cell-type resolution with sample size held exactly fixed,
which is the comparison the main ladder cannot make.

If resolution is what closes the constrained-gene gap, splitting should close
it further. If the gap is unmoved or widens, the closure seen on the main
ladder is coming from donor count, not from resolution.

One asymmetry is handled explicitly: the subtype arm gets several tests per
gene where the pooled arm gets one, so BH is applied within each arm over that
arm's own gene x column tests. Without that the subtype arm would win on
multiplicity alone. End is excluded -- it has no subtypes.

Writes: data/processed/resolution_test.parquet
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from pipeline import config as cfg
from pipeline.provenance import Provenance
from pipeline.s05_detection import DETECTION_LONG, FDR_ALPHA
from pipeline.s06_recovery import N_BOOT, _ci, _weighted_rate

OUT = cfg.DIR_PROCESSED / "resolution_test.parquet"

#: Major class -> its subtypes, as SingleBrain names them. End has none.
SUBTYPES = {
    "Ast": ["Ast1", "Ast2", "Ast3", "Ast4"],
    "Ext": [f"Ext{i}" for i in range(1, 9)],
    "IN": [f"IN{i}" for i in range(1, 8)],
    "MG": ["MG1", "MG2", "MG3", "MG4"],
    "OD": ["OD1", "OD2", "OD3"],
    "OPC": ["OPC1", "OPC2"],
}


def _detect(long: pd.DataFrame, cells: list[str], genes: pd.Index) -> pd.Series:
    """BH within this arm's own tests, then 'detected in any listed cell type'."""
    sub = long[long["cell"].isin(cells)].copy()
    ok = sub["p_bonf"].notna()
    q = pd.Series(np.nan, index=sub.index)
    q[ok] = multipletests(sub.loc[ok, "p_bonf"], method="fdr_bh")[1]
    sub["hit"] = q <= FDR_ALPHA
    any_hit = sub.groupby("gene_id")["hit"].any()
    return any_hit.reindex(genes).fillna(False).astype(bool)


def _gap_with_ci(
    hit: pd.Series,
    cases: pd.Index,
    ctrls: pd.Index,
    weights: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float, float, float, float]:
    case_d = hit.loc[cases].to_numpy()
    ctrl_d = hit.loc[ctrls].to_numpy()
    case_rate = float(case_d.mean())
    ctrl_rate = _weighted_rate(ctrl_d, weights)

    n_c, n_k = len(case_d), len(ctrl_d)
    boots = np.empty(N_BOOT)
    for b in range(N_BOOT):
        ci = rng.integers(0, n_c, n_c)
        wi = rng.integers(0, n_k, n_k)
        boots[b] = _weighted_rate(ctrl_d[wi], weights[wi]) - case_d[ci].mean()
    lo, hi = _ci(boots)
    return case_rate, ctrl_rate, ctrl_rate - case_rate, lo, hi


def main() -> None:
    prov = Provenance("s08_resolution_test")
    rng = np.random.default_rng(cfg.RANDOM_SEED)

    long = pd.read_parquet(DETECTION_LONG)
    gene_sets = pd.read_parquet(cfg.GENE_SETS_TABLE)
    gene_sets = gene_sets[gene_sets["set"].isin(["constrained", "control"])]
    cases = gene_sets.index[gene_sets["set"] == "constrained"]
    ctrls = gene_sets.index[gene_sets["set"] == "control"]
    weights = gene_sets.loc[ctrls, "weight"].to_numpy(dtype=float)
    genes = gene_sets.index

    rows = []
    for major, subs in SUBTYPES.items():
        present = [c for c in subs if c in set(long["cell"])]
        if major not in set(long["cell"]) or not present:
            continue

        for arm, cells in (("pooled", [major]), ("split", present)):
            hit = _detect(long, cells, genes)
            cr, kr, gap, lo, hi = _gap_with_ci(hit, cases, ctrls, weights, rng)
            rows.append(
                {
                    "major": major,
                    "arm": arm,
                    "n_columns": len(cells),
                    "constrained_rate": cr,
                    "control_rate": kr,
                    "gap": gap,
                    "gap_lo": lo,
                    "gap_hi": hi,
                }
            )

    out = pd.DataFrame(rows)
    out.to_parquet(OUT, index=False)

    print("Splitting a cell class into its subtypes, donors held fixed")
    print(f"  {'class':<6}{'cols':>5}{'arm':>9}{'constr':>9}{'ctrl':>8}"
          f"{'gap':>9}{'95% CI':>18}")
    deltas = []
    for major, g in out.groupby("major", sort=False):
        for _, r in g.iterrows():
            print(
                f"  {r['major']:<6}{r['n_columns']:>5}{r['arm']:>9}"
                f"{r['constrained_rate']:>9.3f}{r['control_rate']:>8.3f}"
                f"{r['gap']:>9.3f}   [{r['gap_lo']:.3f}, {r['gap_hi']:.3f}]"
            )
        p = g[g["arm"] == "pooled"]["gap"].iloc[0]
        s = g[g["arm"] == "split"]["gap"].iloc[0]
        deltas.append(s - p)
        print(f"  {'':<6}{'':>5}{'change':>9}{'':>17}{s - p:>+9.3f}"
              f"   {'(split narrows)' if s < p else '(split widens)'}")

    mean_delta = float(np.mean(deltas))
    n_narrow = sum(d < 0 for d in deltas)
    verdict = (
        "splitting narrows the gap -- resolution is doing work"
        if mean_delta < -0.01
        else "splitting does NOT narrow the gap -- resolution is not the "
        "active ingredient"
    )
    print(f"\n  mean change in gap on splitting: {mean_delta:+.3f}")
    print(f"  classes where splitting narrowed the gap: {n_narrow} of {len(deltas)}")
    print(f"  => {verdict}")

    prov.record(
        "within-SingleBrain resolution test",
        len(deltas),
        n_narrow,
        detail=(
            f"mean change in constrained-gene gap on splitting a class into "
            f"subtypes: {mean_delta:+.3f}; narrowed in {n_narrow} of "
            f"{len(deltas)} classes. Donors, pipeline and normalisation are "
            f"identical between arms."
        ),
    )
    prov.save(cfg.DIR_LOGS / "s08_resolution_test.json")
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
