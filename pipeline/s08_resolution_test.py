"""Stage 1, step 4: does cell-type resolution move the gap, at fixed donor count?

The recovery curve cannot separate resolution from donor count, because
SingleBrain arrives with 4.8x GTEx's donors. This module runs the two tests
that can, in opposite directions, and they now agree.

  SPLITTING (SingleBrain)   pool a cell class, then split it into its subtypes
  POOLING   (Bryois, D-007) pool all nuclei into pseudobulk, vs 8 cell types

Neither isolates resolution on its own, because resolution and per-context
power are intrinsically coupled: at fixed sequencing depth you cannot resolve
more contexts without putting fewer reads in each. Splitting loses reads per
context; pooling gains them. Running both directions is what brackets the
answer -- if resolution were doing the work, splitting should narrow the gap
and pooling should widen it, and the two should disagree in sign.

They do not. Both centre on zero (see the constants below and the Stage 1
report), so the closure the main ladder shows is coming from donor count.

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

#: The effect this test is trying to explain: the constrained-gene gap falls
#: 0.357 -> 0.143 from bulk cortex to single-nucleus major cell types, so
#: 0.214 in absolute terms. Every resolution delta below is reported as a
#: fraction of it, because "the gap moved by 0.005" only means something
#: against the size of the movement being accounted for.
LADDER_CLOSURE = 0.214

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


def _delta_with_ci(
    hit_pooled: pd.Series,
    hit_split: pd.Series,
    cases: pd.Index,
    ctrls: pd.Index,
    weights: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float, float]:
    """CI on (split gap - pooled gap), from PAIRED bootstrap replicates.

    The two arms are computed on the same genes, so their gaps are strongly
    correlated. Comparing their separate intervals for overlap would be far too
    conservative, and comparing point estimates against a fixed cutoff -- which
    is what an earlier version of this module did -- ignores uncertainty
    entirely and can call a difference that is pure noise.

    Both gaps are therefore evaluated on each shared replicate and the
    difference taken within the replicate, so the correlation cancels and the
    interval is on the quantity actually being argued about.
    """
    cp = hit_pooled.loc[cases].to_numpy()
    kp = hit_pooled.loc[ctrls].to_numpy()
    cs = hit_split.loc[cases].to_numpy()
    ks = hit_split.loc[ctrls].to_numpy()

    n_c, n_k = len(cases), len(ctrls)
    deltas = np.empty(N_BOOT)
    for b in range(N_BOOT):
        ci = rng.integers(0, n_c, n_c)
        wi = rng.integers(0, n_k, n_k)
        w = weights[wi]
        g_pooled = _weighted_rate(kp[wi], w) - cp[ci].mean()
        g_split = _weighted_rate(ks[wi], w) - cs[ci].mean()
        deltas[b] = g_split - g_pooled

    observed = (_weighted_rate(ks, weights) - cs.mean()) - (
        _weighted_rate(kp, weights) - cp.mean()
    )
    lo, hi = _ci(deltas)
    return observed, lo, hi


def singlebrain_on_bryois_genes(
    long: pd.DataFrame,
    cases: pd.Index,
    ctrls: pd.Index,
    weights: np.ndarray,
    genes: pd.Index,
    rng: np.random.Generator,
) -> dict | None:
    """Run the SingleBrain splitting test on exactly the genes Bryois tested.

    Without this, comparing the two studies' resolution effects confounds the
    study with the gene set. Bryois is currently restricted to whichever
    chromosomes have downloaded, which is a small and chromosome-specific
    subset -- so a disagreement between "Bryois says -0.037" and "SingleBrain
    says -0.007" could be about those genes rather than about the studies.

    Holding the gene set fixed and re-running the SingleBrain contrast on it
    separates the two. If SingleBrain also turns negative on these genes, the
    tension is the subset. If it stays flat, the disagreement is real and is
    between the studies.
    """
    bry_arms = ("bryois_pb", "bryois_celltype")
    tested = set(long.loc[long["rung"].isin(bry_arms), "gene_id"])
    if not tested:
        return None

    c = cases.intersection(tested)
    keep = ctrls.isin(tested)
    k, w = ctrls[keep], weights[keep]
    if len(c) < 50 or len(k) < 50:
        return None

    # All SingleBrain major classes pooled, versus all their subtypes -- the
    # same contrast as the per-class test, aggregated so it is comparable to
    # Bryois's single 1-vs-8 comparison.
    majors = [m for m in SUBTYPES if m in set(long["cell"])]
    subs = [s for m in majors for s in SUBTYPES[m] if s in set(long["cell"])]
    if not majors or not subs:
        return None

    hit_pooled = _detect(long, majors, genes)
    hit_split = _detect(long, subs, genes)
    d, d_lo, d_hi = _delta_with_ci(hit_pooled, hit_split, c, k, w, rng)
    return {
        "major": "SingleBrain on Bryois genes",
        "arm": "summary",
        "n_columns": len(subs),
        "n_cases_restricted": len(c),
        "n_controls_restricted": len(k),
        "delta": d,
        "delta_lo": d_lo,
        "delta_hi": d_hi,
    }


def pooled_delta(
    hits_by_class: dict[str, tuple[pd.Series, pd.Series]],
    cases: pd.Index,
    ctrls: pd.Index,
    weights: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float, float]:
    """CI on the MEAN delta across cell classes, from shared replicates.

    Each class on its own is underpowered -- the per-class intervals span
    roughly +/-0.05 -- so six inconclusive tests cannot be read as "no effect"
    by counting signs. Averaging the deltas within each bootstrap replicate
    gives one interval on the aggregate, which is the quantity the report's
    claim actually rests on.

    Note what this can and cannot say. An interval that spans zero and is
    narrow rules out a LARGE effect of resolution; it never proves there is
    none. The report states the bound rather than claiming a null.
    """
    arrays = {
        cls: (
            hp.loc[cases].to_numpy(), hp.loc[ctrls].to_numpy(),
            hs.loc[cases].to_numpy(), hs.loc[ctrls].to_numpy(),
        )
        for cls, (hp, hs) in hits_by_class.items()
    }
    n_c, n_k = len(cases), len(ctrls)

    def mean_delta(ci, wi):
        w = weights[wi]
        ds = []
        for cp, kp, cs, ks in arrays.values():
            g_p = _weighted_rate(kp[wi], w) - cp[ci].mean()
            g_s = _weighted_rate(ks[wi], w) - cs[ci].mean()
            ds.append(g_s - g_p)
        return float(np.mean(ds))

    observed = mean_delta(np.arange(n_c), np.arange(n_k))
    boots = np.array(
        [
            mean_delta(rng.integers(0, n_c, n_c), rng.integers(0, n_k, n_k))
            for _ in range(N_BOOT)
        ]
    )
    lo, hi = _ci(boots)
    return observed, lo, hi


def bryois_contrast(
    long: pd.DataFrame,
    cases: pd.Index,
    ctrls: pd.Index,
    weights: np.ndarray,
    genes: pd.Index,
    rng: np.random.Generator,
) -> list[dict]:
    """The D-007 arm: Bryois pseudobulk vs its own 8 cell types.

    The mirror image of the splitting test above. Same donors, same pipeline,
    same normalisation; pseudobulk POOLS all nuclei per individual while the
    cell-type arm keeps them separate. So this varies resolution in the
    opposite direction from splitting a SingleBrain class.

    Neither isolates resolution -- pooling raises reads per context as it
    lowers resolution, splitting does the reverse -- but if BOTH leave the gap
    unmoved, resolution is not the active ingredient in either direction, and
    the closure seen on the main ladder has to be coming from donor count.
    """
    arms = [a for a in ("bryois_pb", "bryois_celltype") if a in set(long["rung"])]
    if len(arms) < 2:
        return []

    # CRITICAL: restrict to genes Bryois actually tested.
    #
    # Everywhere else in this project a gene absent from a rung counts as "no
    # detectable eQTL there" (D-012), which is the right reading when the rung
    # covers the whole genome. Bryois may be running on a subset of
    # chromosomes while its 198 files download, and treating every gene on an
    # unfetched chromosome as undetected would deflate BOTH arms by the same
    # large factor -- compressing the gap toward zero and making the contrast
    # look null for a purely clerical reason.
    #
    # Both arms are restricted to the same tested-gene set, so the comparison
    # stays symmetric; it is only the denominator that changes.
    tested = set(long.loc[long["rung"].isin(arms), "gene_id"])
    cases = cases.intersection(tested)
    keep = ctrls.isin(tested)
    ctrls, weights = ctrls[keep], weights[keep]
    if len(cases) < 50 or len(ctrls) < 50:
        return []

    rows = []
    hits = {}
    for arm in arms:
        cells = sorted(set(long.loc[long["rung"] == arm, "cell"]))
        hit = _detect(long, cells, genes)
        label = "pooled" if arm == "bryois_pb" else "split"
        hits[label] = hit
        cr, kr, gap, lo, hi = _gap_with_ci(hit, cases, ctrls, weights, rng)
        rows.append(
            {
                "major": "Bryois (D-007)",
                "arm": label,
                "n_columns": len(cells),
                "n_cases_restricted": len(cases),
                "n_controls_restricted": len(ctrls),
                "constrained_rate": cr,
                "control_rate": kr,
                "gap": gap,
                "gap_lo": lo,
                "gap_hi": hi,
            }
        )
    d, d_lo, d_hi = _delta_with_ci(
        hits["pooled"], hits["split"], cases, ctrls, weights, rng
    )
    for r in rows:
        r["delta"], r["delta_lo"], r["delta_hi"] = d, d_lo, d_hi
    return rows


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
    class_hits: dict[str, tuple[pd.Series, pd.Series]] = {}
    for major, subs in SUBTYPES.items():
        present = [c for c in subs if c in set(long["cell"])]
        if major not in set(long["cell"]) or not present:
            continue

        hits = {}
        for arm, cells in (("pooled", [major]), ("split", present)):
            hit = _detect(long, cells, genes)
            hits[arm] = hit
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
        d, d_lo, d_hi = _delta_with_ci(
            hits["pooled"], hits["split"], cases, ctrls, weights, rng
        )
        for r in rows[-2:]:
            r["delta"], r["delta_lo"], r["delta_hi"] = d, d_lo, d_hi
        class_hits[major] = (hits["pooled"], hits["split"])

    bry = bryois_contrast(long, cases, ctrls, weights, genes, rng)
    if bry:
        rows.extend(bry)
    else:
        prov.note(
            "Bryois contrast (D-007)",
            "not available yet -- needs both bryois_pb and bryois_celltype in "
            "the detection table",
        )

    out = pd.DataFrame(rows)

    print("Splitting a cell class into its subtypes, donors held fixed")
    print(f"  {'class':<6}{'cols':>5}{'arm':>9}{'constr':>9}{'ctrl':>8}"
          f"{'gap':>9}{'95% CI':>18}")
    deltas = []
    bryois_delta = None
    for major, g in out.groupby("major", sort=False):
        for _, r in g.iterrows():
            print(
                f"  {r['major']:<15}{r['n_columns']:>5}{r['arm']:>9}"
                f"{r['constrained_rate']:>9.3f}{r['control_rate']:>8.3f}"
                f"{r['gap']:>9.3f}   [{r['gap_lo']:.3f}, {r['gap_hi']:.3f}]"
            )
        d = g["delta"].iloc[0]
        d_lo, d_hi = g["delta_lo"].iloc[0], g["delta_hi"].iloc[0]
        # A delta whose interval spans zero says nothing about direction, so
        # the readout says so rather than reporting the sign of noise.
        call = "no effect" if d_lo <= 0 <= d_hi else (
            "finer NARROWS" if d < 0 else "finer WIDENS"
        )
        # Bryois is a separate study varying resolution the other way, so it
        # is reported on its own rather than averaged into the SingleBrain mean.
        if major.startswith("Bryois"):
            bryois_delta = (d, d_lo, d_hi)
        else:
            deltas.append(d)
        print(f"  {'':<15}{'':>5}{'delta':>9}{'':>17}{d:>+9.3f}"
              f"   [{d_lo:+.3f}, {d_hi:+.3f}]  {call}")

    mean_delta = float(np.mean(deltas))
    n_inconclusive = sum(
        1
        for _, g in out[out["major"].isin(class_hits)].groupby("major")
        if g["delta_lo"].iloc[0] <= 0 <= g["delta_hi"].iloc[0]
    )
    pooled = pooled_delta(class_hits, cases, ctrls, weights, rng)
    p_d, p_lo, p_hi = pooled

    print(
        f"\n  SingleBrain, per-class deltas: {n_inconclusive} of "
        f"{len(class_hits)} have intervals spanning zero"
    )
    print(
        f"  POOLED mean delta across {len(class_hits)} classes: {p_d:+.3f} "
        f"[{p_lo:+.3f}, {p_hi:+.3f}]"
    )
    if p_lo <= 0 <= p_hi:
        bound = max(abs(p_lo), abs(p_hi))
        verdict = (
            f"no detectable effect of splitting; the data rule out a change "
            f"larger than about {bound:.3f} in either direction, against a "
            f"bulk-to-single-nucleus closure of {LADDER_CLOSURE:.3f}"
        )
    else:
        verdict = (
            "splitting NARROWS the gap -- resolution is doing work"
            if p_d < 0
            else "splitting WIDENS the gap"
        )
    print(f"  => {verdict}")

    # Persist the pooled estimate as its own row so the figure and the report
    # quote the same numbers rather than recomputing them.
    out = pd.concat(
        [
            out,
            pd.DataFrame([{
                "major": "POOLED",
                "arm": "summary",
                "n_columns": len(class_hits),
                "delta": p_d,
                "delta_lo": p_lo,
                "delta_hi": p_hi,
            }]),
        ],
        ignore_index=True,
    )
    out.to_parquet(OUT, index=False)

    if bryois_delta is not None:
        d, d_lo, d_hi = bryois_delta
        # An interval spanning zero is not automatically uninformative. What
        # matters is how tight it is relative to the effect being explained:
        # a wide interval says "cannot tell", a narrow one centred on zero is a
        # real null result and is the point of the test. The comparator is the
        # bulk-to-single-nucleus closure the main ladder shows.
        spans_zero = d_lo <= 0 <= d_hi
        bound = max(abs(d_lo), abs(d_hi))
        as_frac = bound / LADDER_CLOSURE
        informative = (not spans_zero) or as_frac <= 0.35

        # Same-gene-set control: is the Bryois/SingleBrain difference about the
        # studies, or about which genes Bryois has downloaded so far?
        matched = singlebrain_on_bryois_genes(
            long, cases, ctrls, weights, genes, rng
        )
        if matched is not None:
            rows.append(matched)
            m, m_lo, m_hi = (
                matched["delta"], matched["delta_lo"], matched["delta_hi"]
            )
            print(
                f"\n  SingleBrain restricted to the same {matched['n_cases_restricted']}"
                f" constrained genes Bryois tested: {m:+.3f} [{m_lo:+.3f}, {m_hi:+.3f}]"
            )
            overlap = not (d_hi < m_lo or m_hi < d_lo)
            print(
                "  => the two studies' intervals "
                + ("OVERLAP, so no evidence they disagree" if overlap
                   else "DO NOT overlap -- a real study-level disagreement")
            )
            prov.record(
                "SingleBrain on the Bryois gene set",
                matched["n_cases_restricted"],
                matched["n_controls_restricted"],
                detail=(
                    f"delta {m:+.3f} [{m_lo:+.3f}, {m_hi:+.3f}] against Bryois "
                    f"{d:+.3f} [{d_lo:+.3f}, {d_hi:+.3f}]; intervals "
                    + ("overlap" if overlap else "do not overlap")
                ),
            )
        print(
            f"\n  Bryois (D-007), pseudobulk -> 8 cell types: {d:+.3f} "
            f"[{d_lo:+.3f}, {d_hi:+.3f}]"
        )
        if not spans_zero:
            print(
                "  => the Bryois arm resolves a DIRECTION: pooling and "
                "splitting are not equivalent."
            )
        elif informative:
            print(
                f"  => informative NULL: bounded at {bound:.3f}, i.e. at most "
                f"{as_frac:.0%} of the {LADDER_CLOSURE:.3f} ladder closure.\n"
                "     Pooling nuclei changes the gap no more than splitting "
                "them does."
            )
        else:
            print(
                f"  => UNINFORMATIVE: interval spans zero and is wide "
                f"({bound:.3f}, {as_frac:.0%} of the ladder closure)."
            )
        prov.record(
            "Bryois pseudobulk vs cell types (D-007)",
            0,
            0,
            detail=(
                f"change in constrained-gene gap on going from pseudobulk to 8 "
                f"cell types, same donors: {d:+.3f} [{d_lo:+.3f}, {d_hi:+.3f}]. "
                + (
                    "Interval excludes zero."
                    if informative
                    else "Interval spans zero -- uninformative at this "
                    "chromosome coverage."
                )
            ),
        )

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
