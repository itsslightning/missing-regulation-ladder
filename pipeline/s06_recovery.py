"""Stage 1, step 2: the recovery curve, with confidence intervals.

For each rung, what fraction of constrained genes have a detectable cis-eQTL,
against the same fraction for matched unconstrained controls? The difference
between those two curves is the "gap", and how the gap behaves across rungs is
the whole experiment:

  gap roughly flat across rungs      -> leans H1, selection
  gap narrows toward zero            -> leans H2, power/resolution
  gap narrows partially              -> the expected answer; report the split

WHY THE CONTROL ARM NEEDS DIFFERENT STATISTICS
----------------------------------------------
D-002 matched with replacement, so 1,294 unique control genes serve 2,767
matched case-slots. A control used seven times is one observation, not seven,
and a naive binomial interval over 2,767 control slots would be far too narrow.

So the two arms are treated differently and deliberately:

  constrained   2,767 independent genes -> Wilson interval
  control       frequency-weighted rate, CI from a cluster bootstrap that
                resamples the 1,294 unique control GENES (carrying their
                weights), not the slots

The gap and its interval come from resampling both arms jointly, which is what
makes the gap interval honest about the reuse.

All bootstraps draw from a generator seeded with RANDOM_SEED, so intervals are
reproducible to the digit.

Writes: data/processed/recovery_by_rung.parquet
        data/processed/recovery_by_loeuf.parquet
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint

from pipeline import config as cfg
from pipeline.provenance import Provenance
from pipeline.s01_gene_universe import GENE_UNIVERSE
from pipeline.s04_schema_gene_sets import SCHEMA_SETS
from pipeline.s05_detection import DETECTION_BY_RUNG

RECOVERY_BY_LOEUF = cfg.DIR_PROCESSED / "recovery_by_loeuf.parquet"
#: The same, holding expression roughly fixed. The raw decile view is confounded
#: badly enough to invert the result, so this is what should be read.
RECOVERY_BY_LOEUF_EXPR = cfg.DIR_PROCESSED / "recovery_by_loeuf_expr.parquet"

#: Bootstrap replicates for the control arm and the gap. 2,000 is enough for a
#: stable 95% percentile interval and cheap at this table size.
N_BOOT = 2000

#: Rung display order. This is the x-axis of the curve, so it is fixed here
#: rather than inferred from whatever happens to be in the table.
RUNG_ORDER = ["gtex_cortex", "bulk_brain", "sn_major", "sn_subtype"]
RUNG_LABELS = {
    "gtex_cortex": "GTEx cortex\n(bulk tissue)",
    "bulk_brain": "PsychENCODE\n(bulk brain)",
    "sn_major": "SingleBrain\n7 major types",
    "sn_subtype": "SingleBrain\n28 subtypes",
    "bryois_pb": "Bryois\npseudobulk",
    "bryois_celltype": "Bryois\n8 cell types",
}

#: Donor counts, carried so the curve can be plotted against N as well as rung
#: index. The report names cross-study power as the main technical risk, so a
#: recovery curve without these numbers beside it is not interpretable.
RUNG_N_DONORS = {
    "gtex_cortex": 205,
    "bulk_brain": 1387,
    "sn_major": 983,
    "sn_subtype": 983,
    "bryois_pb": 192,
    "bryois_celltype": 192,
}

#: Which detection column each arm reads. The primary is D-001 as recorded;
#: the others are the sensitivity arms it was chosen over.
ARMS = {
    "uniform": "detected",
    "native": "detected_native",
    "effect": "detected_effect",
}


def _wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
    return float(lo), float(hi)


def _weighted_rate(detected: np.ndarray, weights: np.ndarray) -> float:
    total = weights.sum()
    return float((detected * weights).sum() / total) if total else np.nan


def _boot_control(
    detected: np.ndarray, weights: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """Cluster bootstrap over unique control genes, carrying their weights."""
    n = len(detected)
    idx = rng.integers(0, n, size=(N_BOOT, n))
    d = detected[idx]
    w = weights[idx]
    return (d * w).sum(axis=1) / w.sum(axis=1)


def _boot_cases(detected: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n = len(detected)
    idx = rng.integers(0, n, size=(N_BOOT, n))
    return detected[idx].mean(axis=1)


def _ci(samples: np.ndarray) -> tuple[float, float]:
    return tuple(np.percentile(samples, [2.5, 97.5]).astype(float))


def compute(
    gene_sets: pd.DataFrame,
    detection: pd.DataFrame,
    arm: str,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Recovery, gap and intervals for every rung under one detection arm."""
    col = ARMS[arm]
    det = detection.pivot(index="gene_id", columns="rung", values=col)

    cases = gene_sets[gene_sets["set"] == "constrained"]
    ctrls = gene_sets[gene_sets["set"] == "control"]

    rows = []
    for rung in [r for r in RUNG_ORDER + list(RUNG_LABELS) if r in det.columns]:
        if any(r["rung"] == rung for r in rows):
            continue
        # A gene absent from a rung's table was never tested there. Under D-012
        # that counts as "no detectable eQTL at this rung", which is the
        # honest reading: a gene with no eQTL findable in a cell type is
        # missing regulation there whether or not it passed an expression
        # filter.
        hit = det[rung].reindex(gene_sets.index).fillna(False).astype(bool)

        case_d = hit.loc[cases.index].to_numpy()
        ctrl_d = hit.loc[ctrls.index].to_numpy()
        ctrl_w = ctrls["weight"].to_numpy(dtype=float)

        case_rate = float(case_d.mean())
        ctrl_rate = _weighted_rate(ctrl_d, ctrl_w)
        case_lo, case_hi = _wilson(int(case_d.sum()), len(case_d))

        boot_ctrl = _boot_control(ctrl_d, ctrl_w, rng)
        boot_case = _boot_cases(case_d, rng)
        ctrl_lo, ctrl_hi = _ci(boot_ctrl)
        gap_boot = boot_ctrl - boot_case
        gap_lo, gap_hi = _ci(gap_boot)

        rows.append(
            {
                "rung": rung,
                "arm": arm,
                "n_donors": RUNG_N_DONORS.get(rung),
                "n_constrained": len(case_d),
                "n_control_genes": len(ctrl_d),
                "constrained_rate": case_rate,
                "constrained_lo": case_lo,
                "constrained_hi": case_hi,
                "control_rate": ctrl_rate,
                "control_lo": ctrl_lo,
                "control_hi": ctrl_hi,
                "gap": ctrl_rate - case_rate,
                "gap_lo": gap_lo,
                "gap_hi": gap_hi,
                # Two-sided bootstrap p for "gap != 0", used only for the
                # small family of rung-level contrasts (D-003).
                "gap_p": float(2 * min((gap_boot <= 0).mean(), (gap_boot >= 0).mean())),
            }
        )
    return pd.DataFrame(rows)


def by_loeuf(
    universe: pd.DataFrame, detection: pd.DataFrame, arm: str
) -> pd.DataFrame:
    """Recovery per LOEUF decile per rung, on the full universe.

    Stratification is on the whole universe rather than the matched sets: the
    matched design answers "constrained vs comparable unconstrained", while
    the decile view answers "how does recovery vary with constraint", and
    forcing the second through the matched sets would throw away eight deciles.
    """
    col = ARMS[arm]
    det = detection.pivot(index="gene_id", columns="rung", values=col)
    rows = []
    for rung in [r for r in RUNG_ORDER + list(RUNG_LABELS) if r in det.columns]:
        if any(r["rung"] == rung for r in rows):
            continue
        hit = det[rung].reindex(universe.index).fillna(False).astype(bool)
        for decile, idx in universe.groupby("loeuf_decile").groups.items():
            d = hit.loc[idx]
            k, n = int(d.sum()), len(d)
            lo, hi = _wilson(k, n)
            rows.append(
                {
                    "rung": rung,
                    "arm": arm,
                    "loeuf_decile": int(decile),
                    "n": n,
                    "k": k,
                    "rate": k / n if n else np.nan,
                    "lo": lo,
                    "hi": hi,
                }
            )
    return pd.DataFrame(rows)


def gap_closure(
    gene_sets: pd.DataFrame,
    detection: pd.DataFrame,
    arm: str,
    rng: np.random.Generator,
    from_rung: str = "gtex_cortex",
    to_rung: str = "sn_major",
) -> dict:
    """The headline decomposition: what fraction of the bulk gap closes?

    closure = (gap_from - gap_to) / gap_from

    Both gaps are computed on the SAME bootstrap replicate so the ratio's
    interval accounts for the correlation between them -- the two rungs share
    the same genes, so computing the intervals separately and dividing would
    overstate the uncertainty considerably.

    A closure of 1.0 means the gap vanished (pure H2, power/resolution). A
    closure of 0.0 means it did not move at all (pure H1, selection). The
    expected answer is in between, and the point of this project is the
    interval around it, not the point estimate.
    """
    col = ARMS[arm]
    det = detection.pivot(index="gene_id", columns="rung", values=col)
    if from_rung not in det.columns or to_rung not in det.columns:
        return {}

    cases = gene_sets[gene_sets["set"] == "constrained"]
    ctrls = gene_sets[gene_sets["set"] == "control"]
    ctrl_w = ctrls["weight"].to_numpy(dtype=float)

    hits = {}
    for rung in (from_rung, to_rung):
        h = det[rung].reindex(gene_sets.index).fillna(False).astype(bool)
        hits[rung] = (
            h.loc[cases.index].to_numpy(),
            h.loc[ctrls.index].to_numpy(),
        )

    def gap_of(case_d, ctrl_d, ci, wi):
        return _weighted_rate(ctrl_d[wi], ctrl_w[wi]) - case_d[ci].mean()

    n_case, n_ctrl = len(cases), len(ctrls)
    closures = np.empty(N_BOOT)
    for b in range(N_BOOT):
        ci = rng.integers(0, n_case, n_case)
        wi = rng.integers(0, n_ctrl, n_ctrl)
        g_from = gap_of(*hits[from_rung], ci, wi)
        g_to = gap_of(*hits[to_rung], ci, wi)
        closures[b] = (g_from - g_to) / g_from if g_from else np.nan

    g_from = gap_of(*hits[from_rung], np.arange(n_case), np.arange(n_ctrl))
    g_to = gap_of(*hits[to_rung], np.arange(n_case), np.arange(n_ctrl))
    lo, hi = _ci(closures[~np.isnan(closures)])
    return {
        "arm": arm,
        "from_rung": from_rung,
        "to_rung": to_rung,
        "gap_from": g_from,
        "gap_to": g_to,
        "closure": (g_from - g_to) / g_from if g_from else np.nan,
        "closure_lo": lo,
        "closure_hi": hi,
    }


def by_loeuf_within_expression(
    universe: pd.DataFrame, detection: pd.DataFrame, arm: str
) -> pd.DataFrame:
    """Recovery per LOEUF decile *within* expression tertiles.

    The raw decile view is badly confounded and inverts the matched result, so
    it must not be shown on its own. Across LOEUF deciles 0 -> 9 the median
    cortex TPM falls from 12.71 to 0.18, coding exons from 17.9 to 3.4, and the
    number of brain tissues expressed from 12.1 to 4.5. The least-constrained
    decile is barely expressed, so it has few detectable eQTLs for reasons that
    have nothing to do with selection -- and at the single-nucleus rungs that
    confound is strong enough to make the MOST constrained genes look like the
    BEST recovered.

    Holding expression roughly fixed and then looking across constraint is what
    the decile view has to do to say anything. It is the unmatched analogue of
    what D-002's matching does, and the two should agree.
    """
    col = ARMS[arm]
    det = detection.pivot(index="gene_id", columns="rung", values=col)
    u = universe.copy()
    u["expr_tertile"] = pd.qcut(
        u["tpm_reference"].rank(method="first"), 3, labels=["low", "mid", "high"]
    )

    rows = []
    seen = set()
    for rung in [r for r in RUNG_ORDER + list(RUNG_LABELS) if r in det.columns]:
        if rung in seen:
            continue
        seen.add(rung)
        hit = det[rung].reindex(u.index).fillna(False).astype(bool)
        for (tertile, decile), idx in u.groupby(
            ["expr_tertile", "loeuf_decile"], observed=True
        ).groups.items():
            d = hit.loc[idx]
            k, n = int(d.sum()), len(d)
            if n < 30:  # too few to plot an interval worth reading
                continue
            lo, hi = _wilson(k, n)
            rows.append(
                {
                    "rung": rung,
                    "arm": arm,
                    "expr_tertile": str(tertile),
                    "loeuf_decile": int(decile),
                    "n": n,
                    "k": k,
                    "rate": k / n,
                    "lo": lo,
                    "hi": hi,
                }
            )
    return pd.DataFrame(rows)


def schema_overlay(
    schema: pd.DataFrame, detection: pd.DataFrame, arm: str
) -> pd.DataFrame:
    """The 32 published SCHEMA genes, as a small sharp overlay on the curve."""
    col = ARMS[arm]
    det = detection.pivot(index="gene_id", columns="rung", values=col)
    genes = schema.index[schema["is_schema_published_fdr"]]
    rows = []
    for rung in [r for r in RUNG_ORDER + list(RUNG_LABELS) if r in det.columns]:
        if any(r["rung"] == rung for r in rows):
            continue
        d = det[rung].reindex(genes).fillna(False).astype(bool)
        k, n = int(d.sum()), len(d)
        lo, hi = _wilson(k, n)
        rows.append(
            {
                "rung": rung,
                "arm": arm,
                "n_schema": n,
                "k": k,
                "rate": k / n if n else np.nan,
                "lo": lo,
                "hi": hi,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    prov = Provenance("s06_recovery")
    rng = np.random.default_rng(cfg.RANDOM_SEED)

    universe = pd.read_parquet(GENE_UNIVERSE)
    gene_sets = pd.read_parquet(cfg.GENE_SETS_TABLE)
    gene_sets = gene_sets[gene_sets["set"].isin(["constrained", "control"])]
    schema = pd.read_parquet(SCHEMA_SETS)
    detection = pd.read_parquet(DETECTION_BY_RUNG)

    prov.note(
        "sets entering the curve",
        f"{int((gene_sets['set'] == 'constrained').sum()):,} constrained, "
        f"{int((gene_sets['set'] == 'control').sum()):,} unique controls "
        f"(weights sum to {gene_sets.loc[gene_sets['set'] == 'control', 'weight'].sum():,.0f})",
    )

    curves, loeuf, loeuf_expr, schema_rows = [], [], [], []
    for arm in ARMS:
        curves.append(compute(gene_sets, detection, arm, rng))
        loeuf.append(by_loeuf(universe, detection, arm))
        loeuf_expr.append(by_loeuf_within_expression(universe, detection, arm))
        schema_rows.append(schema_overlay(schema, detection, arm))

    pd.concat(loeuf_expr, ignore_index=True).to_parquet(
        RECOVERY_BY_LOEUF_EXPR, index=False
    )

    curve = pd.concat(curves, ignore_index=True)
    # D-003: the rung-level gap contrasts are their own small family.
    for arm, grp in curve.groupby("arm"):
        curve.loc[grp.index, "gap_q"] = multipletests(
            grp["gap_p"].clip(lower=1e-12), method="fdr_bh"
        )[1]

    curve.to_parquet(cfg.RECOVERY_TABLE, index=False)
    pd.concat(loeuf, ignore_index=True).to_parquet(RECOVERY_BY_LOEUF, index=False)
    pd.concat(schema_rows, ignore_index=True).to_parquet(
        cfg.DIR_PROCESSED / "recovery_schema.parquet", index=False
    )

    # -- report -------------------------------------------------------------
    for arm in ARMS:
        sub = curve[curve["arm"] == arm]
        print(f"\n=== {arm} ===")
        print(
            f"  {'rung':<14}{'N':>6}{'constrained':>22}{'control':>22}{'gap':>22}"
        )
        for _, r in sub.iterrows():
            print(
                f"  {r['rung']:<14}{r['n_donors'] or 0:>6}"
                f"{r['constrained_rate']:>10.3f} [{r['constrained_lo']:.3f},{r['constrained_hi']:.3f}]"
                f"{r['control_rate']:>10.3f} [{r['control_lo']:.3f},{r['control_hi']:.3f}]"
                f"{r['gap']:>10.3f} [{r['gap_lo']:.3f},{r['gap_hi']:.3f}]"
            )

    closures = [
        c
        for arm in ARMS
        if (c := gap_closure(gene_sets, detection, arm, rng))
    ]
    if closures:
        cl = pd.DataFrame(closures)
        cl.to_parquet(cfg.DIR_PROCESSED / "gap_closure.parquet", index=False)
        print("\n=== headline: fraction of the bulk-cortex gap closed at "
              "single-nucleus major-cell-type resolution ===")
        for _, r in cl.iterrows():
            print(
                f"  {r['arm']:<8} gap {r['gap_from']:.3f} -> {r['gap_to']:.3f}   "
                f"closure {r['closure']:6.1%}  "
                f"[{r['closure_lo']:.1%}, {r['closure_hi']:.1%}]"
            )
            prov.record(
                f"gap closure, {r['arm']} arm",
                0,
                0,
                detail=(
                    f"gap {r['gap_from']:.3f} -> {r['gap_to']:.3f}; "
                    f"closure {r['closure']:.1%} "
                    f"[{r['closure_lo']:.1%}, {r['closure_hi']:.1%}]"
                ),
            )

    prov.save(cfg.DIR_LOGS / "s06_recovery.json")
    print(f"\nWrote {cfg.RECOVERY_TABLE}")


if __name__ == "__main__":
    main()
