"""Stage 1, step 5: does the headline survive the choices it rests on?

The Stage 1 result is one number, 59.9% of the constrained-gene gap closes
between bulk cortex and single-nucleus major cell types, and it sits on top of
four choices that could each have gone another way. Two of those were promised
a robustness check in DECISIONS.md and had not been run. This runs them all.

  1. CONSTRAINT SOURCE (D-005). gnomAD v2.1.1 LOEUF is primary because it is the
     vintage Mostafavi and SCHEMA used. v4.1 recomputes constraint on ~730k
     exomes and reorders genes. If the conclusion moves under v4.1 that is a
     finding, not a bug.
  2. CONSTRAINT METRIC. LOEUF < 0.35 is primary; pLI >= 0.9 is the other
     conventional cut and selects a slightly different 3,025 genes.
  3. MULTIPLE TESTING (D-003). bh_within_rung is primary. by_across_all was
     rejected as the calling rule on coherence grounds and listed as an
     unrun sensitivity arm. It is run here.
  4. SCHEMA RELEASE (D-009). The published 32-gene set is primary; the browser
     release's 50 genes are the sensitivity arm.

Every variant rebuilds the matched control set from scratch under the same
D-002 rules and the same seed, because a different case set needs its own
controls, reusing the primary controls would silently mismatch.

Writes: data/processed/robustness.parquet
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from pipeline import config as cfg
from pipeline.provenance import Provenance
from pipeline.s01_gene_universe import GENE_UNIVERSE, strip_version
from pipeline.s03_control_matching import (
    CALIPER_SD,
    CONTROL_LOEUF_MIN,
    _greedy_caliper_match,
    prepare,
)
from pipeline.s04_schema_gene_sets import SCHEMA_SETS
from pipeline.s05_detection import DETECTION_LONG, FDR_ALPHA
from pipeline.s06_recovery import N_BOOT, _ci, _weighted_rate

OUT = cfg.DIR_PROCESSED / "robustness.parquet"

FROM_RUNG, TO_RUNG = "gtex_cortex", "sn_major"
MATCH_COLS = ["log_tpm", "n_coding_exons"]


def load_v4_constraint() -> pd.DataFrame:
    """gnomAD v4.1 LOEUF, one row per gene.

    v4.1 ships one row per transcript with MANE Select and canonical flags.
    MANE Select is preferred where present: it is the transcript the field
    now treats as the gene's representative, falling back to canonical, then
    to the lowest LOEUF, which is how the gnomAD browser presents a gene.
    """
    path = (
        cfg.DIR_RAW / "gnomad_constraint_v4" / "gnomad.v4.1.constraint_metrics.tsv"
    )
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=["gene", "gene_id", "canonical", "mane_select",
                 "lof.oe_ci.upper", "lof.pLI"],
        low_memory=False,
    )
    df = df.rename(columns={"lof.oe_ci.upper": "loeuf_v4", "lof.pLI": "pli_v4"})
    df = df.dropna(subset=["loeuf_v4", "gene_id"])
    df["gene_id"] = strip_version(df["gene_id"])
    df["rank"] = np.where(df["mane_select"] == True, 0,  # noqa: E712
                          np.where(df["canonical"] == True, 1, 2))  # noqa: E712
    df = df.sort_values(["rank", "loeuf_v4"]).drop_duplicates("gene_id", keep="first")
    return df.set_index("gene_id")[["loeuf_v4", "pli_v4"]]


def build_matched(
    df: pd.DataFrame, is_case: pd.Series, is_pool: pd.Series, seed: int
) -> tuple[pd.Index, pd.Index, np.ndarray]:
    """Rebuild cases, controls and weights under D-002 for an arbitrary case set."""
    rng = np.random.default_rng(seed)
    cases = df[is_case]
    pool = df[is_pool]
    pairs = _greedy_caliper_match(cases, pool, MATCH_COLS, CALIPER_SD, rng)
    use = pairs["control_gene"].value_counts()
    return pd.Index(pairs["case_gene"].unique()), pd.Index(use.index), use.to_numpy(float)


def gap_and_closure(
    hits: dict[str, pd.Series],
    cases: pd.Index,
    ctrls: pd.Index,
    weights: np.ndarray,
    rng: np.random.Generator,
) -> dict:
    """Gap at both rungs plus the closure between them, on shared replicates."""
    arrays = {
        r: (h.reindex(cases).fillna(False).to_numpy(bool),
            h.reindex(ctrls).fillna(False).to_numpy(bool))
        for r, h in hits.items()
    }

    def gap(r, ci, wi):
        c, k = arrays[r]
        return _weighted_rate(k[wi], weights[wi]) - c[ci].mean()

    full_c, full_k = np.arange(len(cases)), np.arange(len(ctrls))
    g_from = gap(FROM_RUNG, full_c, full_k)
    g_to = gap(TO_RUNG, full_c, full_k)

    closures = np.empty(N_BOOT)
    gaps_to = np.empty(N_BOOT)
    for b in range(N_BOOT):
        ci = rng.integers(0, len(cases), len(cases))
        wi = rng.integers(0, len(ctrls), len(ctrls))
        gf, gt = gap(FROM_RUNG, ci, wi), gap(TO_RUNG, ci, wi)
        closures[b] = (gf - gt) / gf if gf else np.nan
        gaps_to[b] = gt
    lo, hi = _ci(closures[~np.isnan(closures)])
    glo, ghi = _ci(gaps_to)
    return {
        "n_cases": len(cases),
        "n_controls": len(ctrls),
        "gap_from": g_from,
        "gap_to": g_to,
        "gap_to_lo": glo,
        "gap_to_hi": ghi,
        "closure": (g_from - g_to) / g_from if g_from else np.nan,
        "closure_lo": lo,
        "closure_hi": hi,
    }


def detection_hits(long: pd.DataFrame, correction: str) -> dict[str, pd.Series]:
    """Per-rung 'detected in any cell type', under a given D-003 correction."""
    d = long.copy()
    ok = d["p_bonf"].notna()
    if correction == "bh_within_rung":
        q = pd.Series(np.nan, index=d.index)
        for _, idx in d[ok].groupby("rung").groups.items():
            q[idx] = multipletests(d.loc[idx, "p_bonf"], method="fdr_bh")[1]
    elif correction == "by_across_all":
        q = pd.Series(np.nan, index=d.index)
        q[ok] = multipletests(d.loc[ok, "p_bonf"], method="fdr_by")[1]
    else:
        raise ValueError(correction)
    d["hit"] = q <= FDR_ALPHA
    return {
        rung: g.groupby("gene_id")["hit"].any()
        for rung, g in d.groupby("rung")
        if rung in (FROM_RUNG, TO_RUNG)
    }


def main() -> None:
    prov = Provenance("s09_robustness")
    universe = pd.read_parquet(GENE_UNIVERSE)
    long = pd.read_parquet(DETECTION_LONG)
    schema = pd.read_parquet(SCHEMA_SETS)

    v4 = load_v4_constraint()
    universe = universe.join(v4, how="left")
    prov.record(
        "universe x gnomAD v4.1 constraint",
        len(universe),
        int(universe["loeuf_v4"].notna().sum()),
        detail="genes without a v4.1 LOEUF cannot enter the v4.1 variant",
    )
    df = prepare(universe)

    hits_primary = detection_hits(long, "bh_within_rung")
    hits_by = detection_hits(long, "by_across_all")

    variants = [
        ("constraint source", "gnomAD v2.1.1 LOEUF < 0.35  (primary)",
         df["loeuf"] < 0.35, df["loeuf"] >= CONTROL_LOEUF_MIN, hits_primary),
        ("constraint source", "gnomAD v4.1 LOEUF < 0.35",
         df["loeuf_v4"] < 0.35, df["loeuf_v4"] >= CONTROL_LOEUF_MIN, hits_primary),
        ("constraint metric", "gnomAD v2.1.1 pLI >= 0.9",
         df["pli"] >= 0.9, df["pli"] < 0.1, hits_primary),
        ("multiple testing", "BH within rung  (primary)",
         df["loeuf"] < 0.35, df["loeuf"] >= CONTROL_LOEUF_MIN, hits_primary),
        ("multiple testing", "BY across the whole grid",
         df["loeuf"] < 0.35, df["loeuf"] >= CONTROL_LOEUF_MIN, hits_by),
    ]

    rows = []
    for family, label, is_case, is_pool, hits in variants:
        is_case = is_case.fillna(False)
        is_pool = is_pool.fillna(False)
        cases, ctrls, w = build_matched(df, is_case, is_pool, cfg.RANDOM_SEED)
        rng = np.random.default_rng(cfg.RANDOM_SEED)
        res = gap_and_closure(hits, cases, ctrls, w, rng)
        rows.append({"family": family, "variant": label, **res})

    # SCHEMA sets are not matched: they are a named gene list, so the
    # comparison is the raw recovery of that list, not a gap.
    for label, col in (
        ("published, FDR < 0.05  (primary)", "is_schema_published_fdr"),
        ("browser 2026-08-21, local BH FDR < 0.05", "is_schema_browser_fdr"),
    ):
        genes = schema.index[schema[col].fillna(False)]
        rec = {}
        for rung, h in hits_primary.items():
            d = h.reindex(genes).fillna(False)
            rec[rung] = float(d.mean())
        rows.append({
            "family": "SCHEMA release",
            "variant": label,
            "n_cases": len(genes),
            "n_controls": np.nan,
            "gap_from": np.nan,
            "gap_to": np.nan,
            "gap_to_lo": np.nan,
            "gap_to_hi": np.nan,
            "closure": np.nan,
            "closure_lo": np.nan,
            "closure_hi": np.nan,
            "schema_rate_from": rec.get(FROM_RUNG),
            "schema_rate_to": rec.get(TO_RUNG),
        })

    out = pd.DataFrame(rows)
    out.to_parquet(OUT, index=False)

    print(f"Headline under each variant  ({FROM_RUNG} -> {TO_RUNG})\n")
    for family, g in out.groupby("family", sort=False):
        print(f"  {family}")
        for _, r in g.iterrows():
            if pd.notna(r["closure"]):
                print(
                    f"    {r['variant']:<42} n={int(r['n_cases']):>5}  "
                    f"gap {r['gap_from']:.3f} -> {r['gap_to']:.3f}  "
                    f"closure {r['closure']:5.1%} "
                    f"[{r['closure_lo']:.1%}, {r['closure_hi']:.1%}]"
                )
            else:
                print(
                    f"    {r['variant']:<42} n={int(r['n_cases']):>5}  "
                    f"recovery {r['schema_rate_from']:.3f} -> "
                    f"{r['schema_rate_to']:.3f}"
                )
        print()

    closures = out["closure"].dropna()
    prov.record(
        "robustness of the headline closure",
        len(closures),
        int(((closures > 0.3) & (closures < 0.9)).sum()),
        detail=(
            f"closure ranges {closures.min():.1%} to {closures.max():.1%} across "
            f"{len(closures)} variants of the constraint source, constraint "
            f"metric and multiple-testing rule"
        ),
    )
    prov.save(cfg.DIR_LOGS / "s09_robustness.json")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
