"""Stage 1, step 7: what actually closes the gap, donors, granularity, or assay?

This module exists because the answer changed once rung 2 was repaired, and it
changed in a way that overturned the previous conclusion.

Three candidate explanations for the ~60% closure between bulk and
single-nucleus, and this project can now test all three directly:

  DONOR COUNT       more samples -> more power -> weaker eQTLs found
  GRANULARITY       finer cell-type labels -> context-specific eQTLs revealed
  ASSAY             measuring nuclei rather than tissue

The design separates them because the five arms cross donor count and assay
almost orthogonally:

    arm                assay                       donors
    GTEx cortex        bulk tissue RNA-seq            205
    PsychENCODE        bulk tissue RNA-seq          1,387
    Bryois pseudobulk  snRNA-seq, nuclei pooled       192
    Bryois cell types  snRNA-seq, 8 types             192
    SingleBrain        snRNA-seq, 7 types             983

Two bulk studies differing 6.8x in donors, three single-nucleus arms differing
5.1x in donors. If donor count drove the gap, it would track N within each
assay class. If assay drove it, the two classes would separate and each class
would be internally flat.

Everything is computed on ONE gene set: the genes Bryois tested, so the
comparison is not confounded by which genes each study happened to cover.
Without that restriction the bulk arms are scored on ~17,000 genes and the
Bryois arms on ~14,700, and the difference in denominators would contaminate
the very contrast being drawn.

Writes: data/processed/assay_contrast.parquet
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline.provenance import Provenance
from pipeline.s05_detection import DETECTION_BY_RUNG, DETECTION_LONG
from pipeline.s06_recovery import N_BOOT, _ci, _weighted_rate

OUT = cfg.DIR_PROCESSED / "assay_contrast.parquet"

BRYOIS_ARMS_RUNGS = ("bryois_pb", "bryois_celltype")

#: (rung, assay class, human label, donor count). Order is the story: bulk
#: first, ascending N, then single-nucleus, ascending N.
ARMS = [
    ("gtex_cortex", "bulk", "GTEx cortex", 205),
    ("bulk_brain", "bulk", "PsychENCODE", 1387),
    ("bryois_pb", "single-nucleus", "Bryois pseudobulk", 192),
    ("bryois_celltype", "single-nucleus", "Bryois, 8 cell types", 192),
    ("sn_major", "single-nucleus", "SingleBrain, 7 cell types", 983),
]


def main() -> None:
    prov = Provenance("s11_assay_contrast")

    gene_sets = pd.read_parquet(cfg.GENE_SETS_TABLE)
    gene_sets = gene_sets[gene_sets["set"].isin(["constrained", "control"])]
    det = pd.read_parquet(DETECTION_BY_RUNG)
    long = pd.read_parquet(DETECTION_LONG)

    tested = set(long.loc[long["rung"].isin(BRYOIS_ARMS_RUNGS), "gene_id"])
    cases = gene_sets.index[gene_sets["set"] == "constrained"].intersection(tested)
    ctrl_mask = gene_sets.index.isin(tested) & (gene_sets["set"] == "control").values
    ctrls = gene_sets.index[ctrl_mask]
    weights = gene_sets.loc[ctrls, "weight"].to_numpy(dtype=float)

    prov.record(
        "restrict every arm to one common gene set",
        int((gene_sets["set"] == "constrained").sum()),
        len(cases),
        detail=(
            "The genes Bryois tested. Without this the bulk arms score on "
            "~17,000 genes and the single-nucleus arms on ~14,700, and the "
            "denominator difference would contaminate the assay contrast."
        ),
    )

    piv = det.pivot(index="gene_id", columns="rung", values="detected")
    rng = np.random.default_rng(cfg.RANDOM_SEED)

    rows = []
    for rung, assay, label, n_donors in ARMS:
        if rung not in piv.columns:
            continue
        hit = piv[rung].reindex(gene_sets.index).fillna(False).astype(bool)
        cd = hit.loc[cases].to_numpy()
        kd = hit.loc[ctrls].to_numpy()
        gap = _weighted_rate(kd, weights) - cd.mean()

        boots = np.empty(N_BOOT)
        for i in range(N_BOOT):
            ci = rng.integers(0, len(cd), len(cd))
            wi = rng.integers(0, len(kd), len(kd))
            boots[i] = _weighted_rate(kd[wi], weights[wi]) - cd[ci].mean()
        lo, hi = _ci(boots)

        rows.append(
            {
                "rung": rung,
                "assay": assay,
                "label": label,
                "n_donors": n_donors,
                "constrained_rate": float(cd.mean()),
                "control_rate": _weighted_rate(kd, weights),
                "gap": gap,
                "gap_lo": lo,
                "gap_hi": hi,
                "n_cases": len(cases),
                "n_controls": len(ctrls),
            }
        )

    out = pd.DataFrame(rows)
    out.to_parquet(OUT, index=False)

    print(f"All arms on the same {len(cases):,} constrained / {len(ctrls):,} "
          f"control genes\n")
    print(f"  {'arm':<26}{'assay':<17}{'donors':>8}{'gap':>8}   95% CI")
    for _, r in out.iterrows():
        print(
            f"  {r['label']:<26}{r['assay']:<17}{r['n_donors']:>8}"
            f"{r['gap']:>8.3f}   [{r['gap_lo']:.3f}, {r['gap_hi']:.3f}]"
        )

    bulk = out[out["assay"] == "bulk"]
    sn = out[out["assay"] == "single-nucleus"]
    if bulk.empty or sn.empty:
        prov.save(cfg.DIR_LOGS / "s11_assay_contrast.json")
        return

    # Within-class spread against between-class separation. If donor count
    # drove the gap, the within-class spread would be large (each class spans
    # a 5-7x range of N) and the classes would overlap.
    bulk_spread = float(bulk["gap"].max() - bulk["gap"].min())
    sn_spread = float(sn["gap"].max() - sn["gap"].min())
    separation = float(bulk["gap"].min() - sn["gap"].max())
    n_ratio_bulk = bulk["n_donors"].max() / bulk["n_donors"].min()
    n_ratio_sn = sn["n_donors"].max() / sn["n_donors"].min()

    print(
        f"\n  bulk arms:            spread {bulk_spread:.3f} across a "
        f"{n_ratio_bulk:.1f}x range of donors"
    )
    print(
        f"  single-nucleus arms:  spread {sn_spread:.3f} across a "
        f"{n_ratio_sn:.1f}x range of donors"
    )
    print(f"  gap between classes:  {separation:.3f}")
    verdict = (
        "ASSAY, not donor count: the classes separate and neither is ordered "
        "by N"
        if separation > max(bulk_spread, sn_spread)
        else "inconclusive: within-class spread is comparable to the "
        "between-class separation"
    )
    print(f"  => {verdict}")

    prov.record(
        "assay vs donor count",
        len(out),
        len(sn),
        detail=(
            f"bulk gap spread {bulk_spread:.3f} over {n_ratio_bulk:.1f}x donors; "
            f"single-nucleus spread {sn_spread:.3f} over {n_ratio_sn:.1f}x "
            f"donors; between-class separation {separation:.3f}. {verdict}"
        ),
    )
    prov.save(cfg.DIR_LOGS / "s11_assay_contrast.json")
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
