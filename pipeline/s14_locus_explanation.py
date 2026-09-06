"""Stage 2, step 3: how many GWAS loci gain an eQTL explanation, per rung?

This is the Stage 2 headline the project brief asks for. The recovery curve in
Stage 1 counted genes; this counts LOCI, which is the unit a GWAS actually
delivers and the unit a drug programme starts from.

DEFINING A LOCUS
----------------
Distance-based clumping on the PGC3 genome-wide significant variants: take all
variants at p < 5e-8, sort by position, and merge any within MERGE_KB of each
other into one locus. No LD reference is used.

That is a deliberate simplification and it is stated rather than buried.
LD-based clumping is more standard, but it needs a reference panel matched to
the GWAS ancestry, and the result here is a COUNT COMPARED ACROSS RUNGS, the
same locus definition applies to every rung, so a locus that is slightly too
wide or too narrow affects all rungs identically. The absolute number of loci
should not be quoted against a published locus count; the ratio between rungs
is what this measures.

The MHC is excluded. It is a single enormous LD block whose association signal
spans megabases, so distance clumping would merge it into one locus containing
hundreds of genes and it would dominate any per-locus count. Excluding it is
conventional and is reported.

ATTRIBUTING A GENE TO A LOCUS
-----------------------------
An SMR result is a statement about a specific variant: the gene's top cis
instrument. So a gene is attributed to the locus containing THAT variant, not
to the locus nearest its transcription start site. This is the honest mapping:
it is the variant the causal claim runs through.

Writes: data/processed/locus_explanation.parquet
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline.provenance import Provenance
from pipeline.s13_smr import OUT as SMR_RESULTS
from pipeline.s13_smr import load_pgc3

OUT = cfg.DIR_PROCESSED / "locus_explanation.parquet"

GWAS_SIG = 5e-8
#: Genome-wide significant variants within this distance are one locus.
MERGE_KB = 250

#: MHC on GRCh37. One LD block spanning megabases; distance clumping would
#: collapse it into a single locus holding hundreds of genes and it would
#: dominate every count.
MHC = ("6", 25_000_000, 34_000_000)

RUNG_ORDER = ["gtex_cortex", "bulk_brain", "sn_major", "sn_subtype"]


def define_loci(gwas: pd.DataFrame, prov: Provenance) -> pd.DataFrame:
    sig = gwas[gwas["gwas_p"] < GWAS_SIG].copy()
    prov.record(
        f"PGC3 variants at p < {GWAS_SIG:g}",
        len(gwas),
        len(sig),
        detail="the raw material for locus definition",
    )

    sig["chrom"] = sig["chrom"].astype(str)
    in_mhc = (
        (sig["chrom"] == MHC[0])
        & (sig["pos"] >= MHC[1])
        & (sig["pos"] <= MHC[2])
    )
    prov.record(
        "genome-wide significant variants outside the MHC",
        len(sig),
        int((~in_mhc).sum()),
        detail=(
            f"MHC chr{MHC[0]}:{MHC[1]:,}-{MHC[2]:,} excluded: one LD block "
            "spanning megabases, which distance clumping would collapse into a "
            "single locus containing hundreds of genes."
        ),
    )
    sig = sig[~in_mhc].sort_values(["chrom", "pos"])

    loci, lid = [], 0
    for chrom, g in sig.groupby("chrom", sort=False):
        start = prev = None
        best_p, best_rs = None, None
        for pos, p, rs in zip(g["pos"], g["gwas_p"], g["rsid"]):
            if start is None:
                start = prev = pos
                best_p, best_rs = p, rs
                continue
            if pos - prev <= MERGE_KB * 1000:
                prev = pos
                if p < best_p:
                    best_p, best_rs = p, rs
            else:
                loci.append((lid, chrom, start, prev, best_rs, best_p))
                lid += 1
                start = prev = pos
                best_p, best_rs = p, rs
        if start is not None:
            loci.append((lid, chrom, start, prev, best_rs, best_p))
            lid += 1

    out = pd.DataFrame(
        loci, columns=["locus", "chrom", "start", "end", "lead_rsid", "lead_p"]
    )
    prov.record(
        "loci after distance clumping",
        len(sig),
        len(out),
        detail=f"variants merged within {MERGE_KB} kb; no LD reference used",
    )
    return out


def assign(smr: pd.DataFrame, loci: pd.DataFrame) -> pd.DataFrame:
    """Attribute each SMR result to the locus containing its instrument."""
    smr = smr.copy()
    smr["chrom"] = smr["chrom"].astype(str)
    out = []
    for chrom, g in smr.groupby("chrom", sort=False):
        L = loci[loci["chrom"] == chrom]
        if L.empty:
            g = g.assign(locus=pd.NA)
            out.append(g)
            continue
        starts = L["start"].to_numpy()
        ends = L["end"].to_numpy()
        ids = L["locus"].to_numpy()
        pos = g["pos"].to_numpy()
        idx = np.searchsorted(starts, pos, side="right") - 1
        ok = (idx >= 0) & (pos <= np.where(idx >= 0, ends[idx], -1))
        g = g.assign(locus=np.where(ok, ids[idx], pd.NA))
        out.append(g)
    return pd.concat(out, ignore_index=True)


def main() -> None:
    prov = Provenance("s14_locus_explanation")

    gwas = load_pgc3()
    loci = define_loci(gwas, prov)
    smr = pd.read_parquet(SMR_RESULTS)

    smr = assign(smr, loci)
    sig = smr[smr["smr_sig"] & smr["locus"].notna()]

    rows = []
    seen_rungs = [r for r in RUNG_ORDER if r in set(smr["rung"])]
    for rung in seen_rungs:
        g = sig[sig["rung"] == rung]
        rows.append(
            {
                "rung": rung,
                "loci_total": len(loci),
                "loci_explained": g["locus"].nunique(),
                "genes": g["gene_id"].nunique(),
                "gene_locus_pairs": len(g.drop_duplicates(["gene_id", "locus"])),
            }
        )
    res = pd.DataFrame(rows)
    res["frac_explained"] = res["loci_explained"] / res["loci_total"]

    # Which loci are explained ONLY at single-nucleus resolution, the
    # question the brief actually asks.
    bulk = set(sig.loc[sig["rung"].isin(["gtex_cortex", "bulk_brain"]), "locus"])
    sn = set(sig.loc[sig["rung"].isin(["sn_major", "sn_subtype"]), "locus"])
    gained = sn - bulk
    lost = bulk - sn

    res.to_parquet(OUT, index=False)

    print(f"GWAS loci (p < {GWAS_SIG:g}, MHC excluded, {MERGE_KB} kb clumping): "
          f"{len(loci)}\n")
    print(f"  {'rung':<14}{'loci explained':>16}{'of total':>10}{'genes':>8}")
    for _, r in res.iterrows():
        print(
            f"  {r['rung']:<14}{r['loci_explained']:>16}{r['frac_explained']:>10.1%}"
            f"{r['genes']:>8}"
        )

    print(
        f"\n  loci explained in bulk:            {len(bulk)}"
        f"\n  loci explained at single-nucleus:  {len(sn)}"
        f"\n  GAINED only at single-nucleus:     {len(gained)}"
        f"\n  explained in bulk but not sn:      {len(lost)}"
    )
    prov.record(
        "GWAS loci gaining an eQTL explanation at single-nucleus resolution",
        len(bulk),
        len(sn),
        detail=(
            f"{len(gained)} loci explained only in single-nucleus data, "
            f"{len(lost)} only in bulk. This is the Stage 2 headline."
        ),
    )

    prov.save(cfg.DIR_LOGS / "s14_locus_explanation.json")
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
