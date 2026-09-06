"""Stage 2, step 4: are the newly-explained genes druggable?

The brief asks: for genes that gain a causal eQTL link only at high resolution,
check Open Targets for tractability, and flag any connection back to
CHRM4/muscarinic mechanisms.

The framing matters more than the lookup. A gene that is druggable AND newly
visible is interesting precisely because it was invisible to the assay the
field has been using, it is a candidate that bulk-tissue eQTL screens would
have missed. A gene that is druggable and was always visible is not news.
So genes are split three ways before any tractability is attached:

    bulk_and_sn   causal link visible in bulk tissue too
    sn_only       causal link ONLY in single-nucleus data   <- the interesting set
    bulk_only     visible in bulk but not single-nucleus

METHOD
------
Open Targets GraphQL, batched and cached, following the same route and for the
same reason as the sibling scz-target-prioritization repo: the bulk download is
one to two gigabytes of which this reads a single nested field, against a few
hundred genes that batch into a handful of requests.

Open Targets returns 28 boolean tractability buckets grouped by modality and
ordered strongest-evidence-first within each group. The tier for a modality is
the first bucket that is true; a gene with no true bucket is recorded as null
rather than as a worst-place tier, because "no evidence of tractability" and
"tractable, weakly" are different claims.

Writes: data/processed/druggability.parquet
"""

from __future__ import annotations

import json
import time

import pandas as pd
import requests

from pipeline import config as cfg
from pipeline.provenance import Provenance
from pipeline.s01_gene_universe import GENE_UNIVERSE
from pipeline.s04_schema_gene_sets import SCHEMA_SETS
from pipeline.s13_smr import OUT as SMR_RESULTS

OUT = cfg.DIR_PROCESSED / "druggability.parquet"
CACHE = cfg.DIR_INTERIM / "opentargets_tractability.json"

API = "https://api.platform.opentargets.org/api/v4/graphql"
BATCH_SIZE = 200

QUERY = """
query($ids:[String!]!){
  targets(ensemblIds:$ids){
    id
    approvedSymbol
    tractability { modality value label }
  }
}
"""

#: Small-molecule buckets, strongest evidence first. Verified ordering, same as
#: the sibling repo uses.
SM_ORDER = [
    "Approved Drug",
    "Advanced Clinical",
    "Phase 1 Clinical",
    "Structure with Ligand",
    "High-Quality Ligand",
    "High-Quality Pocket",
    "Med-Quality Pocket",
    "Druggable Family",
]
AB_ORDER = [
    "Approved Drug",
    "Advanced Clinical",
    "Phase 1 Clinical",
    "UniProt loc high conf",
    "GO CC high conf",
    "UniProt loc med conf",
    "UniProt SigP or TMHMM",
    "GO CC med conf",
    "Human Protein Atlas loc",
]

#: The CHRM4/Cobenfy callback. Muscarinic receptors plus the machinery a
#: muscarinic story would run through.
MUSCARINIC = {
    "CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5",
    "CHAT", "ACHE", "BCHE", "SLC5A7", "SLC18A3",
}


def fetch_tractability(ensembl_ids: list[str]) -> dict:
    """Batched Open Targets lookup, cached so a rerun is free and offline."""
    cache = {}
    if CACHE.exists():
        cache = json.loads(CACHE.read_text(encoding="utf-8"))

    todo = [g for g in ensembl_ids if g not in cache]
    for i in range(0, len(todo), BATCH_SIZE):
        batch = todo[i : i + BATCH_SIZE]
        try:
            r = requests.post(
                API, json={"query": QUERY, "variables": {"ids": batch}}, timeout=90
            )
            r.raise_for_status()
            payload = r.json()
        except (requests.RequestException, ValueError) as exc:
            print(f"    batch {i // BATCH_SIZE}: {type(exc).__name__}: {exc}")
            continue
        for t in (payload.get("data") or {}).get("targets") or []:
            cache[t["id"]] = {
                "symbol": t.get("approvedSymbol"),
                "tractability": [
                    x for x in (t.get("tractability") or []) if x.get("value")
                ],
            }
        # Genes the API knows nothing about are cached as empty, so a rerun
        # does not ask again.
        for g in batch:
            cache.setdefault(g, {"symbol": None, "tractability": []})
        time.sleep(0.2)

    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(cache, indent=1), encoding="utf-8")
    return cache


def tier(buckets: list[dict], modality: str, order: list[str]) -> str | None:
    labels = {b["label"] for b in buckets if b.get("modality") == modality}
    for label in order:
        if label in labels:
            return label
    return None


def main() -> None:
    prov = Provenance("s15_druggability")

    smr = pd.read_parquet(SMR_RESULTS)
    sig = smr[smr["smr_sig"]]
    universe = pd.read_parquet(GENE_UNIVERSE)
    schema = pd.read_parquet(SCHEMA_SETS)

    bulk_rungs = {"gtex_cortex", "bulk_brain"}
    sn_rungs = {"sn_major", "sn_subtype"}
    in_bulk = set(sig.loc[sig["rung"].isin(bulk_rungs), "gene_id"])
    in_sn = set(sig.loc[sig["rung"].isin(sn_rungs), "gene_id"])

    groups = {
        "sn_only": sorted(in_sn - in_bulk),
        "bulk_and_sn": sorted(in_sn & in_bulk),
        "bulk_only": sorted(in_bulk - in_sn),
    }
    for name, genes in groups.items():
        prov.note(f"SMR genes, {name}", f"{len(genes):,}")

    all_genes = sorted(in_bulk | in_sn)
    print(f"Querying Open Targets for {len(all_genes):,} SMR-significant genes")
    cache = fetch_tractability(all_genes)

    rows = []
    for name, genes in groups.items():
        for g in genes:
            entry = cache.get(g) or {}
            buckets = entry.get("tractability") or []
            symbol = entry.get("symbol") or universe["symbol"].get(g)
            rows.append(
                {
                    "gene_id": g,
                    "symbol": symbol,
                    "group": name,
                    "sm_tier": tier(buckets, "SM", SM_ORDER),
                    "ab_tier": tier(buckets, "AB", AB_ORDER),
                    "is_schema": bool(
                        schema["is_schema_published_fdr"].get(g, False)
                    ),
                    "is_muscarinic": (symbol or "") in MUSCARINIC,
                }
            )

    out = pd.DataFrame(rows)
    out["druggable_sm"] = out["sm_tier"].notna()
    out["clinically_advanced"] = out["sm_tier"].isin(
        ["Approved Drug", "Advanced Clinical", "Phase 1 Clinical"]
    )
    out.to_parquet(OUT, index=False)

    print(f"\n  {'group':<14}{'genes':>7}{'SM tractable':>14}"
          f"{'clinical':>10}{'SCHEMA':>8}")
    for name in ("sn_only", "bulk_and_sn", "bulk_only"):
        g = out[out["group"] == name]
        if g.empty:
            continue
        print(
            f"  {name:<14}{len(g):>7}{int(g['druggable_sm'].sum()):>14}"
            f"{int(g['clinically_advanced'].sum()):>10}"
            f"{int(g['is_schema'].sum()):>8}"
        )
        prov.record(
            f"{name}: SMR genes -> small-molecule tractable",
            len(g),
            int(g["druggable_sm"].sum()),
            detail=(
                f"{int(g['clinically_advanced'].sum())} at Phase 1 or beyond; "
                f"{int(g['is_schema'].sum())} are published SCHEMA genes"
            ),
        )

    top = out[
        (out["group"] == "sn_only") & out["clinically_advanced"]
    ].sort_values("symbol")
    if not top.empty:
        print(
            f"\n  Clinically advanced targets visible ONLY in single-nucleus "
            f"data ({len(top)}):"
        )
        for _, r in top.iterrows():
            mark = "  <- muscarinic" if r["is_muscarinic"] else ""
            schema_mark = " [SCHEMA]" if r["is_schema"] else ""
            print(f"    {r['symbol']:<12} {r['sm_tier']}{schema_mark}{mark}")

    musc = out[out["is_muscarinic"]]
    if musc.empty:
        prov.note(
            "CHRM4 / muscarinic callback",
            "no muscarinic-pathway gene reaches SMR significance at any rung",
        )
        print("\n  Muscarinic pathway: no gene reaches SMR significance.")
    else:
        print("\n  Muscarinic pathway genes with SMR evidence:")
        for _, r in musc.iterrows():
            print(f"    {r['symbol']:<10} {r['group']:<12} {r['sm_tier']}")
        prov.note(
            "CHRM4 / muscarinic callback",
            ", ".join(f"{r['symbol']} ({r['group']})" for _, r in musc.iterrows()),
        )

    prov.save(cfg.DIR_LOGS / "s15_druggability.json")
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
