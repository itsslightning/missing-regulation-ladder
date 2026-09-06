"""Stage 1, step 6: which SCHEMA genes switch on, and where.

Section 2 of the Stage 1 report says 15 SCHEMA genes acquire a detectable
cis-eQTL between bulk cortex and single-nucleus resolution, and calls that the
Stage 2 target list. This names them, and says which cell types each one
switches on in.

That second part is the useful bit. A SCHEMA gene that becomes visible only in
microglia is a different therapeutic proposition from one visible in every
neuronal subtype, and Stage 2's colocalization work needs to know which
cell-type summary statistics to pull, the SingleBrain full-association files
are 4-10 GB each, so pulling them per gene of interest rather than in bulk is
the whole reason D-006 deferred them.

Output feeds the Stage 3 dashboard's per-cell-type heatmap directly.

Writes: data/processed/schema_switch.parquet
        docs/schema_switch.md
"""

from __future__ import annotations

import pandas as pd

from pipeline import config as cfg
from pipeline.provenance import Provenance
from pipeline.s02_audit_rungs import SB_MAJOR
from pipeline.s04_schema_gene_sets import SCHEMA_SETS
from pipeline.s05_detection import DETECTION_BY_RUNG, DETECTION_LONG

OUT = cfg.DIR_PROCESSED / "schema_switch.parquet"
DOC = cfg.DIR_DOCS / "schema_switch.md"

BULK_RUNGS = ["gtex_cortex", "bulk_brain"]
SN_RUNG = "sn_major"


def main() -> None:
    prov = Provenance("s10_schema_switch")

    schema = pd.read_parquet(SCHEMA_SETS)
    by_rung = pd.read_parquet(DETECTION_BY_RUNG)
    long = pd.read_parquet(DETECTION_LONG)

    genes = schema.index[schema["is_schema_published_fdr"].fillna(False)]
    det = by_rung.pivot(index="gene_id", columns="rung", values="detected")

    tbl = pd.DataFrame(index=genes)
    tbl["symbol"] = schema.loc[genes, "symbol"]
    tbl["q_schema"] = schema.loc[genes, "schema_q_browser"]

    for rung in BULK_RUNGS + [SN_RUNG, "sn_subtype"]:
        tbl[rung] = (
            det[rung].reindex(genes).fillna(False).astype(bool)
            if rung in det.columns
            else False
        )

    # `bulk` is the union of the bulk rungs, both of which are now scored
    # from full association tables rather than significant-only files.
    tbl["bulk_any"] = tbl[[r for r in BULK_RUNGS if r in tbl]].any(axis=1)
    tbl["switches_on"] = (~tbl["bulk_any"]) & tbl[SN_RUNG]
    tbl["already_visible"] = tbl["bulk_any"] & tbl[SN_RUNG]
    tbl["never_seen"] = (~tbl["bulk_any"]) & (~tbl[SN_RUNG])

    # Which major cell types each gene is detected in.
    sn = long[(long["rung"] == SN_RUNG) & long["detected"]]
    per_cell = (
        sn[sn["gene_id"].isin(genes)]
        .groupby("gene_id")["cell"]
        .apply(lambda s: ",".join(sorted(set(s), key=lambda c: SB_MAJOR.index(c))))
    )
    tbl["cell_types"] = per_cell.reindex(genes).fillna("")
    tbl["n_cell_types"] = tbl["cell_types"].str.count(",").where(
        tbl["cell_types"] != "", -1
    ) + 1

    tbl = tbl.sort_values(
        ["switches_on", "n_cell_types"], ascending=[False, False]
    )
    tbl.to_parquet(OUT)

    n_switch = int(tbl["switches_on"].sum())
    n_already = int(tbl["already_visible"].sum())
    n_never = int(tbl["never_seen"].sum())

    prov.record(
        "SCHEMA genes gaining an eQTL at single-nucleus resolution",
        len(tbl),
        n_switch,
        detail=(
            f"{n_switch} switch on, {n_already} already visible in bulk, "
            f"{n_never} still undetected at single-nucleus resolution"
        ),
    )

    lines = [
        "# SCHEMA genes across the ladder\n",
        f"The {len(tbl)} published SCHEMA genes (Singh et al. 2022, FDR < 0.05), "
        "scored under the uniform detection rule (D-001).\n",
        f"- **{n_switch} switch on**: no detectable cis-eQTL in bulk, one at "
        "single-nucleus major-cell-type resolution. This is the Stage 2 "
        "colocalization target list.",
        f"- **{n_already} were already visible** in bulk.",
        f"- **{n_never} remain undetected** even at single-nucleus resolution.\n",
        "| Gene | Bulk cortex | Bulk brain | sn major | sn subtype | Cell types (sn major) |",
        "|---|---|---|---|---|---|",
    ]
    tick = {True: "yes", False: "—"}
    for _, r in tbl.iterrows():
        mark = " **←switches on**" if r["switches_on"] else ""
        lines.append(
            f"| {r['symbol']}{mark} | {tick[bool(r['gtex_cortex'])]} | "
            f"{tick[bool(r['bulk_brain'])]} | {tick[bool(r[SN_RUNG])]} | "
            f"{tick[bool(r['sn_subtype'])]} | {r['cell_types'] or '—'} |"
        )
    DOC.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    # ASCII for the console: the markdown file is UTF-8, but a Windows console
    # is cp1252 and turns an em dash into a replacement character.
    ascii_tick = {True: "yes", False: "no"}
    print(f"SCHEMA genes: {n_switch} switch on, {n_already} already visible in "
          f"bulk, {n_never} never seen\n")
    print(f"  {'gene':<10}{'bulk':>6}{'sn':>5}  cell types")
    for _, r in tbl.iterrows():
        if r["switches_on"]:
            print(f"  {r['symbol']:<10}{ascii_tick[bool(r['bulk_any'])]:>6}"
                  f"{ascii_tick[bool(r[SN_RUNG])]:>5}  {r['cell_types']}")

    prov.save(cfg.DIR_LOGS / "s10_schema_switch.json")
    print(f"\nWrote {OUT}\n      {DOC}")


if __name__ == "__main__":
    main()
