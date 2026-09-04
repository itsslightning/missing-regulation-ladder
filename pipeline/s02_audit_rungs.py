"""Stage 0, step 2: audit what each rung actually ships.

The report's data plan assumes four comparable rungs. This script checks that
assumption against the files rather than against the plan, and writes
docs/stage0_audit.md.

What it is looking for, in order of how much damage each does if missed:

1. Whether a rung supplies a *denominator*. Recovery is a fraction. A file of
   significant associations only -- which is what PsychENCODE's Bonferroni
   release is -- gives a numerator and nothing else. Scoring it against the
   fixed gene universe is possible; scoring it against "genes tested" is not,
   because that set is not in the file.

2. Whether the significance statistic means the same thing at each rung. GTEx's
   qval is a Storey q on a beta-approximated permutation p-value; SingleBrain's
   qval is a Storey q on a within-gene variant-corrected p-value. Those are
   comparable. Bryois ships nominal p-values with no per-gene correction at
   all, which is not.

3. Whether the outcome measure saturates. If nearly every tested gene is an
   eGene at a rung, the recovery fraction has no headroom left to show a gap,
   and the curve will close for reasons that have nothing to do with selection.
   This is reported per cell type so it cannot be averaged away.

Nothing here decides anything. It measures, and prints what the measurements
imply for the decisions still open in DECISIONS.md.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd

from pipeline.config import DIR_DOCS, DIR_RAW, DIR_RESTRICTED
from pipeline.s01_gene_universe import GENE_UNIVERSE, strip_version

AUDIT_DOC = DIR_DOCS / "stage0_audit.md"

#: SingleBrain's 7 major cell types. Everything else in the record is either a
#: numbered subtype or MiGA3, which is a different meta-analysis again.
SB_MAJOR = ("Ast", "End", "Ext", "IN", "MG", "OD", "OPC")
SB_EXCLUDED = ("MiGA3",)


def _universe() -> set[str]:
    if not GENE_UNIVERSE.exists():
        raise FileNotFoundError(
            f"{GENE_UNIVERSE} missing -- run pipeline.s01_gene_universe first."
        )
    return set(pd.read_parquet(GENE_UNIVERSE).index)


# ---------------------------------------------------------------------------
# Rungs 3 and 4 -- SingleBrain
# ---------------------------------------------------------------------------


def audit_singlebrain(universe: set[str]) -> pd.DataFrame:
    rows = []
    for path in sorted((DIR_RAW / "singlebrain").glob("*_top_assoc.tsv.gz")):
        cell = path.name.split("_")[0]
        with gzip.open(path, "rt") as fh:
            df = pd.read_csv(
                fh,
                sep="\t",
                usecols=["feature", "qval", "Fixed_bonf", "Fixed_P"],
                low_memory=False,
            )
        df["gene_id"] = strip_version(df["feature"])
        in_universe = df["gene_id"].isin(universe)

        if cell in SB_EXCLUDED:
            tier = "excluded"
        elif cell in SB_MAJOR:
            tier = "major (rung 3)"
        else:
            tier = "subtype (rung 4)"

        rows.append(
            {
                "rung": tier,
                "cell": cell,
                "genes_tested": len(df),
                "in_universe": int(in_universe.sum()),
                # Two candidate definitions of "detectable", both computed so
                # the size of the D-001 choice is visible rather than argued.
                "egenes_qval": int((df["qval"] <= 0.05).sum()),
                "egenes_bonf": int((df["Fixed_bonf"] <= 0.05).sum()),
            }
        )
    out = pd.DataFrame(rows)
    out["frac_qval_of_tested"] = out["egenes_qval"] / out["genes_tested"]
    out["frac_bonf_of_tested"] = out["egenes_bonf"] / out["genes_tested"]
    out["frac_qval_of_universe"] = out["egenes_qval"] / len(universe)
    return out


# ---------------------------------------------------------------------------
# Rung 1 -- GTEx
# ---------------------------------------------------------------------------


def audit_gtex(universe: set[str]) -> pd.DataFrame:
    rows = []
    for path in sorted((DIR_RAW / "gtex_v10").glob("*.eGenes.txt.gz")):
        tissue = path.name.split(".v10.")[0]
        with gzip.open(path, "rt") as fh:
            df = pd.read_csv(
                fh, sep="\t", usecols=["gene_id", "qval", "pval_beta"], low_memory=False
            )
        df["gene_id_bare"] = strip_version(df["gene_id"])
        rows.append(
            {
                "rung": "bulk tissue (rung 1)",
                "cell": tissue,
                "genes_tested": len(df),
                "in_universe": int(df["gene_id_bare"].isin(universe).sum()),
                "egenes_qval": int((df["qval"] <= 0.05).sum()),
                "egenes_bonf": pd.NA,
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["frac_qval_of_tested"] = out["egenes_qval"] / out["genes_tested"]
    out["frac_bonf_of_tested"] = pd.NA
    out["frac_qval_of_universe"] = out["egenes_qval"] / len(universe)
    return out


# ---------------------------------------------------------------------------
# Rung 2 -- PsychENCODE
# ---------------------------------------------------------------------------


def audit_psychencode(universe: set[str]) -> dict:
    path = DIR_RESTRICTED / "psychencode" / "DER-08b_hg38_eQTL.bonferroni.txt"
    if not path.exists():
        return {"available": False}
    df = pd.read_csv(path, sep="\t", usecols=["gene_id", "FDR", "top_SNP"])
    genes = set(strip_version(df["gene_id"]).unique())
    return {
        "available": True,
        "rows": len(df),
        "unique_genes": len(genes),
        "in_universe": len(genes & universe),
        "has_denominator": False,
        "note": (
            "Significant pairs only. Every gene in this file is an eGene, so the "
            "file cannot supply the set of genes tested-but-not-significant."
        ),
    }


# ---------------------------------------------------------------------------
# Power contrast -- Bryois
# ---------------------------------------------------------------------------


def audit_bryois(universe: set[str]) -> list[dict]:
    """Bryois files are headerless, space-separated, one per chromosome.

    Only the chr1 files are read here: the point of the audit is the format and
    the gene-ID shape, and reading 200 files to establish that would be waste.
    Gene IDs are a compound `SYMBOL_ENSG`, unversioned.
    """
    out = []
    for name, label in (("Astrocytes.1.gz", "cell type"), ("pb.1.gz", "pseudobulk")):
        path = DIR_RAW / "bryois" / name
        if not path.exists():
            continue
        try:
            df = pd.read_csv(
                path,
                sep=r"\s+",
                header=None,
                names=["gene", "snp", "tss_distance", "p_nominal", "beta"],
                engine="python",
            )
        except EOFError:
            # A truncated download, not a format problem. Worth surfacing
            # loudly rather than skipping quietly -- a short file would
            # otherwise show up as a genuinely smaller gene count.
            print(f"  WARNING: {path.name} is truncated; re-download it.")
            continue
        # `AL627309.1_ENSG00000238009` -> ENSG00000238009
        ensg = df["gene"].str.extract(r"(ENSG\d+)", expand=False)
        genes = set(ensg.dropna().unique())
        out.append(
            {
                "file": name,
                "arm": label,
                "rows_chr1": len(df),
                "genes_chr1": len(genes),
                "in_universe_chr1": len(genes & universe),
                "id_parsed_ok": float(ensg.notna().mean()),
                "min_p": float(df["p_nominal"].min()),
                "has_per_gene_correction": False,
            }
        )
    return out


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def _md_table(df: pd.DataFrame, float_cols: tuple[str, ...] = ()) -> str:
    d = df.copy()
    for c in float_cols:
        if c in d:
            d[c] = d[c].map(lambda v: "-" if pd.isna(v) else f"{v:.3f}")
    for c in d.columns:
        if c not in float_cols:
            d[c] = d[c].map(lambda v: "-" if pd.isna(v) else v)
    header = "| " + " | ".join(d.columns) + " |"
    sep = "|" + "|".join("---" for _ in d.columns) + "|"
    body = "\n".join("| " + " | ".join(str(v) for v in r) + " |" for r in d.values)
    return "\n".join([header, sep, body])


def main() -> None:
    universe = _universe()
    print(f"Gene universe: {len(universe):,} genes\n")

    sb = audit_singlebrain(universe)
    gtex = audit_gtex(universe)
    pec = audit_psychencode(universe)
    bry = audit_bryois(universe)

    fl = ("frac_qval_of_tested", "frac_bonf_of_tested", "frac_qval_of_universe")
    parts = ["# Stage 0 audit: what each rung actually ships\n"]
    parts.append(f"Gene universe: **{len(universe):,}** genes "
                 "(protein-coding, gnomAD v2.1.1 LOEUF, GTEx brain expression).\n")

    parts.append("\n## Rung 1 — GTEx v10 brain\n")
    if gtex.empty:
        parts.append("_Not yet extracted._\n")
    else:
        parts.append(_md_table(gtex, fl) + "\n")

    parts.append("\n## Rung 2 — PsychENCODE\n")
    if not pec["available"]:
        parts.append("_Not downloaded._\n")
    else:
        parts.append(
            f"- rows: {pec['rows']:,}\n"
            f"- unique genes: {pec['unique_genes']:,} "
            f"({pec['in_universe']:,} in universe)\n"
            f"- **supplies a denominator: no.** {pec['note']}\n"
        )

    parts.append("\n## Rungs 3–4 — SingleBrain\n")
    parts.append(_md_table(sb, fl) + "\n")

    parts.append("\n## Power contrast — Bryois\n")
    if bry:
        parts.append(_md_table(pd.DataFrame(bry), ("id_parsed_ok", "min_p")) + "\n")
    else:
        parts.append("_Not downloaded._\n")

    AUDIT_DOC.parent.mkdir(parents=True, exist_ok=True)
    AUDIT_DOC.write_text("\n".join(parts), encoding="utf-8", newline="\n")

    # Console summary: the two numbers that decide whether the design survives.
    print(sb.groupby("rung")[["genes_tested", "egenes_qval", "egenes_bonf"]].sum())
    print("\nSaturation check (fraction of TESTED genes that are eGenes):")
    for rung, grp in sb[sb["rung"] != "excluded"].groupby("rung"):
        print(
            f"  {rung:20s} qval   min={grp['frac_qval_of_tested'].min():.3f} "
            f"max={grp['frac_qval_of_tested'].max():.3f}"
        )
        print(
            f"  {'':20s} bonf   min={grp['frac_bonf_of_tested'].min():.3f} "
            f"max={grp['frac_bonf_of_tested'].max():.3f}"
        )
    print(f"\nWrote {AUDIT_DOC}")


if __name__ == "__main__":
    main()
