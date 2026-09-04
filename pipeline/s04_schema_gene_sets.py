"""Stage 0, step 4: define the SCHEMA gene sets.

D-009 was settled as "both, published primary": the peer-reviewed Singh et al.
2022 set drives the headline result, and the 2026-08-21 browser release --
roughly 3.6x the cases -- runs alongside as a sensitivity analysis. This module
builds both and keeps them clearly labelled, because the two are not
interchangeable and a figure that silently mixed them would be indefensible.

The two releases differ in more than sample size:

  Singh 2022 (published)   18,324 genes, keyed by GENE SYMBOL, ships `Q meta`
                           (the FDR the paper's thresholds are defined on)
  Browser 2026-08-21       18,573 genes, keyed by ENSG, ships p-values ONLY

So the published set can use the paper's own significance calls, while any set
drawn from the browser release requires computing an FDR here -- a judgment the
published table does not force. That asymmetry is the main reason the published
set is primary.

Joining the published table costs a symbol -> ENSG hop, because it carries no
Ensembl IDs. gnomAD supplies the bridge (it has both), and the loss is reported
rather than assumed: gene symbols drift between GENCODE releases far more than
Ensembl IDs do.

Writes: data/processed/schema_gene_sets.parquet
"""

from __future__ import annotations

import gzip

import pandas as pd
from statsmodels.stats.multitest import multipletests

from pipeline import downloads, sources
from pipeline.config import DIR_PROCESSED
from pipeline.s01_gene_universe import GENE_UNIVERSE, strip_version

SCHEMA_SETS = DIR_PROCESSED / "schema_gene_sets.parquet"

#: The paper's own two thresholds. Singh et al. report 10 genes at exome-wide
#: significance and 32 at FDR < 0.05; both are carried so downstream code can
#: choose the stringency without recomputing anything.
FDR_THRESHOLD = 0.05
EXOME_WIDE_P = 2.14e-6

#: Sheet name in the Singh et al. supplementary workbook.
S5_SHEET = "Table S5 - Gene Results"


def load_published() -> pd.DataFrame:
    """Singh et al. 2022 Supplementary Table 5, joined to ENSG through gnomAD."""
    path = (
        sources.SCHEMA_PUBLISHED.directory / "singh2022_supplementary_tables.xlsx"
    )
    if not path.exists():
        path = downloads.fetch(
            sources.SCHEMA_PUBLISHED, filename="singh2022_supplementary_tables.xlsx"
        )

    df = pd.read_excel(path, sheet_name=S5_SHEET)
    df = df.rename(
        columns={
            "Gene Symbol": "symbol",
            "P meta": "schema_p_meta",
            "Q meta": "schema_q_meta",
        }
    )
    df = df[["symbol", "schema_p_meta", "schema_q_meta"]].dropna(subset=["symbol"])
    df["symbol"] = df["symbol"].astype("string").str.strip()
    df = df.drop_duplicates("symbol", keep="first")
    print(f"  published (Singh 2022): {len(df):,} genes with symbol")
    return df


def load_browser() -> pd.DataFrame:
    """Browser release 2026-08-21, keyed by ENSG, FDR computed here.

    The release ships no q-value, so an FDR has to be produced locally. It is
    Benjamini-Hochberg on the combined case-control-plus-de-novo p-value, which
    is the column corresponding to the published `P meta`. This is a judgment
    the published table does not require, and is the reason this set is the
    sensitivity arm rather than the primary one.
    """
    src = sources.SCHEMA_GENES
    path = src.directory / "SCHEMA_gene_results.tsv.bgz"
    if not path.exists():
        path = downloads.fetch(src, filename="SCHEMA_gene_results.tsv.bgz")

    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", low_memory=False)
    df = df[df["group"] == "meta"]
    df = df.rename(columns={"case_control_plus_de_novo_p_value": "browser_p_meta"})
    df["gene_id"] = strip_version(df["gene_id"])
    df = df[["gene_id", "gene_symbol", "browser_p_meta"]].dropna(
        subset=["browser_p_meta"]
    )
    df = df.drop_duplicates("gene_id", keep="first")

    df["browser_q_meta"] = multipletests(df["browser_p_meta"], method="fdr_bh")[1]
    print(f"  browser (2026-08-21): {len(df):,} genes, BH-FDR computed locally")
    return df.set_index("gene_id")


def build() -> pd.DataFrame:
    universe = pd.read_parquet(GENE_UNIVERSE)
    print("Loading SCHEMA releases")
    published = load_published()
    browser = load_browser()

    # Symbol -> ENSG bridge. gnomAD is already the constraint source, so using
    # it here keeps the universe and the gene set on one identifier authority.
    bridge = universe.reset_index()[["gene_id", "symbol"]].dropna()
    bridge["symbol"] = bridge["symbol"].astype("string").str.strip()
    bridge = bridge.drop_duplicates("symbol", keep="first")

    pub = published.merge(bridge, on="symbol", how="inner").set_index("gene_id")
    print(
        f"\n  symbol -> ENSG: {len(pub):,} of {len(published):,} published genes "
        f"mapped ({len(published) - len(pub):,} lost to symbol drift)"
    )

    out = pd.DataFrame(index=universe.index)
    out["symbol"] = universe["symbol"]
    out["schema_p_published"] = pub["schema_p_meta"]
    out["schema_q_published"] = pub["schema_q_meta"]
    out["schema_p_browser"] = browser["browser_p_meta"]
    out["schema_q_browser"] = browser["browser_q_meta"]

    # PRIMARY: the published FDR < 0.05 set.
    out["is_schema_published_fdr"] = out["schema_q_published"] < FDR_THRESHOLD
    out["is_schema_published_exome_wide"] = out["schema_p_published"] < EXOME_WIDE_P
    # SENSITIVITY: the browser release under a locally computed BH FDR.
    out["is_schema_browser_fdr"] = out["schema_q_browser"] < FDR_THRESHOLD

    out.to_parquet(SCHEMA_SETS)
    return out


def main() -> None:
    out = build()

    print("\nSCHEMA gene sets within the 18,481-gene universe:")
    print(
        f"  PRIMARY   published, FDR < {FDR_THRESHOLD}:        "
        f"{out['is_schema_published_fdr'].sum():>5,}"
    )
    print(
        f"            published, exome-wide p < {EXOME_WIDE_P:g}: "
        f"{out['is_schema_published_exome_wide'].sum():>5,}"
    )
    print(
        f"  SENSITIVITY browser, local BH FDR < {FDR_THRESHOLD}:  "
        f"{out['is_schema_browser_fdr'].sum():>5,}"
    )

    both = (out["is_schema_published_fdr"] & out["is_schema_browser_fdr"]).sum()
    only_pub = (out["is_schema_published_fdr"] & ~out["is_schema_browser_fdr"]).sum()
    only_bro = (~out["is_schema_published_fdr"] & out["is_schema_browser_fdr"]).sum()
    print(
        f"\n  overlap: {both} in both, {only_pub} published-only, "
        f"{only_bro} browser-only"
    )

    top = out[out["is_schema_published_fdr"]].nsmallest(12, "schema_q_published")
    print("\n  most significant published SCHEMA genes:")
    for gid, row in top.iterrows():
        print(f"    {row['symbol']:<12} q={row['schema_q_published']:.3g}")

    print(f"\nWrote {SCHEMA_SETS}")


if __name__ == "__main__":
    main()
