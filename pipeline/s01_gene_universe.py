"""Stage 0, step 1: build the gene universe every rung is scored against.

The recovery curve is a fraction, and a fraction needs a denominator that does
not move between rungs. If each rung supplied its own denominator -- the genes
that study happened to test -- then a cell type expressing fewer genes would
score a higher recovery fraction purely by testing fewer genes, and constrained
genes, which are more broadly expressed, would be systematically advantaged.
So the denominator is fixed here, once, and every rung is scored against it.

The universe is the set of genes with:
  - a gnomAD v2.1.1 LOEUF value (this is what "constrained" is defined from), and
  - protein-coding biotype, and
  - a GTEx v10 brain expression value (needed as a matching covariate).

Three joins have to survive gene-ID drift, because no two sources here use the
same GENCODE vintage:

  gnomAD v2.1.1   ENSG, unversioned, GENCODE v19 era
  SCHEMA          ENSG, unversioned, GRCh38
  GTEx v10        ENSG, versioned    (e.g. ENSG00000227232.5)
  SingleBrain     ENSG, versioned,   GENCODE v38
  PsychENCODE     ENSG, versioned,   GENCODE v19 era
  Bryois          SYMBOL_ENSG compound, unversioned

Stripping the version suffix and joining on bare ENSG is the only key that
reaches all six. That is not free -- genes retired or merged between GENCODE
v19 and v38 drop out -- so this script reports the loss at every join rather
than letting it disappear silently.

Writes: data/processed/gene_universe.parquet
"""

from __future__ import annotations

import gzip

import pandas as pd

from pipeline import provenance, sources
from pipeline.config import DIR_PROCESSED, DIR_RAW

GENE_UNIVERSE = DIR_PROCESSED / "gene_universe.parquet"

#: The 13 GTEx v10 brain tissues, in the order the portal lists them.
GTEX_BRAIN_TISSUES = (
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
)

#: Cortex is the reference tissue for the expression covariate: it is the
#: anatomical match to PsychENCODE prefrontal cortex and to SingleBrain's
#: neocortical nuclei. Which tissue defines *rung 1* is a separate and still
#: open question (DECISIONS.md D-010); this constant only chooses where the
#: matching covariate is read from.
EXPRESSION_REFERENCE_TISSUE = "Brain_Cortex"


def strip_version(series: pd.Series) -> pd.Series:
    """ENSG00000227232.5 -> ENSG00000227232. Idempotent on unversioned IDs."""
    return series.astype("string").str.split(".", n=1).str[0]


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------


def load_constraint() -> pd.DataFrame:
    """gnomAD v2.1.1 LoF constraint, one row per gene.

    The file is one row per *transcript* -- the canonical transcript per gene --
    so it is already gene-level, but a handful of genes appear twice where two
    transcripts were scored. The lower LOEUF (more constrained) is kept, which
    matches how gnomAD's own browser presents a gene.

    Note that `brain_expression` is dropped rather than used: it is NA for all
    19,704 rows in this release, despite being the obvious candidate for the
    expression covariate. GTEx median TPM supplies that instead.
    """
    src = sources.GNOMAD_CONSTRAINT
    path = src.directory / "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
    if not path.exists():
        path = provenance.fetch(src)

    keep = {
        "gene": "symbol",
        "gene_id": "gene_id",
        "oe_lof_upper": "loeuf",
        "oe_lof_upper_bin": "loeuf_decile",
        "pLI": "pli",
        "num_coding_exons": "n_coding_exons",
        "cds_length": "cds_length",
        "gene_length": "gene_length",
        "gene_type": "gene_type",
        "chromosome": "chrom",
    }
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", usecols=list(keep), low_memory=False)
    df = df.rename(columns=keep)

    before = len(df)
    df = df[df["gene_type"] == "protein_coding"]
    df = df.dropna(subset=["loeuf"])
    df = df.sort_values("loeuf").drop_duplicates("gene_id", keep="first")
    df["gene_id"] = strip_version(df["gene_id"])

    print(f"  constraint: {before} rows -> {len(df)} protein-coding genes with LOEUF")
    return df.set_index("gene_id")


def load_brain_expression() -> pd.DataFrame:
    """GTEx v10 median TPM across the 13 brain tissues.

    Two covariates come out of this: the reference-tissue TPM used for
    matching, and the number of brain tissues in which the gene is detectably
    expressed (TPM >= 1), which is a crude breadth-of-expression measure and a
    useful sanity check on the matching.
    """
    src = sources.GTEX_V10
    path = src.directory / "GTEx_v10_gene_median_tpm.gct.gz"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing. It is fetched separately from the eQTL archive; "
            "see the Stage 0 report."
        )

    # .gct format: two header lines before the real header row.
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", skiprows=2, low_memory=False)

    df["gene_id"] = strip_version(df["Name"])
    brain = [t for t in GTEX_BRAIN_TISSUES if t in df.columns]
    out = pd.DataFrame(
        {
            "gene_id": df["gene_id"],
            "tpm_reference": df[EXPRESSION_REFERENCE_TISSUE],
            "n_brain_tissues_expressed": (df[brain] >= 1).sum(axis=1),
        }
    )
    # The GCT carries PAR_Y duplicates of X-linked genes; keep the first.
    out = out.drop_duplicates("gene_id", keep="first")
    print(f"  expression: {len(out)} genes, {len(brain)} brain tissues")
    return out.set_index("gene_id")


def load_schema() -> pd.DataFrame:
    """SCHEMA per-gene exome association statistics.

    Deliberately does NOT define "a SCHEMA gene". The browser release ships no
    q-value column, so any gene set requires choosing a p-value column and a
    correction -- and choosing between this release and the published Singh et
    al. 2022 set (DECISIONS.md D-009). Both are conclusion-shaping, so this
    loader carries the statistics through and stops there.
    """
    src = sources.SCHEMA_GENES
    path = src.directory / "SCHEMA_gene_results.tsv.bgz"
    if not path.exists():
        path = provenance.fetch(src, filename="SCHEMA_gene_results.tsv.bgz")

    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", low_memory=False)

    df = df[df["group"] == "meta"]
    keep = {
        "gene_id": "gene_id",
        "gene_symbol": "schema_symbol",
        "ptv_p_value": "schema_ptv_p",
        "schema_case_control_p_value": "schema_cc_p",
        "case_control_plus_de_novo_p_value": "schema_cc_dn_p",
        "n_cases": "schema_n_cases",
        "n_controls": "schema_n_controls",
    }
    df = df[list(keep)].rename(columns=keep)
    df["gene_id"] = strip_version(df["gene_id"])
    df = df.drop_duplicates("gene_id", keep="first")
    n_cases = df["schema_n_cases"].dropna().max()
    print(f"  schema: {len(df)} genes, max n_cases={n_cases:,.0f}")
    return df.set_index("gene_id")


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build() -> pd.DataFrame:
    print("Loading annotation sources")
    constraint = load_constraint()
    expression = load_brain_expression()
    schema = load_schema()

    print("\nJoining on unversioned ENSG")
    universe = constraint.join(expression, how="inner")
    print(
        f"  constraint x expression: {len(constraint)} -> {len(universe)} "
        f"({len(constraint) - len(universe)} lost to GENCODE drift)"
    )

    universe = universe.join(schema, how="left")
    matched = universe["schema_cc_p"].notna().sum()
    print(
        f"  + schema (left join): {matched} of {len(universe)} universe genes "
        f"carry SCHEMA statistics ({len(universe) - matched} without)"
    )

    # Constraint labels. These are gnomAD's and SCHEMA's own conventional cut
    # points, not a threshold invented here -- see config.py.
    from pipeline.config import LOEUF_CONSTRAINED_MAX, PLI_CONSTRAINED_MIN

    universe["is_constrained_loeuf"] = universe["loeuf"] < LOEUF_CONSTRAINED_MAX
    universe["is_constrained_pli"] = universe["pli"] >= PLI_CONSTRAINED_MIN

    universe = universe.sort_values("loeuf")
    DIR_PROCESSED.mkdir(parents=True, exist_ok=True)
    universe.to_parquet(GENE_UNIVERSE)
    return universe


def main() -> None:
    universe = build()

    print(f"\nGene universe: {len(universe):,} genes -> {GENE_UNIVERSE}")
    print(
        f"  constrained (LOEUF < 0.35):  {universe['is_constrained_loeuf'].sum():,}\n"
        f"  constrained (pLI >= 0.9):    {universe['is_constrained_pli'].sum():,}"
    )
    print("\n  LOEUF decile sizes:")
    counts = universe["loeuf_decile"].value_counts().sort_index()
    for decile, n in counts.items():
        print(f"    decile {int(decile)}: {n:,}")

    from pipeline import decisions

    print("\nDecision status:")
    print(decisions.summary())


if __name__ == "__main__":
    main()
