"""Stage 1, step 1: decide, for every gene at every rung, whether an eQTL is detectable.

This is the y-axis of the recovery curve, so it is the file where the project
is most easily fooled. D-001 governs it.

THE UNIFORM STATISTIC
---------------------
Each rung ships a different significance column, and they are not on a common
evidential scale. GTEx and SingleBrain both ship Storey q-values, and a Storey
q depends on pi0 estimated within that study -- so q <= 0.05 is a *more*
permissive bar in a better-powered study. Using each study's native call would
therefore loosen the threshold exactly where power is highest, manufacturing
the recovery the power account predicts.

So detection is recomputed identically at every rung, in two steps:

  1. Per-gene Bonferroni over the cis variants tested for that gene. This is
     computable everywhere from what each source ships:
       GTEx         pval_nominal x num_var
       SingleBrain  Fixed_bonf (already exactly this construction)
       PsychENCODE  nominal_pval x number_of_SNPs_tested
       Bryois       min nominal p x variants tested for that gene
  2. Benjamini-Hochberg across genes within the rung (D-003).

Bonferroni over cis variants is conservative -- it ignores LD, so it is
stricter than the permutation p-value GTEx would use. That is accepted
deliberately: it is conservative *by the same construction at every rung*,
which is what a cross-rung comparison needs. Absolute eGene counts will be
below every published figure. That is the price of comparability and is stated
rather than hidden.

WHAT COUNTS AS A RUNG
---------------------
Rungs 3 and 4 contain many cell types, and rung 4 contains four times as many
as rung 3. A gene is called detected at a rung if it is detected in ANY cell
type of that rung -- which is the natural reading of "does the regulatory
signal appear at this resolution" -- but that gives rung 4 more chances purely
by having more columns. So BH is applied across ALL gene x cell-type tests
within the rung, not per cell type. A rung with 28 cell types is therefore
penalised for its 28 opportunities, and the extra detections it keeps are real
rather than multiplicity.

SENSITIVITY ARMS
----------------
D-001 and D-003 were delegated rather than chosen by the project owner, so the
rejected options run alongside:
  native      each study's own q <= 0.05 (what the papers report)
  effect      uniform significance PLUS a minimum |beta|
  by_grid     Benjamini-Yekutieli across the whole gene x rung grid

Writes: data/processed/detection_long.parquet  (gene x rung x cell)
        data/processed/detection_by_rung.parquet  (gene x rung, any-cell-type)
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from pipeline import config as cfg
from pipeline import decisions
from pipeline.provenance import Provenance
from pipeline.s01_gene_universe import GENE_UNIVERSE, strip_version
from pipeline.s02_audit_rungs import SB_EXCLUDED, SB_MAJOR, _read_gz

DETECTION_LONG = cfg.DIR_PROCESSED / "detection_long.parquet"
DETECTION_BY_RUNG = cfg.DIR_PROCESSED / "detection_by_rung.parquet"

#: The across-gene FDR level. 0.05 throughout, matching the convention every
#: source study uses, so only the *construction* differs from theirs, not the
#: nominal level.
FDR_ALPHA = 0.05

#: Sensitivity arm: minimum |beta| for the effect-size floor. GTEx slopes are
#: on inverse-normal-transformed expression and SingleBrain betas on
#: scaled/quantile-normalised expression, so both are roughly in SD units of
#: normalised expression. 0.1 SD is a small but non-trivial effect. This arm is
#: indicative, not definitive -- see D-001.
MIN_ABS_BETA = 0.1


def _bonferroni(p: pd.Series, n_var: pd.Series) -> pd.Series:
    """Per-gene Bonferroni over cis variants, capped at 1."""
    return np.minimum(1.0, p.astype(float) * n_var.astype(float))


# ---------------------------------------------------------------------------
# Per-source loaders. Each returns a long frame with one row per gene x cell.
# ---------------------------------------------------------------------------


def load_gtex(tissue: str = "Brain_Cortex") -> pd.DataFrame:
    """Rung 1. D-010 fixed this to Brain_Cortex alone."""
    path = (
        cfg.DIR_RAW / "gtex_v10" / f"{tissue}.v10.eGenes.txt.gz"
    )
    df = _read_gz(
        path,
        sep="\t",
        usecols=["gene_id", "biotype", "num_var", "pval_nominal", "qval", "slope"],
        low_memory=False,
    )
    if df is None:
        raise FileNotFoundError(f"{path} missing or incomplete")
    df["gene_id"] = strip_version(df["gene_id"])
    return pd.DataFrame(
        {
            "gene_id": df["gene_id"],
            "rung": "gtex_cortex",
            "cell": tissue,
            "p_bonf": _bonferroni(df["pval_nominal"], df["num_var"]),
            "q_native": df["qval"],
            "beta": df["slope"].abs(),
        }
    )


def load_singlebrain() -> pd.DataFrame:
    """Rungs 3 and 4.

    `Fixed_bonf` is already the per-gene Bonferroni over cis variants, so no
    recomputation is needed -- which is also the check that the construction
    used for the other rungs is the right one.

    MiGA3 is excluded: it is a different meta-analysis that additionally
    integrates isoMiGA, tests 21,059 genes against ~12,000 for the others, and
    is not a cell type of the same kind.
    """
    frames = []
    for path in sorted((cfg.DIR_RAW / "singlebrain").glob("*_top_assoc.tsv.gz")):
        cell = path.name.split("_")[0]
        if cell in SB_EXCLUDED:
            continue
        df = _read_gz(
            path,
            sep="\t",
            usecols=["feature", "Fixed_bonf", "qval", "fixed_beta"],
            low_memory=False,
        )
        if df is None:
            continue
        frames.append(
            pd.DataFrame(
                {
                    "gene_id": strip_version(df["feature"]),
                    "rung": "sn_major" if cell in SB_MAJOR else "sn_subtype",
                    "cell": cell,
                    "p_bonf": df["Fixed_bonf"].clip(upper=1.0),
                    "q_native": df["qval"],
                    "beta": df["fixed_beta"].abs(),
                }
            )
        )
    return pd.concat(frames, ignore_index=True)


def load_psychencode() -> pd.DataFrame | None:
    """Rung 2.

    Prefers the full association file, which contains every tested gene and so
    supplies a real denominator. Falls back to the Bonferroni-filtered release,
    which contains only significant genes -- usable as a numerator but flagged,
    because absence from that file cannot be distinguished from "tested and
    null".
    """
    full = cfg.DIR_RESTRICTED / "psychencode" / "Full_hg19_cis-eQTL.txt.gz"
    if full.exists() and full.stat().st_size > 3_000_000_000:
        return _load_psychencode_full(full)

    sig = cfg.DIR_RESTRICTED / "psychencode" / "DER-08b_hg38_eQTL.bonferroni.txt"
    if not sig.exists():
        return None
    df = pd.read_csv(
        sig,
        sep="\t",
        usecols=[
            "gene_id",
            "nominal_pval",
            "number_of_SNPs_tested",
            "regression_slope",
            "top_SNP",
            "FDR",
        ],
    )
    df = df[df["top_SNP"] == 1]
    df["gene_id"] = strip_version(df["gene_id"])
    out = pd.DataFrame(
        {
            "gene_id": df["gene_id"],
            "rung": "bulk_brain",
            "cell": "PsychENCODE_PFC",
            "p_bonf": _bonferroni(df["nominal_pval"], df["number_of_SNPs_tested"]),
            "q_native": df["FDR"],
            "beta": df["regression_slope"].abs(),
        }
    )
    out.attrs["partial_denominator"] = True
    return out


def _load_psychencode_full(path: Path) -> pd.DataFrame:
    """Reduce the 3.3 GB full association file to one row per gene.

    Only three columns are needed and the file does not fit comfortably in
    memory, so it is streamed in chunks and reduced as it goes: the minimum
    nominal p-value per gene, and the number of variants tested for that gene.
    Coordinates are ignored entirely, which is why the hg19 build of this file
    costs nothing -- gene-level detection needs no liftOver.
    """
    best: dict[str, float] = {}
    counts: dict[str, int] = {}
    slopes: dict[str, float] = {}

    reader = pd.read_csv(
        path,
        sep="\t",
        usecols=["gene_id", "nominal_pval", "regression_slope"],
        chunksize=5_000_000,
        low_memory=False,
    )
    for chunk in reader:
        chunk["gene_id"] = strip_version(chunk["gene_id"])
        grp = chunk.groupby("gene_id", sort=False)
        for gid, n in grp.size().items():
            counts[gid] = counts.get(gid, 0) + int(n)
        idx = grp["nominal_pval"].idxmin()
        top = chunk.loc[idx]
        for gid, p, b in zip(top["gene_id"], top["nominal_pval"], top["regression_slope"]):
            if gid not in best or p < best[gid]:
                best[gid] = float(p)
                slopes[gid] = abs(float(b))

    genes = list(best)
    df = pd.DataFrame(
        {
            "gene_id": genes,
            "rung": "bulk_brain",
            "cell": "PsychENCODE_PFC",
            "p_bonf": [min(1.0, best[g] * counts[g]) for g in genes],
            "q_native": np.nan,  # recomputed below; the full file ships no FDR
            "beta": [slopes[g] for g in genes],
        }
    )
    df["q_native"] = multipletests(
        pd.Series([best[g] for g in genes]), method="fdr_bh"
    )[1]
    return df


#: Bryois arms: the 8 cell types, plus pseudobulk from the same donors.
BRYOIS_CELLS = (
    "Astrocytes",
    "Endothelial.cells",
    "Excitatory.neurons",
    "Inhibitory.neurons",
    "Microglia",
    "OPCs...COPs",
    "Oligodendrocytes",
    "Pericytes",
)
BRYOIS_ARMS = BRYOIS_CELLS + ("pb",)


def bryois_complete_chromosomes() -> list[int]:
    """Chromosomes for which every Bryois arm is present on disk.

    The comparison is pseudobulk against cell types, so it is only valid on
    chromosomes where BOTH arms have data -- otherwise one arm would be scored
    on a gene set the other never saw. Restricting to complete chromosomes
    keeps the contrast internally valid while the remaining files download,
    which matters because there are 198 of them.
    """
    d = cfg.DIR_RAW / "bryois"
    return [
        c
        for c in range(1, 23)
        if all((d / f"{cell}.{c}.gz").exists() for cell in BRYOIS_ARMS)
    ]


def load_bryois(min_chromosomes: int = 2) -> pd.DataFrame | None:
    """The D-007 power-control arm: pseudobulk vs cell types, same donors.

    Headerless, space-separated, one file per cell type per chromosome, nominal
    p-values only. Reduced the same way as PsychENCODE: min p and variant count
    per gene, then the same Bonferroni.

    Runs on whatever chromosomes are complete across all nine arms, so a
    partial download still yields a valid -- if lower-powered -- contrast. The
    chromosomes used are recorded on the frame so the report can state them.
    """
    d = cfg.DIR_RAW / "bryois"
    chroms = bryois_complete_chromosomes()
    if len(chroms) < min_chromosomes:
        return None

    files = [d / f"{cell}.{c}.gz" for c in chroms for cell in BRYOIS_ARMS]

    acc: dict[tuple[str, str], dict[str, list]] = {}
    for path in files:
        stem = path.name[: -len(".gz")]
        cell, _, _ = stem.rpartition(".")
        cell = cell or stem
        arm = "bryois_pb" if cell == "pb" else "bryois_celltype"
        df = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            names=["gene", "snp", "tss_distance", "p_nominal", "beta"],
            engine="c",
        )
        df["gene_id"] = df["gene"].str.extract(r"(ENSG\d+)", expand=False)
        df = df.dropna(subset=["gene_id"])
        grp = df.groupby("gene_id")
        agg = grp.agg(
            p_min=("p_nominal", "min"),
            n_var=("p_nominal", "size"),
        )
        top = df.loc[grp["p_nominal"].idxmin()].set_index("gene_id")["beta"].abs()
        agg["beta"] = top
        key = (arm, cell)
        acc.setdefault(key, []).append(agg.reset_index())

    frames = []
    for (arm, cell), parts in acc.items():
        a = pd.concat(parts, ignore_index=True)
        a = a.groupby("gene_id", as_index=False).agg(
            p_min=("p_min", "min"), n_var=("n_var", "sum"), beta=("beta", "max")
        )
        frames.append(
            pd.DataFrame(
                {
                    "gene_id": a["gene_id"],
                    "rung": arm,
                    "cell": cell,
                    "p_bonf": _bonferroni(a["p_min"], a["n_var"]),
                    "q_native": np.nan,
                    "beta": a["beta"],
                }
            )
        )
    out = pd.concat(frames, ignore_index=True)
    out.attrs["chromosomes"] = chroms
    return out


# ---------------------------------------------------------------------------
# Detection calling
# ---------------------------------------------------------------------------


def call_detection(long: pd.DataFrame) -> pd.DataFrame:
    """Apply D-001 and D-003 plus the sensitivity arms, per rung."""
    if decisions.DETECTABLE_EQTL.value != "uniform_recomputed_fdr":
        raise NotImplementedError(
            f"D-001 is {decisions.DETECTABLE_EQTL.value!r}; only "
            "'uniform_recomputed_fdr' is implemented as primary."
        )
    if decisions.MULTIPLE_TESTING.value != "bh_within_rung":
        raise NotImplementedError(
            f"D-003 is {decisions.MULTIPLE_TESTING.value!r}; only "
            "'bh_within_rung' is implemented as primary."
        )

    out = []
    for rung, grp in long.groupby("rung", sort=False):
        g = grp.copy()
        ok = g["p_bonf"].notna()
        # PRIMARY: BH across every gene x cell-type test within the rung, so a
        # rung with 28 cell types pays for its 28 opportunities.
        q = pd.Series(np.nan, index=g.index)
        q[ok] = multipletests(g.loc[ok, "p_bonf"], method="fdr_bh")[1]
        g["q_uniform"] = q
        g["detected"] = g["q_uniform"] <= FDR_ALPHA

        # SENSITIVITY: each study's own call.
        g["detected_native"] = g["q_native"] <= FDR_ALPHA

        # SENSITIVITY: uniform significance plus an effect-size floor.
        g["detected_effect"] = g["detected"] & (g["beta"] >= MIN_ABS_BETA)
        out.append(g)
    return pd.concat(out, ignore_index=True)


def collapse_to_rung(called: pd.DataFrame) -> pd.DataFrame:
    """A gene is detected at a rung if detected in any of its cell types."""
    agg = called.groupby(["gene_id", "rung"], as_index=False).agg(
        detected=("detected", "any"),
        detected_native=("detected_native", "any"),
        detected_effect=("detected_effect", "any"),
        n_cells_tested=("cell", "nunique"),
        n_cells_detected=("detected", "sum"),
        best_q=("q_uniform", "min"),
        max_beta=("beta", "max"),
    )
    return agg


def main() -> None:
    prov = Provenance("s05_detection")
    universe = set(pd.read_parquet(GENE_UNIVERSE).index)

    frames = [load_gtex(), load_singlebrain()]

    pec = load_psychencode()
    if pec is None:
        prov.note("rung 2 (PsychENCODE)", "NOT AVAILABLE -- no file on disk")
    else:
        if pec.attrs.get("partial_denominator"):
            prov.note(
                "rung 2 (PsychENCODE)",
                "using the significant-only release: absence from the file "
                "cannot be distinguished from tested-and-null, so rung 2 "
                "detection is an upper bound on the numerator with an assumed "
                "denominator.",
            )
        frames.append(pec)

    bry = load_bryois()
    if bry is None:
        prov.note(
            "Bryois arm (D-007)",
            "NOT AVAILABLE -- fewer than 2 chromosomes complete across all "
            "nine arms",
        )
    else:
        chroms = bry.attrs.get("chromosomes", [])
        prov.note(
            "Bryois arm (D-007)",
            f"pseudobulk vs 8 cell types on {len(chroms)} complete "
            f"chromosome(s): {chroms}. Both arms are restricted to the same "
            f"chromosomes, so the contrast is valid; it gains power as the "
            f"remaining files arrive.",
        )
        frames.append(bry)

    long = pd.concat(frames, ignore_index=True)
    before = len(long)
    long = long[long["gene_id"].isin(universe)]
    prov.record(
        "restrict all rungs to the fixed gene universe",
        before,
        len(long),
        detail=(
            "Tests on genes outside the 18,481-gene universe are dropped. The "
            "universe is the denominator (D-012), so a test on a gene it does "
            "not contain has nowhere to land."
        ),
    )

    called = call_detection(long)
    by_rung = collapse_to_rung(called)

    called.to_parquet(DETECTION_LONG, index=False)
    by_rung.to_parquet(DETECTION_BY_RUNG, index=False)

    print("\nDetected genes per rung (of 18,481-gene universe):")
    print(f"  {'rung':<18}{'cells':>6}{'tested':>9}{'uniform':>9}{'native':>9}{'effect':>9}")
    for rung, g in by_rung.groupby("rung", sort=False):
        cells = called.loc[called["rung"] == rung, "cell"].nunique()
        print(
            f"  {rung:<18}{cells:>6}{len(g):>9,}{int(g['detected'].sum()):>9,}"
            f"{int(g['detected_native'].sum()):>9,}"
            f"{int(g['detected_effect'].sum()):>9,}"
        )
        prov.record(
            f"rung {rung}: tested -> detected (uniform)",
            len(g),
            int(g["detected"].sum()),
            detail=f"{cells} cell type(s); native call would give "
            f"{int(g['detected_native'].sum()):,}",
        )

    prov.save(cfg.DIR_LOGS / "s05_detection.json")
    print(f"\nWrote {DETECTION_LONG}\n      {DETECTION_BY_RUNG}")


if __name__ == "__main__":
    main()
