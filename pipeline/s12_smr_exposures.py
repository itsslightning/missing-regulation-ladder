"""Stage 2, step 1: the eQTL side of SMR, harmonised to one variant key.

SMR needs, per gene per rung, the top cis-eQTL variant and its effect and
standard error. This builds that table for every arm, keyed on rsID, and does
so without any GWAS input, so it runs before PGC3 access arrives and is the
thing PGC3 is joined onto when it does.

WHY rsID
--------
The join key for Stage 2 is the VARIANT, not the gene, and the arms do not
agree on how to name one:

    SingleBrain   variant_id IS an rsID          (plus chr/pos on GRCh38)
    GTEx v10      chr_pos_ref_alt_b38, plus rs_id_dbSNP155_GRCh38p13
    Bryois        rsID native                    (positions on GRCh37)
    PsychENCODE   "chr:pos" on hg19, no rsID

rsID reaches three of the four directly and is the only key that does. It also
sidesteps a genome-build problem that would otherwise need liftOver: the arms
span GRCh37 and GRCh38, and PGC3 is on GRCh37. An rsID means the same variant
in both builds.

PsychENCODE needs a position -> rsID lookup, and Bryois's own snp_pos.txt.gz
happens to be exactly that: rsID with GRCh37 coordinates, which is the build
PsychENCODE's full file uses. Using it costs nothing extra and is recorded as
the bridge.

STANDARD ERRORS
---------------
GTEx and SingleBrain ship one. PsychENCODE and Bryois ship an effect and a
nominal p-value, so it is recovered as se = |beta| / |z|, with z the two-sided
normal quantile of p. That is exact for a Wald test, which is what all four
studies report, and it is the same reconstruction each time.

Writes: data/processed/smr_exposures.parquet
"""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
from scipy.stats import norm

from pipeline import config as cfg
from pipeline import sources
from pipeline.provenance import Provenance
from pipeline.s01_gene_universe import GENE_UNIVERSE, strip_version
from pipeline.s02_audit_rungs import SB_EXCLUDED, SB_MAJOR, _read_gz

OUT = cfg.DIR_PROCESSED / "smr_exposures.parquet"

#: Below this p-value the normal quantile saturates in double precision and the
#: reconstructed SE would be garbage. Such variants are kept but flagged, since
#: an eQTL that extreme is significant under any method and its exact SE does
#: not change an SMR call.
P_UNDERFLOW = 1e-300


def se_from_beta_p(beta: pd.Series, p: pd.Series) -> pd.Series:
    """Recover a Wald standard error from an effect size and a two-sided p."""
    p = p.astype(float).clip(lower=P_UNDERFLOW, upper=1.0)
    z = np.abs(norm.isf(p / 2.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.abs(beta.astype(float)) / z
    return pd.Series(np.where(np.isfinite(se) & (z > 0), se, np.nan), index=beta.index)


def load_gtex(tissue: str = "Brain_Cortex") -> pd.DataFrame:
    df = _read_gz(
        cfg.DIR_RAW / "gtex_v10" / f"{tissue}.v10.eGenes.txt.gz",
        sep="\t",
        usecols=[
            "gene_id", "rs_id_dbSNP155_GRCh38p13", "variant_id", "ref", "alt",
            "slope", "slope_se", "pval_nominal", "num_var",
        ],
        low_memory=False,
    )
    if df is None:
        return pd.DataFrame()
    return pd.DataFrame(
        {
            "gene_id": strip_version(df["gene_id"]),
            "rung": "gtex_cortex",
            "cell": tissue,
            "rsid": df["rs_id_dbSNP155_GRCh38p13"].astype("string"),
            # GTEx reports `slope` with respect to the ALT allele.
            "effect_allele": df["alt"].astype("string").str.upper(),
            "other_allele": df["ref"].astype("string").str.upper(),
            "beta": df["slope"].astype(float),
            "se": df["slope_se"].astype(float),
            "p_nominal": df["pval_nominal"].astype(float),
            "n_var": df["num_var"].astype("float"),
            "se_source": "reported",
        }
    )


def load_singlebrain() -> pd.DataFrame:
    frames = []
    for path in sorted((cfg.DIR_RAW / "singlebrain").glob("*_top_assoc.tsv.gz")):
        cell = path.name.split("_")[0]
        if cell in SB_EXCLUDED:
            continue
        df = _read_gz(
            path,
            sep="\t",
            usecols=["feature", "variant_id", "ref", "alt", "Allele",
                     "fixed_beta", "fixed_sd", "Fixed_P", "Fixed_bonf"],
            low_memory=False,
        )
        if df is None:
            continue
        # Fixed_bonf = Fixed_P x n_var, so the variant count is recoverable.
        with np.errstate(divide="ignore", invalid="ignore"):
            n_var = df["Fixed_bonf"].astype(float) / df["Fixed_P"].astype(float)
        # `Allele` is the tested (effect) allele; the other is whichever of
        # ref/alt it is not.
        eff = df["Allele"].astype("string").str.upper()
        ref = df["ref"].astype("string").str.upper()
        alt = df["alt"].astype("string").str.upper()
        oth = ref.where(eff != ref, alt)
        frames.append(
            pd.DataFrame(
                {
                    "gene_id": strip_version(df["feature"]),
                    "rung": "sn_major" if cell in SB_MAJOR else "sn_subtype",
                    "cell": cell,
                    "rsid": df["variant_id"].astype("string"),
                    "effect_allele": eff,
                    "other_allele": oth,
                    "beta": df["fixed_beta"].astype(float),
                    "se": df["fixed_sd"].astype(float),
                    "p_nominal": df["Fixed_P"].astype(float),
                    "n_var": n_var,
                    "se_source": "reported",
                }
            )
        )
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


VARIANT_MAP = cfg.DIR_PROCESSED / "variant_map.parquet"


def load_variant_map() -> pd.DataFrame:
    """rsID, both genome builds, and the allele pair: the Stage 2 keystone.

    Bryois's snp_pos.txt.gz turns out to carry far more than positions:

        SNP  SNP_id_hg38  SNP_id_hg19  effect_allele  other_allele  MAF

    That makes it a complete harmonisation table for this project. Three things
    fall out of it at once:

    1. **The PsychENCODE bridge.** Its full file gives chr:pos on hg19 only,
       and this maps hg19 positions to rsIDs.
    2. **Build independence.** The arms span GRCh37 and GRCh38 and PGC3 is
       GRCh37; having both columns keyed to one rsID removes any need for
       liftOver.
    3. **Allele alignment, which SMR cannot be correct without.** An SMR
       estimate is an eQTL effect propagated through a GWAS effect at the same
       variant. If the two studies report effects relative to opposite alleles
       and that is not detected, the estimate's SIGN flips, turning a gene
       whose increased expression raises risk into one that appears protective.
       That is a silent, plausible-looking error, so the allele pair is carried
       from here rather than assumed.
    """
    path = cfg.DIR_RAW / "bryois" / "snp_pos.txt.gz"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(
        path,
        sep="\t",
        dtype="string",
        engine="c",
    )
    df = df.rename(
        columns={
            "SNP": "rsid",
            "SNP_id_hg38": "pos_hg38",
            "SNP_id_hg19": "pos_hg19",
            "effect_allele": "effect_allele",
            "other_allele": "other_allele",
            "MAF": "maf",
        }
    )
    keep = [
        c
        for c in ("rsid", "pos_hg38", "pos_hg19", "effect_allele",
                  "other_allele", "maf")
        if c in df.columns
    ]
    df = df[keep].dropna(subset=["rsid"]).drop_duplicates("rsid", keep="first")
    df["maf"] = pd.to_numeric(df["maf"], errors="coerce")
    return df


def load_psychencode(vmap: pd.DataFrame) -> pd.DataFrame:
    """Rung 2 exposures, via the variant map.

    The full PsychENCODE file identifies variants as "chr:pos" on hg19 with no
    rsID and no alleles, so on its own it cannot enter an SMR analysis keyed on
    rsID and harmonised on alleles. The variant map supplies both.

    Streamed in chunks like s05's reducer, but keeping the SNP identity of the
    best variant rather than only its p-value. SMR needs to know WHICH
    variant the instrument is.
    """
    from pipeline.s05_detection import PSYCHENCODE_FULL_COLUMNS

    path = cfg.DIR_RESTRICTED / "psychencode" / sources.PSYCHENCODE_FULL_FILE
    if not path.exists() or path.stat().st_size != sources.PSYCHENCODE_FULL_BYTES:
        return pd.DataFrame()

    best: dict[str, tuple[float, float, str, int]] = {}
    reader = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=PSYCHENCODE_FULL_COLUMNS,
        usecols=["gene_id", "SNP_id", "number_of_SNPs_tested",
                 "nominal_pval", "regression_slope"],
        dtype={"gene_id": "string", "SNP_id": "string",
               "number_of_SNPs_tested": "int32"},
        chunksize=5_000_000,
        engine="c",
    )
    for chunk in reader:
        chunk["gene_id"] = strip_version(chunk["gene_id"])
        top = chunk.loc[chunk.groupby("gene_id", sort=False)["nominal_pval"].idxmin()]
        for gid, p, b, snp, n in zip(
            top["gene_id"], top["nominal_pval"], top["regression_slope"],
            top["SNP_id"], top["number_of_SNPs_tested"],
        ):
            cur = best.get(gid)
            if cur is None or p < cur[0]:
                best[gid] = (float(p), float(b), str(snp), int(n))

    if not best:
        return pd.DataFrame()
    genes = list(best)
    df = pd.DataFrame(
        {
            "gene_id": genes,
            "snp_hg19": ["chr" + best[g][2] for g in genes],
            "p_nominal": [best[g][0] for g in genes],
            "beta": [best[g][1] for g in genes],
            "n_var": [float(best[g][3]) for g in genes],
        }
    )
    # "3:28130472" -> "chr3:28130472" to match the map's SNP_id_hg19 form.
    m = vmap[["rsid", "pos_hg19", "effect_allele", "other_allele"]].rename(
        columns={"pos_hg19": "snp_hg19"}
    )
    df = df.merge(m, on="snp_hg19", how="inner")
    df["se"] = se_from_beta_p(df["beta"], df["p_nominal"])
    return pd.DataFrame(
        {
            "gene_id": df["gene_id"],
            "rung": "bulk_brain",
            "cell": "PsychENCODE_PFC",
            "rsid": df["rsid"].astype("string"),
            "effect_allele": df["effect_allele"].str.upper(),
            "other_allele": df["other_allele"].str.upper(),
            "beta": df["beta"],
            "se": df["se"],
            "p_nominal": df["p_nominal"],
            "n_var": df["n_var"],
            "se_source": "derived_from_p",
        }
    )


def main() -> None:
    prov = Provenance("s12_smr_exposures")
    universe = set(pd.read_parquet(GENE_UNIVERSE).index)

    vmap = load_variant_map()
    frames = [load_gtex(), load_singlebrain()]
    if not vmap.empty:
        pec = load_psychencode(vmap)
        if not pec.empty:
            frames.append(pec)
    frames = [f for f in frames if not f.empty]
    if not frames:
        raise RuntimeError("no exposure arms could be built")
    expo = pd.concat(frames, ignore_index=True)

    before = len(expo)
    expo = expo[expo["gene_id"].isin(universe)]
    prov.record(
        "restrict exposures to the fixed gene universe",
        before,
        len(expo),
        detail="same denominator as the recovery curve (D-012)",
    )

    have_rsid = expo["rsid"].notna() & expo["rsid"].str.startswith("rs", na=False)
    prov.record(
        "exposures carrying a usable rsID",
        len(expo),
        int(have_rsid.sum()),
        detail=(
            "rsID is the Stage 2 join key: it reaches SingleBrain, GTEx and "
            "Bryois directly and is build-independent, which matters because "
            "the arms span GRCh37 and GRCh38 and PGC3 is GRCh37."
        ),
        dropped=sorted(expo.loc[~have_rsid, "gene_id"].unique())[:8],
    )
    expo = expo[have_rsid]

    usable_se = expo["se"].notna() & (expo["se"] > 0)
    prov.record(
        "exposures with a usable standard error",
        len(expo),
        int(usable_se.sum()),
        detail="SMR needs beta and se; rows without one cannot be tested",
    )
    expo = expo[usable_se]

    if vmap.empty:
        prov.note("variant map", "snp_pos.txt.gz not available")
    else:
        vmap.to_parquet(VARIANT_MAP, index=False)
        covered = expo["rsid"].isin(set(vmap["rsid"]))
        prov.record(
            "exposure variants covered by the variant map",
            len(expo),
            int(covered.sum()),
            detail=(
                f"{len(vmap):,} variants with rsID, hg19 and hg38 positions, "
                "and the allele pair. The allele pair is the part SMR cannot "
                "be correct without: if the eQTL and GWAS report effects "
                "against opposite alleles and that is missed, the SMR estimate "
                "flips sign and a risk-increasing gene reads as protective. "
                "Uncovered variants are not lost. PGC3 supplies its own "
                "alleles, but they lose the independent cross-check."
            ),
        )

    expo.to_parquet(OUT, index=False)

    print("SMR exposures, one row per gene per cell type:\n")
    print(f"  {'rung':<14}{'cells':>6}{'genes':>9}{'rows':>9}{'median |z|':>12}")
    for rung, g in expo.groupby("rung", sort=False):
        z = (g["beta"].abs() / g["se"]).median()
        print(
            f"  {rung:<14}{g['cell'].nunique():>6}{g['gene_id'].nunique():>9,}"
            f"{len(g):>9,}{z:>12.2f}"
        )
    print(f"\n  total variants referenced: {expo['rsid'].nunique():,}")

    prov.save(cfg.DIR_LOGS / "s12_smr_exposures.json")
    print(f"\nWrote {OUT}")
    print(
        "\nStage 2 is now blocked only on the PGC3 outcome side: per-variant "
        "beta, se and effect allele keyed by rsID."
    )


if __name__ == "__main__":
    main()
