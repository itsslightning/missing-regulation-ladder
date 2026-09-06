"""Stage 2, step 2: SMR of SCZ risk on brain gene expression, per rung.

D-004 chose SMR over coloc because it is the only method computable identically
at every rung -- see DECISIONS.md for why coloc is not, and for the cost of
that choice.

THE STATISTIC
-------------
SMR (Zhu et al. 2016) propagates a variant's effect on expression through the
same variant's effect on the trait:

    b_xy = b_gwas / b_eqtl

and tests it with

    T_SMR = (z_x^2 * z_y^2) / (z_x^2 + z_y^2)  ~  chi-square, 1 df

where z_x and z_y are the eQTL and GWAS z-statistics. The chi-square form is
used rather than the delta-method variance of the ratio because it is exact and
does not degrade when b_eqtl is near zero.

WHAT SMR DOES NOT DO, AND WHY IT IS STILL THE RIGHT CHOICE HERE
--------------------------------------------------------------
It cannot distinguish a shared causal variant from linkage between two distinct
causal variants, so counts are inflated relative to true colocalization. That is
accepted deliberately: the bias has the same construction at every rung, and
Stage 2 asks a CROSS-RUNG question -- how many more loci gain an eQTL
explanation as the assay changes -- not an absolute count comparable to a
published coloc figure. HEIDI, the usual remedy, needs regional data that only
two arms have; it runs separately.

THREE THINGS THAT SILENTLY BREAK SMR, ALL HANDLED HERE
------------------------------------------------------
1. **Allele orientation.** If the eQTL and GWAS report effects against opposite
   alleles and it is not detected, b_xy flips sign -- a risk-increasing gene
   reads as protective. Every variant is harmonised on its allele pair and
   dropped if the pair cannot be reconciled.
2. **Strand-ambiguous variants.** A/T and C/G variants look identical on either
   strand, so allele matching alone cannot orient them. They are excluded from
   the primary and counted, rather than silently kept.
3. **Weak instruments.** b_xy is a ratio, so a near-zero denominator explodes
   it. SMR's convention is to test only genes with an eQTL at p < 5e-8, which
   is applied here and reported per rung -- the number of TESTABLE genes is
   itself part of the answer.

Writes: data/processed/smr_results.parquet
"""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

from pipeline import config as cfg
from pipeline import sources
from pipeline.provenance import Provenance
from pipeline.s12_smr_exposures import OUT as EXPOSURES

OUT = cfg.DIR_PROCESSED / "smr_results.parquet"

#: SMR's conventional instrument threshold. A gene needs a cis-eQTL at least
#: this strong before the ratio is stable enough to test.
P_EQTL_INSTRUMENT = 5e-8

#: FDR level for SMR significance, matching D-003's within-rung BH.
FDR_ALPHA = 0.05

#: Strand-ambiguous pairs: indistinguishable after a strand flip.
PALINDROMIC = {frozenset("AT"), frozenset("CG")}

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def load_pgc3() -> pd.DataFrame:
    """PGC3 SCZ wave 3, European autosomes.

    'PGCsumstatsVCFv1.0': 73 '##' metadata lines carrying UTF-8 text, then a
    tab-separated header. A1 is the effect allele for BETA.
    """
    path = sources.PGC3_SCZ.directory / (
        "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
    )
    n_meta = 0
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            n_meta += 1

    df = pd.read_csv(
        path,
        sep="\t",
        skiprows=n_meta,
        usecols=["CHROM", "POS", "ID", "A1", "A2", "BETA", "SE", "PVAL"],
        dtype={"ID": "string", "A1": "string", "A2": "string"},
        encoding="utf-8",
        engine="c",
    )
    df = df.rename(
        columns={
            "CHROM": "chrom",
            "POS": "pos",
            "ID": "rsid",
            "A1": "gwas_effect_allele",
            "A2": "gwas_other_allele",
            "BETA": "gwas_beta",
            "SE": "gwas_se",
            "PVAL": "gwas_p",
        }
    )
    for c in ("gwas_effect_allele", "gwas_other_allele"):
        df[c] = df[c].str.upper()
    return df.dropna(subset=["gwas_beta", "gwas_se"]).drop_duplicates(
        "rsid", keep="first"
    )


def harmonise(df: pd.DataFrame) -> pd.DataFrame:
    """Orient the GWAS effect to the eQTL effect allele.

    Returns the frame with `gwas_beta_aligned` and a `harmonisation` label, and
    with unresolvable or strand-ambiguous variants marked rather than removed,
    so the counts can be reported.
    """
    e1, e2 = df["effect_allele"], df["other_allele"]
    g1, g2 = df["gwas_effect_allele"], df["gwas_other_allele"]

    same = (e1 == g1) & (e2 == g2)
    flipped = (e1 == g2) & (e2 == g1)

    comp = lambda s: s.map(COMPLEMENT).astype("string")  # noqa: E731
    same_strandflip = (e1 == comp(g1)) & (e2 == comp(g2))
    flip_strandflip = (e1 == comp(g2)) & (e2 == comp(g1))

    palindromic = [
        frozenset({a, b}) in PALINDROMIC if isinstance(a, str) and isinstance(b, str)
        else False
        for a, b in zip(df["effect_allele"], df["other_allele"])
    ]
    df = df.assign(palindromic=palindromic)

    label = pd.Series("unresolved", index=df.index, dtype="object")
    label[same] = "same"
    label[flipped] = "flipped"
    label[same_strandflip & ~same & ~flipped] = "strand"
    label[flip_strandflip & ~same & ~flipped] = "strand_flipped"

    sign = pd.Series(np.nan, index=df.index, dtype=float)
    sign[label == "same"] = 1.0
    sign[label == "flipped"] = -1.0
    sign[label == "strand"] = 1.0
    sign[label == "strand_flipped"] = -1.0

    df["harmonisation"] = label
    df["gwas_beta_aligned"] = df["gwas_beta"] * sign
    return df


def smr(df: pd.DataFrame) -> pd.DataFrame:
    """SMR chi-square and the causal estimate b_xy."""
    z_x = df["beta"] / df["se"]
    z_y = df["gwas_beta_aligned"] / df["gwas_se"]
    zx2, zy2 = z_x**2, z_y**2
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (zx2 * zy2) / (zx2 + zy2)
        b_xy = df["gwas_beta_aligned"] / df["beta"]
        # Delta method, treating the two studies as independent.
        se_xy = np.abs(b_xy) * np.sqrt(
            (df["gwas_se"] / df["gwas_beta_aligned"]) ** 2
            + (df["se"] / df["beta"]) ** 2
        )
    df = df.assign(
        smr_chi2=t,
        smr_p=chi2.sf(t, df=1),
        b_xy=b_xy,
        se_xy=se_xy,
    )
    return df


def main() -> None:
    prov = Provenance("s13_smr")

    expo = pd.read_parquet(EXPOSURES)
    gwas = load_pgc3()
    prov.note(
        "PGC3 outcome loaded",
        f"{len(gwas):,} variants, European autosomes, 52,017 cases / 75,889 "
        "controls (Trubetskoy et al. 2022, public v3 release)",
    )

    before = len(expo)
    expo = expo[expo["p_nominal"] < P_EQTL_INSTRUMENT]
    prov.record(
        f"exposures passing the SMR instrument threshold (p < {P_EQTL_INSTRUMENT:g})",
        before,
        len(expo),
        detail=(
            "SMR estimates a ratio, so a weak instrument makes it explode. The "
            "number of testable genes differs by rung and is part of the "
            "Stage 2 answer, not a nuisance."
        ),
    )

    merged = expo.merge(gwas, on="rsid", how="inner")
    prov.record(
        "exposures matched to a PGC3 variant by rsID",
        len(expo),
        len(merged),
        detail=(
            "rsID join. Unmatched variants are typically absent from the GWAS "
            "imputation panel."
        ),
    )

    merged = harmonise(merged)
    counts = merged["harmonisation"].value_counts().to_dict()
    prov.note(
        "allele harmonisation",
        ", ".join(f"{k}={v:,}" for k, v in sorted(counts.items())),
    )

    usable = merged["harmonisation"] != "unresolved"
    prov.record(
        "variants with a reconcilable allele pair",
        len(merged),
        int(usable.sum()),
        detail=(
            "An unreconciled pair cannot be oriented, and guessing would flip "
            "the sign of b_xy for an unknown subset."
        ),
    )
    merged = merged[usable]

    n_palindromic = int(merged["palindromic"].sum())
    merged = merged[~merged["palindromic"]]
    prov.record(
        "variants surviving the strand-ambiguity filter",
        len(merged) + n_palindromic,
        len(merged),
        detail=(
            f"{n_palindromic:,} A/T or C/G variants excluded: they are "
            "identical on either strand, so allele matching alone cannot "
            "orient them."
        ),
    )

    res = smr(merged)
    res = res[res["smr_p"].notna()]

    # D-003: BH across all gene x cell-type tests within a rung.
    res["smr_q"] = np.nan
    for rung, idx in res.groupby("rung").groups.items():
        res.loc[idx, "smr_q"] = multipletests(
            res.loc[idx, "smr_p"], method="fdr_bh"
        )[1]
    res["smr_sig"] = res["smr_q"] <= FDR_ALPHA

    res.to_parquet(OUT, index=False)

    print("\nSMR results per rung (BH within rung, D-003):\n")
    print(f"  {'rung':<14}{'tests':>9}{'genes':>9}{'sig tests':>11}"
          f"{'sig genes':>11}")
    for rung, g in res.groupby("rung", sort=False):
        sig = g[g["smr_sig"]]
        print(
            f"  {rung:<14}{len(g):>9,}{g['gene_id'].nunique():>9,}"
            f"{len(sig):>11,}{sig['gene_id'].nunique():>11,}"
        )
        prov.record(
            f"SMR {rung}: tested -> significant genes",
            g["gene_id"].nunique(),
            sig["gene_id"].nunique(),
            detail=f"{len(g):,} gene x cell tests, {len(sig):,} significant",
        )

    prov.save(cfg.DIR_LOGS / "s13_smr.json")
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
