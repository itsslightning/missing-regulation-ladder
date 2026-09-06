"""Gene-ID harmonisation is the join every rung depends on.

No two sources here share a GENCODE vintage, so bare unversioned ENSG is the
only key that reaches all six. A silent failure in this function would not
crash anything, it would just quietly drop genes from the denominator and
shift the recovery curve.
"""

from __future__ import annotations

import pandas as pd

from pipeline.s01_gene_universe import strip_version


def test_strips_a_version_suffix() -> None:
    s = pd.Series(["ENSG00000227232.5", "ENSG00000108510.12"])
    assert list(strip_version(s)) == ["ENSG00000227232", "ENSG00000108510"]


def test_is_idempotent_on_unversioned_ids() -> None:
    """gnomAD and SCHEMA ship unversioned IDs; they must pass through intact."""
    s = pd.Series(["ENSG00000164190"])
    once = strip_version(s)
    assert list(strip_version(once)) == ["ENSG00000164190"]


def test_par_y_collapses_onto_the_x_copy() -> None:
    """GTEx carries PAR_Y duplicates of pseudoautosomal genes, whose IDs look
    like ENSG00000182378.14_PAR_Y. Splitting on the first dot discards the
    PAR_Y marker, so the Y copy collapses onto the X copy's ID.

    That is why every loader deduplicates after stripping. If it did not, a
    PAR gene would appear twice in the universe and be counted twice in the
    denominator."""
    s = pd.Series(["ENSG00000182378.14_PAR_Y", "ENSG00000182378.14"])
    assert list(strip_version(s)) == ["ENSG00000182378", "ENSG00000182378"]


def test_missing_values_survive() -> None:
    s = pd.Series(["ENSG00000227232.5", None])
    out = strip_version(s)
    assert out.iloc[0] == "ENSG00000227232"
    assert pd.isna(out.iloc[1])
