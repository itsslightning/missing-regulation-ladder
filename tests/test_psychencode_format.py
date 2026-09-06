"""The PsychENCODE full-association file's format is load-bearing and undocumented.

Two releases of the same data disagree about their own format:

  DER-08b_hg38_eQTL.bonferroni.txt   headed, TAB-separated, 15 columns
  Full_hg19_cis-eQTL.txt.gz          NO header, WHITESPACE-separated, 14 columns

The loader was originally written against the first and pointed at the second,
which would have failed only after a 3.3 GB download finished. These tests pin
the layout so that regression is caught in half a second instead.
"""

from __future__ import annotations

import gzip

import pandas as pd
import pytest

from pipeline.s05_detection import (
    PSYCHENCODE_FULL_COLUMNS,
    _load_psychencode_full,
)

#: Two real lines from the file, same gene, verbatim.
SAMPLE = (
    "ENSG00000215004.3 chr3 29128887 29128887 + 4558 -998415 3:28130472 "
    "chr3 28130472 28130472 0.935933 0.00173423 0\n"
    "ENSG00000215004.3 chr3 29128887 29128887 + 4558 -996771 3:28132116 "
    "chr3 28132116 28132116 0.445269 0.0164145 0\n"
)


def test_column_layout_is_fourteen_fields() -> None:
    assert len(PSYCHENCODE_FULL_COLUMNS) == 14
    first = SAMPLE.splitlines()[0].split()
    assert len(first) == len(PSYCHENCODE_FULL_COLUMNS)


def test_the_columns_the_loader_needs_are_present() -> None:
    for needed in (
        "gene_id",
        "number_of_SNPs_tested",
        "nominal_pval",
        "regression_slope",
    ):
        assert needed in PSYCHENCODE_FULL_COLUMNS


@pytest.fixture
def sample_file(tmp_path):
    p = tmp_path / "Full_hg19_cis-eQTL.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(SAMPLE)
    return p


def test_reduces_to_one_row_per_gene(sample_file) -> None:
    df = _load_psychencode_full(sample_file)
    assert len(df) == 1
    assert df["gene_id"].iloc[0] == "ENSG00000215004"


def test_takes_the_best_variant_and_the_files_own_variant_count(
    sample_file,
) -> None:
    """p_bonf must be min(p) x the study's number_of_SNPs_tested, capped at 1.

    Counting rows instead would give 2 here rather than 4558, the bug this
    pins. With the real count the correction saturates, which is correct: one
    weak variant among 4,558 tested is not evidence.
    """
    df = _load_psychencode_full(sample_file)
    assert df["p_bonf"].iloc[0] == 1.0  # 0.445269 * 4558, capped
    assert df["beta"].iloc[0] == pytest.approx(0.0164145)


def test_a_tab_separated_headed_read_would_not_work() -> None:
    """Guards the assumption that broke: this file is not a headed TSV."""
    parsed = pd.read_csv(
        pd.io.common.StringIO(SAMPLE), sep="\t", header=None, nrows=1
    )
    # Tab-separated parsing collapses the whole record into a single column.
    assert parsed.shape[1] == 1
