"""Data access for the dashboard.

The app computes nothing. It reads the tables the pipeline already wrote, so
what is on screen is exactly what the reproducible pipeline produced and the
deployed app needs no data pull at runtime. Same contract as the sibling
scz-target-prioritization dashboard.

Every number the app shows therefore traces to a parquet file in
data/processed/, which traces to a stage script, which traces to a decision in
DECISIONS.md.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
LOGS = ROOT / "logs"
DOCS = ROOT / "docs"

#: Ladder order and display labels. The rung order IS the experiment, so it is
#: fixed here rather than inferred from whatever a table happens to contain.
RUNG_ORDER = ["gtex_cortex", "bulk_brain", "sn_major", "sn_subtype"]

RUNG_LABEL = {
    "gtex_cortex": "GTEx cortex",
    "bulk_brain": "PsychENCODE",
    "sn_major": "SingleBrain, 7 types",
    "sn_subtype": "SingleBrain, 28 subtypes",
    "bryois_pb": "Bryois pseudobulk",
    "bryois_celltype": "Bryois, 8 cell types",
}

RUNG_ASSAY = {
    "gtex_cortex": "bulk tissue",
    "bulk_brain": "bulk tissue",
    "sn_major": "single-nucleus",
    "sn_subtype": "single-nucleus",
    "bryois_pb": "single-nucleus",
    "bryois_celltype": "single-nucleus",
}

RUNG_DONORS = {
    "gtex_cortex": 205,
    "bulk_brain": 1387,
    "sn_major": 983,
    "sn_subtype": 983,
    "bryois_pb": 192,
    "bryois_celltype": 192,
}

#: Detection arms. The primary is D-001 as recorded; the others are the
#: options it was chosen over, shipped so the choice can be inspected.
ARM_LABEL = {
    "uniform": "Uniform Bonferroni + BH  (D-001, primary)",
    "native": "Each study's own q ≤ 0.05  (sensitivity)",
    "effect": "Uniform + |β| ≥ 0.1  (sensitivity)",
}

#: The same three arms named for a narrow sidebar control, where the full
#: label wraps to three lines and stops being readable.
ARM_SHORT = {
    "uniform": "Uniform Bonferroni + BH",
    "native": "Each study's own q ≤ 0.05",
    "effect": "Uniform + |β| ≥ 0.1",
}

#: Which arm D-001 actually selected. Everything else is a shipped sensitivity.
PRIMARY_ARM = "uniform"

#: The sibling repo's Streamlit primaryColor, plus a warm contrast.
TEAL = "#3d7a6f"
RUST = "#b5651d"
GREY = "#6b6b6b"


def _read(name: str) -> pd.DataFrame:
    path = PROCESSED / name
    if not path.exists():
        return pd.DataFrame()
    return pd.read_parquet(path)


@st.cache_data(show_spinner=False)
def recovery() -> pd.DataFrame:
    return _read("recovery_by_rung.parquet")


@st.cache_data(show_spinner=False)
def gap_closure() -> pd.DataFrame:
    return _read("gap_closure.parquet")


@st.cache_data(show_spinner=False)
def assay_contrast() -> pd.DataFrame:
    return _read("assay_contrast.parquet")


@st.cache_data(show_spinner=False)
def resolution_test() -> pd.DataFrame:
    return _read("resolution_test.parquet")


@st.cache_data(show_spinner=False)
def robustness() -> pd.DataFrame:
    return _read("robustness.parquet")


@st.cache_data(show_spinner=False)
def by_loeuf() -> pd.DataFrame:
    return _read("recovery_by_loeuf.parquet")


@st.cache_data(show_spinner=False)
def by_loeuf_expr() -> pd.DataFrame:
    return _read("recovery_by_loeuf_expr.parquet")


@st.cache_data(show_spinner=False)
def by_pli() -> pd.DataFrame:
    return _read("recovery_by_pli.parquet")


@st.cache_data(show_spinner=False)
def schema_switch() -> pd.DataFrame:
    return _read("schema_switch.parquet")


@st.cache_data(show_spinner=False)
def schema_recovery() -> pd.DataFrame:
    return _read("recovery_schema.parquet")


@st.cache_data(show_spinner=False)
def detection_long() -> pd.DataFrame:
    return _read("detection_long.parquet")


@st.cache_data(show_spinner=False)
def smr_results() -> pd.DataFrame:
    return _read("smr_results.parquet")


@st.cache_data(show_spinner=False)
def locus_explanation() -> pd.DataFrame:
    return _read("locus_explanation.parquet")


@st.cache_data(show_spinner=False)
def druggability() -> pd.DataFrame:
    return _read("druggability.parquet")


@st.cache_data(show_spinner=False)
def gene_universe() -> pd.DataFrame:
    return _read("gene_universe.parquet")


@st.cache_data(show_spinner=False)
def decisions() -> list[dict]:
    """The decision log, newest entry per key."""
    path = LOGS / "decisions.jsonl"
    if not path.exists():
        return []
    latest: dict[str, dict] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            rec = json.loads(line)
            latest[rec["key"]] = rec
    return [latest[k] for k in sorted(latest)]


@st.cache_data(show_spinner=False)
def downloads() -> pd.DataFrame:
    """Provenance: which bytes each run read."""
    path = LOGS / "downloads.json"
    if not path.exists():
        return pd.DataFrame()
    log = json.loads(path.read_text(encoding="utf-8"))
    rows = [e for files in log.values() for e in files.values()]
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    return df[
        ["source_name", "file", "version", "licence", "redistributable",
         "bytes", "downloaded_utc", "sha256", "citation"]
    ].sort_values(["source_name", "file"])


@st.cache_data(show_spinner=False)
def doc(name: str) -> str:
    path = DOCS / name
    return path.read_text(encoding="utf-8") if path.exists() else ""


def ordered_rungs(df: pd.DataFrame, col: str = "rung") -> list[str]:
    present = set(df[col])
    return [r for r in RUNG_ORDER if r in present]
