"""The Missing Regulation Ladder: dashboard entry point.

Reads the tables the pipeline wrote and presents them. It computes nothing, so
what is on screen and what the reproducible pipeline produced cannot drift
apart -- the same contract as the sibling scz-target-prioritization dashboard.

Two design requirements are load-bearing rather than decorative:

1. The DETECTION ARM is a visible toggle, not a hidden default. D-001 was
   delegated rather than chosen by the project owner, and the rejected options
   move the headline from 53% to 88%. A reader who cannot see that number move
   cannot judge the claim.

2. The ASSAY, not the rung index, is the organising variable. The project's
   result is that bulk tissue and single-nucleus separate while donor count
   does not matter, so the pages lead with that rather than with the ladder
   the project originally set out to climb.

Run:  streamlit run app/streamlit_app.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))

st.set_page_config(
    page_title="The Missing Regulation Ladder",
    page_icon=":material/stairs:",
    layout="wide",
)

if "arm" not in st.session_state:
    st.session_state.arm = "uniform"

page = st.navigation(
    [
        st.Page(
            "app_pages/overview.py",
            title="Recovery curve",
            icon=":material/stairs:",
            default=True,
        ),
        st.Page(
            "app_pages/assay.py",
            title="What closes the gap",
            icon=":material/biotech:",
        ),
        st.Page(
            "app_pages/schema_genes.py",
            title="SCHEMA genes",
            icon=":material/genetics:",
        ),
        st.Page(
            "app_pages/causal.py",
            title="Causal & druggable",
            icon=":material/medication:",
        ),
        st.Page(
            "app_pages/methodology.py",
            title="Methods",
            icon=":material/science:",
        ),
    ],
    position="top",
)

page.run()
