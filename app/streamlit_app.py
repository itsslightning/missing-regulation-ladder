"""The Missing Regulation Ladder: dashboard entry point.

Reads the tables the pipeline wrote and presents them. It computes nothing, so
what is on screen and what the reproducible pipeline produced cannot drift
apart, the same contract as the sibling scz-target-prioritization dashboard.

Three design requirements are load-bearing rather than decorative:

1. The DETECTION ARM is a visible control, not a hidden default. D-001 was
   delegated rather than chosen by the project owner, and the rejected options
   move the headline from 53% to 88%. A reader who cannot see that number move
   cannot judge the claim. It sits in the sidebar because two pages depend on
   it, and a control that governs more than the page you are on has no business
   being buried in one of them.

2. The ASSAY, not the rung index, is the organising variable. The project's
   result is that bulk tissue and single-nucleus separate while donor count
   does not matter, so the pages lead with that rather than with the ladder the
   project originally set out to climb.

3. Every page opens with what it found, then shows the evidence, then the
   caveats. Reading order follows importance rather than the order the analysis
   happened to run in.

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
    initial_sidebar_state="expanded",
)

from lib import data as D  # noqa: E402  (needs the sys.path line above)
from lib import ui as U  # noqa: E402

if "arm" not in st.session_state:
    st.session_state.arm = D.PRIMARY_ARM

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
            title="Causal and druggable",
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

with st.sidebar:
    st.markdown("### Missing regulation ladder")
    st.caption(
        "Do schizophrenia risk genes lack cis-eQTLs because selection removed "
        "them, or because we were not looking closely enough?"
    )
    st.space("small")

    arm = U.arm_control()
    if arm == D.PRIMARY_ARM:
        st.caption(":material/check_circle: The rule D-001 selected.")
    else:
        st.badge("Sensitivity arm", icon=":material/science:", color="orange")
        st.caption(
            "A rejected option, shipped so the choice can be inspected. "
            "Numbers on the recovery and SCHEMA pages change; the assay, "
            "causal and methods pages do not depend on this."
        )

    st.space("medium")
    st.caption(
        "[Source and data](https://github.com/itsslightning/missing-regulation-ladder)"
    )

page.run()
