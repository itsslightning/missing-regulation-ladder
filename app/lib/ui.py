"""Presentation helpers shared by the dashboard pages.

Keeping chart furniture and empty states here means every page gets the same
margins, legend placement and hover behaviour without repeating a dozen-line
``update_layout`` call, and a change to the house style happens once.

Colour carries meaning in this app and the meaning is not the same everywhere:
on the recovery pages the pair separates constrained genes from their matched
controls, while on the assay pages it separates bulk tissue from single-nucleus.
The static figures in ``docs/`` use the same two colours the same two ways, so
they are not remapped here. Where the risk of misreading is real, the assay
charts also vary marker SHAPE, which makes the encoding redundant rather than
relying on hue alone.
"""

from __future__ import annotations

from typing import Any

import plotly.graph_objects as go
import streamlit as st

#: Re-exported so pages import one module for presentation. The values live in
#: data.py and match pipeline/s07_figures.py, which is what keeps the dashboard
#: charts and the committed static figures the same colours.
from lib.data import GREY, RUST, TEAL  # noqa: F401

#: Redundant encoding for assay class, so hue is never the only cue.
ASSAY_SYMBOL = {"bulk tissue": "circle", "single-nucleus": "diamond"}


def is_dark() -> bool:
    """True when the viewer is in dark mode, for the few colours Plotly's
    Streamlit template does not derive on its own."""
    try:
        return st.context.theme.type == "dark"
    except Exception:
        return False


def muted() -> str:
    """A neutral that reads as 'absent' against the current background."""
    return "#2b3532" if is_dark() else "#eef1f0"


def layout(fig: go.Figure, height: int = 420, **kwargs: Any) -> go.Figure:
    """Apply the house chart style, then whatever the caller overrides."""
    fig.update_layout(
        height=height,
        # Room above for a horizontal legend and to the sides for point
        # labels, both of which get clipped by a flush margin.
        margin=dict(t=44, b=12, l=12, r=24),
        legend=dict(
            orientation="h", y=1.02, x=0,
            yanchor="bottom", xanchor="left",
            bgcolor="rgba(0,0,0,0)",
        ),
        hoverlabel=dict(font_size=13),
    )
    fig.update_layout(**kwargs)
    return fig


def plot(fig: go.Figure, height: int = 420, **kwargs: Any) -> None:
    st.plotly_chart(layout(fig, height, **kwargs), width="stretch")


def table(df, height: int | None = None, **kwargs: Any) -> None:
    # height is only forwarded when set: Streamlit rejects None, and its own
    # default ("content") is what an unbounded table should use.
    if height is not None:
        kwargs["height"] = height
    st.dataframe(df, hide_index=True, width="stretch", **kwargs)


def card(title: str, help: str | None = None):
    """A bordered card with a label. Used instead of a horizontal rule: the
    border does the grouping a divider used to do, and does it better."""
    box = st.container(border=True)
    box.markdown(f"**{title}**", help=help)
    return box


def missing(what: str) -> None:
    """The empty state when a pipeline output has not been generated yet.

    Calm and actionable rather than a red error box, and the command is in a
    code block so it can be copied instead of retyped.
    """
    with st.container(border=True):
        st.subheader(f":material/database_off: {what} is not available yet")
        st.markdown("This page reads a table the pipeline writes. Generate it with:")
        st.code("python scripts/run_all.py", language="bash")
        st.caption(
            "The dashboard computes nothing itself, so a missing table means "
            "the stage that produces it has not been run in this checkout."
        )
    st.stop()


def arm_control(label: str = "Detection rule") -> str:
    """The D-001 detection-rule selector.

    It lives in the sidebar rather than on one page because two pages depend on
    it. A reader looking at the SCHEMA page previously had no way to see which
    rule was in force, let alone change it.
    """
    from lib import data as D

    if "arm" not in st.session_state:
        st.session_state.arm = "uniform"
    st.radio(
        label,
        options=list(D.ARM_LABEL),
        format_func=lambda a: D.ARM_SHORT[a],
        key="arm",
        help=(
            "D-001 was delegated to the assistant rather than chosen by the "
            "project owner, so the rejected options ship alongside the primary "
            "one. The headline moves from 53% to 88% between them, which is "
            "why the study-native rule was not made primary."
        ),
    )
    return st.session_state.arm or "uniform"
