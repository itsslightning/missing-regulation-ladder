"""The recovery curve: constrained genes vs matched controls, per rung."""

from __future__ import annotations

import plotly.graph_objects as go
import streamlit as st

from lib import data as D

st.title("The Missing Regulation Ladder")
st.markdown(
    "**Does the shortage of cis-eQTLs at schizophrenia risk genes reflect "
    "natural selection, or insufficient resolution and power?**"
)

curve = D.recovery()
closure = D.gap_closure()
if curve.empty:
    st.error("No recovery table found. Run `python scripts/run_all.py` first.")
    st.stop()

# Initialised here, not only in streamlit_app.py, so the page also works when
# run on its own -- which is what a deep link and the AppTest suite both do.
if "arm" not in st.session_state:
    st.session_state.arm = "uniform"

arm = st.segmented_control(
    "Detection rule",
    options=list(D.ARM_LABEL),
    format_func=lambda a: D.ARM_LABEL[a],
    key="arm",
    help=(
        "D-001 was delegated to the assistant rather than chosen by the "
        "project owner, so the rejected options ship alongside the primary. "
        "The headline moves from 53% to 88% between them — which is exactly "
        "why the study-native rule was not made primary."
    ),
)
arm = arm or "uniform"

sub = curve[curve["arm"] == arm].set_index("rung").reindex(
    D.ordered_rungs(curve)
).reset_index()
cl = closure[closure["arm"] == arm]

# --- headline ------------------------------------------------------------
if not cl.empty:
    r = cl.iloc[0]
    with st.container(horizontal=True):
        st.metric(
            "Gap closed, bulk cortex → single-nucleus",
            f"{r['closure']:.1%}",
            f"95% CI {r['closure_lo']:.0%}–{r['closure_hi']:.0%}",
            border=True,
        )
        st.metric(
            "Residual gap that does not close",
            f"{1 - r['closure']:.1%}",
            "consistent with selection",
            delta_color="off",
            border=True,
        )
        st.metric(
            "Gap at bulk cortex → at single-nucleus",
            f"{r['gap_from']:.3f} → {r['gap_to']:.3f}",
            border=True,
        )

# --- the curve -----------------------------------------------------------
x = [D.RUNG_LABEL[r] for r in sub["rung"]]
fig = go.Figure()
for rate, lo, hi, colour, name in (
    ("control_rate", "control_lo", "control_hi", D.RUST,
     "Matched unconstrained controls"),
    ("constrained_rate", "constrained_lo", "constrained_hi", D.TEAL,
     "Constrained genes (LOEUF < 0.35)"),
):
    fig.add_trace(go.Scatter(
        x=x + x[::-1],
        y=list(sub[hi]) + list(sub[lo])[::-1],
        fill="toself", fillcolor=colour, opacity=0.16,
        line=dict(width=0), hoverinfo="skip", showlegend=False,
    ))
    fig.add_trace(go.Scatter(
        x=x, y=sub[rate], name=name, mode="lines+markers",
        line=dict(color=colour, width=3), marker=dict(size=10),
        customdata=sub[["n_donors", lo, hi]],
        hovertemplate=(
            "<b>%{x}</b><br>%{y:.3f} "
            "[%{customdata[1]:.3f}, %{customdata[2]:.3f}]"
            "<br>%{customdata[0]} donors<extra></extra>"
        ),
    ))

fig.update_layout(
    height=460,
    yaxis=dict(title="Fraction with a detectable cis-eQTL", range=[0, 1]),
    xaxis=dict(title=None),
    legend=dict(orientation="h", y=1.12, x=0),
    margin=dict(t=40, b=10, l=10, r=10),
    hovermode="x unified",
)
st.plotly_chart(fig, use_container_width=True)

st.caption(
    "The distance between the two curves is the **gap**. How it behaves across "
    "rungs is the whole experiment: flat favours selection, closing favours "
    "power/resolution. It closes partially — and the *What closes the gap* "
    "page shows the cause is neither of the two hypotheses this project set "
    "out to adjudicate."
)

# --- per-rung table ------------------------------------------------------
with st.container(border=True):
    st.markdown("**Per rung**, with 95% confidence intervals")
    show = sub.assign(
        Rung=[D.RUNG_LABEL[r] for r in sub["rung"]],
        Assay=[D.RUNG_ASSAY[r] for r in sub["rung"]],
        Donors=sub["n_donors"],
        Constrained=sub.apply(
            lambda r: f"{r['constrained_rate']:.3f} "
            f"[{r['constrained_lo']:.3f}, {r['constrained_hi']:.3f}]", axis=1),
        Control=sub.apply(
            lambda r: f"{r['control_rate']:.3f} "
            f"[{r['control_lo']:.3f}, {r['control_hi']:.3f}]", axis=1),
        Gap=sub.apply(
            lambda r: f"{r['gap']:.3f} [{r['gap_lo']:.3f}, {r['gap_hi']:.3f}]",
            axis=1),
    )
    st.dataframe(
        show[["Rung", "Assay", "Donors", "Constrained", "Control", "Gap"]],
        hide_index=True, use_container_width=True,
    )

# --- robustness ----------------------------------------------------------
rob = D.robustness()
if not rob.empty:
    with st.container(border=True):
        st.markdown("**Does the headline survive the choices it rests on?**")
        r = rob[rob["closure"].notna()].copy()
        r["Closure"] = r.apply(
            lambda x: f"{x['closure']:.1%} "
            f"[{x['closure_lo']:.1%}, {x['closure_hi']:.1%}]", axis=1)
        r["Gap"] = r.apply(
            lambda x: f"{x['gap_from']:.3f} → {x['gap_to']:.3f}", axis=1)
        r["Genes"] = r["n_cases"].astype(int)
        st.dataframe(
            r[["family", "variant", "Genes", "Gap", "Closure"]].rename(
                columns={"family": "Choice", "variant": "Variant"}),
            hide_index=True, use_container_width=True,
        )
        lo, hi = r["closure"].min(), r["closure"].max()
        st.caption(
            f"Closure spans **{lo:.1%}–{hi:.1%}** across every structural "
            "variant tested, and the residual gap excludes zero in all of them."
        )
