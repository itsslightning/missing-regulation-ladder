"""The recovery curve: constrained genes vs matched controls, per rung."""

from __future__ import annotations

import plotly.graph_objects as go
import streamlit as st

from lib import data as D
from lib import ui as U

st.title("The missing regulation ladder")
st.markdown(
    "Genes under strong selective constraint carry fewer detectable cis-eQTLs "
    "than comparable genes. Two published accounts explain why: selection has "
    "removed the regulatory variation, or the assays were never resolved "
    "enough to see it. This climbs four rungs of brain expression data to find "
    "out which."
)

curve = D.recovery()
closure = D.gap_closure()
if curve.empty:
    U.missing("The recovery table")

arm = st.session_state.get("arm", D.PRIMARY_ARM)
sub = curve[curve["arm"] == arm].set_index("rung").reindex(
    D.ordered_rungs(curve)
).reset_index()
cl = closure[closure["arm"] == arm]

#: headline ---------------------------------------------------------------
if not cl.empty:
    r = cl.iloc[0]
    with st.container(border=True):
        st.subheader(":material/insights: The gap closes, but only part way")
        st.markdown(
            f"Going from bulk cortex to single-nucleus data closes "
            f"**{r['closure']:.1%}** of the excess, and the remaining "
            f"**{1 - r['closure']:.1%}** does not close at any resolution "
            "tested. Neither published account predicts that, and the "
            "*What closes the gap* page shows the cause is neither of them."
        )
        with st.container(horizontal=True):
            st.metric(
                "Gap closed, bulk cortex to single-nucleus",
                f"{r['closure']:.1%}",
                f"95% CI {r['closure_lo']:.0%}\u2013{r['closure_hi']:.0%}",
                delta_color="off", delta_arrow="off", border=True,
            )
            st.metric(
                "Residual gap that never closes",
                f"{1 - r['closure']:.1%}",
                "consistent with selection",
                delta_color="off", delta_arrow="off", border=True,
            )
            st.metric(
                "Gap, first rung to last",
                f"{r['gap_from']:.3f} \u2192 {r['gap_to']:.3f}",
                delta_color="off", delta_arrow="off", border=True,
            )

#: the curve --------------------------------------------------------------
x = [D.RUNG_LABEL[r] for r in sub["rung"]]
fig = go.Figure()
for rate, lo, hi, colour, name in (
    ("control_rate", "control_lo", "control_hi", U.RUST,
     "Matched controls"),
    ("constrained_rate", "constrained_lo", "constrained_hi", U.TEAL,
     "Constrained genes"),
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

box = st.container(border=True)
with box:
    st.markdown("**Detection rate across the ladder**")
    U.plot(
        fig, height=440,
        yaxis=dict(title="Fraction with a detectable cis-eQTL", range=[0, 1]),
        xaxis=dict(title=None),
        hovermode="x unified",
    )
    st.caption(
        "Teal is the constrained set (LOEUF < 0.35), rust its matched "
        "unconstrained controls. The distance between the two curves is the "
        "gap, and how it behaves across rungs is the whole experiment. Flat "
        "favours selection, closing favours power and resolution. Shaded "
        "bands are 95% confidence intervals."
    )

#: per-rung table ---------------------------------------------------------
with U.card("Every rung, with 95% confidence intervals"):
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
    U.table(
        show[["Rung", "Assay", "Donors", "Constrained", "Control", "Gap"]],
        column_config={
            "Donors": st.column_config.NumberColumn(format="%d"),
            "Assay": st.column_config.TextColumn(width="small"),
        },
    )

#: robustness -------------------------------------------------------------
rob = D.robustness()
if not rob.empty:
    r = rob[rob["closure"].notna()].copy()
    lo, hi = r["closure"].min(), r["closure"].max()
    with st.expander(
        f"Does the headline survive the choices it rests on? "
        f"Closure spans {lo:.1%}\u2013{hi:.1%}",
        icon=":material/rule_settings:",
    ):
        r["Closure"] = r.apply(
            lambda x: f"{x['closure']:.1%} "
            f"[{x['closure_lo']:.1%}, {x['closure_hi']:.1%}]", axis=1)
        r["Gap"] = r.apply(
            lambda x: f"{x['gap_from']:.3f} \u2192 {x['gap_to']:.3f}", axis=1)
        r["Genes"] = r["n_cases"].astype(int)
        U.table(
            r[["family", "variant", "Genes", "Gap", "Closure"]].rename(
                columns={"family": "Choice", "variant": "Variant"}),
        )
        st.caption(
            "Every structural variant tested, including the two the detection "
            "rule was chosen over. The residual gap excludes zero in all of "
            "them, which is the part of the result that does not depend on "
            "the choices."
        )
