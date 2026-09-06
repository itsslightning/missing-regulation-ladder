"""Which SCHEMA genes switch on, and in which cell types."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from lib import data as D

st.header("SCHEMA genes across the ladder")
st.markdown(
    "The 32 genes with exome-wide schizophrenia evidence "
    "(Singh et al. 2022, FDR < 0.05), and where a cis-eQTL becomes detectable."
)

switch = D.schema_switch()
long = D.detection_long()
if switch.empty:
    st.error("No SCHEMA switch table. Run `python scripts/run_all.py` first.")
    st.stop()

n_switch = int(switch["switches_on"].sum())
n_bulk = int(switch["already_visible"].sum())
n_never = int(switch["never_seen"].sum())

with st.container(horizontal=True):
    st.metric("Switch on at single-nucleus", n_switch,
              "invisible in bulk", delta_color="off", border=True)
    st.metric("Already visible in bulk", n_bulk, border=True)
    st.metric("Never detected anywhere", n_never,
              "the residual", delta_color="off", border=True)

# --- the heatmap ---------------------------------------------------------
st.subheader("Per-cell-type detection")

SB_MAJOR = ["Ast", "End", "Ext", "IN", "MG", "OD", "OPC"]
genes = switch.index.tolist()
symbols = switch["symbol"].tolist()

sn = long[(long["rung"] == "sn_major") & long["detected"]]
hit = {(g, c) for g, c in zip(sn["gene_id"], sn["cell"])}

z = np.array([[1 if (g, c) in hit else 0 for c in SB_MAJOR] for g in genes])
group = np.where(
    switch["switches_on"], "switches on",
    np.where(switch["already_visible"], "visible in bulk", "never detected"),
)

fig = go.Figure(go.Heatmap(
    z=z,
    x=SB_MAJOR,
    y=symbols,
    colorscale=[[0, "#f2f2f2"], [1, D.TEAL]],
    showscale=False,
    xgap=2, ygap=2,
    customdata=np.array([[g] * len(SB_MAJOR) for g in group]),
    hovertemplate=(
        "<b>%{y}</b> in %{x}<br>%{customdata}<extra></extra>"
    ),
))
fig.update_layout(
    height=max(420, 22 * len(genes)),
    xaxis=dict(title="SingleBrain major cell type", side="top"),
    yaxis=dict(title=None, autorange="reversed"),
    margin=dict(t=60, b=10, l=10, r=10),
)
st.plotly_chart(fig, use_container_width=True)
st.caption(
    "Filled = a detectable cis-eQTL in that cell type under the uniform rule "
    "(D-001). Genes are ordered with those that switch on first."
)

# --- the table -----------------------------------------------------------
tick = {True: "●", False: "—"}
tbl = switch.assign(
    Gene=switch["symbol"],
    Group=group,
    **{
        "Bulk cortex": switch["gtex_cortex"].map(tick),
        "Bulk brain": switch["bulk_brain"].map(tick),
        "sn major": switch["sn_major"].map(tick),
        "sn subtype": switch["sn_subtype"].map(tick),
        "Cell types": switch["cell_types"].replace("", "—"),
    },
)
with st.container(border=True):
    st.markdown("**All 32 genes**")
    st.dataframe(
        tbl[["Gene", "Group", "Bulk cortex", "Bulk brain", "sn major",
             "sn subtype", "Cell types"]],
        hide_index=True, use_container_width=True,
    )

never = switch[switch["never_seen"]]["symbol"].tolist()
if never:
    st.warning(
        "**No detectable cis-eQTL in any assay at any resolution:** "
        + ", ".join(never)
        + ". This is the residual missing regulation — the group most "
        "relevant to the selection hypothesis, and where CHRM4 also sits.",
        icon=":material/priority_high:",
    )

# --- recovery curve for the SCHEMA set ------------------------------------
rec = D.schema_recovery()
if not rec.empty:
    arm = st.session_state.get("arm", "uniform")
    s = rec[rec["arm"] == arm].set_index("rung").reindex(
        D.ordered_rungs(rec)).reset_index()
    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(
        x=[D.RUNG_LABEL[r] for r in s["rung"]], y=s["rate"],
        mode="lines+markers", line=dict(color=D.GREY, width=3, dash="dash"),
        marker=dict(size=10, symbol="square"),
        error_y=dict(
            type="data", symmetric=False,
            array=s["hi"] - s["rate"], arrayminus=s["rate"] - s["lo"],
            color=D.GREY, thickness=1.5, width=6),
        name="SCHEMA genes",
        hovertemplate="<b>%{x}</b><br>%{y:.3f}<extra></extra>",
    ))
    fig2.update_layout(
        height=340,
        yaxis=dict(title="Fraction with a detectable cis-eQTL", range=[0, 1]),
        margin=dict(t=30, b=10, l=10, r=10), showlegend=False,
    )
    st.subheader("Recovery for the SCHEMA set")
    st.plotly_chart(fig2, use_container_width=True)
    st.caption(
        f"Only {int(s['n_schema'].iloc[0])} genes, so the intervals are wide. "
        "The LOEUF-constrained set (2,937 genes) is the statistical workhorse; "
        "SCHEMA is the sharper, smaller overlay."
    )
