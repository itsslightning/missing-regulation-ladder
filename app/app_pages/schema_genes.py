"""Which SCHEMA genes switch on, and in which cell types."""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from lib import data as D
from lib import ui as U

st.header("SCHEMA genes across the ladder")
st.markdown(
    "The 32 genes with exome-wide schizophrenia evidence (Singh et al. 2022, "
    "FDR < 0.05). The constrained set is the statistical workhorse; these are "
    "the genes the question was really about, so it is worth asking which of "
    "them a cis-eQTL ever becomes visible for."
)

switch = D.schema_switch()
long = D.detection_long()
if switch.empty:
    U.missing("The SCHEMA switch table")

n_switch = int(switch["switches_on"].sum())
n_bulk = int(switch["already_visible"].sum())
n_never = int(switch["never_seen"].sum())

with st.container(horizontal=True):
    st.metric("Switch on at single-nucleus", n_switch,
              "invisible in bulk",
              delta_color="off", delta_arrow="off", border=True)
    st.metric("Already visible in bulk", n_bulk, border=True)
    st.metric("Never detected anywhere", n_never,
              "the residual",
              delta_color="off", delta_arrow="off", border=True)

never = switch[switch["never_seen"]]["symbol"].tolist()
if never:
    with st.container(border=True):
        st.markdown(
            f"**{len(never)} genes have no detectable cis-eQTL in any assay at "
            "any resolution tested**"
        )
        st.markdown(" ".join(f"`{g}`" for g in never))
        st.caption(
            "This is the residual missing regulation: the group the selection "
            "hypothesis speaks to most directly, and where CHRM4 also sits."
        )

#: the heatmap ------------------------------------------------------------
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
    colorscale=[[0, U.muted()], [1, U.TEAL]],
    showscale=False,
    xgap=3, ygap=3,
    customdata=np.array([[g] * len(SB_MAJOR) for g in group]),
    hovertemplate="<b>%{y}</b> in %{x}<br>%{customdata}<extra></extra>",
))

with st.container(border=True):
    st.markdown("**Detection by SingleBrain major cell type**")
    U.plot(
        fig,
        height=max(420, 22 * len(genes)),
        xaxis=dict(title=None, side="top"),
        yaxis=dict(title=None, autorange="reversed"),
        margin=dict(t=40, b=8, l=8, r=8),
    )
    st.caption(
        "Filled means a detectable cis-eQTL in that cell type under the "
        "uniform rule (D-001). Genes are ordered with those that switch on "
        "first, so the empty rows at the bottom are the residual."
    )

#: the table --------------------------------------------------------------
with U.card("All 32 genes"):
    tbl = switch.assign(
        Gene=switch["symbol"],
        Group=group,
        **{
            "Bulk cortex": switch["gtex_cortex"],
            "Bulk brain": switch["bulk_brain"],
            "sn major": switch["sn_major"],
            "sn subtype": switch["sn_subtype"],
            "Cell types": switch["cell_types"].replace("", "—"),
        },
    )
    # Booleans rather than tick glyphs: they sort, they filter, and the column
    # reads the same way in both themes.
    ticks = ["Bulk cortex", "Bulk brain", "sn major", "sn subtype"]
    U.table(
        tbl[["Gene", "Group"] + ticks + ["Cell types"]],
        column_config={
            c: st.column_config.CheckboxColumn(width="small") for c in ticks
        },
    )
    st.caption(
        "Sort by any rung column to group the genes, or use the search icon "
        "in the table toolbar to find one."
    )

#: recovery curve for the SCHEMA set --------------------------------------
rec = D.schema_recovery()
if not rec.empty:
    arm = st.session_state.get("arm", D.PRIMARY_ARM)
    s = rec[rec["arm"] == arm].set_index("rung").reindex(
        D.ordered_rungs(rec)).reset_index()
    fig2 = go.Figure()
    fig2.add_trace(go.Scatter(
        x=[D.RUNG_LABEL[r] for r in s["rung"]], y=s["rate"],
        mode="lines+markers", line=dict(color=U.GREY, width=3, dash="dash"),
        marker=dict(size=11, symbol="square"),
        error_y=dict(
            type="data", symmetric=False,
            array=s["hi"] - s["rate"], arrayminus=s["rate"] - s["lo"],
            color=U.GREY, thickness=1.5, width=6),
        name="SCHEMA genes",
        hovertemplate="<b>%{x}</b><br>%{y:.3f}<extra></extra>",
    ))
    with st.container(border=True):
        st.markdown("**Recovery for the SCHEMA set**")
        U.plot(
            fig2, height=340,
            yaxis=dict(title="Fraction with a detectable cis-eQTL",
                       range=[0, 1]),
            showlegend=False,
        )
        st.caption(
            f"Only {int(s['n_schema'].iloc[0])} genes, so the intervals are "
            "wide. The LOEUF-constrained set of 2,937 genes carries the "
            "statistical weight; this is the sharper, smaller overlay."
        )
