"""What closes the gap: assay, donor count, or cell-type granularity?"""

from __future__ import annotations

import plotly.graph_objects as go
import streamlit as st

from lib import data as D

st.header("What closes the gap")
st.markdown(
    "Three explanations were on the table. All three are directly testable "
    "here, and **two are excluded**."
)

assay = D.assay_contrast()
res = D.resolution_test()

if assay.empty:
    st.error("No assay contrast table. Run `python scripts/run_all.py` first.")
    st.stop()

# --- verdict table -------------------------------------------------------
with st.container(border=True):
    st.markdown("**The three candidate explanations**")
    st.markdown(
        """
| Explanation | Test | Verdict |
|---|---|---|
| **Donor count** — the power account | ×6.8 donors within bulk tissue | :red[**excluded**] — gap unchanged |
| **Cell-type resolution** | splitting *and* pooling, two studies | :red[**excluded**] — bounded ≤10% |
| **Assay** — nuclei vs whole tissue | bulk vs snRNA-seq, matched genes | :green[**survives**] |
"""
    )

# --- the assay figure ----------------------------------------------------
fig = go.Figure()
for label, colour, name in (
    ("bulk", D.RUST, "Bulk tissue RNA-seq"),
    ("single-nucleus", D.TEAL, "Single-nucleus RNA-seq"),
):
    d = assay[assay["assay"] == label].sort_values("n_donors")
    # Bryois contributes two arms at the same donor count, so a single
    # textposition would stack their labels on top of each other.
    positions = {
        "bryois_pb": "bottom left",
        "bryois_celltype": "top left",
        "sn_major": "bottom center",
    }
    fig.add_trace(go.Scatter(
        x=d["n_donors"], y=d["gap"], name=name, mode="lines+markers+text",
        text=d["label"],
        textposition=[positions.get(r, "top center") for r in d["rung"]],
        textfont=dict(size=10, color=D.GREY),
        line=dict(color=colour, width=3),
        marker=dict(size=12),
        error_y=dict(
            type="data", symmetric=False,
            array=d["gap_hi"] - d["gap"], arrayminus=d["gap"] - d["gap_lo"],
            color=colour, thickness=1.5, width=6,
        ),
        hovertemplate=(
            "<b>%{text}</b><br>gap %{y:.3f}<br>%{x} donors<extra></extra>"
        ),
    ))

fig.update_layout(
    height=470,
    xaxis=dict(title="Donors (log scale)", type="log"),
    yaxis=dict(title="Constrained-gene gap", range=[0, 0.47]),
    legend=dict(orientation="h", y=1.12, x=0),
    margin=dict(t=40, b=10, l=10, r=10),
)
st.plotly_chart(fig, use_container_width=True)

bulk = assay[assay["assay"] == "bulk"]
sn = assay[assay["assay"] == "single-nucleus"]
n_cases = int(assay["n_cases"].iloc[0])
n_ctrl = int(assay["n_controls"].iloc[0])

with st.container(horizontal=True):
    st.metric(
        "Bulk spread, across 6.8× donors",
        f"{bulk['gap'].max() - bulk['gap'].min():.3f}",
        "two studies, same answer", delta_color="off", border=True,
    )
    st.metric(
        "Single-nucleus spread, across 5.1× donors",
        f"{sn['gap'].max() - sn['gap'].min():.3f}",
        "three arms, same answer", delta_color="off", border=True,
    )
    st.metric(
        "Separation between assays",
        f"{bulk['gap'].min() - sn['gap'].max():.3f}",
        "six times the larger within-class spread",
        delta_color="off", border=True,
    )

st.caption(
    f"All arms scored on the same {n_cases:,} constrained and {n_ctrl:,} "
    "matched control genes, so no arm is advantaged by which genes it covers. "
    "The sharpest single comparison is **Bryois pseudobulk** — single-nucleus "
    "data with all nuclei pooled, the *smallest* study here at 192 donors — "
    "showing a gap of 0.144 against **PsychENCODE** bulk tissue at 1,387 "
    "donors showing 0.367. Seven times the donors, more than twice the gap."
)

# --- resolution test -----------------------------------------------------
if not res.empty:
    st.subheader("Why resolution is excluded")
    pooled = res[res["major"] == "POOLED"]
    bry = res[res["major"].str.startswith("Bryois") & (res["arm"] == "split")]

    cols = st.columns(2, border=True)
    with cols[0]:
        st.markdown("**Splitting** — SingleBrain classes into their subtypes")
        if not pooled.empty:
            p = pooled.iloc[0]
            st.metric(
                "Change in gap (pooled across 6 classes)",
                f"{p['delta']:+.3f}",
                f"95% CI [{p['delta_lo']:+.3f}, {p['delta_hi']:+.3f}]",
                delta_color="off",
            )
        st.caption(
            "Same donors, same nuclei, same pipeline — only the grouping "
            "changes. Every per-class interval spans zero."
        )
    with cols[1]:
        st.markdown("**Pooling** — Bryois pseudobulk vs its own 8 cell types")
        if not bry.empty:
            b = bry.iloc[0]
            st.metric(
                "Change in gap",
                f"{b['delta']:+.3f}",
                f"95% CI [{b['delta_lo']:+.3f}, {b['delta_hi']:+.3f}]",
                delta_color="off",
            )
        st.caption(
            "The same donors again, varying resolution in the **opposite** "
            "direction. A matched-gene-set control confirms the two studies "
            "do not disagree."
        )

    st.info(
        "Resolution and per-context power are **intrinsically coupled** in "
        "single-cell data: at fixed sequencing depth you cannot resolve more "
        "contexts without fewer reads in each. Neither test isolates "
        "resolution alone — but they vary it in opposite directions, and both "
        "land on zero. That brackets the answer.",
        icon=":material/info:",
    )

# --- per-arm table -------------------------------------------------------
with st.container(border=True):
    st.markdown("**Every arm, one gene set**")
    show = assay.assign(
        Arm=assay["label"],
        Assay=assay["assay"],
        Donors=assay["n_donors"],
        Constrained=assay["constrained_rate"].map("{:.3f}".format),
        Control=assay["control_rate"].map("{:.3f}".format),
        Gap=assay.apply(
            lambda r: f"{r['gap']:.3f} [{r['gap_lo']:.3f}, {r['gap_hi']:.3f}]",
            axis=1),
    )
    st.dataframe(
        show[["Arm", "Assay", "Donors", "Constrained", "Control", "Gap"]],
        hide_index=True, use_container_width=True,
    )
