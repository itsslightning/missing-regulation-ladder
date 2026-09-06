"""What closes the gap: assay, donor count, or cell-type granularity?"""

from __future__ import annotations

import plotly.graph_objects as go
import streamlit as st

from lib import data as D
from lib import ui as U

st.header("What closes the gap")
st.markdown(
    "Three things change as you climb the ladder: the number of donors, how "
    "finely cells are split, and whether the assay reads whole tissue or "
    "isolated nuclei. Each one is separable here, and two of them are ruled "
    "out."
)

assay = D.assay_contrast()
res = D.resolution_test()

if assay.empty:
    U.missing("The assay contrast table")

#: verdicts ---------------------------------------------------------------
verdicts = [
    (
        "Donor count",
        "excluded", "red",
        "The power account. Tested across a 6.8-fold range of donors within "
        "bulk tissue, where nothing else changes.",
        "The gap does not move.",
    ),
    (
        "Cell-type resolution",
        "excluded", "red",
        "Tested twice, by splitting classes into subtypes and by pooling "
        "subtypes back together, in two independent studies.",
        "Bounded at 10% of the closure, in both directions.",
    ),
    (
        "Assay",
        "survives", "green",
        "Whole tissue against isolated nuclei, scored on one shared gene set "
        "so no arm is helped by which genes it happens to cover.",
        "The only variable that separates the arms.",
    ),
]
cols = st.columns(3, border=True)
for col, (name, verdict, colour, how, result) in zip(cols, verdicts):
    with col:
        st.markdown(f"**{name}**")
        st.badge(
            verdict,
            icon=":material/close:" if colour == "red" else ":material/check:",
            color=colour,
        )
        st.caption(how)
        st.markdown(result)

#: the assay figure -------------------------------------------------------
fig = go.Figure()
for label, colour, name in (
    ("bulk", U.RUST, "Bulk tissue"),
    ("single-nucleus", U.TEAL, "Single-nucleus"),
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
        textfont=dict(size=10, color=U.GREY),
        line=dict(color=colour, width=3),
        # Shape as well as hue, so the two classes stay distinguishable
        # without relying on colour vision.
        marker=dict(
            size=13,
            symbol="circle" if label == "bulk" else "diamond",
        ),
        error_y=dict(
            type="data", symmetric=False,
            array=d["gap_hi"] - d["gap"], arrayminus=d["gap"] - d["gap_lo"],
            color=colour, thickness=1.5, width=6,
        ),
        hovertemplate=(
            "<b>%{text}</b><br>gap %{y:.3f}<br>%{x} donors<extra></extra>"
        ),
    ))

bulk = assay[assay["assay"] == "bulk"]
sn = assay[assay["assay"] == "single-nucleus"]
n_cases = int(assay["n_cases"].iloc[0])
n_ctrl = int(assay["n_controls"].iloc[0])

with st.container(border=True):
    st.markdown("**Constrained-gene gap against donor count, every arm**")
    U.plot(
        fig, height=470,
        # Padded well past the 192-1,387 donor range: the arm labels sit
        # beside their markers and would otherwise run off the plot. Ticks are
        # named explicitly because a log axis that spans less than a decade
        # falls back to labelling its minor ticks, which reads as 2,3,4...
        xaxis=dict(
            title="Donors (log scale)", type="log", range=[2.12, 3.32],
            tickmode="array",
            tickvals=[200, 300, 500, 700, 1000, 1400],
            ticktext=["200", "300", "500", "700", "1,000", "1,400"],
        ),
        yaxis=dict(title="Constrained-gene gap", range=[0, 0.47]),
    )
    st.caption(
        f"All arms scored on the same {n_cases:,} constrained and {n_ctrl:,} "
        "matched control genes. Circles are bulk tissue, diamonds are "
        "single-nucleus. If donors drove the gap the two series would trend "
        "downward together; instead each series is flat and they sit apart."
    )

with st.container(horizontal=True):
    st.metric(
        "Bulk spread, across 6.8× donors",
        f"{bulk['gap'].max() - bulk['gap'].min():.3f}",
        "two studies, same answer",
        delta_color="off", delta_arrow="off", border=True,
    )
    st.metric(
        "Single-nucleus spread, across 5.1× donors",
        f"{sn['gap'].max() - sn['gap'].min():.3f}",
        "three arms, same answer",
        delta_color="off", delta_arrow="off", border=True,
    )
    st.metric(
        "Separation between the two assays",
        f"{bulk['gap'].min() - sn['gap'].max():.3f}",
        "six times the larger within-class spread",
        delta_color="off", delta_arrow="off", border=True,
    )

st.caption(
    "The sharpest single comparison is **Bryois pseudobulk**: single-nucleus "
    "data with all nuclei pooled, the smallest study here at 192 donors, with "
    "a gap of 0.144, against **PsychENCODE** bulk tissue at 1,387 donors "
    "showing 0.367. Seven times the donors, more than twice the gap."
)

#: resolution test --------------------------------------------------------
if not res.empty:
    st.subheader("Why resolution is excluded")
    st.markdown(
        "Resolution was varied in both directions, in two studies that share "
        "no donors, and neither move shifts the gap."
    )
    pooled = res[res["major"] == "POOLED"]
    bry = res[res["major"].str.startswith("Bryois") & (res["arm"] == "split")]

    cols = st.columns(2, border=True)
    with cols[0]:
        st.markdown("**Splitting**")
        st.caption("SingleBrain classes into their subtypes")
        if not pooled.empty:
            p = pooled.iloc[0]
            st.metric(
                "Change in gap, pooled across 6 classes",
                f"{p['delta']:+.3f}",
                f"95% CI [{p['delta_lo']:+.3f}, {p['delta_hi']:+.3f}]",
                delta_color="off", delta_arrow="off",
            )
        st.caption(
            "Same donors, same nuclei, same pipeline. Only the grouping "
            "changes, and every per-class interval spans zero."
        )
    with cols[1]:
        st.markdown("**Pooling**")
        st.caption("Bryois pseudobulk against its own 8 cell types")
        if not bry.empty:
            b = bry.iloc[0]
            st.metric(
                "Change in gap",
                f"{b['delta']:+.3f}",
                f"95% CI [{b['delta_lo']:+.3f}, {b['delta_hi']:+.3f}]",
                delta_color="off", delta_arrow="off",
            )
        st.caption(
            "The same donors again, varying resolution the opposite way. A "
            "matched-gene-set control confirms the two studies do not disagree."
        )

    with st.expander(
        "What these two tests can and cannot rule out",
        icon=":material/help:",
    ):
        st.markdown(
            "Resolution and per-context power are coupled in single-cell data: "
            "at fixed sequencing depth you cannot resolve more contexts "
            "without putting fewer reads in each. So neither test isolates "
            "resolution on its own. What they do is vary it in opposite "
            "directions, which brackets the answer: if resolution were doing "
            "the work, splitting and pooling could not both land on zero."
        )

#: per-arm table ----------------------------------------------------------
with U.card("Every arm, one gene set"):
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
    U.table(
        show[["Arm", "Assay", "Donors", "Constrained", "Control", "Gap"]],
        column_config={"Donors": st.column_config.NumberColumn(format="%d")},
    )
