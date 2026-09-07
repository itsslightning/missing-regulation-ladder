"""Why the matched control set exists: the expression confound.

This page answers the first objection anyone raises, which is that constrained
genes might simply be lowly expressed. It is the strongest single argument for
D-002, and it is worth showing rather than asserting, because the unadjusted
view does not merely weaken the result: it reverses it.

The decile summary in the middle panel is a median over a committed column of
gene_universe.parquet. That is presentation-level aggregation of a table the
pipeline already wrote, not a re-derivation of any result, so the app still
computes no analysis of its own.
"""

from __future__ import annotations

import plotly.graph_objects as go
import streamlit as st

from lib import data as D
from lib import ui as U

st.header("Why the controls matter")
st.markdown(
    "Constrained genes are also, on average, the most highly expressed genes "
    "in the genome. Any comparison that ignores that is measuring expression "
    "as much as constraint. This page shows what happens when you ignore it, "
    "which is the reason every headline on this dashboard is scored against "
    "matched controls (D-002) rather than against the genome."
)

raw = D.by_loeuf()
adj = D.by_loeuf_expr()
if raw.empty or adj.empty:
    U.missing("The LOEUF decile tables")

arm = st.session_state.get("arm", D.PRIMARY_ARM)
RUNGS = ["gtex_cortex", "bulk_brain", "sn_major", "sn_subtype"]
STYLE = {
    "gtex_cortex": (U.RUST, "solid"),
    "bulk_brain": (U.RUST, "dash"),
    "sn_major": (U.TEAL, "solid"),
    "sn_subtype": (U.TEAL, "dash"),
}

#: The two decile charts sit in half-width columns, where the full rung labels
#: run out of legend and get clipped mid-word.
SHORT = {
    "gtex_cortex": "GTEx",
    "bulk_brain": "PsychENCODE",
    "sn_major": "SB 7 types",
    "sn_subtype": "SB 28 subtypes",
}


def _curves(df, title_hint: str) -> go.Figure:
    fig = go.Figure()
    for rung in RUNGS:
        d = df[df["rung"] == rung].sort_values("loeuf_decile")
        if d.empty:
            continue
        colour, dash = STYLE[rung]
        fig.add_trace(go.Scatter(
            x=d["loeuf_decile"], y=d["rate"],
            name=SHORT[rung], mode="lines+markers",
            line=dict(color=colour, width=2.5, dash=dash),
            marker=dict(size=7),
            hovertemplate=(
                f"<b>{D.RUNG_LABEL[rung]}</b><br>decile %{{x}}"
                f"<br>{title_hint} %{{y:.3f}}<extra></extra>"
            ),
        ))
    return fig


#: the inversion, stated as a number ---------------------------------------
def _spread(df) -> dict[str, float]:
    out = {}
    for rung in RUNGS:
        d = df[df["rung"] == rung].set_index("loeuf_decile")["rate"]
        if {0, 9}.issubset(d.index):
            out[rung] = float(d.loc[0] - d.loc[9])
    return out

raw_u = raw[(raw["arm"] == arm) & (raw["metric"] == "loeuf")]
adj_u = adj[(adj["arm"] == arm) & (adj["expr_tertile"] == "high")]
raw_gap = _spread(raw_u)
adj_gap = _spread(adj_u)

if raw_gap and adj_gap:
    with st.container(border=True):
        st.subheader(":material/swap_vert: Unadjusted, the result comes out backwards")
        st.markdown(
            "Read the most-constrained decile against the least-constrained "
            "one. A **negative** number means constrained genes are recovered "
            "*worse*, which is the effect this project is about. A "
            "**positive** number means they look *better* recovered, which is "
            "the opposite."
        )
        with st.container(horizontal=True):
            st.metric(
                "SingleBrain, 7 types — unadjusted",
                f"{raw_gap['sn_major']:+.3f}",
                "constrained look BEST recovered",
                delta_color="off", delta_arrow="off", border=True,
            )
            st.metric(
                "SingleBrain, 7 types — expression held fixed",
                f"{adj_gap['sn_major']:+.3f}",
                "the inversion disappears",
                delta_color="off", delta_arrow="off", border=True,
            )
            st.metric(
                "GTEx cortex — expression held fixed",
                f"{adj_gap['gtex_cortex']:+.3f}",
                "expected direction, seven times clearer",
                delta_color="off", delta_arrow="off", border=True,
            )

#: raw vs adjusted, side by side -------------------------------------------
cols = st.columns(2, border=True)
with cols[0]:
    st.markdown("**Unadjusted**")
    st.caption("Every gene, binned by LOEUF decile. Decile 0 is most constrained.")
    U.plot(
        _curves(raw_u, "recovery"), height=380,
        xaxis=dict(title="LOEUF decile (0 = most constrained)",
                   dtick=1, tickmode="linear"),
        yaxis=dict(title="Fraction with a detectable cis-eQTL", range=[0, 0.8]),
    )
    st.caption(
        "At the single-nucleus rungs the curve runs downhill from decile 0, so "
        "the most constrained genes appear to be the best served by the data. "
        "Taken at face value this says there is no missing regulation at all."
    )
with cols[1]:
    st.markdown("**Expression held fixed**")
    st.caption("The same genes, restricted to the high-expression tertile.")
    U.plot(
        _curves(adj_u, "recovery"), height=380,
        xaxis=dict(title="LOEUF decile (0 = most constrained)",
                   dtick=1, tickmode="linear"),
        yaxis=dict(title="Fraction with a detectable cis-eQTL", range=[0, 0.8]),
    )
    st.caption(
        "Comparing like with like, the bulk rungs climb steeply: constrained "
        "genes really are worse served. The single-nucleus rungs flatten, "
        "which is the recovery the rest of this dashboard measures."
    )

#: why it happens ----------------------------------------------------------
uni = D.gene_universe()
if not uni.empty and "loeuf_decile" in uni:
    prof = (
        uni.groupby("loeuf_decile")
        .agg(median_tpm=("tpm_reference", "median"),
             mean_exons=("n_coding_exons", "mean"),
             mean_tissues=("n_brain_tissues_expressed", "mean"),
             genes=("symbol", "size"))
        .reset_index()
    )
    fig = go.Figure(go.Scatter(
        x=prof["loeuf_decile"], y=prof["median_tpm"],
        mode="lines+markers", line=dict(color=U.GREY, width=3),
        marker=dict(size=9),
        hovertemplate="decile %{x}<br>median %{y:.2f} TPM<extra></extra>",
    ))
    with st.container(border=True):
        st.markdown("**The cause: expression tracks constraint almost perfectly**")
        U.plot(
            fig, height=320,
            xaxis=dict(title="LOEUF decile (0 = most constrained)",
                       dtick=1, tickmode="linear"),
            yaxis=dict(title="Median cortex TPM (log scale)", type="log"),
            showlegend=False,
        )
        lo = prof.loc[prof["loeuf_decile"] == 0].iloc[0]
        hi = prof.loc[prof["loeuf_decile"] == 9].iloc[0]
        st.caption(
            f"Median cortex expression falls from {lo['median_tpm']:.2f} to "
            f"{hi['median_tpm']:.2f} TPM across the deciles, a "
            f"{lo['median_tpm'] / hi['median_tpm']:.0f}-fold drop. Coding "
            f"exons fall {lo['mean_exons']:.1f} → {hi['mean_exons']:.1f} "
            f"and brain tissues expressed {lo['mean_tissues']:.1f} → "
            f"{hi['mean_tissues']:.1f} (means). The least-constrained decile "
            "is barely transcribed, so it has few detectable eQTLs for reasons "
            "that have nothing to do with selection, and it drags the "
            "unadjusted comparison the wrong way."
        )

    with st.expander("The per-decile numbers", icon=":material/table:"):
        U.table(
            prof.assign(
                Decile=prof["loeuf_decile"].astype(int),
                Genes=prof["genes"].astype(int),
                **{
                    "Median cortex TPM": prof["median_tpm"].round(2),
                    "Mean coding exons": prof["mean_exons"].round(1),
                    "Mean brain tissues": prof["mean_tissues"].round(1),
                },
            )[["Decile", "Genes", "Median cortex TPM", "Mean coding exons",
               "Mean brain tissues"]]
        )

#: pLI, the same story on a different metric --------------------------------
pli = D.by_pli()
if not pli.empty:
    p = pli[pli["arm"] == arm]
    order = ["<0.1", "0.1-0.5", "0.5-0.9", ">=0.9"]
    fig = go.Figure()
    for rung in RUNGS:
        d = p[p["rung"] == rung].set_index("bin_label").reindex(order).reset_index()
        if d["rate"].isna().all():
            continue
        colour, dash = STYLE[rung]
        fig.add_trace(go.Scatter(
            x=d["bin_label"], y=d["rate"], name=SHORT[rung],
            mode="lines+markers", line=dict(color=colour, width=2.5, dash=dash),
            marker=dict(size=8),
            hovertemplate=(
                f"<b>{D.RUNG_LABEL[rung]}</b><br>pLI %{{x}}"
                "<br>recovery %{y:.3f}<extra></extra>"
            ),
        ))
    with st.container(border=True):
        st.markdown("**The same pattern on a different constraint metric**")
        U.plot(
            fig, height=340,
            xaxis=dict(title="pLI bin"),
            yaxis=dict(title="Fraction with a detectable cis-eQTL"),
        )
        st.caption(
            "pLI is a different summary of the same underlying constraint, and "
            "it carries the same expression confound, which is why D-002 "
            "matches on expression rather than on the constraint metric. The "
            "robustness table on the recovery page re-runs the whole ladder on "
            "pLI and lands at 52.2% closure against 59.9%."
        )
