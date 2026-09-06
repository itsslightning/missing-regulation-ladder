"""SMR causal links per rung, and the druggability of newly-visible genes."""

from __future__ import annotations

import plotly.graph_objects as go
import streamlit as st

from lib import data as D

st.header("Causal links and druggability")
st.markdown(
    "SMR of schizophrenia risk on brain gene expression (PGC3 European, "
    "D-004, D-014), and what it means for target discovery."
)

loci = D.locus_explanation()
smr = D.smr_results()
drug = D.druggability()

if loci.empty or smr.empty:
    st.error("No SMR tables. Run `python scripts/run_all.py` first.")
    st.stop()

# --- loci explained ------------------------------------------------------
st.subheader("GWAS loci gaining an eQTL explanation")

order = D.ordered_rungs(loci)
l = loci.set_index("rung").reindex(order).reset_index()
colours = [D.RUST if D.RUNG_ASSAY[r] == "bulk tissue" else D.TEAL
           for r in l["rung"]]

fig = go.Figure(go.Bar(
    x=[D.RUNG_LABEL[r] for r in l["rung"]],
    y=l["loci_explained"],
    marker_color=colours,
    text=[f"{v} ({f:.1%})" for v, f in
          zip(l["loci_explained"], l["frac_explained"])],
    textposition="outside",
    customdata=l[["genes"]].assign(
        donors=[D.RUNG_DONORS[r] for r in l["rung"]],
        assay=[D.RUNG_ASSAY[r] for r in l["rung"]]),
    hovertemplate=(
        "<b>%{x}</b><br>%{y} loci explained"
        "<br>%{customdata[0]} genes"
        "<br>%{customdata[1]} donors, %{customdata[2]}<extra></extra>"
    ),
))
fig.update_layout(
    height=400,
    yaxis=dict(title=f"Loci explained (of {int(l['loci_total'].iloc[0])})"),
    margin=dict(t=30, b=10, l=10, r=10),
)
st.plotly_chart(fig, use_container_width=True)

st.caption(
    "Orange = bulk tissue, teal = single-nucleus. **PsychENCODE has the most "
    "donors on the ladder (1,387) and explains 27 loci; SingleBrain has 29% "
    "fewer donors and explains 59.** Distance-based clumping, MHC excluded — "
    "so the ratio between rungs is what this measures, not the absolute count."
)

with st.container(horizontal=True):
    bulk_loci = l[l["rung"].isin(["gtex_cortex", "bulk_brain"])]["loci_explained"]
    sn_loci = l[l["rung"].isin(["sn_major", "sn_subtype"])]["loci_explained"]
    st.metric("Best bulk rung", f"{int(bulk_loci.max())} loci", border=True)
    st.metric("Best single-nucleus rung", f"{int(sn_loci.max())} loci",
              f"+{int(sn_loci.max() - bulk_loci.max())} with fewer donors",
              border=True)
    rate = smr[smr["smr_sig"]].groupby("rung")["gene_id"].nunique()
    tot = smr.groupby("rung")["gene_id"].nunique()
    st.metric(
        "SMR significance rate",
        f"{(rate / tot).min():.1%}–{(rate / tot).max():.1%}",
        "near-constant across rungs", delta_color="off", border=True,
    )

st.info(
    "The significance **rate** barely moves across rungs, which locates the "
    "gain precisely: single-nucleus data does not make eQTLs more likely to be "
    "causal — it makes **more genes testable at all**.",
    icon=":material/lightbulb:",
)

# --- druggability --------------------------------------------------------
if not drug.empty:
    st.subheader("Druggability of the newly-visible genes")

    summary = drug.groupby("group").agg(
        genes=("gene_id", "nunique"),
        tractable=("druggable_sm", "sum"),
        clinical=("clinically_advanced", "sum"),
        schema=("is_schema", "sum"),
    ).reindex(["sn_only", "bulk_and_sn", "bulk_only"]).dropna(how="all")
    summary["rate"] = summary["tractable"] / summary["genes"]

    label = {
        "sn_only": "Single-nucleus only",
        "bulk_and_sn": "Bulk and single-nucleus",
        "bulk_only": "Bulk only",
    }
    with st.container(border=True):
        st.markdown("**Where the causal link is visible**")
        show = summary.reset_index().assign(
            Group=lambda d: d["group"].map(label),
            Genes=lambda d: d["genes"].astype(int),
            **{
                "Tractable (small molecule)": lambda d: [
                    f"{int(t)} ({r:.1%})" for t, r in zip(d["tractable"], d["rate"])],
                "Phase 1+": lambda d: d["clinical"].astype(int),
                "SCHEMA": lambda d: d["schema"].astype(int),
            },
        )
        st.dataframe(
            show[["Group", "Genes", "Tractable (small molecule)",
                  "Phase 1+", "SCHEMA"]],
            hide_index=True, use_container_width=True,
        )
        st.caption(
            "Tractability rates are near identical across groups — the right "
            "sanity check. Single-nucleus data is not enriching for druggable "
            "genes, it is finding **more genes at the same druggable rate**."
        )

    adv = drug[(drug["group"] == "sn_only") & drug["clinically_advanced"]]
    if not adv.empty:
        with st.container(border=True):
            st.markdown(
                f"**{len(adv)} clinically-advanced targets causally implicated "
                "ONLY in single-nucleus data**"
            )
            st.dataframe(
                adv.sort_values("symbol")[["symbol", "sm_tier", "is_schema"]]
                .rename(columns={"symbol": "Gene", "sm_tier": "Tractability",
                                 "is_schema": "SCHEMA"}),
                hide_index=True, use_container_width=True, height=260,
            )
            st.caption(
                "A bulk-tissue eQTL screen would not have surfaced these as "
                "causal. The calcium-channel family — CACNA1C, CACNA1D, "
                "CACNA1I, CACNB2 — has been a schizophrenia target class for "
                "a decade."
            )

    # --- the CHRM4 callback ---------------------------------------------
    musc = drug[drug["is_muscarinic"]]
    st.subheader("The CHRM4 callback — an informative negative")
    st.markdown(
        "CHRM4 is **in the constrained case set of this very study** "
        "(LOEUF 0.265, decile 0; pLI 0.974), **well expressed in cortex** "
        "(11.14 TPM, so not a low-expression artefact), and **an approved "
        "drug target** since Cobenfy's 2024 approval — yet it has "
        "**no detectable cis-eQTL at any assay or resolution tested here**, "
        "and no muscarinic receptor reaches SMR significance."
    )
    if not musc.empty:
        st.dataframe(
            musc[["symbol", "group", "sm_tier"]].rename(
                columns={"symbol": "Gene", "group": "Visible in",
                         "sm_tier": "Tractability"}),
            hide_index=True, use_container_width=True,
        )
    st.error(
        "**This is the project's thesis in one gene.** An eQTL-based target "
        "discovery pipeline would never have surfaced CHRM4 through "
        "regulatory evidence. It lives in the residual ~40% of the gap that "
        "no assay and no resolution closes — which is the practical cost of "
        "missing regulation.",
        icon=":material/target:",
    )

st.warning(
    "**SMR is not colocalization.** It cannot separate a shared causal variant "
    "from linkage between two distinct causal variants, so these counts are "
    "inflated relative to a true coloc and should not be quoted against "
    "published coloc figures. D-004 accepts this deliberately: the bias has "
    "the same construction at every rung, and the question here is a "
    "cross-rung comparison.",
    icon=":material/warning:",
)
