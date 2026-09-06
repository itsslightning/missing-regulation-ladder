"""SMR causal links per rung, and the druggability of newly-visible genes."""

from __future__ import annotations

import plotly.graph_objects as go
import streamlit as st

from lib import data as D
from lib import ui as U

st.header("Causal links and druggability")
st.markdown(
    "If single-nucleus data only recovers eQTLs that were statistically "
    "marginal, the extra genes will not carry causal signal. This runs SMR of "
    "schizophrenia risk on brain gene expression at every rung (PGC3 European, "
    "D-004, D-014) to see whether they do."
)

loci = D.locus_explanation()
smr = D.smr_results()
drug = D.druggability()

if loci.empty or smr.empty:
    U.missing("The SMR tables")

#: loci explained ---------------------------------------------------------
order = D.ordered_rungs(loci)
l = loci.set_index("rung").reindex(order).reset_index()
colours = [U.RUST if D.RUNG_ASSAY[r] == "bulk tissue" else U.TEAL
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

with st.container(border=True):
    st.markdown("**GWAS loci gaining an eQTL explanation**")
    U.plot(
        fig, height=400,
        yaxis=dict(title=f"Loci explained (of {int(l['loci_total'].iloc[0])})"),
        showlegend=False,
    )
    st.caption(
        "Rust is bulk tissue, teal is single-nucleus. PsychENCODE has the most "
        "donors on the ladder at 1,387 and explains 27 loci; SingleBrain has "
        "29% fewer donors and explains 59. Clumping is distance-based and the "
        "MHC is excluded, so what this measures is the ratio between rungs, "
        "not the absolute count."
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
        "near-constant across rungs",
        delta_color="off", delta_arrow="off", border=True,
    )

st.caption(
    "That last number locates the gain precisely. The significance **rate** "
    "barely moves, so single-nucleus data is not making eQTLs more likely to "
    "be causal. It is making more genes testable at all."
)

#: druggability -----------------------------------------------------------
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
    with U.card("Where the causal link is visible"):
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
        U.table(
            show[["Group", "Genes", "Tractable (small molecule)",
                  "Phase 1+", "SCHEMA"]]
        )
        st.caption(
            "Tractability rates come out near identical across the groups, "
            "which is the sanity check this table exists for. Single-nucleus "
            "data is not enriching for druggable genes, it is finding more "
            "genes at the same druggable rate."
        )

    adv = drug[(drug["group"] == "sn_only") & drug["clinically_advanced"]]
    if not adv.empty:
        with U.card(
            f"{len(adv)} clinically-advanced targets implicated only in "
            "single-nucleus data"
        ):
            U.table(
                adv.sort_values("symbol")[["symbol", "sm_tier", "is_schema"]]
                .rename(columns={"symbol": "Gene", "sm_tier": "Tractability",
                                 "is_schema": "SCHEMA"}),
                height=260,
            )
            st.caption(
                "A bulk-tissue eQTL screen would not have surfaced these as "
                "causal. The calcium-channel family (CACNA1C, CACNA1D, "
                "CACNA1I, CACNB2) has been a schizophrenia target class for a "
                "decade."
            )

    #: the CHRM4 callback -------------------------------------------------
    musc = drug[drug["is_muscarinic"]]
    st.subheader("The CHRM4 callback, an informative negative")
    with st.container(border=True):
        st.badge(
            "No detectable cis-eQTL at any assay or resolution",
            icon=":material/target:", color="orange",
        )
        st.markdown(
            "CHRM4 is in the constrained case set of this very study "
            "(LOEUF 0.265, decile 0; pLI 0.974). It is well expressed in "
            "cortex at 11.14 TPM, so this is not a low-expression artefact. "
            "It has been an approved drug target since Cobenfy's 2024 "
            "approval. And no muscarinic receptor reaches SMR significance "
            "here."
        )
        st.markdown(
            "An eQTL-based target discovery pipeline would never have "
            "surfaced it through regulatory evidence. CHRM4 sits in the "
            "residual part of the gap that no assay and no resolution closes, "
            "which is what missing regulation costs in practice."
        )
        if not musc.empty:
            U.table(
                musc.assign(group=musc["group"].map(label).fillna(musc["group"]))
                [["symbol", "group", "sm_tier"]].rename(
                    columns={"symbol": "Gene", "group": "Visible in",
                             "sm_tier": "Tractability"})
            )
            st.caption(
                "CHRM4 is not in this table, and that is the finding. The "
                "table is built from genes that had an eQTL to test, so a "
                "gene with no detectable eQTL at any rung never enters it."
            )

with st.expander(
    "SMR is not colocalization, and what that means for these counts",
    icon=":material/warning:",
):
    st.markdown(
        "SMR cannot separate a shared causal variant from linkage between two "
        "distinct causal variants, so every count on this page is inflated "
        "relative to a true colocalization and should not be quoted against "
        "published coloc figures. D-004 accepts that deliberately: the bias "
        "has the same construction at every rung, and the question here is a "
        "comparison between rungs rather than an absolute count. The HEIDI "
        "and coloc sensitivity arms are not yet run."
    )
