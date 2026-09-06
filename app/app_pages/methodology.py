"""Decisions, provenance and limitations -- the audit trail behind every number."""

from __future__ import annotations

import streamlit as st

from lib import data as D

st.header("Methods")
st.markdown(
    "Every number in this dashboard traces to a parquet file, which traces to "
    "a stage script, which traces to a recorded decision. This page is that "
    "trail."
)

tabs = st.tabs(
    ["Decisions", "Data provenance", "Limitations", "Detection rule"]
)

# --- decisions -----------------------------------------------------------
with tabs[0]:
    st.markdown(
        "Four choices could each move the headline on their own, so the "
        "pipeline held them as sentinels that **raise on use** until recorded. "
        "Analysis code depending on one could not run and silently default."
    )
    for d in D.decisions():
        with st.container(border=True):
            st.markdown(f"**{d['key']} — {d['question']}**")
            st.markdown(f"Chosen: `{d['chosen']}` — {d['meaning']}")
            with st.expander("Why, and what was rejected"):
                st.markdown(d.get("rationale") or "_no rationale recorded_")
                rejected = d.get("rejected") or {}
                if rejected:
                    st.markdown("**Alternatives not taken:**")
                    for k, v in rejected.items():
                        st.markdown(f"- `{k}` — {v}")
            st.caption(f"Recorded {d['timestamp_utc']}")

# --- provenance ----------------------------------------------------------
with tabs[1]:
    dl = D.downloads()
    if dl.empty:
        st.info("No download log found.")
    else:
        st.markdown(
            f"**{len(dl)} files**, each stamped with the URL it came from, the "
            "UTC date, the byte count and a sha256. eQTL catalogues are "
            "re-released, so *“I downloaded GTEx brain eQTLs”* is not a "
            "reproducible statement — a hash and a date are."
        )
        by_source = dl.groupby(
            ["source_name", "licence", "redistributable"], as_index=False
        ).agg(files=("file", "count"), bytes=("bytes", "sum"))
        by_source["Size"] = (by_source["bytes"] / 1e6).map("{:,.0f} MB".format)
        by_source["Redistributable"] = by_source["redistributable"].map(
            {True: "yes", False: "no — restricted"}
        )
        st.dataframe(
            by_source[["source_name", "files", "Size", "licence",
                       "Redistributable"]].rename(
                columns={"source_name": "Source", "files": "Files",
                         "licence": "Licence"}),
            hide_index=True, use_container_width=True,
        )
        st.caption(
            "Sources marked restricted are written to `data/restricted/`, "
            "which `.gitignore` excludes as a whole directory. Only derived "
            "outputs are published from them, and only where they cannot be "
            "used to reconstruct the source table (D-013)."
        )
        with st.expander("Every file, with hashes"):
            show = dl.copy()
            show["MB"] = (show["bytes"] / 1e6).round(2)
            show["sha256"] = show["sha256"].str.slice(0, 16) + "…"
            st.dataframe(
                show[["source_name", "file", "version", "MB",
                      "downloaded_utc", "sha256"]],
                hide_index=True, use_container_width=True, height=340,
            )

# --- limitations ---------------------------------------------------------
with tabs[2]:
    st.markdown(
        """
### What this project does not claim

- **Not a winner between the two hypotheses.** The output is a decomposition
  with uncertainty. Mostafavi's selection account predicts the gap persists at
  any resolution — it does not. Rosen's power account predicts it closes with
  sample size — it does not do that either. What closes it is the assay, which
  is a third account neither proposed.
- **Not a mechanism.** Bulk and single-nucleus differ in more than assay —
  ancestry, brain region, pipeline, GENCODE vintage. The evidence that this is
  assay rather than a study-level accident is that two bulk studies agree
  exactly despite 6.8× different N and different consortia, and three
  single-nucleus arms agree despite 5.1× different N. A confound would have to
  track assay class across four independent datasets.
- **Not colocalization.** SMR cannot separate a shared causal variant from
  linkage, so Stage 2 counts inflate. The HEIDI and coloc sensitivity arms are
  not yet run.
- **Not settled numbers from preprints.** Rosen et al. 2026 and the Nov-2025
  brain/blood snRNA preprint have not completed peer review.

### Known limitations, stated rather than buried

- **Absolute eGene counts are below every published figure**, deliberately.
  Per-gene Bonferroni over cis variants ignores LD and is stricter than
  permutation. It is applied identically at every rung, which is what a
  cross-rung comparison needs, and it is why these numbers should not be
  compared to a paper's headline eGene count.
- **Resolution and per-context power are intrinsically coupled** in
  single-cell data. No test here isolates resolution alone; two tests vary it
  in opposite directions and bracket it.
- **The unadjusted LOEUF-decile view is misleading and inverts the result.**
  Median cortex TPM falls 12.71 → 0.18 across the deciles, and at
  single-nucleus rungs that confound makes the *most* constrained genes look
  best recovered. This is exactly what the matched control set exists to
  remove.
- **Locus definition is distance-based**, not LD-based, and the MHC is
  excluded.
"""
    )

# --- detection rule ------------------------------------------------------
with tabs[3]:
    st.markdown(
        """
### What counts as a "detectable" cis-eQTL

Each rung ships a different significance column and they are **not on a common
evidential scale**. GTEx and SingleBrain both ship Storey q-values, and a
Storey q depends on π₀ — the estimated fraction of true nulls — computed
*within each study*. A better-powered study has a lower π₀, which makes
q ≤ 0.05 a **more permissive** bar there.

That matters concretely: the threshold loosens exactly where power is highest,
which is the direction that would manufacture the recovery the power account
predicts. Using study-native calls as primary reports **88% gap closure**
against **60%** under the uniform rule.

So detection is recomputed identically at every rung:

1. **Per-gene Bonferroni** over the cis variants tested for that gene —
   computable everywhere from what each source ships.
2. **Benjamini–Hochberg across genes within the rung** (D-003), applied over
   all gene × cell-type tests so a rung with 28 cell types pays for its 28
   opportunities.

Use the toggle on the *Recovery curve* page to see the rejected rules.
"""
    )
    rob = D.robustness()
    if not rob.empty:
        r = rob[rob["closure"].notna()]
        st.metric(
            "Closure across every structural variant tested",
            f"{r['closure'].min():.1%} – {r['closure'].max():.1%}",
            "residual gap excludes zero in all of them",
            delta_color="off", border=True,
        )
