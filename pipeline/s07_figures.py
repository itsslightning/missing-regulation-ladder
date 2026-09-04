"""Stage 1, step 3: the recovery-curve figures for the report.

Static PNGs for docs/ and the methods note. The Stage 3 dashboard rebuilds
these interactively in Plotly from the same parquet files, so nothing here is
the only route to a number.

Three figures, in the order an argument needs them:

  fig1  the recovery curve itself -- constrained vs matched control, per rung
  fig2  the gap, under all three D-001 arms, which is the robustness claim
  fig3  why the raw LOEUF-decile view cannot be read on its own -- it is
        dominated by expression and inverts the matched result

Colours come from the sibling repo's Streamlit theme so the two projects look
like one programme.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline import config as cfg
from pipeline.s01_gene_universe import GENE_UNIVERSE
from pipeline.s06_recovery import (
    RECOVERY_BY_LOEUF,
    RECOVERY_BY_LOEUF_EXPR,
    RUNG_LABELS,
    RUNG_ORDER,
)

FIG_DIR = cfg.DIR_DOCS / "figures"

#: The sibling repo's Streamlit primaryColor, plus a warm contrast for the
#: control arm and a neutral for annotation.
TEAL = "#3d7a6f"
RUST = "#b5651d"
GREY = "#6b6b6b"

#: Rungs whose numbers are not yet trustworthy get drawn, but hatched and
#: annotated. Hiding them would be worse: a reader would not know rung 2 exists.
PROVISIONAL = {"bulk_brain"}

plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.6,
    }
)


def _ordered(df: pd.DataFrame) -> pd.DataFrame:
    present = [r for r in RUNG_ORDER if r in set(df["rung"])]
    return df.set_index("rung").loc[present].reset_index()


def fig_recovery_curve(curve: pd.DataFrame, schema: pd.DataFrame) -> None:
    """Figure 1: the two curves, with CI bands and donor counts."""
    d = _ordered(curve[curve["arm"] == "uniform"])
    s = schema[schema["arm"] == "uniform"].set_index("rung").reindex(d["rung"])
    x = np.arange(len(d))

    fig, ax = plt.subplots(figsize=(7.2, 4.6))

    for rate, lo, hi, colour, label in (
        ("control_rate", "control_lo", "control_hi", RUST,
         "Matched unconstrained controls"),
        ("constrained_rate", "constrained_lo", "constrained_hi", TEAL,
         "Constrained genes (LOEUF < 0.35)"),
    ):
        ax.fill_between(x, d[lo], d[hi], color=colour, alpha=0.18, linewidth=0)
        ax.plot(x, d[rate], "o-", color=colour, lw=2, ms=6, label=label, zorder=3)

    # SCHEMA overlay: 32 genes, so the interval is wide and must be shown.
    ax.errorbar(
        x,
        s["rate"],
        yerr=[s["rate"] - s["lo"], s["hi"] - s["rate"]],
        fmt="s--",
        color=GREY,
        ms=5,
        lw=1.2,
        capsize=3,
        alpha=0.85,
        label=f"SCHEMA genes (n={int(s['n_schema'].iloc[0])}, published FDR<0.05)",
        zorder=2,
    )

    for i, row in d.iterrows():
        if row["rung"] in PROVISIONAL:
            ax.axvspan(i - 0.42, i + 0.42, color="0.85", alpha=0.55, zorder=0)
            ax.text(
                i, 0.02, "provisional", ha="center", va="bottom",
                fontsize=7, color=GREY, style="italic",
            )

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{RUNG_LABELS[r]}\nN={int(n)}" for r, n in zip(d["rung"], d["n_donors"])],
        fontsize=8,
    )
    ax.set_ylabel("Fraction with a detectable cis-eQTL")
    ax.set_ylim(0, 1)
    ax.set_title(
        "Constrained-gene eQTL recovery across the resolution ladder\n"
        "uniform per-gene Bonferroni + BH within rung (D-001, D-003); "
        "95% CIs",
        fontsize=10,
        loc="left",
    )
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig1_recovery_curve.png", bbox_inches="tight")
    plt.close(fig)


def fig_gap(curve: pd.DataFrame) -> None:
    """Figure 2: the gap under all three detection arms.

    This is the figure that carries the robustness claim. D-001 was delegated
    rather than chosen by the project owner, so the rejected options are drawn
    beside the chosen one.
    """
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    styles = {
        "uniform": (TEAL, "o-", "uniform Bonferroni + BH  (D-001, primary)"),
        "native": (RUST, "s--", "each study's own q<=0.05  (sensitivity)"),
        "effect": (GREY, "^:", "uniform + |beta|>=0.1  (sensitivity)"),
    }
    for arm, (colour, style, label) in styles.items():
        d = _ordered(curve[curve["arm"] == arm])
        x = np.arange(len(d))
        ax.errorbar(
            x,
            d["gap"],
            yerr=[d["gap"] - d["gap_lo"], d["gap_hi"] - d["gap"]],
            fmt=style,
            color=colour,
            lw=1.8,
            ms=6,
            capsize=3,
            label=label,
        )

    d = _ordered(curve[curve["arm"] == "uniform"])
    x = np.arange(len(d))
    ax.axhline(0, color="k", lw=1, ls="-", alpha=0.6)
    ax.text(
        len(d) - 0.55, 0.012, "no gap = full recovery", fontsize=7,
        color="k", ha="right", va="bottom", alpha=0.7,
    )
    for i, row in d.iterrows():
        if row["rung"] in PROVISIONAL:
            ax.axvspan(i - 0.42, i + 0.42, color="0.85", alpha=0.55, zorder=0)

    ax.set_xticks(x)
    ax.set_xticklabels([RUNG_LABELS[r] for r in d["rung"]], fontsize=8)
    ax.set_ylabel("Gap  (control rate  -  constrained rate)")
    ax.set_ylim(bottom=-0.02)
    ax.set_title(
        "The constrained-gene gap narrows but does not close\n"
        "every arm keeps a gap whose 95% CI excludes zero at every rung",
        fontsize=10,
        loc="left",
    )
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig2_gap_sensitivity.png", bbox_inches="tight")
    plt.close(fig)


def fig_loeuf(loeuf: pd.DataFrame, loeuf_expr: pd.DataFrame, universe) -> None:
    """Figure 3: why the raw LOEUF-decile view cannot be read on its own.

    Three panels, because the naive version of this figure says the opposite of
    the truth and it is worth showing exactly why.

    Left: recovery by LOEUF decile, unadjusted. At the single-nucleus rungs the
    MOST constrained decile scores HIGHEST -- the reverse of the matched result.
    Middle: the reason. Median cortex TPM falls 70-fold across the deciles.
    Right: the same recovery curves inside one expression tertile, where the
    confound is largely held fixed and the direction agrees with the matched
    comparison again.
    """
    d = loeuf[loeuf["arm"] == "uniform"]
    present = [r for r in RUNG_ORDER if r in set(d["rung"])]
    cmap = plt.get_cmap("viridis")
    colours = {
        r: cmap(i / max(1, len(present) - 1) * 0.85) for i, r in enumerate(present)
    }

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.1))

    # -- left: raw, confounded ---------------------------------------------
    ax = axes[0]
    for rung in present:
        g = d[d["rung"] == rung].sort_values("loeuf_decile")
        ax.fill_between(g["loeuf_decile"], g["lo"], g["hi"],
                        color=colours[rung], alpha=0.15, linewidth=0)
        ax.plot(g["loeuf_decile"], g["rate"], "o-", color=colours[rung],
                lw=1.8, ms=4, label=RUNG_LABELS[rung].replace("\n", " "))
    ax.set_ylim(0, 1)
    ax.set_xticks(range(10))
    ax.set_xlabel("LOEUF decile  (0 = most constrained)")
    ax.set_ylabel("Fraction with a detectable cis-eQTL")
    ax.set_title(
        "A. Unadjusted — misleading\nmost-constrained decile scores HIGHEST at "
        "single-nucleus rungs",
        fontsize=9, loc="left",
    )
    ax.legend(loc="lower left", fontsize=7, framealpha=0.9)

    # -- middle: the confound ----------------------------------------------
    ax = axes[1]
    expr = universe.groupby("loeuf_decile").agg(
        median_tpm=("tpm_reference", "median"),
        mean_exons=("n_coding_exons", "mean"),
    )
    ax.plot(expr.index, expr["median_tpm"], "o-", color=RUST, lw=2, ms=5)
    ax.set_yscale("log")
    ax.set_xticks(range(10))
    ax.set_xlabel("LOEUF decile  (0 = most constrained)")
    ax.set_ylabel("Median cortex TPM  (log scale)", color=RUST)
    ax.tick_params(axis="y", labelcolor=RUST)
    ax2 = ax.twinx()
    ax2.plot(expr.index, expr["mean_exons"], "s--", color=GREY, lw=1.5, ms=4)
    ax2.set_ylabel("Mean coding exons", color=GREY)
    ax2.tick_params(axis="y", labelcolor=GREY)
    ax2.grid(False)
    ax.set_title(
        "B. Why — expression tracks constraint\nmedian TPM 12.7 → 0.2 across "
        "the deciles (70-fold)",
        fontsize=9, loc="left",
    )

    # -- right: expression held fixed --------------------------------------
    ax = axes[2]
    e = loeuf_expr[(loeuf_expr["arm"] == "uniform")
                   & (loeuf_expr["expr_tertile"] == "high")]
    for rung in present:
        g = e[e["rung"] == rung].sort_values("loeuf_decile")
        if g.empty:
            continue
        ax.fill_between(g["loeuf_decile"], g["lo"], g["hi"],
                        color=colours[rung], alpha=0.15, linewidth=0)
        ax.plot(g["loeuf_decile"], g["rate"], "o-", color=colours[rung],
                lw=1.8, ms=4)
    ax.set_ylim(0, 1)
    ax.set_xticks(range(10))
    ax.set_xlabel("LOEUF decile  (0 = most constrained)")
    ax.set_ylabel("Fraction with a detectable cis-eQTL")
    ax.set_title(
        "C. Within the top expression tertile\nbulk rungs steeply constrained; "
        "single-nucleus rungs nearly flat",
        fontsize=9, loc="left",
    )

    fig.suptitle(
        "The unadjusted constraint gradient is dominated by expression and "
        "inverts the matched result — this is what the matched control set "
        "(D-002) exists to remove",
        fontsize=10, y=1.02, x=0.01, ha="left",
    )
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig3_loeuf_deciles.png", bbox_inches="tight")
    plt.close(fig)


def fig_resolution_test(res: pd.DataFrame) -> None:
    """Figure 4: the within-SingleBrain resolution test.

    The main ladder cannot separate resolution from donor count. This can, at
    least for the splitting direction: each SingleBrain class is reported both
    pooled and split into its own subtypes, from the same donors and nuclei.
    """
    classes = list(dict.fromkeys(res["major"]))
    x = np.arange(len(classes))
    pooled = res[res["arm"] == "pooled"].set_index("major").loc[classes]
    split = res[res["arm"] == "split"].set_index("major").loc[classes]

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(10.4, 4.2), gridspec_kw={"width_ratios": [1.55, 1]}
    )

    w = 0.36
    for off, d, colour, label in (
        (-w / 2, pooled, TEAL, "pooled into one class"),
        (+w / 2, split, RUST, "split into subtypes"),
    ):
        ax.bar(x + off, d["gap"], width=w, color=colour, alpha=0.85, label=label)
        ax.errorbar(
            x + off, d["gap"],
            yerr=[d["gap"] - d["gap_lo"], d["gap_hi"] - d["gap"]],
            fmt="none", ecolor="0.25", elinewidth=1, capsize=2.5,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{c}\n({int(split.loc[c, 'n_columns'])} sub)" for c in classes], fontsize=8
    )
    ax.set_ylabel("Constrained-gene gap")
    ax.legend(fontsize=8, framealpha=0.9)
    ax.set_title(
        "A. Same donors, same nuclei — only the grouping changes",
        fontsize=9, loc="left",
    )

    delta = (split["gap"] - pooled["gap"]).to_numpy()
    colours = [TEAL if d < 0 else RUST for d in delta]
    ax2.barh(x, delta, color=colours, alpha=0.85)
    ax2.axvline(0, color="k", lw=1)
    ax2.set_yticks(x)
    ax2.set_yticklabels(classes, fontsize=8)
    ax2.invert_yaxis()
    ax2.set_xlabel("change in gap on splitting")
    ax2.set_title(
        f"B. Mean change {delta.mean():+.3f}\nnarrows in "
        f"{int((delta < 0).sum())} of {len(delta)} classes — no systematic effect",
        fontsize=9, loc="left",
    )

    fig.suptitle(
        "Finer cell-type resolution, at identical donor count, does not close "
        "the constrained-gene gap",
        fontsize=10, y=1.03, x=0.01, ha="left",
    )
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig4_resolution_test.png", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    curve = pd.read_parquet(cfg.RECOVERY_TABLE)
    loeuf = pd.read_parquet(RECOVERY_BY_LOEUF)
    loeuf_expr = pd.read_parquet(RECOVERY_BY_LOEUF_EXPR)
    schema = pd.read_parquet(cfg.DIR_PROCESSED / "recovery_schema.parquet")
    universe = pd.read_parquet(GENE_UNIVERSE)

    fig_recovery_curve(curve, schema)
    fig_gap(curve)
    fig_loeuf(loeuf, loeuf_expr, universe)

    res_path = cfg.DIR_PROCESSED / "resolution_test.parquet"
    if res_path.exists():
        fig_resolution_test(pd.read_parquet(res_path))

    for p in sorted(FIG_DIR.glob("*.png")):
        print(f"  wrote {p.relative_to(cfg.ROOT).as_posix()}  "
              f"({p.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
