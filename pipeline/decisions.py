"""The four choices that can move the headline conclusion, held unset on purpose.

config.py explains why these are not constants. In short: what counts as a
"detectable" eQTL, how tightly control genes are matched, which multiple-testing
correction is used, and the colocalization priors are each capable of deciding
the selection-vs-power question on their own. A default written here would be a
conclusion smuggled in as a config value.

So each one is an `OpenDecision`. It has a question, the alternatives that were
actually on the table, and no value. Reading `.value` before the decision is
recorded raises `UndecidedError` naming the decision and pointing at DECISIONS.md.
Code that depends on one of these therefore cannot run early and quietly pick a
default, it stops.

Recording a decision is a deliberate, auditable act:

    from pipeline import decisions
    decisions.DETECTABLE_EQTL.decide(
        "storey_q_0.05",
        rationale="Matches the threshold GTEx and SingleBrain both ship.",
    )

`decide()` refuses anything not in `alternatives`, so a decision cannot drift
into an option nobody wrote down and weighed. It also refuses to silently
overwrite an existing choice: changing your mind requires `force=True`, which
is logged, because a mid-analysis redefinition of significance is exactly the
kind of thing a methods section has to disclose.

Every call appends to logs/decisions.jsonl, which is the machine-readable twin
of DECISIONS.md. DECISIONS.md remains the human record and is the file the
thesis methods section is written from; this module exists so the code and that
document cannot disagree without someone noticing.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

from pipeline.config import DIR_LOGS

DECISION_LOG = DIR_LOGS / "decisions.jsonl"


class UndecidedError(RuntimeError):
    """Raised when analysis code reaches a choice the human has not yet made."""


@dataclass
class OpenDecision:
    """One methodological choice, its alternatives, and (once made) its value.

    `key` matches the D-00x identifier used in DECISIONS.md so the two records
    can be lined up by eye. `alternatives` maps an option name to the one-line
    statement of what choosing it would mean, written before the data is
    seen, so the menu is not quietly rewritten to fit the result.
    """

    key: str
    question: str
    alternatives: dict[str, str]
    why_it_matters: str
    _value: Any = field(default=None, repr=False)
    _decided: bool = field(default=False, repr=False)
    _rationale: str | None = field(default=None, repr=False)

    @property
    def value(self) -> Any:
        if not self._decided:
            raise UndecidedError(
                f"\n\n  Decision {self.key} has not been made yet.\n"
                f"  Question: {self.question}\n"
                f"  Why it matters: {self.why_it_matters}\n"
                f"  Alternatives on the table:\n"
                + "".join(
                    f"    - {name}: {meaning}\n"
                    for name, meaning in self.alternatives.items()
                )
                + "\n  This is a call for the human running the project, not for\n"
                "  the pipeline. Record it with decisions."
                f"{self.key_as_attr()}.decide(<option>, rationale=...)\n"
                "  and write the corresponding entry in DECISIONS.md.\n"
            )
        return self._value

    @property
    def decided(self) -> bool:
        return self._decided

    @property
    def rationale(self) -> str | None:
        return self._rationale

    def key_as_attr(self) -> str:
        """The module-level name this decision is bound to, for error messages."""
        return _ATTR_BY_KEY.get(self.key, self.key)

    def decide(self, option: str, *, rationale: str, force: bool = False) -> None:
        if option not in self.alternatives:
            raise ValueError(
                f"{self.key}: {option!r} is not one of the alternatives that were "
                f"written down and weighed ({sorted(self.alternatives)}). If it is "
                "genuinely a new option, add it to decisions.py and DECISIONS.md "
                "first, so the record shows it was considered rather than "
                "back-filled."
            )
        if self._decided and not force:
            raise ValueError(
                f"{self.key} is already set to {self._value!r}. Re-deciding a "
                "threshold after seeing results changes what the p-values mean, so "
                "it requires force=True and an entry in DECISIONS.md saying what "
                "prompted the change."
            )
        previous = self._value if self._decided else None
        self._value = option
        self._decided = True
        self._rationale = rationale
        self._append_log(option, rationale, previous)

    def _append_log(self, option: str, rationale: str, previous: Any) -> None:
        DECISION_LOG.parent.mkdir(parents=True, exist_ok=True)
        record = {
            "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "key": self.key,
            "question": self.question,
            "chosen": option,
            "meaning": self.alternatives[option],
            "rejected": {k: v for k, v in self.alternatives.items() if k != option},
            "rationale": rationale,
            "superseded": previous,
        }
        with DECISION_LOG.open("a", encoding="utf-8", newline="\n") as fh:
            fh.write(json.dumps(record) + "\n")


# ---------------------------------------------------------------------------
# The four load-bearing decisions.
#
# Each is registered here with its alternatives stated *before* any result is
# known. The wording is deliberately about what the option means, not which is
# preferable: the argument for a particular choice belongs in DECISIONS.md
# next to the date it was made.
# ---------------------------------------------------------------------------

DETECTABLE_EQTL = OpenDecision(
    key="D-001",
    question=(
        "What counts as a gene having a 'detectable' cis-eQTL at a given rung?"
    ),
    why_it_matters=(
        "This is the y-axis of the recovery curve. A permissive definition makes "
        "every rung look successful and flatters the power account; a strict one "
        "suppresses the small effects that single-nucleus data is supposed to "
        "reveal and flatters the selection account. The rungs also do not ship "
        "the same statistics. GTEx and SingleBrain provide per-gene "
        "permutation q-values, Bryois provides only nominal p-values, so any "
        "definition has to say how that gap is bridged."
    ),
    alternatives={
        "study_native_threshold": (
            "Use each study's own significance call (GTEx qval<=0.05, SingleBrain "
            "qval<=0.05), and for Bryois reproduce an equivalent per-gene "
            "correction locally. Maximises comparability with each paper's "
            "published eGene counts; the correction differs slightly per rung."
        ),
        "uniform_recomputed_fdr": (
            "Ignore study-native calls and recompute one identical per-gene FDR "
            "across every rung from the top-association p-values. Maximises "
            "internal comparability; the numbers will no longer match the "
            "published eGene counts for any individual study."
        ),
        "fixed_nominal_p": (
            "A single fixed nominal p-value threshold on the top variant per gene "
            "at every rung. Simplest and fully transparent, but ignores that the "
            "number of variants tested per gene differs across rungs, so it "
            "quietly favours rungs with denser variant coverage."
        ),
        "effect_size_floor": (
            "Significance plus a minimum absolute effect size, so that 'detected' "
            "means the same biological magnitude at every rung. Directly targets "
            "the power confound, at the cost of discarding real small-effect "
            "eQTLs that are the point of higher resolution."
        ),
    },
)

CONTROL_MATCHING = OpenDecision(
    key="D-002",
    question=(
        "How tightly are unconstrained control genes matched to constrained/SCHEMA "
        "genes, and on which covariates?"
    ),
    why_it_matters=(
        "The control curve is the whole comparison. Constrained genes are longer, "
        "more highly expressed and have more exons than average, and every one of "
        "those also predicts eQTL discovery power. Loose matching leaves a gap "
        "that is really a length/expression artefact; over-tight matching on "
        "covariates that are themselves consequences of constraint can regress "
        "away the very effect being measured."
    ),
    alternatives={
        "caliper_expression_exons": (
            "Nearest-neighbour matching within a caliper on expression decile and "
            "coding-exon count only, per the report's stated caveat. Fewest "
            "assumptions; leaves gene length and TSS density unbalanced."
        ),
        "caliper_plus_length": (
            "As above plus gene length and number of cis variants tested. Closes "
            "the most obvious power confounds; risks conditioning on a mediator, "
            "since regulatory-landscape complexity is part of Mostafavi's proposed "
            "mechanism rather than a nuisance."
        ),
        "propensity_score": (
            "One propensity model over all covariates, matched on the score. Uses "
            "all covariates at once and degrades gracefully; harder to explain in "
            "a viva and hides which covariate is doing the work."
        ),
        "stratified_no_matching": (
            "No matched control set at all, compare across LOEUF deciles and "
            "adjust for covariates in a regression instead. Avoids discarding "
            "genes and makes the adjustment explicit; loses the simple, "
            "interpretable two-curve plot the report asks for."
        ),
    },
)

MULTIPLE_TESTING = OpenDecision(
    key="D-003",
    question=(
        "What multiple-testing correction is applied across genes and across the "
        "resolution rungs?"
    ),
    why_it_matters=(
        "Tests here are strongly dependent in two directions at once: the same "
        "gene is retested at every rung, and neighbouring cell subtypes share "
        "donors and nuclei. Benjamini-Hochberg assumes a dependence structure "
        "that this design does not have, so the choice changes how many genes "
        "count as recovered, which is the headline number."
    ),
    alternatives={
        "bh_within_rung": (
            "Benjamini-Hochberg across genes separately within each rung, with no "
            "correction across rungs. Treats each rung as its own experiment; the "
            "cross-rung comparison is then descriptive rather than tested."
        ),
        "by_across_all": (
            "Benjamini-Yekutieli across the full gene x rung grid. Valid under "
            "arbitrary dependence, which this design has; substantially more "
            "conservative and will reduce apparent recovery at every rung."
        ),
        "bh_across_all": (
            "Benjamini-Hochberg across the full gene x rung grid. More powerful "
            "than BY; its independence/PRDS assumption is not satisfied by "
            "repeated testing of the same gene across correlated cell types."
        ),
        "permutation_null": (
            "Build the null by permuting gene labels within matched strata and "
            "calibrate against that. Makes no parametric dependence assumption "
            "and is the most defensible; materially more compute and more code to "
            "get right."
        ),
    },
)

COLOC_PRIORS = OpenDecision(
    key="D-004",
    question=(
        "What priors and posterior thresholds define a colocalization between the "
        "SCZ GWAS and an eQTL at each rung? (Stage 2)"
    ),
    why_it_matters=(
        "coloc's p12 prior directly sets how readily a shared causal variant is "
        "declared. The default p12=1e-5 is known to be generous, and the number "
        "of loci that 'gain an eQTL explanation' at higher resolution: the "
        "Stage 2 headline, moves with it. Cell-subtype rungs test more "
        "features, compounding the effect."
    ),
    alternatives={
        "coloc_default_priors": (
            "coloc defaults (p1=1e-4, p2=1e-4, p12=1e-5), PP4>=0.8. What most "
            "papers report, so directly comparable; widely argued to over-declare "
            "colocalization."
        ),
        "coloc_conservative_p12": (
            "p12=1e-6 with PP4>=0.8. Guards against the known over-calling; will "
            "reduce absolute counts at every rung, and the cross-rung comparison "
            "is what matters here rather than the absolute number."
        ),
        "coloc_sensitivity_band": (
            "Report the result across a p12 range (1e-6 to 1e-5) as a sensitivity "
            "band rather than a single number. Most honest given the prior is "
            "genuinely uncertain; complicates the single decision-oriented "
            "headline statement the report asks for."
        ),
        "smr_heidi": (
            "SMR with the HEIDI test instead of coloc. Needs only top-SNP "
            "summary statistics, so it works where full associations are too "
            "large to download; HEIDI rejects linkage rather than confirming a "
            "shared variant, so it answers a subtly different question."
        ),
    },
)

_ATTR_BY_KEY = {
    "D-001": "DETECTABLE_EQTL",
    "D-002": "CONTROL_MATCHING",
    "D-003": "MULTIPLE_TESTING",
    "D-004": "COLOC_PRIORS",
}

ALL_DECISIONS = (
    DETECTABLE_EQTL,
    CONTROL_MATCHING,
    MULTIPLE_TESTING,
    COLOC_PRIORS,
)

_BY_KEY = {d.key: d for d in ALL_DECISIONS}


def _replay_log() -> None:
    """Restore decisions already recorded, so they survive across runs.

    Without this the sentinels would reset on every import and a decision made
    last week would look unmade today. The log is the source of truth and is
    replayed in order, so a later `force=True` correction supersedes the
    earlier entry exactly as it did when it was made.

    Restoration deliberately bypasses `decide()`: re-running it would append
    duplicate entries to the log every time the module is imported, and would
    trip the "already set" guard on the second replayed entry for a key.
    """
    if not DECISION_LOG.exists():
        return
    for line in DECISION_LOG.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        decision = _BY_KEY.get(rec.get("key"))
        if decision is None or rec.get("chosen") not in decision.alternatives:
            # An entry for a decision or option that no longer exists means the
            # code moved on from the log. Skipped rather than crashing, but the
            # mismatch is worth surfacing.
            continue
        decision._value = rec["chosen"]
        decision._decided = True
        decision._rationale = rec.get("rationale")


_replay_log()


def outstanding() -> list[OpenDecision]:
    """The decisions still unmade, for the Stage-0 report and the dashboard."""
    return [d for d in ALL_DECISIONS if not d.decided]


def summary() -> str:
    """Human-readable status block, printed at the end of each stage script."""
    lines = []
    for d in ALL_DECISIONS:
        mark = "set" if d.decided else "OPEN"
        detail = f"{d._value}" if d.decided else f"{len(d.alternatives)} alternatives"
        lines.append(f"  [{mark:>4}] {d.key} {d.key_as_attr():<18} {detail}")
    return "\n".join(lines)
