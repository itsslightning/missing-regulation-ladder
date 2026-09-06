"""Record what every pipeline step kept and what it dropped.

Carried over unchanged from the sibling scz-target-prioritization repo, where
silent gene loss was identified as the main correctness risk in a
GWAS-to-anything pipeline. It is at least as much of a risk here: this project
joins six catalogues that share no GENCODE vintage, and every join is a place
where genes can vanish without anything failing.

Every step that changes the number of genes in play calls one of these methods,
and the accumulated record is written to logs/sNN_*.json for the dashboard
methodology tab to read.

The point is not decoration. If a reviewer asks "how many constrained genes
never reached the recovery curve, and why", the answer should already be
written down rather than reconstructed.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import pandas as pd


@dataclass
class Step:
    """One transition, with enough detail to be argued with."""

    label: str
    n_in: int
    n_out: int
    detail: str = ""
    examples_dropped: list[str] = field(default_factory=list)
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def n_dropped(self) -> int:
        return self.n_in - self.n_out

    @property
    def pct_retained(self) -> float:
        return 100.0 * self.n_out / self.n_in if self.n_in else float("nan")

    def as_row(self) -> dict[str, Any]:
        d = asdict(self)
        d["n_dropped"] = self.n_dropped
        d["pct_retained"] = round(self.pct_retained, 2)
        d["examples_dropped"] = ", ".join(self.examples_dropped)
        d.pop("extra")
        return {**d, **self.extra}


class Provenance:
    """Accumulates steps for one pipeline run."""

    #: how many dropped identifiers to keep as illustrative examples
    N_EXAMPLES = 8

    def __init__(self, name: str, verbose: bool = True) -> None:
        self.name = name
        self.verbose = verbose
        self.started = datetime.now(timezone.utc).isoformat(timespec="seconds")
        self.steps: list[Step] = []

    #, recording ----------------------------------------------------------

    def record(
        self,
        label: str,
        n_in: int,
        n_out: int,
        detail: str = "",
        dropped: Iterable[Any] | None = None,
        **extra: Any,
    ) -> Step:
        examples: list[str] = []
        if dropped is not None:
            examples = [str(x) for x in list(dropped)[: self.N_EXAMPLES]]
        step = Step(
            label=label,
            n_in=int(n_in),
            n_out=int(n_out),
            detail=detail,
            examples_dropped=examples,
            extra=extra,
        )
        self.steps.append(step)
        if self.verbose:
            print(self._format(step))
        return step

    def filter(
        self,
        label: str,
        before: pd.DataFrame,
        after: pd.DataFrame,
        key: str,
        detail: str = "",
        **extra: Any,
    ) -> Step:
        """Record a row filter, capturing which keys disappeared."""
        lost = set(before[key]) - set(after[key])
        return self.record(
            label,
            len(before),
            len(after),
            detail=detail,
            dropped=sorted(lost),
            **extra,
        )

    def join(
        self,
        label: str,
        left: pd.DataFrame,
        right: pd.DataFrame,
        merged: pd.DataFrame,
        left_key: str,
        right_key: str | None = None,
        detail: str = "",
        **extra: Any,
    ) -> Step:
        """Record a merge from the left table's point of view.

        Reports how many left-hand keys found a partner, and keeps examples of
        the ones that did not. A join is where genes go missing most quietly,
        so this also records the right-hand table size for context.
        """
        right_key = right_key or left_key
        left_keys = set(left[left_key].dropna())
        right_keys = set(right[right_key].dropna())
        unmatched = left_keys - right_keys
        return self.record(
            label,
            len(left_keys),
            len(left_keys & right_keys),
            detail=detail,
            dropped=sorted(unmatched),
            n_right_keys=len(right_keys),
            n_merged_rows=len(merged),
            **extra,
        )

    def note(self, label: str, detail: str, **extra: Any) -> None:
        """Record something worth knowing that is not a count transition."""
        self.record(label, 0, 0, detail=detail, **extra)

    # -- output -----------------------------------------------------------

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame([s.as_row() for s in self.steps])

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "pipeline": self.name,
            "started_utc": self.started,
            "finished_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "steps": [s.as_row() for s in self.steps],
        }
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        csv_path = path.with_suffix(".csv")
        self.to_frame().to_csv(csv_path, index=False)
        if self.verbose:
            print(f"\nprovenance written to {path} and {csv_path}")
        return path

    def summary(self) -> str:
        counted = [s for s in self.steps if s.n_in]
        if not counted:
            return "no counted steps"
        lines = [f"=== {self.name}: {len(counted)} counted steps ==="]
        for s in counted:
            lines.append(self._format(s))
        return "\n".join(lines)

    @staticmethod
    def _format(step: Step) -> str:
        if not step.n_in:
            return f"  note   {step.label}: {step.detail}"
        head = (
            f"  {step.label:<44s} {step.n_in:>7,} -> {step.n_out:>7,}"
            f"  ({step.pct_retained:5.1f}% kept, {step.n_dropped:,} dropped)"
        )
        if step.detail:
            head += f"\n           {step.detail}"
        if step.examples_dropped:
            head += f"\n           e.g. {', '.join(step.examples_dropped[:5])}"
        return head


def require_columns(df: pd.DataFrame, columns: Sequence[str], where: str) -> None:
    """Fail loudly and early when an upstream file changes shape."""
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError(
            f"{where}: expected column(s) {missing} not found. "
            f"Available: {sorted(df.columns)[:20]}"
        )
