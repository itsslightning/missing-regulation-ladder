"""Record the Stage 0 decisions taken by the project owner on 2026-09-04.

Run once. `decide()` refuses to overwrite an existing choice without
`force=True`, so re-running is a no-op that reports what is already set rather
than quietly re-stamping the log.

Only D-002 lives in decisions.py -- it is one of the four sentinels that gate
analysis code. D-007, D-009 and D-010 are recorded in DECISIONS.md as prose,
because they shape the design rather than a single parameter the pipeline reads.
"""

from __future__ import annotations

from pipeline import decisions


def main() -> None:
    if decisions.CONTROL_MATCHING.decided:
        print(
            f"D-002 already recorded as {decisions.CONTROL_MATCHING.value!r}\n"
            f"  rationale: {decisions.CONTROL_MATCHING.rationale}"
        )
        return

    decisions.CONTROL_MATCHING.decide(
        "caliper_expression_exons",
        rationale=(
            "Chosen by the project owner, 2026-09-04. Matches the report's stated "
            "caveat exactly (expression level and exon count) and carries the "
            "fewest assumptions. Gene length is deliberately NOT matched on: it "
            "is a partial proxy for regulatory-landscape complexity, which is "
            "part of Mostafavi's proposed mechanism rather than a nuisance "
            "covariate, so matching on it would risk conditioning on a mediator "
            "and regressing away the effect being measured. The resulting "
            "residual imbalance in gene length is reported as a known limitation "
            "rather than hidden. Caliper 0.25 SD, the conventional tolerance."
        ),
    )
    print(f"D-002 recorded: {decisions.CONTROL_MATCHING.value}")
    print("\nDecision status:")
    print(decisions.summary())


if __name__ == "__main__":
    main()
