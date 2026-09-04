"""The decision gate is a safety mechanism, so it gets tested like one.

If these tests fail, the project's central claim -- that no conclusion-shaping
choice was made implicitly -- stops being true.
"""

from __future__ import annotations

import pytest

from pipeline import decisions
from pipeline.decisions import OpenDecision, UndecidedError


@pytest.fixture
def sample() -> OpenDecision:
    return OpenDecision(
        key="D-TEST",
        question="Does the gate hold?",
        alternatives={"a": "first option", "b": "second option"},
        why_it_matters="It is the whole point of the module.",
    )


def test_reading_an_undecided_value_raises(sample: OpenDecision) -> None:
    with pytest.raises(UndecidedError) as exc:
        _ = sample.value
    # The error has to be actionable on its own, since it may surface deep in a
    # stack trace with no other context.
    assert "D-TEST" in str(exc.value)
    assert "first option" in str(exc.value)


def test_deciding_makes_the_value_readable(sample: OpenDecision) -> None:
    sample.decide("a", rationale="because")
    assert sample.value == "a"
    assert sample.decided


def test_an_unlisted_option_is_refused(sample: OpenDecision) -> None:
    """Guards against a choice drifting to something nobody wrote down."""
    with pytest.raises(ValueError, match="written down and weighed"):
        sample.decide("c", rationale="improvised")


def test_redeciding_requires_force(sample: OpenDecision) -> None:
    """Changing a threshold mid-analysis must be deliberate and disclosed."""
    sample.decide("a", rationale="first call")
    with pytest.raises(ValueError, match="already set"):
        sample.decide("b", rationale="second thoughts")
    sample.decide("b", rationale="documented change", force=True)
    assert sample.value == "b"


def test_all_four_project_decisions_start_open() -> None:
    """A default silently introduced later would break this."""
    assert len(decisions.ALL_DECISIONS) == 4
    assert {d.key for d in decisions.outstanding()} == {
        "D-001",
        "D-002",
        "D-003",
        "D-004",
    }


def test_every_decision_offers_real_alternatives() -> None:
    """A decision with one option is not a decision."""
    for d in decisions.ALL_DECISIONS:
        assert len(d.alternatives) >= 2, d.key
        assert d.why_it_matters.strip()
        for name, meaning in d.alternatives.items():
            assert meaning.strip(), f"{d.key}:{name} has no stated meaning"
