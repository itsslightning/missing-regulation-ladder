"""Every dashboard page must render without raising.

A Streamlit page that throws shows the user a traceback where a figure should
be, and the server still returns 200, so a health check passes while the app
is broken. These tests run each page headlessly and assert no exception, which
is the only cheap way to know the dashboard actually works.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "app"

PAGES = [
    "overview.py",
    "assay.py",
    "schema_genes.py",
    "causal.py",
    "methodology.py",
]


@pytest.fixture(autouse=True)
def _app_on_path():
    """Pages import `lib.data`, which resolves relative to app/."""
    sys.path.insert(0, str(APP))
    yield
    sys.path.remove(str(APP))


@pytest.mark.parametrize("page", PAGES)
def test_page_renders(page: str) -> None:
    at = AppTest.from_file(str(APP / "app_pages" / page), default_timeout=90)
    at.run()
    assert not at.exception, (
        f"{page} raised: "
        + "; ".join(str(e.value) for e in at.exception)
    )


def test_entry_point_defines_navigation() -> None:
    """The entry point should list every page, so none is orphaned."""
    src = (APP / "streamlit_app.py").read_text(encoding="utf-8")
    for page in PAGES:
        assert page in src, f"{page} is not registered in streamlit_app.py"


def test_arm_control_changes_the_headline() -> None:
    """The D-001 control is load-bearing, not decoration.

    D-001 was delegated rather than chosen by the project owner, and the
    rejected rules move the headline substantially. If the control stopped
    changing anything, the dashboard would be quietly asserting one answer.

    It is driven through the entry point because it lives in the sidebar: two
    pages read it, so it does not belong to any one of them.
    """
    at = AppTest.from_file(str(APP / "streamlit_app.py"), default_timeout=120)
    at.run()
    assert not at.exception
    assert at.radio, "no detection-rule control in the sidebar"

    before = [m.value for m in at.metric]
    at.radio[0].set_value("native").run()
    assert not at.exception
    after = [m.value for m in at.metric]
    assert before != after, "switching the detection rule changed nothing"


@pytest.mark.parametrize("page", PAGES)
def test_page_renders_on_a_sensitivity_arm(page: str) -> None:
    """Pages must render on the rejected detection rules too.

    test_page_renders covers the no-state case, which is what a deep link
    gives. This covers the other path: the sidebar control set away from the
    primary arm, where a page that filters to an arm with no rows would render
    an empty frame or raise on .iloc[0].
    """
    at = AppTest.from_file(str(APP / "app_pages" / page), default_timeout=90)
    at.session_state["arm"] = "native"
    at.run()
    assert not at.exception, (
        f"{page} raised on the native arm: "
        + "; ".join(str(e.value) for e in at.exception)
    )
