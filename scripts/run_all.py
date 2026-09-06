"""Run the whole pipeline in order.

    uv run python scripts/run_all.py            # everything
    uv run python scripts/run_all.py --from s05 # resume from a stage
    uv run python scripts/run_all.py --list     # what would run

Stages are ordered by their sNN prefix and each is a module with a main().
Running them in one process rather than shelling out keeps the decision
sentinels loaded once, so a stage that would trip an unmade decision fails at
that stage rather than after the rest of the work is thrown away.

Stages that depend on data still downloading are skipped with a note rather
than failing the run -- rung 2 and the Bryois arm both behave this way, and a
partial pipeline is the normal state of this project until those land.
"""

from __future__ import annotations

import argparse
import importlib
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

#: In dependency order, which is also sNN order. The name is the module stem.
STAGES = [
    ("s00_record_downloads", "stamp every downloaded file with URL/date/sha256"),
    ("s01_gene_universe", "build the fixed 18,481-gene denominator"),
    ("s02_audit_rungs", "audit what each rung ships"),
    ("record_stage0_decisions", "record D-002"),
    ("record_stage1_decisions", "record D-001 and D-003"),
    ("s03_control_matching", "matched unconstrained control set"),
    ("s04_schema_gene_sets", "SCHEMA published + browser sets"),
    ("s05_detection", "per-gene, per-rung eQTL detection"),
    ("s06_recovery", "recovery curves, gaps, closure, CIs"),
    ("s07_figures", "static figures for the report"),
    ("s08_resolution_test", "within-SingleBrain resolution test"),
    ("s09_robustness", "robustness across constraint and correction choices"),
    ("s10_schema_switch", "which SCHEMA genes switch on, and where"),
    ("s11_assay_contrast", "assay vs donor count vs granularity"),
    ("record_stage2_decisions", "record D-004"),
    ("s12_smr_exposures", "eQTL side of SMR, harmonised on rsID"),
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="start", help="resume from this stage prefix")
    ap.add_argument("--list", action="store_true", help="list stages and exit")
    args = ap.parse_args()

    stages = STAGES
    if args.start:
        idx = next(
            (i for i, (n, _) in enumerate(stages) if n.startswith(args.start)), None
        )
        if idx is None:
            print(f"no stage starting with {args.start!r}")
            return 2
        stages = stages[idx:]

    if args.list:
        for name, what in stages:
            print(f"  {name:<26} {what}")
        return 0

    failed = []
    for name, what in stages:
        print(f"\n{'=' * 78}\n{name}  --  {what}\n{'=' * 78}")
        t0 = time.time()
        try:
            mod = importlib.import_module(f"pipeline.{name}")
            mod.main()
            print(f"\n[ok] {name} in {time.time() - t0:.1f}s")
        except Exception:
            traceback.print_exc()
            failed.append(name)
            print(f"\n[FAILED] {name}")

    print(f"\n{'=' * 78}")
    if failed:
        print(f"{len(failed)} stage(s) failed: {', '.join(failed)}")
        return 1
    print(f"all {len(stages)} stages completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
