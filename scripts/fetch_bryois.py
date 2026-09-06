"""Fetch the Bryois record, resumably, in the order the analysis needs.

    uv run python scripts/fetch_bryois.py            # all 198 files
    uv run python scripts/fetch_bryois.py --chrom 8  # stop after chromosome 8

Replaces the inline shell loop this project started with, which had two flaws
that both cost real time:

  1. It skipped any file that existed and was non-empty, so a file truncated by
     an interrupted run was never repaired -- it just sat there looking done.
     This checks the byte size against the Zenodo manifest instead.
  2. It walked the record cell-type-major. Alphabetically the pseudobulk files
     sort last, so it would have fetched all eight cell types before the arm
     they are compared against. This goes chromosome-major, so every
     chromosome completes across all nine arms together and the D-007 contrast
     becomes valid -- if imprecise -- as early as possible.

Resumes with an HTTP Range request, so an interrupted 45 MB file costs only the
bytes still missing. Verifies each file against the manifest after writing and
reports anything that still does not match rather than moving on quietly.
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import requests

from pipeline.config import DIR_RAW
from pipeline.s05_detection import BRYOIS_ARMS

RECORD = "7276971"
FILE_URL = "https://zenodo.org/api/records/{record}/files/{name}/content"
DEST = DIR_RAW / "bryois"


def fetch_one(name: str, expected: int, timeout: int = 1800) -> tuple[bool, int]:
    """Download or resume one file. Returns (matches_expected, final_size)."""
    path = DEST / name
    have = path.stat().st_size if path.exists() else 0
    if have == expected:
        return True, have
    if have > expected:
        # Longer than upstream says: something is wrong enough that resuming
        # would compound it. Start over rather than append to a bad file.
        path.unlink()
        have = 0

    headers = {"Range": f"bytes={have}-"} if have else {}
    url = FILE_URL.format(record=RECORD, name=name)
    try:
        with requests.get(url, stream=True, timeout=timeout, headers=headers) as r:
            if headers and r.status_code == 200:
                # Server ignored the Range header and is sending the whole
                # file; truncate rather than appending to what is already here.
                have = 0
                path.unlink(missing_ok=True)
            elif headers and r.status_code != 206:
                return False, have
            r.raise_for_status()
            with path.open("ab" if have else "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
    except requests.RequestException as exc:
        print(f"    {name}: {type(exc).__name__}: {exc}")
        return False, path.stat().st_size if path.exists() else 0

    size = path.stat().st_size
    return size == expected, size


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", type=int, default=22, help="stop after this chromosome")
    ap.add_argument(
        "--workers",
        type=int,
        default=4,
        help="concurrent downloads. Zenodo throttles per connection rather "
        "than per client, so a few in parallel is several times faster than "
        "one. Kept deliberately small -- this is someone else's archive.",
    )
    args = ap.parse_args()

    manifest = json.loads((DEST / "manifest.json").read_text(encoding="utf-8"))

    # Chromosome-major, so each chromosome completes across all nine arms
    # together and the D-007 contrast becomes valid as early as possible.
    todo, done = [], 0
    for chrom in range(1, args.chrom + 1):
        for cell in BRYOIS_ARMS:
            name = f"{cell}.{chrom}.gz"
            expected = manifest.get(name)
            if expected is None:
                continue
            path = DEST / name
            if path.exists() and path.stat().st_size == expected:
                done += 1
                continue
            todo.append((chrom, name, expected))

    print(f"{done} already complete, {len(todo)} to fetch, "
          f"{args.workers} workers")

    fixed, failed = 0, []
    if todo:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(fetch_one, name, expected): (chrom, name, expected)
                for chrom, name, expected in todo
            }
            for fut in as_completed(futures):
                chrom, name, expected = futures[fut]
                try:
                    ok, size = fut.result()
                except Exception as exc:  # noqa: BLE001 - report, do not abort
                    ok, size = False, 0
                    print(f"  chr{chrom:<2} {name:<28} {type(exc).__name__}: {exc}")
                if ok:
                    fixed += 1
                    print(f"  chr{chrom:<2} {name:<28} ok {size:,}"
                          f"   ({fixed}/{len(todo)})")
                else:
                    failed.append(name)
                    print(f"  chr{chrom:<2} {name:<28} SHORT {size:,}/{expected:,}")

    print(f"\n{done} already complete, {fixed} fetched, {len(failed)} failed")
    if failed:
        print("failed: " + ", ".join(failed[:12]))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
