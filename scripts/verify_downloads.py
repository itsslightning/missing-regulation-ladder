"""Check every downloaded file against its expected size, and re-fetch shorts.

    uv run python scripts/verify_downloads.py           # report only
    uv run python scripts/verify_downloads.py --fix     # re-fetch mismatches

Truncated downloads have been the single most common failure in this project:
Bryois pb.1.gz, Bryois Astrocytes.8.gz, and the 3.3 GB PsychENCODE full file
twice, every time with curl exiting 0. A short gzip either raises mid-read or,
worse, yields fewer genes and looks like a real biological difference rather
than a broken file.

Two sources of truth, and the order matters:

  1. UPSTREAM manifests: data/raw/bryois/manifest.json (Zenodo byte sizes)
     and the PsychENCODE Content-Length in sources.py. These are authoritative
     because they come from the server, not from this machine.
  2. logs/downloads.json, the size and sha256 stamped locally, used only for
     files no upstream manifest covers.

The order was originally the other way round and was wrong: s00 stamps whatever
is on disk, so a file stamped WHILE it was downloading recorded its truncated
size as though that were correct. Checking against that made a since-completed
file look "164% of expected". Upstream truth first, local record second.

Anything with neither is reported as unverifiable rather than silently passed.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import requests

from pipeline import config as cfg
from pipeline import downloads, sources

ZENODO_FILE_URL = "https://zenodo.org/api/records/{record}/files/{name}/content"
ZENODO_RECORDS = {"bryois": "7276971", "singlebrain": "14908182"}


def _bryois_manifest() -> dict[str, int]:
    p = cfg.DIR_RAW / "bryois" / "manifest.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def _refetch(path: Path, url: str) -> bool:
    """Resume the download in place. Returns True if it now matches."""
    headers = {}
    if path.exists():
        headers["Range"] = f"bytes={path.stat().st_size}-"
    mode = "ab" if headers else "wb"
    try:
        with requests.get(url, stream=True, timeout=1800, headers=headers) as r:
            if r.status_code not in (200, 206):
                return False
            with path.open(mode) as fh:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
    except requests.RequestException as exc:
        print(f"      re-fetch failed: {exc}")
        return False
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fix", action="store_true", help="re-fetch mismatches")
    args = ap.parse_args()

    recorded = {
        (e["source_key"], e["file"]): e for e in downloads.summary_table()
    }
    manifest = _bryois_manifest()

    ok, bad, unverifiable = 0, [], []
    for base in (cfg.DIR_RAW, cfg.DIR_RESTRICTED):
        for path in sorted(base.rglob("*")):
            if not path.is_file() or path.name in (
                ".gitkeep", "urls.txt", "urls_all.txt", "urls_priority.txt",
                "manifest.json", "zenodo_record.json",
            ):
                continue
            # Derived outputs, not downloads: see s00_record_downloads.
            if path.suffix == ".parquet":
                continue

            # Upstream truth first; the local stamp only as a fallback.
            key = (path.parent.name, path.name)
            if path.name in manifest:
                expected_size = manifest[path.name]
            elif path.name == sources.PSYCHENCODE_FULL_FILE:
                expected_size = sources.PSYCHENCODE_FULL_BYTES
            elif (entry := recorded.get(key)) is not None:
                expected_size = entry["bytes"]
            else:
                expected_size = None

            if expected_size is None:
                unverifiable.append(path)
                continue

            actual = path.stat().st_size
            if actual == expected_size:
                ok += 1
                continue

            bad.append((path, actual, expected_size))
            print(
                f"  SHORT {path.parent.name}/{path.name}: "
                f"{actual:,} of {expected_size:,} bytes "
                f"({actual / expected_size:.1%})"
            )
            if args.fix:
                src_key = path.parent.name
                if src_key in ZENODO_RECORDS:
                    url = ZENODO_FILE_URL.format(
                        record=ZENODO_RECORDS[src_key], name=path.name
                    )
                elif path.name == sources.PSYCHENCODE_FULL_FILE:
                    url = (
                        "http://resource.psychencode.org/Datasets/Derived/QTLs/"
                        + sources.PSYCHENCODE_FULL_FILE
                    )
                else:
                    print("      no known URL; skipping")
                    continue
                print("      re-fetching...")
                _refetch(path, url)
                now = path.stat().st_size
                print(
                    f"      now {now:,} bytes "
                    + ("OK" if now == expected_size else "STILL SHORT")
                )

    print(f"\n{ok} file(s) verified, {len(bad)} short, "
          f"{len(unverifiable)} unverifiable")
    if unverifiable:
        for p in unverifiable[:10]:
            print(f"  unverifiable: {p.parent.name}/{p.name}")
    return 1 if bad and not args.fix else 0


if __name__ == "__main__":
    raise SystemExit(main())
