"""Download files and record exactly which bytes were used.

Named `downloads` rather than `provenance` deliberately: in the sibling
scz-target-prioritization repo `provenance.py` means gene-loss accounting, and
that module is carried over here unchanged. This one answers a different
question, not "how many genes survived this step" but "which version of the
upstream file did this run read".

eQTL catalogues are living resources. GTEx has had ten releases, the SCHEMA
browser is currently serving a cohort roughly 3.6x the published one, and
SingleBrain will get follow-ups. "I downloaded GTEx brain eQTLs" is therefore
not a reproducible statement; "sha256 3f9a... fetched 2026-09-04 from <url>" is.

Every fetch writes one entry to logs/downloads.json holding the URL, the UTC
download date, the byte count, the sha256 and the licence recorded in
sources.py. The methods note and the dashboard's data-provenance panel are both
generated from that file, so a stale download cannot be described as a fresh one.

Re-running is cheap and safe: a file already on disk whose sha256 matches its
recorded entry is left alone and not re-fetched. A file on disk whose hash does
*not* match the record is an error rather than an overwrite, because that means the
upstream file changed under a stable URL, which is exactly the failure this
module exists to catch.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import requests

from pipeline.config import DOWNLOAD_LOG, ROOT
from pipeline.sources import Source

#: Streamed in chunks so a 2.5 GB archive never lands in memory.
_CHUNK = 1 << 20


class DownloadError(RuntimeError):
    """An upstream file changed, or a recorded file has gone missing."""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def _load() -> dict:
    if DOWNLOAD_LOG.exists():
        return json.loads(DOWNLOAD_LOG.read_text(encoding="utf-8"))
    return {}


def _save(log: dict) -> None:
    DOWNLOAD_LOG.parent.mkdir(parents=True, exist_ok=True)
    DOWNLOAD_LOG.write_text(
        json.dumps(log, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )


def record(source: Source, path: Path, *, url: str | None = None) -> dict:
    """Stamp a file that is already on disk, and return its record entry.

    Used both by `fetch` and by the manual-acquisition paths. MetaBrain and
    PGC3 arrive through a form and a data-use agreement respectively, so they
    are copied in by hand but still have to be recorded the same way.
    """
    digest = sha256(path)
    entry = {
        "source_key": source.key,
        "source_name": source.name,
        "file": path.name,
        # Stored relative to the repo root so the log is identical on any
        # machine; absolute paths would make two clones look like they used
        # different data.
        "relative_path": path.relative_to(ROOT).as_posix(),
        "url": url or source.url,
        "downloaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "bytes": path.stat().st_size,
        "sha256": digest,
        "licence": source.licence,
        "redistributable": source.redistributable,
        "version": source.version,
        "citation": source.citation,
    }
    log = _load()
    log.setdefault(source.key, {})[path.name] = entry
    _save(log)
    return entry


def already_have(source: Source, path: Path, *, deep: bool = False) -> bool:
    """True if this exact file was fetched before and has not changed on disk.

    A hash mismatch raises rather than returning False: silently re-downloading
    over a changed file would erase the evidence that it changed.

    By default the check is by byte size, which is O(1) per file. Hashing every
    recorded file on every run is O(total bytes), and this project's data
    directory reaches several gigabytes once the Bryois record is complete --
    enough to make a routine re-stamp take minutes. Size catches truncation,
    which is the failure mode that has actually occurred here, repeatedly.

    `deep=True` forces the sha256 comparison, which additionally catches an
    upstream edit that preserved the byte count. scripts/verify_downloads.py is
    the place that wants it.
    """
    if not path.exists():
        return False
    entry = _load().get(source.key, {}).get(path.name)
    if entry is None:
        return False
    if path.stat().st_size != entry["bytes"]:
        raise DownloadError(
            f"{path} is {path.stat().st_size:,} bytes but was recorded as "
            f"{entry['bytes']:,} on {entry['downloaded_utc']}. Either the file "
            "is mid-download, was edited locally, or the upstream release "
            "changed. Resolve deliberately rather than re-stamping over it."
        )
    if not deep:
        return True
    if sha256(path) != entry["sha256"]:
        raise DownloadError(
            f"{path} is on disk but its sha256 does not match the one recorded on "
            f"{entry['downloaded_utc']}. Either the file was edited locally or the "
            "upstream release changed under the same URL. Resolve deliberately: "
            "delete the file to re-fetch, and note the version change in "
            "DECISIONS.md."
        )
    return True


def fetch(
    source: Source,
    *,
    url: str | None = None,
    filename: str | None = None,
    timeout: int = 600,
) -> Path:
    """Download one file for `source` into the directory its licence permits.

    Restricted sources land in data/restricted/ purely by virtue of
    `Source.directory`, so a licence change in sources.py moves the file
    without any caller needing to know.
    """
    url = url or source.url
    filename = filename or url.rstrip("/").split("/")[-1]
    target_dir = source.directory
    target_dir.mkdir(parents=True, exist_ok=True)
    path = target_dir / filename

    if already_have(source, path):
        return path

    # .part first, renamed on success: an interrupted download must never be
    # mistaken for a complete one by a later run.
    part = path.with_suffix(path.suffix + ".part")
    with requests.get(url, stream=True, timeout=timeout) as resp:
        resp.raise_for_status()
        with part.open("wb") as fh:
            for chunk in resp.iter_content(chunk_size=_CHUNK):
                fh.write(chunk)
    part.replace(path)

    record(source, path, url=url)
    return path


def summary_table() -> list[dict]:
    """Flat list of every recorded file, for the methods note and dashboard."""
    rows = []
    for files in _load().values():
        rows.extend(files.values())
    return sorted(rows, key=lambda r: (r["source_key"], r["file"]))


def restricted_files_present() -> list[str]:
    """Recorded files that must not be published. Checked before Stage 3 ships."""
    return [
        f"{r['source_key']}/{r['file']}"
        for r in summary_table()
        if not r["redistributable"]
    ]
