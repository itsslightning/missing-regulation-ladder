"""Stage 0, step 0: stamp every downloaded file into logs/provenance.json.

Downloads during exploration happen with whatever tool is to hand -- curl, a
browser, a colleague's USB stick. That is fine, but it leaves no record, and
"which version of GTEx did this number come from" is unanswerable six weeks
later. This script walks the data directories, matches each file to its entry
in sources.py, and records URL, UTC date, byte count and sha256 for all of them.

Safe to re-run. A file whose hash already matches its record is skipped; a file
whose hash has *changed* raises, because that means a supposedly stable URL
served different bytes and every downstream number is now suspect.

Files it cannot attribute to a known source are listed rather than ignored --
an unattributed file in data/ is either a source nobody registered or a stray,
and both are worth knowing about before Stage 3 licence checks.
"""

from __future__ import annotations

from pathlib import Path

from pipeline import downloads, sources
from pipeline.config import DIR_RAW, DIR_RESTRICTED

#: Files that live directly in data/raw rather than in a per-source
#: subdirectory, mapped to the source they belong to.
#: Every registered source resolves to its own subdirectory through
#: Source.directory, so this map is only a fallback for files dropped straight
#: into data/raw by hand.
LOOSE_FILES: dict[str, sources.Source] = {
    "gnomad.v4.1.constraint_metrics.tsv": sources.GNOMAD_CONSTRAINT_V4,
}

TPM_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/"
    "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
)

#: Per-file URLs for multi-file records, reconstructed from the Zenodo API
#: pattern so each file's provenance points at the file, not just the record.
ZENODO_FILE_URL = "https://zenodo.org/api/records/{record}/files/{name}/content"
ZENODO_RECORDS = {"singlebrain": "14908182", "bryois": "7276971"}


def _url_for(source: sources.Source, path: Path) -> str:
    if path.name == "GTEx_v10_gene_median_tpm.gct.gz":
        return TPM_URL
    if source.key in ZENODO_RECORDS:
        return ZENODO_FILE_URL.format(
            record=ZENODO_RECORDS[source.key], name=path.name
        )
    if source.key == "gtex_v10":
        return source.url  # extracted from the archive; the archive is the source
    return source.url


def main() -> None:
    recorded, skipped, unattributed = 0, 0, []

    for base in (DIR_RAW, DIR_RESTRICTED):
        for path in sorted(base.rglob("*")):
            if not path.is_file() or path.name in (".gitkeep", "urls.txt"):
                continue
            if path.suffix == ".part":
                print(f"  incomplete download left behind: {path}")
                continue

            parent = path.parent.name
            source = sources.SOURCES_BY_KEY.get(parent) or LOOSE_FILES.get(path.name)
            if source is None:
                unattributed.append(path)
                continue

            if downloads.already_have(source, path):
                skipped += 1
                continue

            downloads.record(source, path, url=_url_for(source, path))
            recorded += 1

    print(f"recorded {recorded} file(s), {skipped} already stamped")

    if unattributed:
        print("\nfiles not attributable to a registered source:")
        for p in unattributed:
            print(f"  {p}")

    rows = downloads.summary_table()
    total = sum(r["bytes"] for r in rows)
    print(f"\ndownload log now covers {len(rows)} file(s), {total / 1e9:.2f} GB")

    restricted = downloads.restricted_files_present()
    if restricted:
        print("\nnon-redistributable files on disk (never commit these):")
        for r in restricted:
            print(f"  {r}")


if __name__ == "__main__":
    main()
