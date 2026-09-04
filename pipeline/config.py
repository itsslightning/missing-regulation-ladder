"""Paths, rung definitions, seeds and source provenance for the ladder pipeline.

Everything another script might otherwise hardcode lives here, so the methods
note, the dashboard and the code all quote the same numbers.

ON THE CONSTANTS THAT ARE DELIBERATELY ABSENT
---------------------------------------------
Four choices in this project can move the headline conclusion on their own:
what counts as a "detectable" eQTL, how tightly control genes are matched,
which multiple-testing correction is applied, and the colocalization priors.
Picking a "reasonable default" for any of them silently would mean the answer
was partly decided in this file rather than by the data.

They are therefore held in decisions.py as explicit, unset sentinels that raise
on use until they are set. DECISIONS.md records what was chosen and what the
alternatives were. If you are reading this because an import blew up, that is
the mechanism working.
"""

from __future__ import annotations

from pathlib import Path

# --- repository layout -----------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
DIR_RAW = DATA / "raw"
DIR_INTERIM = DATA / "interim"
DIR_PROCESSED = DATA / "processed"
DIR_LOGS = ROOT / "logs"
DIR_DOCS = ROOT / "docs"

# Anything whose licence forbids redistribution lands here and nowhere else.
# .gitignore excludes the whole directory rather than its contents, so a new
# file cannot be committed by accident just because nobody added a rule for it.
DIR_RESTRICTED = DATA / "restricted"

for _d in (DIR_RAW, DIR_INTERIM, DIR_PROCESSED, DIR_LOGS, DIR_RESTRICTED):
    _d.mkdir(parents=True, exist_ok=True)

# --- reproducibility -------------------------------------------------------

#: Every sampling, permutation or greedy-matching step seeds from this.
#: Changing it changes the control gene set, so it is pinned, not derived.
RANDOM_SEED = 20260903

# --- trait -----------------------------------------------------------------

TRAIT = "SCZ"
TRAIT_LABEL = "Schizophrenia"
GWAS_CITATION = "Trubetskoy et al. 2022, PGC3 schizophrenia (wave 3)"

# --- the resolution ladder -------------------------------------------------
# Frozen at four rungs per the project report. The ordering is the whole
# experiment: rung index is the x-axis of the recovery curve.
#
# n_samples is the donor count the study reports, carried here because the
# report names cross-study power as the main technical risk. Any recovery
# curve that is not also plotted against this number is not interpretable.

RUNGS: dict[str, dict] = {
    "gtex_cortex": {
        "index": 1,
        "label": "GTEx v10 cortex (bulk)",
        "resolution": "bulk, single tissue",
        "assay": "bulk RNA-seq",
        "n_samples": None,  # filled by s01 from the downloaded metadata
        "citation": "GTEx Consortium v10",
    },
    "bulk_brain": {
        "index": 2,
        "label": "PsychENCODE / MetaBrain (bulk brain)",
        "resolution": "bulk, meta-analysed brain",
        "assay": "bulk RNA-seq",
        "n_samples": None,
        "citation": "PsychENCODE / MetaBrain (de Klein et al. 2023)",
    },
    "sn_major": {
        "index": 3,
        "label": "SingleBrain major cell types",
        "resolution": "single-nucleus, 7 neocortical classes",
        "assay": "snRNA-seq",
        "n_samples": 983,
        "citation": "Jang et al. 2026, Nature Genetics 58(4):737-747",
    },
    "sn_subtype": {
        "index": 4,
        "label": "SingleBrain subtypes",
        "resolution": "single-nucleus, 28 subtypes",
        "assay": "snRNA-seq",
        "n_samples": 983,
        "citation": "Jang et al. 2026, Nature Genetics 58(4):737-747",
    },
}

#: Bryois is not a rung. It is a second, independent snRNA brain eQTL study at
#: a much smaller donor count, used as a replication check on rung 3 and as a
#: within-assay power contrast: if recovery tracks donor N rather than
#: resolution, Bryois and SingleBrain should separate despite similar
#: resolution. See DECISIONS.md D-004.
REPLICATION_SETS = ("bryois",)

# --- constraint stratification ---------------------------------------------

#: gnomAD LOEUF deciles are the conventional presentation. The report asks for
#: LOEUF/pLI stratification without fixing the binning, so both are carried:
#: deciles for the curve, and the pLI >= 0.9 / LOEUF < 0.35 cut points for the
#: headline "constrained" label, which is what SCHEMA and gnomAD themselves use.
LOEUF_N_BINS = 10
LOEUF_CONSTRAINED_MAX = 0.35
PLI_CONSTRAINED_MIN = 0.9

# --- output ----------------------------------------------------------------

RECOVERY_TABLE = DIR_PROCESSED / "recovery_by_rung.parquet"
GENE_SETS_TABLE = DIR_PROCESSED / "gene_sets.parquet"

#: Which bytes each run read: URL, UTC date, sha256, licence per file.
#: Written by downloads.py. Distinct from the per-stage logs/sNN_*.json files,
#: which record gene loss through the pipeline (provenance.py).
DOWNLOAD_LOG = DIR_LOGS / "downloads.json"
