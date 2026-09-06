"""Every external dataset the ladder uses: where it comes from, and what its
licence permits.

Two things are recorded here that are usually left implicit and later cost a
day each.

The first is redistribution. `redistributable` is not "is this data public" --
all of it is public -- it is "may this repository host a copy". CC-BY-4.0 says
yes with attribution; a registration form or a data-use agreement says no,
regardless of how freely the file downloads. Sources marked False are written
to data/restricted/, which .gitignore excludes as a whole directory, and only
derived summaries computed from them may be published. Stage 3 re-checks this
list before anything ships.

The second is that eQTL catalogues are re-released. GTEx v10 is not v8,
SingleBrain will get follow-up versions, and the SCHEMA browser is currently
serving a release newer than the published paper. A recovery curve computed
against one snapshot is not reproducible against another, so every download is
stamped with its date and sha256 by provenance.py and the version string is
pinned here rather than discovered at runtime.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from pipeline.config import DIR_RAW, DIR_RESTRICTED


@dataclass(frozen=True)
class Source:
    """One downloadable dataset.

    `rung` is the ladder rung this source supplies, or None for the gene-level
    annotation sources (constraint, gene sets) that all rungs share.
    """

    key: str
    name: str
    url: str
    citation: str
    licence: str
    redistributable: bool
    rung: str | None = None
    version: str | None = None
    notes: str = ""
    #: Multi-file sources (Zenodo records, tar archives) resolve their own file
    #: lists at download time; a single URL is not enough to describe them.
    multifile: bool = False
    filenames: tuple[str, ...] = field(default_factory=tuple)

    @property
    def directory(self):
        """Where a copy of this source is allowed to live."""
        base = DIR_RAW if self.redistributable else DIR_RESTRICTED
        return base / self.key


# ---------------------------------------------------------------------------
# Gene-level annotation, shared by every rung
# ---------------------------------------------------------------------------

GNOMAD_CONSTRAINT = Source(
    key="gnomad_constraint",
    name="gnomAD v2.1.1 LoF constraint metrics by gene",
    url=(
        "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/"
        "constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
    ),
    citation="Karczewski et al. 2020, Nature 581:434-443 (gnomAD v2.1.1)",
    licence="Free of restrictions (gnomAD terms); citation requested",
    redistributable=True,
    version="2.1.1",
    notes=(
        "v2.1.1 rather than v4.1 because LOEUF (oe_lof_upper) from v2.1.1 is the "
        "vintage Mostafavi et al. 2023 and the SCHEMA papers used, and the "
        "decile field oe_lof_upper_bin ships with it. It also carries the exact "
        "covariates the control matching needs -- num_coding_exons, cds_length, "
        "gene_length, brain_expression -- in one table. v4.1 is downloaded "
        "alongside it as a robustness check, not as the primary."
    ),
)

GNOMAD_CONSTRAINT_V4 = Source(
    key="gnomad_constraint_v4",
    name="gnomAD v4.1 constraint metrics",
    url=(
        "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/"
        "constraint/gnomad.v4.1.constraint_metrics.tsv"
    ),
    citation="Chen et al. 2024, Nature 625:92-100 (gnomAD v4.1)",
    licence="Free of restrictions (gnomAD terms); citation requested",
    redistributable=True,
    version="4.1",
    notes=(
        "Robustness check only. v4.1 recomputes LOEUF on ~730k exomes, so gene "
        "rankings shift; if the recovery curve is qualitatively different under "
        "v4.1 that is a finding worth reporting, not a bug to hide."
    ),
)

SCHEMA_GENES = Source(
    key="schema",
    name="SCHEMA exome meta-analysis gene results",
    url=(
        "https://storage.googleapis.com/exome-results-browsers-public/downloads/"
        "2026-08-21/SCHEMA/SCHEMA_gene_results.tsv.bgz"
    ),
    citation=(
        "Singh et al. 2022, Nature 604:509-516; browser release 2026-08-21 "
        "(schema.broadinstitute.org)"
    ),
    licence="SCHEMA browser terms; summary results freely available",
    redistributable=True,
    version="browser-2026-08-21",
    notes=(
        "WARNING, this is a live decision, see the Stage 0 report. The browser "
        "release carries 87,959 cases / 150,587 controls, roughly 3.6x the "
        "24,248 cases of the published Singh et al. 2022 analysis, and it ships "
        "no q-value column -- only per-gene p-values. Using it means the gene "
        "set is larger and better powered but not the one in the citable paper; "
        "using Singh 2022's supplementary table means the published set but "
        "fewer genes. Do not treat this file as 'the SCHEMA gene set' without "
        "settling that."
    ),
)

SCHEMA_PUBLISHED = Source(
    key="schema_published",
    name="SCHEMA published gene results (Singh et al. 2022, Supplementary Table 5)",
    url=(
        "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-04556-w/"
        "MediaObjects/41586_2022_4556_MOESM3_ESM.xlsx"
    ),
    citation="Singh et al. 2022, Nature 604:509-516",
    licence="Springer Nature supplementary material; check terms before rehosting",
    redistributable=False,
    version="Supplementary Table 5",
    notes=(
        "The PRIMARY SCHEMA gene set (D-009). 18,324 genes with `P meta` and, "
        "crucially, `Q meta` -- the FDR the paper's own 10-gene and 32-gene "
        "thresholds are defined on. The browser release ships no q-value, so "
        "using it would mean inventing an FDR; this table does not force that "
        "choice, which is why it is primary.\n\n"
        "Keyed by gene SYMBOL, not Ensembl ID, so it needs a symbol -> ENSG hop "
        "through gnomAD. Held non-redistributable: publisher supplementary "
        "files are not covered by the CC licences the Zenodo records carry."
    ),
)

# ---------------------------------------------------------------------------
# Rung 1 -- bulk, single tissue
# ---------------------------------------------------------------------------

GTEX_V10 = Source(
    key="gtex_v10",
    name="GTEx v10 single-tissue cis-eQTL eGenes (brain tissues)",
    url=(
        "https://storage.googleapis.com/adult-gtex/bulk-qtl/v10/"
        "single-tissue-cis-qtl/GTEx_Analysis_v10_eQTL.tar"
    ),
    citation="GTEx Consortium, v10 release",
    licence="GTEx open-access summary statistics; no restriction on reuse",
    redistributable=True,
    rung="gtex_cortex",
    version="v10",
    multifile=True,
    notes=(
        "Despite the .tar extension the archive is gzip-compressed, so it needs "
        "tar -xz. Only *Brain*.egenes.txt.gz members are extracted. The egenes "
        "files list every gene tested, not only significant ones, which is what "
        "makes a denominator possible -- recovery fraction needs the tested set, "
        "not just the hits."
    ),
)

# ---------------------------------------------------------------------------
# Rung 2 -- bulk brain, meta-analysed
# ---------------------------------------------------------------------------

PSYCHENCODE = Source(
    key="psychencode",
    name="PsychENCODE DER-08 adult prefrontal cortex cis-eQTLs",
    url=(
        "http://resource.psychencode.org/Datasets/Derived/QTLs/"
        "DER-08b_hg38_eQTL.bonferroni.txt"
    ),
    citation="Wang et al. 2018, Science 362:eaat8464 (PsychENCODE)",
    licence="PsychENCODE open resource; check terms before redistributing",
    redistributable=False,
    rung="bulk_brain",
    version="DER-08b (hg38)",
    notes=(
        "Held as restricted until the resource.psychencode.org terms page is "
        "read and quoted in DECISIONS.md -- the data downloads without a click "
        "through, which is not the same as a licence to rehost. Derived counts "
        "are publishable either way."
    ),
)

#: Exact byte count of PsychENCODE's full association file, from the server's
#: Content-Length. A size check rather than a ">3 GB" heuristic because that
#: file has now truncated twice mid-transfer with curl still exiting 0, and a
#: partial gzip would silently produce a wrong rung 2 rather than an error.
PSYCHENCODE_FULL_BYTES = 3_293_218_507
PSYCHENCODE_FULL_FILE = "Full_hg19_cis-eQTL.txt.gz"

METABRAIN = Source(
    key="metabrain",
    name="MetaBrain cortex cis-eQTLs (EUR)",
    url="https://www.metabrain.nl/cis-eqtls.html",
    citation="de Klein et al. 2023, Nature Genetics 55:377-388 (MetaBrain)",
    licence="Access via registration form; redistribution not granted",
    redistributable=False,
    rung="bulk_brain",
    version="2021-07-23",
    notes=(
        "BLOCKED, needs a human. download.metabrain.nl serves no public file "
        "index; access is gated behind a Google Form the maintainers hand out "
        "links through. Nothing here can be automated. PsychENCODE covers rung 2 "
        "on its own if MetaBrain access does not arrive in time."
    ),
)

# ---------------------------------------------------------------------------
# Rungs 3 and 4 -- single-nucleus, major cell types and subtypes
# ---------------------------------------------------------------------------

SINGLEBRAIN = Source(
    key="singlebrain",
    name="SingleBrain single-nucleus cis-eQTL meta-analysis (top associations)",
    url="https://zenodo.org/api/records/14908182",
    citation=(
        "Jang et al. 2026, Nature Genetics 58(4):737-747; Zenodo 10.5281/"
        "zenodo.14908182"
    ),
    licence="CC-BY-4.0",
    redistributable=True,
    rung="sn_major",
    version="Zenodo 14908182 (2025-02-28)",
    multifile=True,
    notes=(
        "Only the 36 *_top_assoc.tsv.gz files are pulled: ~1.6 MB each, one row "
        "per gene with the best variant and a q-value, which is exactly the "
        "recovery-curve input. The matching *_full_assoc.tsv.gz files are 4-10 "
        "GB each (~150 GB for the set) and are needed only for Stage 2 "
        "colocalization; they are fetched per gene of interest, never in bulk. "
        "Phenotypes are on GENCODE v38, which does not match GTEx v10 -- see the "
        "gene-ID harmonisation note in the Stage 0 report."
    ),
)

# ---------------------------------------------------------------------------
# Replication / power contrast -- deliberately not a rung
# ---------------------------------------------------------------------------

BRYOIS = Source(
    key="bryois",
    name="Bryois brain cell-type cis-eQTLs, 8 cell types plus pseudobulk",
    url="https://zenodo.org/api/records/7276971",
    citation=(
        "Bryois et al. 2022, Nature Neuroscience 25:1104-1112; Zenodo 10.5281/"
        "zenodo.7276971"
    ),
    licence="CC-BY-4.0",
    redistributable=True,
    rung=None,
    version="Zenodo 7276971 (Feb 2023 update)",
    multifile=True,
    notes=(
        "The most useful thing in this project that the report did not ask for. "
        "The Feb-2023 update added pb[1-22].gz, a 'tissue-like' pseudobulk "
        "analysis aggregating reads across all nuclei per individual -- the same "
        "donors, pipeline and normalisation as the 8 cell-type files, differing "
        "only in resolution. That is a within-study bulk-vs-cell-type contrast "
        "at fixed N, which is the one comparison in this whole design where "
        "resolution is not confounded with sample size. See the Stage 0 report. "
        "Caveat: nominal p-values only, no per-gene permutation q-values, so it "
        "cannot be scored the way GTEx and SingleBrain are without recomputing a "
        "correction locally (decision D-001)."
    ),
)

# ---------------------------------------------------------------------------
# Stage 2 only
# ---------------------------------------------------------------------------

PGC3_SCZ = Source(
    key="pgc3_scz",
    name="PGC3 schizophrenia GWAS summary statistics (wave 3, European)",
    url="https://ndownloader.figshare.com/files/34517828",
    citation=(
        "Trubetskoy et al. 2022, Nature 604:502-508 (PGC3 SCZ); "
        "figshare 10.6084/m9.figshare.19426775"
    ),
    licence="CC-BY-4.0",
    redistributable=True,
    version="wave3 v3 public release, European autosomes",
    notes=(
        "This is the PGC's own DESIGNATED PUBLIC RELEASE, not the "
        "agreement-gated distribution. The figshare record reports "
        "is_public=True, is_embargoed=False and CC-BY-4.0, and the files are "
        "named '.public.v3' by the depositors. So no data-use agreement is "
        "accepted to obtain it and it may be redistributed with attribution -- "
        "which is why this entry is redistributable=True where an earlier "
        "version of this file assumed otherwise.\n\n"
        "The EUROPEAN subset is used rather than the larger core or primary "
        "releases because SMR assumes the exposure and outcome samples share "
        "an LD structure, and the eQTL arms are European: SingleBrain is "
        "European-ancestry only, Bryois is European, GTEx is predominantly so. "
        "Pairing a multi-ancestry GWAS with European eQTLs would violate that "
        "assumption. See DECISIONS.md D-014.\n\n"
        "Despite the permissive licence the 240 MB file is not committed: "
        "data/raw is gitignored and only derived SMR results are published."
    ),
)


ALL_SOURCES: tuple[Source, ...] = (
    GNOMAD_CONSTRAINT,
    GNOMAD_CONSTRAINT_V4,
    SCHEMA_GENES,
    SCHEMA_PUBLISHED,
    GTEX_V10,
    PSYCHENCODE,
    METABRAIN,
    SINGLEBRAIN,
    BRYOIS,
    PGC3_SCZ,
)

SOURCES_BY_KEY = {s.key: s for s in ALL_SOURCES}


def restricted_sources() -> tuple[Source, ...]:
    """Sources that may not be rehosted -- the Stage 3 pre-publication check."""
    return tuple(s for s in ALL_SOURCES if not s.redistributable)
