# Stage 2 report — causal links and druggability

Date: 2026-09-06. Decisions: [`../DECISIONS.md`](../DECISIONS.md) D-004, D-014.
Tables: `data/processed/{smr_results,locus_explanation,druggability}.parquet`.

---

## 1. Headline

Of **184** genome-wide significant SCZ loci (PGC3 European, MHC excluded,
250 kb clumping), the share with a causal eQTL explanation by SMR:

| Rung | Assay | Donors | Loci explained | Testable genes | SMR-significant |
|---|---|---|---|---|---|
| GTEx cortex | bulk tissue | 205 | 19 (10.3%) | 2,648 | 218 |
| PsychENCODE | bulk tissue | 1,387 | 27 (14.7%) | 4,553 | 316 |
| **SingleBrain, 7 types** | **snRNA-seq** | **983** | **59 (32.1%)** | 7,837 | 574 |
| SingleBrain, 28 subtypes | snRNA-seq | 983 | 56 (30.4%) | 6,820 | 521 |

> **41 loci gain a causal eQTL explanation only in single-nucleus data.**
> 5 are explained in bulk but not single-nucleus. Bulk explains 35 in total;
> single-nucleus explains 71.

Single-nucleus roughly **doubles** the share of schizophrenia loci with a
causal expression explanation, against the best bulk study — while using 29%
*fewer* donors.

## 2. This refines Stage 1 rather than merely repeating it

Stage 1 found that the constrained-gene *gap* is set by assay and is completely
flat with donor count (0.367 at N=205 and at N=1,387). Stage 2's outcome is a
different quantity — absolute locus counts — and there donor count **does**
help:

| Comparison | Donor change | Loci explained |
|---|---|---|
| GTEx → PsychENCODE (within bulk) | ×6.8 | 19 → 27 (**+42%**) |
| PsychENCODE → SingleBrain (across assay) | ×0.71 | 27 → 59 (**+119%**) |

Both are real and they are not in conflict. More donors raise the overall
detection level, so more genes clear the instrument threshold — testable genes
go 2,648 → 4,553 within bulk. What donors do **not** do is close the
constrained-gene deficit specifically. The assay does both.

The SMR significance *rate* is nearly constant across rungs (6.9%–8.2%), which
locates the gain precisely: single-nucleus data does not make eQTLs more likely
to be causal, it makes **more genes testable at all**.

## 3. Validation

The pipeline was checked against known biology before any of the above was
believed.

**C4A is the top-ranked non-MHC-artefact gene, b_xy = +0.189** — higher
expression raises risk. That is the best-established causal direction in
schizophrenia genetics (Sekar et al. 2016) and the analysis recovers it with
the correct sign. CACNA1C, FURIN, KCTD13, ASPHD1, BTN3A2, GATAD2A and GLT8D1
also appear in the top 15.

**Sign balance across significant results is 0.501**, so there is no systematic
orientation error — the thing most likely to be silently wrong in an SMR
analysis.

**A near-universal allele flip was investigated, not assumed.** Harmonisation
reported 57,044 flipped against 193 same, where arbitrary orientation would
give roughly 50/50. The cause is a convention difference: PGC3 defines A1 as
*"SNP reference allele for beta"* while GTEx and SingleBrain report effects
against the **alternate** allele. Because harmonisation matches on actual
allele letters rather than an assumed convention, the flip rate is confirmation
that it is working.

## 4. Druggability of the newly-visible genes

Splitting SMR-significant genes by where the causal link is visible, then
attaching Open Targets tractability:

| Group | Genes | Small-molecule tractable | Phase 1 or beyond | SCHEMA |
|---|---|---|---|---|
| **single-nucleus only** | **516** | **175 (33.9%)** | **24** | **4** |
| bulk and single-nucleus | 198 | 64 (32.3%) | 7 | 0 |
| bulk only | 239 | 74 (31.0%) | 12 | 0 |

Tractability rates are essentially identical across the three groups (~32–34%),
which is the right sanity check: single-nucleus data is not enriching for
druggable genes, it is finding *more* genes at the same druggable rate.

**24 clinically-advanced targets are causally implicated only in single-nucleus
data**, including the calcium-channel family that has been a schizophrenia
target class for a decade:

> **CACNA1C, CACNA1D, CACNA1I, CACNB2**, GRIA4, SCN1A, SCN3A, KCNB1, KCNQ5,
> ERBB4, FGFR1, CRHR1, APP, DPP4, CFTR, and others.

A bulk-tissue eQTL screen would not have surfaced these as causal.

## 5. The CHRM4 callback — an informative negative

The brief asked to flag any connection back to CHRM4 and the muscarinic
mechanism. The honest answer is that **there is none, and the reason is the
point of the project.**

| Gene | Detectable eQTL (rungs) | SMR-testable | SMR-significant |
|---|---|---|---|
| CHRM1 | 0 | 29 | 0 |
| CHRM2 | 2 | 25 | 0 |
| CHRM3 | 2 | 37 | 0 |
| **CHRM4** | **0** | **2** | **0** |
| CHRM5 | 6 | 36 | 0 |
| SLC5A7 | 1 | 1 | 1 |

No muscarinic receptor reaches SMR significance. Only SLC5A7, the choline
transporter, does — in bulk only.

**CHRM4 is this project's thesis in a single gene.** It is:

- **in the constrained case set of this very study** — LOEUF 0.265, decile 0
  (the most constrained), pLI 0.974;
- **well expressed in cortex** (11.14 TPM), so its invisibility is not a
  detection failure from low expression;
- **an approved drug target** — the M1/M4 agonist xanomeline-trospium
  (Cobenfy) was approved in 2024;
- **and it has no detectable cis-eQTL in brain at any assay or resolution
  tested here.**

A genetics-led target discovery pipeline built on eQTLs — including the sibling
[scz-target-prioritization](https://github.com/itsslightning/scz-target-prioritization)
repo, which surfaced CHRM4 by a route that did not depend on eQTLs — would
never have found it through regulatory evidence. CHRM4 lives in the residual
~40% of the constrained-gene gap that **no assay and no resolution closes**.

That is the practical cost of missing regulation, and it argues that eQTL-based
target discovery has a systematic blind spot precisely at the dosage-sensitive
genes most likely to matter.

## 6. Caveats

- **SMR is not colocalization.** It cannot separate a shared causal variant
  from linkage between two distinct causal variants, so absolute counts are
  inflated relative to a true coloc. D-004 accepts this deliberately: the bias
  has the same construction at every rung and Stage 2 asks a cross-rung
  question. **These counts should not be quoted against published coloc
  figures.** The HEIDI and coloc sensitivity arms on the two arms with regional
  data are not yet run.
- **Locus definition is distance-based**, not LD-based, because an
  ancestry-matched LD reference is not in hand. The same definition applies to
  every rung, so the ratio is what this measures — not the absolute locus count.
- **The MHC is excluded** as a single megabase-spanning LD block that distance
  clumping would collapse into one locus holding hundreds of genes.
- **rung 2 uses a reconstructed standard error** (from beta and nominal p),
  exact for a Wald test but a reconstruction nonetheless.
- **Ancestry is matched but not perfectly** — PGC3 European against
  European-ancestry eQTLs (D-014), which is the right pairing for SMR's LD
  assumption, but GTEx is only predominantly European.
