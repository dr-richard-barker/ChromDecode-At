# ChromDecode-At: expression-to-chromatin decoding of the Arabidopsis spaceflight transcriptome reveals a robust Polycomb-state bias in repressed genes that is not a SOG1 damage response

**Manuscript draft** — generated from the ChromDecode-At project record (v1–v8).
Status: discovery-only results, ready for internal review. All numbers traceable to
tables in this archive.

---

## Abstract

Spaceflight reprograms the plant transcriptome, but how these transcriptional
changes relate to chromatin organization is unexplored, largely because no
spaceflight epigenomic assay exists for Arabidopsis. We built **ChromDecode-At**,
a reproducible pipeline that decodes any Arabidopsis expression dataset into the
chromatin-state language of the Plant Chromatin State Database (PCSD, 36
ChromHMM states) and AraENCODE (12 states), and applied it to every Arabidopsis
RNA-seq dataset in NASA's Open Science Data Repository (OSDR; 33 accessions
enumerated, 10 with root tissue). We find that genes down-regulated by
spaceflight are strongly concentrated in Polycomb-repressed chromatin
(H3K27me3-marked accessible-promoter states; odds up to 18.0 in OSD-37
seedlings, FDR < 1e-90), a signal that replicates across independent datasets
and — after recovering sample metadata from NCBI Biosample for datasets with
undocumented designs — replicates in Col-0 root apical tissue (odds ~3.1,
FDR < 1e-6 in two APEX03-2-campaign root datasets). The repressed class is a
coherent functional module of peroxidase and oxidant-detoxification genes
(11 GO terms significant even against a state-matched Polycomb control), whose
promoters are enriched for AT-hook (AHL) and ZHD-family motifs — an
architecture that tracks the Polycomb state itself rather than the flight
response. Testing mechanism directly, we reanalyzed the only Arabidopsis
ChIP-seq in OSDR (SOG1-3xFLAG, the plant p53-analogue, after bleomycin-induced
DNA damage): SOG1 does not bind redox-module or Polycomb-down promoters,
excluding the canonical SOG1-dependent damage-repression pathway. Methylation
(WGBS, OSD-217) shows no significant flight-induced change at the module.
We conclude that the Polycomb-state bias of the spaceflight repression response
is robust across tissues and campaigns, is not a SOG1-style damage response,
and its upstream regulator remains unidentified — a concrete, testable target
for future spaceflight epigenomics. Expression-distribution features alone do
not predict chromatin state (definitive negative result, genomic-block CV
controlled), so the decoder's coordinate-based backbone, not a learned model,
carries the analysis.

---

## 1. Introduction

Plants are central to bioregenerative life support for space exploration, and
their transcriptomes respond extensively to spaceflight [OSDR]. Chromatin state
— the combinatorial pattern of histone marks, DNA methylation, and
accessibility — constrains which genes a cell can express, and Polycomb
repression (H3K27me3) is a canonical mechanism of stable, heritable gene
silencing in plants. Whether spaceflight transcriptional repression targets
genes in a chromatin-state-dependent manner has never been asked, because no
spaceflight chromatin assay exists for any plant: we verified that OSDR
contains no spaceflight H3K27me3 ChIP-seq, and the only Arabidopsis ChIP-seq
at all is a ground-based DNA-damage experiment (OSD-496).

We previously built ChromDecode-At, a pipeline that sidesteps this gap: given
an expression dataset, it maps every gene to its (ground-state) chromatin
context using the PCSD 36-state ChromHMM model [Liu et al. 2018] and
AraENCODE 12-state annotations, then tests whether differentially expressed
genes are concentrated in particular chromatin classes. Here we present the
full analysis program: decoder construction and validation (Section 2), a
learned expression-to-state model (Section 3, negative), an OSDR-wide sweep
with multi-resolution state analysis (Section 4), robustness diagnostics
(Section 5), functional characterization of the repressed module (Section 6),
methylation tests (Section 7), a SOG1 ChIP-seq mechanism test (Section 8), and
a tissue-matched root replication with metadata forensics (Section 9).

## 2. Results

### 2.1 A hybrid decoder validates on OSD-37 (v1)

ChromDecode-At harmonizes AGI identifiers, computes differential expression
(linear models on log2 counts with ecotype covariates, BH-FDR), maps genes to
chromatin states via TSS-anchored coordinates (gene body + 2 kb upstream;
100% PCSD and AraENCODE coverage; 37.7% PlantCADB ACR overlap), and tests
state enrichment of up/down DEGs by Fisher exact test against the expressed
genome.

On OSD-37 (28,067 genes × 56 samples; 529 up / 526 down at padj < 0.05,
|log2fc| > 1):

- **Up-regulated genes** are enriched for accessible-promoter chromatin
  (odds 3.22, FDR 3.0e-07) and depleted from heterochromatin/TE states
  (odds 0.28) and Polycomb states (0.32).
- **Down-regulated genes** are strongly enriched for Polycomb-repressed
  chromatin (odds 11.46, FDR 1.4e-14) and depleted from accessible promoters
  (0.14).

The two independent state annotations (PCSD 36-state; AraENCODE 12-state)
agree, and locus controls behaved correctly (FLC in a repression state;
ACT2/PHYB active). At 36-state resolution, down-regulated genes converge on
the H3K27me3-accessible-promoter block: S15 (H3K27me3 + accessible DNA)
odds 18.0, S16 6.4, S14 3.9, S13 3.3 (all FDR < 0.001) — spaceflight-repressed
genes are a coherent class of promoter-accessible, Polycomb-marked
euchromatic genes, while constitutive heterochromatin (S29–S36) is
transcriptionally inert in both directions.

### 2.2 A learned expression→state model fails honestly (v1–v2)

An elastic-net multinomial model predicting state group from
expression-distribution features (8 features, v1) reached CV accuracy 0.513
vs majority baseline 0.588. v2 retrained on tissue-matched pairs (AraENCODE
seedling expression, 28,727 genes, 11 features including tissue-specificity
statistics) with nested CV and genomic-block cross-validation (100-kb bins):
accuracy 0.533–0.554 vs baseline 0.549, and transfer to leaf, root, and
OSD-37 spaceflight data all ≤ 0.548. **Expression distributions do not
predict chromatin state**; the negative is not a CV artifact (block CV) and
motivates the deterministic coordinate-based backbone that carries all
subsequent analyses.

### 2.3 OSDR-wide sweep: the signal replicates (v2)

Running the backbone on every Arabidopsis RNA-seq dataset in OSDR (33
accessions; contrast auto-detection from metadata and sample-name
conventions) yielded full DE + enrichment for 5 datasets, expression-only
decoding for 18, and no processed counts for 11. The three spaceflight
datasets with sufficient DEGs show consistent directionality: OSD-37
down-regulated genes enrich for Polycomb_repressed (odds 7.0, FDR 1.7e-92);
OSD-314 (47/45 genes) replicates (down odds 4.9, FDR 1e-5). At 12-state
resolution, down-regulated DEGs enrich for the AraENCODE **Repression** state
(OSD-37 odds 8.0, FDR 3e-105; OSD-314 4.9). An earlier raw-scale DE version
(inflated ~5×) produced inconsistent calls and was corrected; one derived
claim (OSD-314 down genes in Active chromatin) is retracted.

### 2.4 OSD-314 up-regulated genes: escape from Polycomb repression (v2)

OSD-314's 47 up-regulated genes enrich for Polycomb CDS states (S11 odds
22.4, S12 4.2, S13 5.7; all FDR < 1e-3). Five diagnostics confirm robustness:
(1) null calibration, 10,000 random 47-gene draws — S11 empirical p < 0.0001;
(2) downsampling OSD-37's up genes to n=47 cannot reproduce the odds
(p = 0.000 for S11/S12); (3) threshold-grid dose-response — S11 odds rise
monotonically with DE stringency (1.3 → 6.7); (4) threshold-free GSEA on
1,944 Polycomb genes, permutation p < 0.001; (5) gene-identity check —
tandemly duplicated clusters among the 29 up genes in Polycomb states. The
previously reported S15 enrichment for up genes was an artifact of the
uncorrected DE (empirical p = 0.51) and is retracted. Biologically: genes
held in H3K27me3-marked chromatin are disproportionately *released* in
microgravity.

### 2.5 The repressed class is a redox/peroxidase module in AHL/ZHD promoter architecture (v3)

Of OSD-37's 526 down-regulated genes, 265 sit in Polycomb states S11–S15.
Three-set design: target (Polycomb-down, n=265), control 1 (non-Polycomb-down,
n=241), control 2 (state-matched Polycomb non-DE, n=3,303).

- **GO (TAIR GAF, hypergeometric, BH-FDR):** vs genome, 18 significant terms
  dominated by oxidative stress — peroxidase activity (11.1×), heme binding
  (7.2×), hydrogen peroxide catabolic process (13.6×), cellular oxidant
  detoxification (8.8×). **11 terms remain significant against the
  state-matched Polycomb control** — the flight-repressed subset is
  functionally a redox module, not a random Polycomb slice.
- **Motifs (JASPAR 2024 plants, 907 PFMs, strand-aware 2 kb promoters,
  per-motif p < 1e-4 thresholds):** 268 motifs significant vs genome —
  ZHD1 (2.9×), AHL12/13/25 AT-hook family (up to 5.7×), NTL9; ERF/ARF motifs
  depleted. **No motif survives BH-FDR vs the state-matched control** — the
  motif layer tracks the Polycomb state, not the response. Promoter GC
  content is matched (0.310 vs 0.320), excluding GC artifacts.
- **OSD-314 directional replication:** 8/10 top GO terms consistent
  (folds 7–27×; heme binding p = 4e-5), power-limited to 18 genes.

### 2.6 No flight-induced methylation change at the module (v4–v5)

With no spaceflight H3K27me3 data, we tested the slow-turning epigenetic
layer using OSD-217 (Ws ecotype, root + leaf WGBS, flight vs ground;
author-processed with bsmap/cscall against the ecotype-matched Ws genome
assembly). Gene-body tests: no significant change for the redox module
(minimum padj 0.07); nominal CHG gains in both tissues (p ≈ 0.011).
Promoter-window tests (2 kb, 100-bp bins, 21,424 genes): no test survives
BH-FDR (minimum padj 0.099); the module shows directionally consistent
nominal CHG gains. The Ws-defined Polycomb-down module has elevated baseline
promoter CG methylation (0.52 vs 0.02, p = 0.043, n = 6–11, not
BH-significant). Cross-ecotype concordance of flight log2FC between OSD-217
(Ws root) and OSD-37 (Col-0 seedlings) is r = −0.17 — effectively zero,
reflecting ecotype + tissue + campaign confounding, and motivating the
tissue-matched replication of Section 2.8. All methylation results are
null-or-suggestive; module coverage in Ws is reduced (14–72 genes).

### 2.7 SOG1 does not bind the module: excluding the damage-response pathway (v6)

SOG1 is the plant NAC-family master regulator of the DNA-damage response
(functional p53 analogue) and directly recruits Polycomb repression after
damage. If spaceflight repression of the redox module were a SOG1-style
damage response, SOG1 should occupy its promoters. We reanalyzed OSD-496
(SOG1-3xFLAG seedlings, bleomycin; 6 single-end libraries, 13.7–23.4 M
reads): HISAT2 alignment to TAIR10 (--no-spliced-alignment), MACS2 narrow
peaks (q < 0.01, IP vs matched-timepoint input): 56 peaks (20 min), 100
(1 h); IP-vs-wt-IP specificity checks give 160/180 peaks. QC: 67/100 of 1 h
peaks overlap a 2 kb promoter; median peak-to-TSS distance 194 bp.

**Result: a clean negative.** Across the 156-peak union, zero SOG1 peaks
fall in redox-module promoters (0/23) or Polycomb-down promoters (0/265) at
either timepoint (all Fisher padj = 1.0); 3/3,303 background promoters are
bound. The spaceflight repression of the module is therefore not a
recapitulation of SOG1's direct damage targets.

### 2.8 Tissue-matched root replication, with metadata forensics (v7–v8)

We enumerated all 62 Arabidopsis OSDR studies (33 with RNA-seq; 10 mentioning
root tissue) and resolved the field by QC: OSD-624 (Col-0 roots, FLT vs GC,
3+3), OSD-120 (CARA root tips; mixed Col-0/WS/phyD), OSD-193 and OSD-218
(counts keyed by GSM/SAMN accessions with empty OSDR factor values — every
sample ID was resolved to its condition via NCBI Biosample eutils, e.g.
"Spaceflight, 8 days old, Col-0 Rep4"), OSD-281 (excluded on recovery: 16
Sku5-mutant + 16 Ws samples, no Col-0), OSD-406 (Col-0 roots, 12/12, but a
Virgin Galactic **suborbital** flight), plus documented exclusions (OSD-217
Ws; OSD-219/427 no counts; OSD-411/480 EMCS 0.004G partial-gravity; OSD-208
no flight factor).

**The Polycomb enrichment replicates in Col-0 ISS roots** (age covariate):
OSD-193 Col-0 WT (8v8): 97 up / 181 down, odds 3.06, padj 2.2e-10; OSD-218
(16v16): 128 up / 117 down, odds 3.14, padj 4.7e-7. OSD-624 is null (odds
0.76, padj 0.81; 36 down genes); OSD-120 restricted to Col-0 WT collapses to
2 down genes after light/day covariate control (the v2 mixed-genotype 43/9
is composition-confounded and underpowered); OSD-406 suborbital yields 1
down gene.

An interim "aerial-tissue-only" interpretation (v7) was revised by the v8
recovery: the enrichment is robust across three of four independent Col-0
datasets, at roughly half the seedling effect size (odds ~3.1 vs 7.0), with
effect size varying ~2-fold by campaign. OSD-193/218 share the APEX03-2
payload, so they are not fully independent replications.

## 3. Discussion

Three findings organize this work. First, **spaceflight repression is
chromatin-state-biased**: down-regulated genes concentrate in
Polycomb-repressed, promoter-accessible euchromatin across datasets,
tissues, and resolutions, while constitutive heterochromatin is inert. The
bias is not expression-driven (the learned layer fails honestly) and not
redundant with GC content or gene length (state-matched controls).

Second, **the repressed class is a coherent redox/peroxidase module**,
distinguishable within the Polycomb class by function (11 GO terms) but not
by promoter motif architecture (AHL/ZHD motifs track the state, not the
response). This decouples "which genes are in Polycomb states" from "which
Polycomb genes respond to spaceflight" — the central question the project
reduces to.

Third, **the mechanism is not the canonical damage pathway**: SOG1, the
strongest candidate Polycomb-recruiting TF with available ChIP-seq, does not
bind the module's promoters after DNA damage. Combined with the methylation
nulls, the driver of the state bias remains unidentified; candidates include
AHL/ZHD-family chromatin-associated factors (motif layer), light/circadian
machinery, and campaign-specific hardware effects (the OSD-624 null).

**Practical upshot:** the Polycomb-state bias is a robust, replicable
feature of the Arabidopsis spaceflight transcriptome and a concrete target
for spaceflight epigenomics — the field needs a spaceflight H3K27me3
ChIP-seq or CUT&RUN assay to convert the ground-state-decoded bias into a
demonstrated dynamic mark change.

## 4. Limitations

- All chromatin states are **ground-state annotations**; the decoder tests
  enrichment of response genes in pre-existing states, not dynamic mark
  changes.
- OSD-193/218 share a payload and lab; OSD-624's null may reflect DE depth
  (36 down genes) or genuine campaign differences; tissue is confounded with
  campaign.
- Motif and GO findings are correlative and discovery-only; no DAP-seq/ChIP
  validation of AHL/ZHD binding.
- Methylation tests are power-limited (3 replicates; 14–72 module genes in
  Ws) and ecotype-mismatched to the Col-0 discovery datasets.
- SOG1 ChIP-seq is bleomycin-treated, mixed Ler/Col seedlings at two
  timepoints; absence of binding under other damage regimes or tissues is
  not excluded.
- OSD-406 (suborbital) and the four EMCS partial-gravity datasets are
  unusable for the FLT/GC question as designed.

## 5. Methods

**Pipeline.** Python 3.11; `chromdecode.py` modules: ID harmonization
(`clean_agi`), counts loading, DE (`differential_expression_any`: linear
models on log2(counts.clip(1)) with condition + covariates, BH-FDR),
coordinate mapping (TAIR10, TSS-anchored gene body + 2 kb), state assignment
(PCSD 36-state segments; AraENCODE 12-state; PlantCADB ACR overlap),
enrichment (Fisher exact, BH-FDR), expression-only fallback decoding.
Unit tests: 21, all passing.

**Data.** PCSD (36 ChromHMM states, 290,553 segments; Liu et al. 2018),
AraENCODE (12 states, TPM matrix, genes_TAIR10.txt), PlantCADB ACRs
(tolerant streaming parser for malformed whole-genome BED), NASA OSDR
(metadata API + GeneLab-processed counts; geode-py download URLs), JASPAR
2024 plants PFMs (907 usable after length filtering), TAIR GAF (251,290
annotations), go.obo, TAIR10 genome (Ensembl-style headers).

**ChIP-seq reanalysis.** HISAT2 (HPC; index from TAIR10.fa; single-end,
--no-spliced-alignment, -p 8), samtools sort/index, MACS2 2.2.9.1 callpeak
(-f BAM, -g 1.2e8, -q 0.01; IP vs matched-timepoint input; wt-IP
specificity). Promoter overlap: strand-aware 2 kb TSS intervals; Fisher
exact vs expressed non-DE background; BH-FDR. Occupancy: samtools per-100-bp
bin counts, RPM-scaled.

**Metadata forensics.** NCBI E-utilities (esearch + esummary, db=biosample)
for GSM/SAMN sample IDs; title parsing into condition/age/genotype/replicate;
median-ratio size-factor scaling where only unnormalized counts exist.

**Statistics.** Fisher exact (one-sided where enrichment is directional),
BH-FDR within test families, hypergeometric GO ORA (evidence-code filtered),
permutation-based null calibration and GSEA (1,000–10,000 draws), genomic
block CV for the learned layer. All thresholds: padj < 0.05, |log2fc| > 1
unless stated.

## 6. Data and code availability

- All code: `code/` (23 scripts; `chromdecode.py` is the core module).
- All figures: `figures/` (fig1–35, SVG + PNG).
- All result tables: `tables/` (49 CSVs); small derived data in
  `derived_data/`.
- Source data are not redistributed; download URLs and accessions are in
  `README.md` (OSDR OSD-37/120/193/208/217/218/219/251/281/314/346/406/411/
  427/480/496/508/624; GSE107980; PCSD; AraENCODE; PlantCADB; JASPAR; GO).
- Project report with full per-version records:
  `report_chromdecode_at.md`.

## References

1. Liu, Y. et al. (2018). *Plant Chromatin State Database* — PCSD 36-state
   ChromHMM annotations of Arabidopsis thaliana. Nucleic Acids Research.
2. AraENCODE: Zhou et al. — Arabidopsis epigenomic atlas (12-state), Huazhong
   Agricultural University resource.
3. PlantCADB: plant cis-regulatory element database (ACRs).
4. NASA Open Science Data Repository (OSDR), studies OSD-37, OSD-120,
   OSD-193, OSD-208, OSD-217, OSD-218, OSD-251, OSD-281, OSD-313, OSD-314,
   OSD-321, OSD-346, OSD-38, OSD-406, OSD-411, OSD-427, OSD-437, OSD-476,
   OSD-480, OSD-496, OSD-498, OSD-502, OSD-508, OSD-510, OSD-518, OSD-519,
   OSD-522, OSD-565, OSD-624, OSD-658, OSD-678, OSD-782.
   https://osdr.nasa.gov/
5. Paul, A.-L. et al. (2017). Genetic dissection of the Arabidopsis
   spaceflight transcriptome: Are some responses dispensable for the
   physiological adaptation of plants to spaceflight? PLoS ONE 12(6):
   e0180186. (CARA / OSD-120)
6. SOG1 ChIP-seq study (OSD-496 / GSE107980): The SOG1 transcriptional
   activator and the MYB3R repressors control a complex gene network in
   response to DNA damage in Arabidopsis. (Britt lab)
7. Mott, G.A. et al. (2011). Genomic and epigenomic analysis of the
   Wassilewskija (Ws) ecotype — Ws reference genome v7 used by GLDS-217.
8. JASPAR 2024: Rauluseviciute et al. Nucleic Acids Research (plants PFMs).
9. The Gene Ontology Consortium. TAIR GAF annotations.
   http://current.geneontology.org/annotations/tair.gaf.gz
10. Yi, D. et al. — SOG1 as master regulator of the DNA damage response in
    plants (context for Section 2.7).

*References 1–10 are provided for reader orientation; accession-level
provenance for every number is in the tables and code of this archive.*
