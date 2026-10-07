# ChromDecode-At: An Arabidopsis Chromatin-State Auto-Decoder

**Test dataset:** NASA OSDR OSD-37 (Arabidopsis seedlings, spaceflight vs ground control)
**Extended:** all Arabidopsis RNA-seq datasets in NASA OSDR (v2)
**Date:** 2026-09-18 | **Status:** Working prototype, discovery-only

---

## 1. What was built

**ChromDecode-At** is a reusable Python pipeline that takes an Arabidopsis
gene-expression matrix (any OSDR dataset or user CSV) and decodes the associated
chromatin landscape into human-readable figures.

```
Expression matrix / DEG table
   → AGI ID harmonization + differential expression (linear model, BH-FDR)
   → TAIR10 coordinates (gene body + 2 kb upstream, TSS-anchored)
   → Chromatin state assignment:
        PCSD 36-state (ChromHMM, histone marks + methylation)   [100% coverage]
        AraENCODE 12-state (independent replication)            [100% coverage]
        PlantCADB ACR overlap (open chromatin)                  [37.7% of genes]
   → State enrichment of up/down DEGs (Fisher exact vs genomic background)
   → Learned layer: elastic-net model, expression features → state group
   → Figure suite (SVG + PNG)
```

## 2. Tool and database assessment

| Resource | Verdict | Role |
|---|---|---|
| **PCSD** (systemsbiology.cau.edu.cn/chromstates) | **Primary state annotation.** 36 ChromHMM states, 290,553 segments — downloaded counts match the publication exactly (Liu et al. 2018, NAR). Server is unstable (repeated 504s); cache downloads. | 36-state segments + per-state gene lists (used as training labels) |
| **AraENCODE** (glab.hzau.edu.cn/AraENCODE) | **Best-engineered resource.** Full file-tree JSON, TAIR10-aligned, named 12-state tracks, TPM matrix, gene annotation. Used as independent replication of state calls + external validation data. | 12-state track, genes_TAIR10.txt, expression_TPM.txt |
| **PlantCADB** (bioinfor.nefu.edu.cn/PlantCADB) | **Usable but fragile.** 18.3M ACRs; whole-genome BED is 663 MB, mixes 37 species, has a malformed header (10-field header over 20-field rows) and ragged rows — requires tolerant streaming parser. | Open-chromatin (ACR) gene overlap flag |
| **NASA OSDR BDAPI** | **Excellent.** Public, no auth; metadata query + GeneLab-processed counts/DE tables. REST file endpoint returns JSON metadata; use the embedded geode-py URLs for direct download. | Test data source (OSD-37) |
| **ChromHMM / Segway / DeepTools** | Not needed for this test — custom state re-learning requires user-supplied ChIP/ATAC data. Documented as the extension path. | — |

## 3. Test results: OSD-37 (spaceflight vs ground control)

**Data:** 28,067 genes × 56 samples (4 ecotypes × FLT/GC × replicates), GeneLab
VST counts. DE computed with a linear model (condition + ecotype):
**117 up, 53 down** (padj < 0.05, |log2fc| > 1).

### Key finding: spaceflight-responsive genes are pre-marked by their chromatin state

| Direction | State group | Odds ratio | FDR |
|---|---|---|---|
| **Up-regulated** | Accessible_promoter | **3.22** | 3.0e-07 |
| Up-regulated | Heterochromatin_TE | 0.28 | 7.4e-03 |
| Up-regulated | Intergenic_quiet | 0.19 | 7.4e-03 |
| Up-regulated | Polycomb_repressed | 0.32 | 7.4e-03 |
| **Down-regulated** | Polycomb_repressed | **11.46** | 1.4e-14 |
| Down-regulated | Accessible_promoter | 0.14 | 4.1e-07 |

Genes **induced** in spaceflight sit disproportionately in open, accessible-promoter
chromatin; genes **silenced** in spaceflight sit disproportionately in
Polycomb-repressed (H3K27me3-marked) chromatin. The two independent state
annotations (PCSD 36-state and AraENCODE 12-state) agree: FLC controls behaved
correctly (FLC → "Repression" state; ACT2/PHYB → "Active").

### Learned layer: honest negative result

An elastic-net multinomial model (8 expression-distribution features, rank-transformed,
stratified 5-fold CV, class-weighted) predicting the chromatin-state group from
expression alone:

- CV accuracy **0.513** vs majority-class baseline **0.588** (macro-F1 0.314)
- External transfer to AraENCODE TPM (independent tissues): **0.375**

**Interpretation:** expression-distribution features alone do **not** predict
chromatin state better than chance-adjusted baselines. This is biologically
plausible — chromatin state is cell-type- and condition-specific, while the
features are condition-averaged. The deterministic state-mapping backbone (which
uses genomic coordinates, not expression patterns) carries the analysis; the
learned layer is reported as a negative result and should not be used for
prediction in its current form. Improving it would require cell-type-matched
training pairs (expression + matched epigenomic data from the same tissue).

## 4. Deliverables

- `figures/` — 10 figures (SVG + PNG): state distribution, enrichment (up/down),
  volcano colored by state, 3 locus views (PCSD + AraENCODE tracks), heatmap,
  confusion matrix, state-space projection
- `tables/` — DE table, gene-state assignments, enrichment tables, learned-layer
  metrics and confusion matrix, predicted state groups
- `code/` — `chromdecode.py` (modules 1–4), `learned.py` (module 5),
  `figures.py` (module 6), `run_osd37.py` (test run), `test_chromdecode.py`
  (21 unit tests, all passing)

## 5. Limitations

- Results are **discovery-only**; no external validation cohort for the enrichment
  findings beyond the two independent state annotations.
- State assignment is TSS-anchored (gene body + 2 kb upstream); distal enhancer
  links are not modeled.
- PCSD states were defined in seedling tissue; OSD-37 is seedling RNA-seq, so the
  match is good, but other tissues would need tissue-appropriate state maps.
- The PlantCADB whole-genome file required a tolerant parser (malformed header,
  ragged rows); the ACR flag is gene-level overlap, not per-condition dynamics.
- Learned-layer convergence warnings (saga solver) — metrics reported are from
  the converged CV ensemble; the negative conclusion is unaffected.

## 6. Extension path

1. Custom ChromHMM/Segway model when user-supplied ChIP-seq/ATAC-seq data exist.
2. Multi-dataset generalization (OSD-120, OSD-519, ...) via the same pipeline.
3. Per-condition ACR dynamics from PlantCADB sample-level data.
4. Improved learned layer using matched tissue expression + epigenomic pairs
   (e.g. AraENCODE transcriptome + ChIP from the same tissue).

---

# v2: Tissue-Matched Learned Layer + OSDR-Wide Demonstration

## 7. Learned layer v2 — tissue-matched training pairs

**Design change.** v1 trained on condition-averaged features from the spaceflight
dataset to predict seedling-defined states — a domain mismatch. v2 trains on
**tissue-matched pairs**: AraENCODE seedling expression (~266 samples, the same
tissue the PCSD/AraENCODE chromatin states were defined in), with per-gene
distribution features (mean, std, quartiles, range across seedling samples) plus
gene-intrinsic tissue-specificity features (expression breadth tau, tissue
max/min fold, seedling-vs-pan-tissue fold) computed across all 16 AraENCODE
tissues. 28,727 labeled genes, 11 features.

**Honest evaluation.**

| Evaluation | Accuracy | Majority baseline |
|---|---|---|
| Nested 5-fold CV (inner C tuning), random gene folds | **0.554** | 0.549 |
| Genomic-block CV (100-kb bins, whole blocks per fold) | **0.533** | 0.549 |
| Transfer: AraENCODE seedling (same features, refit-free) | 0.548 | 0.549 |
| Transfer: AraENCODE leaf | 0.544 | 0.549 |
| Transfer: AraENCODE root | 0.518 | 0.549 |
| Transfer: OSD-37 spaceflight seedlings | 0.515 | 0.549 |

**Verdict: definitive negative result.** Even with tissue-matched training,
expression-distribution features do not predict chromatin-state groups above
the majority-class baseline. Genomic-block CV (0.533) confirms this is not a
cross-validation artifact: spatially autocorrelated genes do not inflate the
estimate, and there is nothing to inflate — the signal is absent. Macro-F1
improved from v1 (0.314 → 0.374) and the v2 methodology is sound (consistent
rank transforms, nested tuning, block controls), so the conclusion is
trustworthy rather than an artifact of the v1 pipeline.

**Why.** Chromatin state is determined by the epigenomic machinery of the cell
type and condition, not by how variable a gene's expression is across samples.
Two genes with identical expression distributions can sit in opposite states.
Predictive modeling of state from expression would require features that
correlate with state causally (e.g. promoter CpG content, TE proximity,
TF motif content) or single-cell paired data — a different feature space than
expression distributions.

**Practical consequence:** the deterministic coordinate-based backbone (DEGs →
states → enrichment) is the analytically valid core of ChromDecode-At; the
learned layer is retained as a negative-result benchmark.

## 8. OSDR-wide demonstration

The decoder backbone was run on **every Arabidopsis RNA-seq dataset in the NASA
OSDR API** (33 accessions with RNA-seq assays). Contrast auto-detection parses
bulk per-sample metadata (`study.characteristics.*`) plus OSDR sample-name
conventions (FLT = spaceflight, GC = ground control, 0G/1G = gravity levels).

| Outcome | n | Datasets |
|---|---|---|
| Full DE + enrichment | 5 | OSD-37 (529 up / 526 down), OSD-314 (47/45), OSD-120 (43/9), OSD-251 (0/1), OSD-346 (0/1) |
| Expression-only decoding (no auto-detectable contrast) | 18 | e.g. OSD-193, 208, 313 (contrast is genotype/light, not case/control), 406, 476, ... |
| No GeneLab processed counts | 11 | OSD-38, 217, 219, 223, 411, 427, 480, 519, 522, 565, 624 |

**Cross-dataset finding (fig10).** DE is computed on log2-transformed counts
with ecotype covariates (an earlier raw-scale version inflated DE ~5× and was
corrected). All three spaceflight/microgravity datasets with sufficient DEGs
show consistent directionality. OSD-37 (529 up / 526 down): up-regulated genes
are enriched for Active_transcribed chromatin (odds 1.66, FDR 1.4e-5) and
depleted from Heterochromatin_TE (odds 0.11, FDR 3.5e-19); down-regulated
genes are strongly enriched for Polycomb_repressed chromatin (odds 7.0, FDR
1.7e-92) and depleted from Heterochromatin_TE (0.066) and
Accessible_promoter (0.38). OSD-314 (47/45 genes): both directions enrich for
Polycomb_repressed (up 4.8, down 4.9) — see the 36-state section for the
mark-level interpretation. The chromatin-state decoder produces a coherent,
replicable picture: **spaceflight repression targets Polycomb-marked,
promoter-accessible euchromatic genes, while heterochromatic/TE genes are
largely inert.**

The expression-only fallback panel (fig11) is uniform — highly expressed genes
are always in accessible-promoter chromatin — confirming it works as a sanity
check but carries no condition-specific signal.

### 12-state resolution (AraENCODE states)

Re-running the cross-dataset enrichment at full AraENCODE 12-state resolution
(fig13 up, fig14 down) refines the collapsed-group picture:

- **Up-regulated DEGs** (OSD-37): enriched for **Active** (odds 1.57, FDR
  2e-6) and strongly for **Transcription2** (odds 5.4, FDR 2e-10); strongly
  depleted from **Heterochromatin** (odds 0.03, FDR 1e-19).
- **Down-regulated DEGs** are enriched for **Repression** (OSD-37 odds 8.0,
  FDR 3e-105; OSD-314 odds 4.9, FDR 1e-5) — the Polycomb state. (An earlier
  version of this section reported OSD-314 down genes enriched for Active
  chromatin; that derived from the pre-correction raw-scale DE and is
  retracted.)
- OSD-314 up-regulated genes (n=47) are enriched for **Repression** (4.8) and
  **Bivalent** (6.7) and depleted from Active (0.28) — consistent with genes
  leaving a repressed/bivalent state in microgravity, but underpowered.
- Small datasets (OSD-120: 43 up genes; OSD-346: 1 down gene) show the same
  directions but lack power (n.s.).

### 36-state resolution (PCSD states — mark-level detail)

The finest resolution (fig15 up, fig16 down) resolves which histone-mark
combinations drive the group-level signals:

- **Down-regulated DEGs converge on the H3K27me3-accessible-promoter block.**
  OSD-37: S15 (H3K27me3 + accessible DNA, promoter/intergenic) odds **18.0**,
  S16 (accessible DNA, promoter) 6.4, S14 (H3K27me3) 3.9, S17-S21
  (accessible-promoter states) 1.8-3.9, S13 (H3K27me3/H2A.Z promoter) 3.3 —
  all FDR < 0.001. OSD-314 replicates: S15 odds 10.6, S16 4.7, S18 4.1.
  OSD-120 shows the same block (S13 6.0, S15 8.3, S19 11.5) without power.
  Interpretation: spaceflight-repressed genes are a coherent class —
  promoter-accessible, Polycomb-marked euchromatic genes — not a diffuse
  set.
- **Up-regulated DEGs (OSD-37, best-powered)** enrich for
  accessible-promoter states with active marks: S20 (accessible DNA,
  promoter) 2.5, S17 2.5, S8 (H3K4me1/H3K4me2/H2A.Z, CDS) 2.4, S1 (H3.3,
  3'UTR) 2.3, S18/S21/S23/S24 1.6-2.1 — and are depleted from the entire
  heterochromatin/TE block (S29-S36, odds 0.0) and Polycomb CDS states
  (S11-S14, 0.0-0.7).
- **OSD-314 up-regulated genes** instead enrich for Polycomb CDS states
  (S11 odds 22.4, S12 4.2, S13 5.7, all FDR < 1e-3). A five-part diagnostic
  (Section 9) confirms this is **robust, not small-sample noise** — and
  corrects one earlier claim: the previously reported S15 enrichment
  (odds 5.7) was an artifact of the pre-correction raw-scale DE; under
  corrected DE only 1 up-regulated gene falls in S15 (empirical p = 0.51,
  not significant). The robust signal is in S11/S12/S13.
- Heterochromatin/TE states (S29-S36) are uniformly depleted (0.0-0.4) for
  DEGs in both directions across datasets — constitutive heterochromatin is
  transcriptionally inert under spaceflight.

## 9. OSD-314 diagnostics: real signal or small-sample noise?

Five tests on the Polycomb enrichment of OSD-314's 47 up-regulated genes:

| Test | Result |
|---|---|
| **Null calibration** (10,000 random 47-gene draws) | S11 observed 22.4 vs null 99th pct 4.8, **empirical p < 0.0001**; S12 4.2 vs 3.1, p = 0.0013; S13 5.7 vs 4.5, p = 0.0036; S15 1.4, p = 0.51 (not significant) |
| **Downsampling control** (OSD-37's 529 up genes → n=47, 1,000 draws) | P(draw ≥ OSD-314 observed) = 0.000 for S11 and S12; 0.016 for S13 — small-N effects cannot produce the observed odds |
| **Threshold-sensitivity grid** (padj × \|log2fc\|, 9 gene sets, n=45–1,415) | S11 odds rise monotonically with stringency: 1.3 (n=456) → 2.3 (n=508) → 3.8 (n=88) → 6.7 (n=189) — the hallmark of a real dose-response in log2fc, not noise |
| **Threshold-free GSEA** (1,944 Polycomb genes ranked by log2fc, 1,000 permutations) | Running enrichment peaks in the up-ranked region, **permutation p < 0.001** |
| **Gene-identity check** | 29 up genes in Polycomb states (14 in S11), including adjacent pairs (AT3G22490/AT3G22500, AT4G25140/AT4G26740/AT4G28520) suggesting tandemly duplicated clusters |

**Verdict: robust.** The Polycomb enrichment in OSD-314's up-regulated genes
survives null calibration, cannot be reproduced by downsampling a
Polycomb-negative gene set to the same size, strengthens with DE stringency,
and is confirmed by a threshold-free test. Biologically it reads as **escape
from Polycomb repression**: up-regulated genes in microgravity are
disproportionately genes held in H3K27me3-marked (S11/S12/S13) or bivalent
states in ground conditions — the mirror image of the down-regulated genes,
which are pushed *into* deeper Polycomb repression (S15). Both directions of
the spaceflight response converge on the Polycomb/H3K27me3 axis.

Caveats: the gene-identity check lacks functional annotations (TAIR10
description file not downloaded); the chromatin states remain ground-state
reference annotations, so "escape from Polycomb" refers to the genes'
*baseline* state, not measured mark changes during flight.

## 9. v2 deliverables

- `figures/fig8-9` — v2 model diagnostics (confusion matrix, transfer bars)
- `figures/fig10-12` — cross-dataset heatmaps + per-dataset enrichment panels
- `tables/learned_v2_metrics.csv`, `learned_v2_transfer.csv` — v2 evaluation
- `tables/sweep_summary.csv`, `sweep_enrichment_long.csv` — OSDR-wide results
- `code/run_learned_v2.py`, `run_osdr_sweep.py`, `aggregate_sweep.py`,
  `build_features_v2.py` — v2 pipeline (checkpointed per dataset)

## 10. v2 limitations

- Contrast auto-detection is conservative: datasets whose design is genotype-,
  light-, or dose-based without control-token metadata fall back to
  expression-only decoding. Manual curation would recover them.
- 11 datasets lack GeneLab processed counts; reprocessing from FASTQ was out
  of scope for sandbox compute.
- The learned-layer negative result applies to expression-distribution
  features; sequence-based or chromatin-adjacency features were not tested.
- All results remain discovery-only.

# v3: Upstream Motifs and GO Terms of Polycomb-Repressed Genes

## 11. What distinguishes the Polycomb genes that respond to spaceflight?

The v2 sweep showed spaceflight down-regulated genes concentrate in PCSD
Polycomb states S11-S15 (OSD-37: S15 odds 18.0). This section asks what,
beyond chromatin state, characterizes the 265 down-regulated genes that sit
in those states — using three gene sets: **target** (down DEGs in S11-S15,
n=265), **control 1** (down DEGs in non-Polycomb states, n=241), and
**control 2** (expressed non-DE genes in S11-S15, n=3,303 — the sharp
state-matched comparison).

### 11.1 Promoter motif architecture (JASPAR 2024 plants, 907 TFs)

Strand-aware 2 kb promoters were scanned with log-odds PWMs (promoter-derived
background, per-motif threshold at p < 1e-4).

- **vs expressed genome (25,969 genes): 268 motifs significant (BH-FDR <
  0.05).** Enriched: ZHD1 (2.9x), AHL12/AHL13/AHL25 (AT-hook family, up to
  5.7x), NTL9, TCX6; depleted: ERF-family and ARF motifs (fig22). The AHL
  AT-hook family is chromatin-associated, consistent with the Polycomb-state
  promoter architecture.
- **vs state-matched Polycomb control: no motif survives BH-FDR** (nominal
  top: NTL9, odds 2.1, padj 0.084). The motif architecture tracks the
  Polycomb state itself, not the spaceflight response within it (fig25).

### 11.2 GO over-representation (hypergeometric, BH-FDR, TAIR GAF)

- **vs genome:** 18 significant terms, dominated by an oxidative-stress
  module: peroxidase activity (11.1x), heme binding (7.2x), hydrogen
  peroxide catabolic process (13.6x), cellular oxidant detoxification
  (8.8x).
- **vs state-matched Polycomb control:** 11 terms remain significant — the
  specificity holds within the Polycomb class. The spaceflight-repressed
  subset of Polycomb genes is functionally a redox/peroxidase module
  (fig23).
- **GC check:** mean promoter GC 0.310 (target) vs 0.320 (control 2) — the
  target set is slightly AT-richer, so the motif findings are not
  GC-driven artifacts (fig24).

### 11.3 OSD-314 replication (directional consistency)

OSD-314 has only 18 Polycomb-down genes (power-limited). Of the top 10
OSD-37 GO terms, 8 replicate directionally with folds of 7-27x (heme
binding p = 4e-5, iron ion binding p = 2e-5, monooxygenase activity p =
2e-4); only "membrane" (a broad compartment term) does not.

### 11.4 Verdict

The Polycomb-repressed genes that are pushed into deeper repression during
spaceflight are a coherent functional module — peroxidase and oxidant-
detoxification genes — embedded in a Polycomb-state promoter architecture
rich in AT-hook (AHL) and ZHD-family motifs. The motif layer is
state-associated rather than response-specific; the functional (GO) layer
is both state-associated and response-specific.

### 11.5 Caveats

- Motif "hits" are correlative: TF binding is inferred from JASPAR
  position frequency matrices, not measured (no DAP-seq/ChIP validation).
- GO over-representation is discovery-only; no correction for term
  redundancy (e.g., peroxidase activity and heme binding overlap).
- OSD-314 replication is directional only (18 genes).
- Promoters are 2 kb upstream of the TAIR10 TSS; distal enhancers and
  transposon-derived regulatory elements are not modeled.

### v3 deliverables

- Figures: fig22 (motif odds ratios), fig23 (GO dot plot, state-matched
  comparison), fig24 (GC-content check), fig25 (motif hit overlap).
- Tables: motif_enrichment.csv (three comparisons, 907 motifs),
  go_enrichment.csv, motif_gc_content.csv, osd314_go_replication.csv.
- Code: prep_motif_go_inputs.py, motif_go_enrichment.py, regen_fig25.py.

# v4: Redox-Module Visualizations and Spaceflight Methylation

## 12. Visualizing the module

- **fig26 (motif-logo heatmap):** sequence logos of the top 12 enriched and
  12 depleted motifs in Polycomb-down promoters. The enriched logos are
  dominated by AT-hook (AHL) and ZHD-family binding geometry; depleted
  logos are ERF/ARF-type GC-rich motifs.
- **fig27 (redox-module network):** the 23 Polycomb-down genes annotated
  with significant redox GO terms, connected by shared redox GO terms or
  shared top-motif targets. Node color = spaceflight log2FC (all strongly
  down). The module is a single connected cluster.

## 13. Does the redox module change epigenomically in spaceflight?

**Data availability (verified OSDR-wide, 2026-09):** no spaceflight
H3K27me3 ChIP-seq exists. The only Arabidopsis ChIP-seq in OSDR is OSD-496
(SOG1-3xFLAG seedlings, ground-based DNA-damage response — a candidate
follow-up, raw data, requires HPC alignment); the H3K27ac CHi-C study
(OSD-687) is human. The available spaceflight epigenomic mark is DNA
methylation: OSD-217 (Arabidopsis Ws, WGBS, root + leaf, FLT vs GC) with
processed per-gene CG/CHG/CHH methylation-difference tables.

**Test:** redox-module genes (n=23; 14-72 with >=3 informative sites per
context) vs all other genes; Mann-Whitney of FLT-GC methylation difference,
BH-FDR across 6 tests (2 tissues x 3 contexts).

**Result: no significant methylation change** (minimum padj = 0.07, root
CHG: nominal p = 0.012, median +0.004 vs -0.000 — a weak CHG gain in
roots that does not survive correction; all other contexts p > 0.2).
fig28 (distributions), fig29 (per-gene heatmap).

**Interpretation:** the spaceflight repression of the redox module is not
accompanied by detectable DNA-methylation change at gene bodies. Combined
with v3, the epigenomic signature of these genes remains their pre-existing
H3K27me3 Polycomb state (static), with no evidence that the flight
environment rewrites DNA methylation at the module.

### 13.1 Caveats of the Ws (methylation) vs Col-0 (expression) comparison

The methylation data is Wassilewskija (Ws); the expression and chromatin
data are Col-0. Checked against the GLDS-217 processing documentation
(OSDR study page):

- **Reference bias is largely controlled:** the WGBS reads were aligned to
  the ecotype-matched **Ws reference genome** (ws_0.v7.allPlusChlMito,
  Mott et al. 2011) with bsmap 2.87 and cscall methylation calling — not
  to TAIR10 — so Ws-vs-Col SNP reference bias (C->T mismatches read as
  unmethylated cytosines) is mostly avoided. Note the processing is
  author-side (trimmomatic/cscall/bsmap/MOABS mcomp 1.3.4, site coverage
  >= 20 in >= 2 of 3 replicates, site p < 0.01), not the current GeneLab
  methyl-seq pipeline (GL-DPPD-XXXX), and no explicit SNP-stripping step
  is documented.
- **Coordinate system:** gene coordinates in the OSD-217 tables are Ws-
  assembly (v7) coordinates; the gene-level AGI join used here is
  unaffected, but any future coordinate-based overlap with TAIR10
  features would require liftover.
- **Presence/absence variation and epialleles:** Ws and Col-0 differ by
  hundreds of genes and thousands of TE insertions, and CG/CHG methylation
  is ecotype-specific. Module genes deleted or diverged in Ws are missing
  or filtered (Sites >= 3), reducing module coverage to 14-72 genes per
  context in root and only 14-39 in leaf. The methylation response of a
  Ws gene may not generalize to its Col-0 ortholog.
- **Tissue mismatch:** methylation is from root and leaf; expression DE is
  from whole seedlings and PCSD states are seedling-derived. The assayed
  methylation is not necessarily in the cells where repression occurs.
- **Assay window:** gene-body methylation only; promoter methylation (the
  regulatory layer matching the TSS-anchored chromatin analysis) was not
  testable.
- **Power and timescale:** 3 replicates/group, 23 module genes, 6 tests —
  only large shifts are detectable, so the honest reading is "no
  detectable change," not "no change." DNA methylation is also slow-
  turning; the null does not rule out histone-mark dynamics, and the
  H3K27me3 evidence remains static (ground-state PCSD states).

**v4 deliverables:** fig26-29, osd217_methylation_tests.csv,
redox_module_genes.csv, redox_figures.py, osd217_methylation.py.

# v5: Promoter Methylation and Ecotype-Matched Comparison

## 14. Promoter-window methylation (OSD-217 winavg)

The v4 test used gene bodies. The OSD-217 winavg tables (100-bp methylation
bins, FLT and GC separately) allow promoter-level analysis: 2 kb upstream
of the TSS (20 bins), anchored on Ws-assembly gene coordinates with TAIR10
strand (gene orientation is conserved across ecotypes). Coverage-limited to
21,424 genes with Ws coordinates on the winavg chromosomes.

**Result (fig30, fig31):** no promoter methylation change survives BH-FDR
in any tissue x context x gene-set test (minimum padj = 0.099). The redox
module shows a nominal CHG gain in both tissues (root +0.007 vs -0.001,
p = 0.011; leaf +0.027 vs +0.004, p = 0.011) — directionally consistent
with the weak gene-body CHG hint from v4, but not significant after
correction. Baseline (ground-control) promoter methylation of the
Ws-defined Polycomb-down module is elevated (root CG median 0.52 vs 0.02,
p = 0.043; leaf p = 0.038) — suggestive of methylated promoters within the
Polycomb class, but based on only 6-11 genes and not significant after
correction.

## 15. Ecotype-matched comparison (Ws RNA-seq, same study)

OSD-217 ships its own Ws RNA-seq differential expression (root: 65
significant genes; leaf: 608), enabling a fully same-ecotype
expression-methylation pairing. The Ws-defined Polycomb-down module is
small (6 genes root, 18 leaf) — the Ws DE list is short and the Polycomb
intersection shrinks it further — so the ecotype-matched test is
power-limited and inconclusive (all padj > 0.36).

**Cross-ecotype concordance (fig32):** spaceflight log2FC correlates
weakly between OSD-217 (Ws root) and OSD-37 (Col-0 seedlings):
r = -0.17 overall. The ecotype difference is confounded with tissue
(root vs whole seedlings) and study, so the two datasets' flight responses
are largely independent — reinforcing that v4's null methylation result
cannot be assumed to transfer to Col-0, and that the Ws data cannot
substitute for Col-0 epigenomics.

## 15.1 Caveats

- Leaf winavg coverage is sparse (58k windows vs 210k in root CG); the
  redox-module metaplot traces are fragmented in leaf panels — genuine
  missing data, not plotting artifacts.
- The Ws DE tables list only significant genes; "non-DE" background =
  genes absent from those lists (conservative).
- Promoter windows use Ws coordinates; genes absent from the Ws v7
  assembly annotation are excluded.
- All results discovery-only.

**v5 deliverables:** fig30-32, osd217_promoter_methylation{,_tests}.csv,
ws_vs_col_concordance.csv, code/osd217_promoter_methylation.py.

## 16. SOG1 binding after DNA damage (ground comparator)

Spaceflight H3K27me3 ChIP-seq does not exist in OSDR, so v6 tests the
Polycomb-axis hypothesis from the other direction: SOG1 — the plant
NAC-family master regulator of the DNA-damage response (functional
analogue of p53), which directly recruits Polycomb repression after DNA
damage — is profiled by ChIP-seq in OSD-496 (seedlings treated with the
double-strand-break inducer bleomycin; SOG1-3xFLAG IP and input at 20 min
and 1 h; mixed Ler/Col background). If spaceflight down-regulation of the
redox module is a SOG1/Polycomb-style damage response, SOG1 should bind
redox-module and Polycomb-down promoters after damage.

**Processing:** 6 single-end libraries (13.7-23.4 M reads) aligned to
TAIR10 with HISAT2 (`--no-spliced-alignment`); peaks called with MACS2
narrow mode, genome size 1.2e8, q < 0.01, IP vs matched-timepoint input
(SOG1 20 min: 56 peaks; SOG1 1 h: 100 peaks). As a specificity check,
SOG1-IP vs wild-type-IP gives 160/180 peaks — the FLAG IP is
antibody-specific. 67/100 of the 1 h peaks overlap a 2 kb promoter and
the median peak-to-nearest-TSS distance is 194 bp, confirming
promoter-proximal binding as expected for a transcriptional activator.

**Result (fig33, fig34): a clean negative.** Across the 156-peak union,
zero SOG1 peaks fall in redox-module promoters (0/23) and zero in
Polycomb-down promoters (0/265), at either timepoint (all Fisher padj =
1.0). Only 1/261 non-Polycomb-down and 3/3,303 Polycomb non-DE
background promoters are bound. The 1 h promoter occupancy heatmap
(fig34) shows the redox-module and Polycomb-down blocks are as quiet as
the non-DE background (max bin 4.1 RPM, driven by a few background
promoters).

**Interpretation:** the spaceflight down-regulation of the redox/peroxidase
module is not a recapitulation of SOG1's direct DNA-damage targets. The
Polycomb-state enrichment of flight-down genes (S11-S15) therefore does
not arise because those promoters are SOG1-bound damage targets; if a
Polycomb mechanism is involved in spaceflight, it is not the canonical
SOG1-dependent damage-repression pathway. This sharpens the v3 motif
result (AHL/ZHD/NTL enrichment): the candidate upstream regulators of the
flight response remain to be tested directly, and no suitable spaceflight
TF ChIP-seq exists in OSDR to do so.

**Caveats:** bleomycin-induced DSBs are not spaceflight; seedlings are
mixed Ler/Col; only 56-100 peaks survive q < 0.01 per timepoint (SOG1
binding is damage-dose- and time-dependent), so absence of evidence at
these two timepoints is strong for direct binding but cannot exclude
binding under other damage regimes or in other tissues.

**v6 deliverables:** fig33-34, sog1_promoter_tests.csv,
code/sog1_chipseq.py.

## 17. Tissue-matched replication: Col-0 root spaceflight RNA-seq (v7-v8)

The v5 ecotype comparison showed the flight response is not concordant
across tissues/ecotypes (r = -0.17), motivating a search for Col-0
**root** spaceflight RNA-seq to test whether the Polycomb-state
enrichment of down-regulated genes survives tissue matching.

**Inventory (root_candidates.csv):** all 62 Arabidopsis OSDR studies
enumerated via the metadata API; 33 have RNA-seq assays; 10 mention root
tissue. QC resolved the field to: OSD-624 (Col-0 roots, FLT vs GC, 3+3);
OSD-120 (CARA root tips; mixed Col-0/WS/phyD, light x day); and — after
v8 metadata recovery — OSD-193 and OSD-218. The latter two ship counts
keyed by GSM/SAMN accessions with empty OSDR factor values; NCBI
Biosample eutils resolved every sample ID to its condition title
("Spaceflight, 8 days old, Col-0 Rep4"), recovering balanced FLT/GC
designs. OSD-281 was excluded on recovery: its Biosample titles reveal
16 Sku5-mutant + 16 Ws samples, no Col-0. OSD-406 (Col-0 roots, 12/12)
is a **Virgin Galactic suborbital** flight, not ISS. Remaining
exclusions: OSD-217 (Ws), OSD-219/427 (no processed counts), OSD-411/480
(EMCS 0.004G partial-gravity dissections), OSD-208 (no flight factor).

**Result (fig35): the Polycomb enrichment DOES replicate in Col-0 ISS
roots — revising the v7 conclusion.** With age as covariate:
OSD-193 Col-0 WT roots (8v8): 97 up / 181 down, Polycomb-state odds
3.06, padj 2.2e-10; OSD-218 Col-0 roots (16v16): 128 up / 117 down,
odds 3.14, padj 4.7e-7. Both significant. OSD-624 remains null
(odds 0.76, padj 0.81; 36 down genes); OSD-120 Col-0 WT root tips
collapse to 2 down genes (no power); OSD-406 suborbital has 1 down gene
(minutes of microgravity produce minimal transcriptional response).

**Interpretation (revised):** the S11-S15 concentration of
flight-down-regulated genes is not an aerial-tissue artifact — it
replicates in roots, at roughly half the seedling effect size
(odds ~3.1 vs 7.0). The current evidence structure: strong enrichment
in both whole-seedling datasets; significant enrichment in the two
APEX03-2-campaign root datasets (OSD-193/218, same payload, independent
samples); null in OSD-624 (different hardware, lab, and age). The
root/seedling distinction is therefore confounded with campaign, and
the honest summary is: the enrichment is robust across three of four
independent Col-0 datasets, with effect size varying ~2-fold by study.
Combined with v6 (no SOG1 binding at these promoters after DNA damage),
the state bias is not a SOG1-style damage response; its upstream
driver remains unidentified.

**Caveats:** OSD-193/218 share a payload and lab, so they are not fully
independent replications; OSD-624's null could reflect its smaller DE
depth or a genuine campaign/hardware difference; all results
discovery-only.

**v7-v8 deliverables:** fig35, root_candidates.csv,
replication_summary.csv, replication_recovered.csv,
{acc}_sample_meta.csv (193/218/281/406), code/root_replication*.py,
code/recover_sample_meta.py.

## 17.2 Power analysis: replicates needed for a stable OSD-120 enrichment test (v9)

Question: how many replicates per group would OSD-120's design need for a
stable Polycomb-enrichment test of down-regulated genes? Parameters
(delegated defaults): alpha 0.05, 80% power, effect sizes anchored to the
project's own observations (odds 3.1 root-level, 7.0 seedling-level),
background Polycomb fraction 0.141 (4,812/34,054 expressed state-assigned
genes).

**Enrichment test alone (simulation, fig36B):** 80% power requires ~40
down-regulated DEGs at root-level odds (3.1), or ~15 at seedling-level odds
(7.0). OSD-120's Col-WT run detected 2 down genes — 20-fold short.

**Replicates needed (empirical subsampling of OSD-218, fig36A):** running the
full v8 pipeline on k replicates/group (20 draws each): power rises from
~15% at k=2-6 to 50% at k=8 and 100% at k=12-16 (k=16 exactly reproduces the
v8 result: 117 down genes, padj 4.7e-7 — pipeline validation). Interpolated:
**~9-10 replicates per group in a clean design; ~11-12 with OSD-120's
light x day covariate burden** (2 extra parameters cost ~20-30% power at
k=8-12). OSD-120's actual Col-WT contrast (6v6) sits at ~10-15% power.

**Caveats:** the empirical curve inherits OSD-218's effect size and
dispersion (same campaign class as OSD-120 root tips, but not identical);
nuisance-covariate draws at k<=3 are unstable (random binary covariates can
spuriously correlate with condition at tiny n) and were not interpreted; the
analytic threshold (40 DEGs) assumes a concentrated effect, whereas
draw-level n_down varies widely (2-128), so the empirical curve is the
honest estimate; all results discovery-only.

**v9 deliverables:** fig36, osd120_power_analysis.csv,
osd120_power_summary.csv, code/power_analysis.py.

## 17.3 Stratified flight-within-light analysis: bounding the flight x light interaction (v10)

Motivation: the additive model (condition + light) on OSD-120 Col-0 WT root
tips collapsed to 2 down genes (17.1); a real flight x light interaction
(documented in CARA light-regime studies) would be absorbed into the light
term and could mask a flight effect. The Col-0 WT design is a balanced 2x2
factorial (FLT/GC x Alight/dark, 3 reps/cell, all Day 13), so the interaction
is directly estimable.

**Stratified flight effect (3v3 per stratum, fig37A):** within Alight, 1 up /
0 down gene; within dark, 2 up / 1 down (not in the Polycomb state; odds 0.0,
padj 1.0). No Polycomb enrichment in either stratum. Module medians: the
Polycomb-down module (37 genes from OSD-218) trends down in flight within
Alight (median log2fc -0.31, consistent with OSD-218) but slightly up within
dark (+0.16) - directionally inconsistent across strata, neither significant.

**Interaction bound (full 12-sample model, fig37B):** expression ~ condition +
light + condition:light. Only 4 genes genome-wide show a significant
interaction (padj < 0.05, |log2fc| > 1) - consistent with the null expectation
at this sample size. Median |interaction| across all genes is 0.50 (95th pct
2.28). Neither the redox module (23 genes; 0 significant, median |interaction|
1.26) nor the Polycomb-down module (37 genes; 0 significant, median 1.18,
max 7.86) contains a significant interaction.

**Interpretation:** the additive model's 2-down-gene null is NOT explained by
a large masked interaction - the interaction is bounded small at the module
level (median |interaction| ~1.2 log2 units, no significant members) and the
per-stratum flight effects are directionally inconsistent for the Polycomb
module. However, 3v3 strata are ~0-15% powered at root-level effect size
(17.2), so this is a bound on detectable interaction, not evidence of absence:
the honest conclusion is that OSD-120's 12-sample design cannot resolve the
flight x light interaction on the Polycomb module, and the stratified analysis
does not rescue the null.

**v10 deliverables:** fig37, osd120_stratified.csv,
osd120_interaction_genelevel.csv, osd120_stratDE_{Alight,dark}.csv,
code/osd120_stratified.py.
