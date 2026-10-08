# Supplementary tables

All tables are CSV (UTF-8, comma-separated). File names are the ones the analysis
code writes, so a re-run overwrites them in place under `results/`. The
S-numbers below are the ones the manuscript cites. Statistical conventions,
unless a row says otherwise: BH-FDR within each test family; DEG thresholds
padj < 0.05 and |log2FC| > 1; "Polycomb" = PCSD states S11–S15.

Small derived inputs reused across analyses (the gene→state map, sample
metadata recovered from NCBI BioSample, the redox-module gene list) are in
[`../data/derived/`](../data/derived/) (Table D1–D7).

| Table | File | Content | Written by | Manuscript |
|---|---|---|---|---|
| **S1** | `OSD-37_DE_FLT_vs_GC.csv` | v1 OSD-37 differential expression (GeneLab VST counts, condition + ecotype); 117 up / 53 down | `run_osd37.py` | §2.1 |
| **S2** | `enrichment_up.csv` | v1 OSD-37 state-group enrichment, up-regulated genes (Fisher) | `run_osd37.py` | §2.1 |
| **S3** | `enrichment_down.csv` | v1 OSD-37 state-group enrichment, down-regulated genes | `run_osd37.py` | §2.1 |
| **S4** | `learned_metrics.csv` | v1 learned layer: CV accuracy, macro-F1, majority baseline | `run_osd37.py` / `learned.py` | §2.2 |
| **S5** | `learned_confusion_matrix.csv` | v1 learned-layer confusion matrix | `run_osd37.py` | §2.2 |
| **S6** | `learned_external_validation.csv` | v1 transfer accuracy on AraENCODE TPM | `run_osd37.py` | §2.2 |
| **S7** | `learned_v2_metrics.csv` | v2 learned layer: nested CV and genomic-block CV | `run_learned_v2.py` | §2.2 |
| **S8** | `learned_v2_confusion.csv` | v2 confusion matrix | `run_learned_v2.py` | §2.2 |
| **S9** | `learned_v2_transfer.csv` | v2 transfer to AraENCODE seedling/leaf/root and OSD-37 | `run_learned_v2.py` | §2.2 |
| **S10** | `sweep_summary.csv` | OSDR-wide sweep: status, contrast, n up/down per accession (33) | `run_osdr_sweep.py` | §2.3 |
| **S11** | `sweep_enrichment12_long.csv` | Sweep enrichment at AraENCODE 12-state resolution | `heatmap_12state.py` | §2.3 |
| **S12** | `sweep_enrichment36_long.csv` | Sweep enrichment at PCSD 36-state resolution | `heatmap_36state.py` | §2.3 |
| **S13** | `osd314_null_calibration.csv` | OSD-314 up genes: 10,000 random 47-gene draws per state | `diagnose_osd314.py` | §2.4 |
| **S14** | `osd314_downsampling_control.csv` | OSD-37 up genes downsampled to n = 47 (1,000 draws) | `diagnose_osd314.py` | §2.4 |
| **S15** | `osd314_threshold_grid.csv` | Odds across a padj × \|log2FC\| threshold grid | `diagnose_osd314.py` | §2.4 |
| **S16** | `osd314_gsea.csv` | Threshold-free GSEA on 1,944 Polycomb genes | `diagnose_osd314.py` | §2.4 |
| **S17** | `osd314_up_polycomb_genes.csv` | The 29 OSD-314 up genes in Polycomb states | `diagnose_osd314.py` | §2.4 |
| **S18** | `motif_enrichment.csv` | JASPAR 2024 plant motif enrichment, three comparisons × 907 motifs | `motif_go_enrichment.py` | §2.5 |
| **S19** | `go_enrichment.csv` | GO over-representation (TAIR GAF), vs genome and vs state-matched control | `motif_go_enrichment.py` | §2.5 |
| **S20** | `motif_gc_content.csv` | Mean promoter GC content per gene set | `motif_go_enrichment.py` | §2.5 |
| **S21** | `osd314_go_replication.csv` | OSD-314 directional replication of top OSD-37 GO terms | `motif_go_enrichment.py` | §2.5 |
| **S22** | `osd217_methylation_tests.csv` | OSD-217 gene-body methylation, redox module vs background (2 tissues × 3 contexts) | `osd217_methylation.py` | §2.6 |
| **S23** | `osd217_promoter_methylation.csv` | Per-gene 2 kb promoter methylation (100-bp bins), FLT and GC | `osd217_promoter_methylation.py` | §2.6 |
| **S24** | `osd217_promoter_methylation_tests.csv` | Promoter-window tests per tissue × context × gene set | `osd217_promoter_methylation.py` | §2.6 |
| **S25** | `ws_vs_col_concordance.csv` | Flight log2FC concordance, OSD-217 (Ws root) vs OSD-37 (Col-0) | `osd217_promoter_methylation.py` | §2.6 |
| **S26** | `sog1_promoter_tests.csv` | SOG1 peak overlap with gene-set promoters (20 min, 1 h, union, vs-WT) | `sog1_chipseq.py` | §2.7 |
| **S27** | `root_candidates.csv` | Inventory of OSDR root RNA-seq candidates and QC decisions | not in archive (interactive) | §2.8 |
| **S28** | `replication_summary.csv` | First-pass root replication run (v7) | `root_replication.py` | §2.8 |
| **S29** | `replication_recovered.csv` | Root replication after BioSample metadata recovery (OSD-193/218/406) | `root_replication_recovered.py` | §2.8 |
| **S30** | `replication_unnorm.csv` | OSD-624 run from unnormalized counts (median-ratio scaling); stored as a single serialized record | `root_replication_unnorm.py` | §2.8 |
| **S31** | `OSD-193_DE.csv` | OSD-193 Col-0 WT root DE (8 v 8, age covariate) | `root_replication_recovered.py` | §2.8 |
| **S32** | `OSD-193_enrichment_down.csv` | OSD-193 down-gene state-group enrichment | `root_replication_recovered.py` | §2.8 |
| **S33** | `OSD-218_DE.csv` | OSD-218 Col-0 root DE (16 v 16) | `root_replication_recovered.py` | §2.8 |
| **S34** | `OSD-218_enrichment_down.csv` | OSD-218 down-gene enrichment | `root_replication_recovered.py` | §2.8 |
| **S35** | `OSD-624_DE.csv` | OSD-624 Col-0 root DE (3 v 3) | `root_replication_unnorm.py` | §2.8 |
| **S36** | `OSD-624_enrichment_up.csv` | OSD-624 up-gene enrichment | `root_replication_unnorm.py` | §2.8 |
| **S37** | `OSD-624_enrichment_down.csv` | OSD-624 down-gene enrichment | `root_replication_unnorm.py` | §2.8 |
| **S38** | `OSD-120_ColWT_DE.csv` | OSD-120 Col-0 WT root tips DE (light + day covariates) | not in archive (interactive) | §2.8 |
| **S39** | `OSD-120_ColWT_enrichment_down.csv` | OSD-120 Col-0 WT down-gene enrichment | not in archive (interactive) | §2.8 |
| **S40** | `OSD-406_DE.csv` | OSD-406 suborbital Col-0 root DE (12 v 12) | `root_replication_recovered.py` | §2.8 |
| **S41** | `OSD-406_enrichment_down.csv` | OSD-406 down-gene enrichment | `root_replication_recovered.py` | §2.8 |
| **S42** | `osd120_power_analysis.csv` | Power analysis, per-draw results (empirical subsampling of OSD-218 + simulation) | `power_analysis.py` | §2.9 |
| **S43** | `osd120_power_summary.csv` | Power by replicates per group (k) and nuisance covariates | `power_analysis.py` | §2.9 |
| **S44** | `osd120_stratified.csv` | OSD-120 flight-within-light strata, interaction summary, module bounds | `osd120_stratified.py` | §2.10 |
| **S45** | `osd120_interaction_genelevel.csv` | Gene-level condition × light interaction (12-sample model) | `osd120_stratified.py` | §2.10 |
| **S46** | `osd120_stratDE_Alight.csv` | Flight DE within the Alight stratum (3 v 3) | `osd120_stratified.py` | §2.10 |
| **S47** | `osd120_stratDE_dark.csv` | Flight DE within the dark stratum (3 v 3) | `osd120_stratified.py` | §2.10 |
| **S48** | `qc_genotype_rerun.csv` | QC re-run of OSD-218/406 by genotype (rows 1–2 reproduce v8 exactly) | `qc_genotype_rerun.py` | §2.8, Methods |
| **S49** | `OSD-218_Col0_only_DE.csv` | OSD-218 Col-0-only root DE (8 v 8, age covariate) | `qc_genotype_rerun.py` | §2.8 |
| **S50** | `qc_sweep_rerun_check.csv` | Sweep re-run vs Tables S10–S12 (counts, max odds difference) | `qc_sweep_rerun.py` | §2.3, Methods |
| **S51** | `sweep_enrichment5_corrected.csv` | Corrected 5-group sweep enrichment (replaces `superseded/sweep_enrichment_long.csv`) | `qc_sweep_rerun.py` | §2.3, §2.8 |
| **S52** | `OSD-37_DE_v2_log2.csv` | Corrected (v2, log2-count) OSD-37 DE: 529 up / 526 down | `run_osdr_sweep.py` (via `qc_sweep_rerun.py`) | §2.3 |
| **S53** | `qc_osd314_sample_design.csv` | OSD-314 gravity × light design, and how the sweep coded each sample | `qc_sweep_rerun.py` | §2.4 |
| **S54** | `qc_osd314_design_rerun.csv` | OSD-314 corrected contrasts (0g or 0.3g vs 1g, ± light) and their enrichments | `qc_sweep_rerun.py` | §2.4 |
| **S55** | `OSD-314_0g_vs_1g_light_DE.csv` | OSD-314 DE, 0g vs 1g + light | `qc_sweep_rerun.py` | §2.4 |
| **S56** | `OSD-314_03g_vs_1g_light_DE.csv` | OSD-314 DE, 0.3g vs 1g + light | `qc_sweep_rerun.py` | §2.4 |
| **S57** | `OSD-314_DE_sweep_miscoded.csv` | OSD-314 DE as run by the sweep (0g + 0.3g pooled vs 1g); kept for traceability | `run_osdr_sweep.py` | §2.4 |
| **S58** | `qc_osd314_null_calibration_corrected.csv` | OSD-314 null calibration (10,000 draws) for the published and both corrected contrasts | `qc_osd314_diagnostics.py` | §2.4 |
| **S59** | `qc_osd314_downsampling_control_corrected.csv` | OSD-37 downsampling control; not applicable when n_up > 529 (stated per row) | `qc_osd314_diagnostics.py` | §2.4 |
| **S60** | `qc_osd314_threshold_grid_corrected.csv` | Threshold grid (padj × \|log2FC\|) per contrast, light covariate where applicable | `qc_osd314_diagnostics.py` | §2.4 |
| **S61** | `qc_osd314_gsea_corrected.csv` | Threshold-free running enrichment, both tails, 1,000 permutations | `qc_osd314_diagnostics.py` | §2.4 |
| **S62** | `qc_osd314_up_polycomb_genes_corrected.csv` | Up-regulated genes in Polycomb states, per contrast | `qc_osd314_diagnostics.py` | §2.4 |
| **S63** | `qc_osd314_go_replication_corrected.csv` | OSD-314 replication of the top 10 OSD-37 GO terms, per contrast, BH within contrast | `qc_osd314_diagnostics.py` | §2.5 |
| **S64** | `qc_osd314_diagnostics_reproduction_check.csv` | Published contrast re-run vs S13/S15/S16/S19/S21 (differences ≤ 4e-15; identical counts) | `qc_osd314_diagnostics.py` | Methods |
| **S65** | `qc_power_summary_col0.csv` | Empirical power by k: pooled Col-0 + WS (reproduces S43) and Col-0 only (8 v 8, k ≤ 8) | `qc_col0_power_bound.py` | §2.9 |
| **S66** | `qc_power_draws_col0.csv` | Per-draw results for the Col-0-only power curve | `qc_col0_power_bound.py` | §2.9 |
| **S67** | `qc_osd120_module_bound_col0.csv` | OSD-120 module bound with redox, pooled, Col-0-only and WS-only OSD-218 Polycomb-down modules | `qc_col0_power_bound.py` | §2.10 |
| **S68** | `qc_col0_power_bound_reproduction_check.csv` | Re-run vs S42/S43/S44 (differences ≤ 3e-15) and the k = 8 vs S49 check | `qc_col0_power_bound.py` | Methods |

## Derived data (`../data/derived/`)

| Table | File | Content | Written by |
|---|---|---|---|
| **D1** | `gene_state_assignments.csv` | Every TAIR10 gene → PCSD 36-state, AraENCODE 12-state, state group, ACR flag | `run_osd37.py` |
| **D2** | `predicted_state_groups_OSD-37.csv` | v1 learned-layer predictions per gene | `run_osd37.py` |
| **D3** | `redox_module_genes.csv` | The 23 Polycomb-down redox/peroxidase module genes | `redox_figures.py` |
| **D4–D7** | `OSD-{193,218,281,406}_sample_meta.csv` | Sample conditions recovered from NCBI BioSample titles | `recover_sample_meta.py` |

## Superseded (`superseded/`)

`sweep_enrichment_long.csv` is **not cited and should not be used**. It holds the
5-group sweep enrichment from the pre-correction raw-scale DE: OSD-37 has
≈2,100 up genes in it rather than the corrected 529, and OSD-37 down-gene
Polycomb odds of 1.70. **It is replaced by S51**, which `qc_sweep_rerun.py` regenerated
from the public counts. The table is kept only so the archive is complete.

## Pre-submission audit (2026-10)

Tables S48–S68 come from the QC re-runs described in the manuscript Methods. They also
change how some earlier tables should be read:

- **S12/S13–S17/S21 (OSD-314):** computed on the sweep's contrast, which pooled the 0.3g
  (Mars) arm with 0g. Corrected contrasts are in S54–S56, and the diagnostics and GO replication re-run
  on them are in S58–S63 (reproduction check S64).
- **S29/S33/S34 (OSD-218):** the v8 run pooled 16 Col-0 and 16 WS samples. Genotype-specific
  re-runs are in S48/S49.
- **S40/S41 (OSD-406):** a pooled run across Col-0/WS/sku5 and two suborbital rockets.
- **S35–S37 (OSD-624):** a Virgin Galactic suborbital flight, not an orbital dataset.
- **S42–S43 (power) and S44 module rows:** the pool mixes OSD-218 Col-0 and WS. The Col-0-only re-runs are in S65–S67 (check S68).

## Known gaps
- Peak files (MACS2 narrowPeak) for the SOG1 ChIP-seq were not retained. See
  [`../code/hpc/sog1_alignment.md`](../code/hpc/sog1_alignment.md) for how to regenerate them.
