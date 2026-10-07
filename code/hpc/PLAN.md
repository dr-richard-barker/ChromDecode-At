# Plan: Stratified flight-within-light analysis for OSD-120 Col-0 WT — bounding the flight × light interaction

## Summary
The additive model (condition + light) on OSD-120 Col-0 WT root tips collapsed to
2 down genes; the flagged concern was that flight × light interaction (real per the
CARA paper) gets absorbed into the light term. The design is a balanced 2×2
factorial (FLT/GC × Alight/dark, 3 reps/cell, all Day 13), so the interaction is
directly estimable. Plan: stratified DE + enrichment per light regime, plus a
gene-level interaction model, to bound how much the additive model masked.
Delegated defaults throughout (alpha 0.05 BH, |log2fc| > 1, same pipeline).

## Step 1 — Stratified DE and enrichment
- Within Alight (3v3) and within dark (3v3) separately:
  `cd.differential_expression_any` (condition only; no covariates — strata are
  balanced), padj < 0.05, |log2fc| > 1.
- Polycomb-state enrichment of down genes per stratum (odds, padj) against the
  project state assignments (background fraction 0.141).
- Power context: k=3 per group is ~0-15% powered at root-level odds (v9 power
  analysis) — per-stratum results are bounds, not estimates.

## Step 2 — Interaction bound
- Full 12-sample model: expression ~ condition + light + condition:light.
  Gene-level interaction term: count genes with significant interaction
  (padj < 0.05) and report the |interaction log2fc| distribution (quantiles).
- Module-specific bound: max |interaction| among the 23 redox-module genes and
  the 265 Polycomb-down genes; per-stratum log2FC agreement for module genes
  (are they down in both strata, one, or neither?).
- Enrichment-level bound: difference in Polycomb odds between strata with
  Fisher-based CIs; a shared-direction statement (is flight effect directionally
  consistent across strata?).

## Step 3 — Deliverables
- `results/root_replication/osd120_stratified.csv` (per-stratum DE/enrichment +
  interaction summary rows)
- fig37: per-stratum results (n_down + Polycomb odds per stratum; interaction
  log2fc distribution) — SVG+PNG, media QA
- Report subsection 17.3 (full-file rewrite in /workspace, `rm -f` + `cp` to
  /mnt/results); copy deliverables to chromdecode_at folder; refresh zenodo
  manifest

## Compute / execution target
Sandbox worker-0; two 3v3 DE runs + one 12-sample interaction model on cached
counts — seconds to minutes. No HPC.

## Assumptions (delegated defaults)
- Day is constant (all Day 13) — no day term needed; light is the only stratum
  variable.
- Per-stratum tests are one-sided enrichment, BH within stratum as in the
  pipeline; interaction p-values BH-corrected across genes.
- 3v3 strata are severely underpowered individually; results are interpreted as
  bounds on the interaction, not discovery claims.
