# ChromDecode-At

**Expression-to-chromatin decoding of the *Arabidopsis* spaceflight transcriptome.**

ChromDecode-At maps every *Arabidopsis thaliana* gene to its ground-state chromatin context and asks
whether the genes a treatment changes are concentrated in particular chromatin classes. The context comes from
PCSD 36-state ChromHMM, AraENCODE 12-state and PlantCADB accessible regions. Here it is applied to the
*Arabidopsis* RNA-seq datasets in NASA's Open Science Data Repository (OSDR).

**Project page:** <https://dr-richard-barker.github.io/ChromDecode-At/> ·
**Manuscript:** [`manuscript/`](manuscript/) · **Tables:** [`supplementary_tables/`](supplementary_tables/README.md)

## Findings

Every number below is in a table under `supplementary_tables/`, and the S-numbers are listed in its README.

- **Flight-repressed genes sit in Polycomb chromatin.** In OSD-37 seedlings (529 up / 526 down), down-regulated
  genes are enriched for AraENCODE Repression (odds 8.0, FDR 3e-105) and for PCSD S15, H3K27me3 + accessible DNA
  (odds 18.0) (S11, S12, S51).
- **The bias replicates in ISS-grown roots, in two ecotypes.** The odds are 3.06 for OSD-193 Col-0, and 3.12 for
  Col-0 and 3.74 for WS in OSD-218 (S29, S48).
- **The repressed class is a redox/peroxidase module.** 18 GO terms are significant against the genome and 11
  against a state-matched Polycomb control. Its promoters are AHL/ZHD-rich, but no motif survives the state-matched
  comparison (S18, S19).
- **Not a SOG1 damage response.** No SOG1 ChIP-seq peak falls in a redox-module promoter (0/23) or a
  Polycomb-down promoter (0/265) after bleomycin (OSD-496; S26).
- **No detectable methylation change** at the module in OSD-217 WGBS (minimum padj 0.07 gene body, 0.099
  promoter; S22, S24).
- **Design guidance.** Empirical subsampling suggests about 9–12 replicates per group for a stable root-level
  test (S42, S43).
- **Negative result.** Expression-distribution features do not predict chromatin state (block-CV accuracy
  0.533 vs baseline 0.549; S7).

### Pre-submission audit (2026-10)

Every contrast was re-derived from the OSDR study records and re-run from the public counts with the pinned
environment. The re-run reproduced the archived numbers exactly (S50) and found four sample-coding or labelling
errors, now corrected in the manuscript (Methods; fig38):

| Dataset | Problem | Effect of the correction |
|---|---|---|
| OSD-314 | The Mars-gravity (0.3g) arm was pooled into "0G" | Same direction, much smaller effect: S11 up-gene odds 22.4 → about 2.0 (S54) |
| OSD-218 | 16 WS samples were labelled Col-0 | Polycomb odds unchanged (Col-0 3.12, WS 3.74), so the result now covers two ecotypes (S48) |
| OSD-406 | Pooled Col-0/WS/*sku5* across two suborbital rockets | With covariates, odds 3.17, padj 0.0074, at discovery level only (S48) |
| OSD-624 | Described as orbital; it is a Virgin Galactic suborbital flight | Its null no longer counts against the ISS replications |

A stale 5-group table was replaced (S51), and a star-misalignment bug in `code/aggregate_sweep.py` was fixed. The
figure audit and its regeneration log are in [`manuscript/FIGURE_AUDIT.md`](manuscript/FIGURE_AUDIT.md).

## Layout

```
code/                  analysis scripts; chromdecode.py is the core module
code/hpc/              SOG1 ChIP-seq alignment provenance (Biomni HPC job log + commands)
supplementary_tables/  Tables S1–S57 (CSV) + README index; superseded/ holds one retired table
data/derived/          Tables D1–D7: gene→state map, redox-module genes, recovered sample metadata
figures/               fig1–fig38, SVG + PNG
manuscript/            LaTeX manuscript + supplement, references.bib (Crossref-verified), Makefile
reports/               full per-version project record (v1–v10) and the original v1–v8 draft
docs/                  GitHub Pages site
DATA_SOURCES.md        every external input, its source, and where the code expects it
```

## Reproducing

```bash
python3.11 -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
export CHROMDECODE_ROOT=$PWD          # scripts read data/ and write results/ + figures/ here
python code/qc_sweep_rerun.py         # re-runs the OSDR sweep from public counts and checks it against S10–S12
python code/qc_genotype_rerun.py      # genotype-correct root replication
```

Both QC scripts need only `data/derived/` and network access to OSDR. The full pipeline also needs the PCSD,
AraENCODE, PlantCADB, JASPAR and GO inputs listed in [`DATA_SOURCES.md`](DATA_SOURCES.md). The SOG1 ChIP-seq
needs HISAT2/samtools/MACS2 ([`code/hpc/sog1_alignment.md`](code/hpc/sog1_alignment.md)).

Build the manuscript with `make -C manuscript` (latexmk).

## Licence and citation

Code: MIT ([`LICENSE`](LICENSE)). Tables, figures, derived data and text: CC-BY 4.0
([`LICENSE-DATA`](LICENSE-DATA)). Source data are not redistributed. Cite via [`CITATION.cff`](CITATION.cff); a
Zenodo DOI will be added on first release.

Author: Richard J. Barker ([ORCID 0000-0001-5681-9857](https://orcid.org/0000-0001-5681-9857)),
Center of Space Exploration (CoSE).
