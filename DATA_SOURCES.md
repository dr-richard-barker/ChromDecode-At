# Data sources

ChromDecode-At redistributes **no source data**. Everything the scripts read
comes from public resources and lives under `$CHROMDECODE_ROOT/data/`. If
`CHROMDECODE_ROOT` isn't set, it defaults to the repository root. This page
lists each input, where it comes from, and the path the code expects.

"Fetched by code" means a script downloads the file itself. "Manual" means it
has to be downloaded before the scripts are run. The analysis environment
cached these files, so no downloader for them is in the archive.

## Chromatin-state references (manual)

| Resource | Expected path | Notes |
|---|---|---|
| PCSD 36-state ChromHMM segments + state descriptions (Liu et al. 2018) | `data/pcsd/` (segment files) and `data/pcsd/state_descriptions.tsv` | systemsbiology.cau.edu.cn/chromstates. The server is unstable, so cache the downloads. 290,553 segments |
| AraENCODE 12-state seedling track | `data/araencode/seedling_12_dense.redefinecolor.bed.lst.gz` (+ `.tbl`) | glab.hzau.edu.cn/AraENCODE |
| AraENCODE gene annotation | `data/araencode/genes_TAIR10.txt` | TSS-anchored coordinates for every gene |
| AraENCODE expression matrix | `data/araencode/expression_TPM.txt` | Used for the learned layer (v1 transfer, v2 training) |
| TAIR10 functional descriptions | `data/araencode/TAIR10_functional_descriptions_20130831.txt` | Gene-identity annotation |
| PlantCADB accessible chromatin regions | `data/plantcadb/Whole_ACR_AssociatedGenes.bed.gz` | bioinfor.nefu.edu.cn/PlantCADB. The 663 MB whole-genome file has a malformed header, and `chromdecode.assign_acr_overlap` parses it tolerantly |

## Motif and GO inputs (fetched by `code/prep_motif_go_inputs.py`)

| Resource | Expected path | Source |
|---|---|---|
| JASPAR 2024 plant PFMs (907 usable) | `data/motif_go/pfms/*.jaspar` | `https://jaspar.elixir.no/api/v1/matrix` |
| TAIR10 genome | `data/motif_go/TAIR10.fa.gz` | Ensembl Plants release 57 |
| GO annotations | `data/motif_go/tair.gaf.gz`, `goa_arabidopsis.gaf` | Gene Ontology / QuickGO |
| GO ontology | `data/motif_go/go.obo` | `https://purl.obolibrary.org/obo/go.obo` |

## NASA OSDR expression data (fetched by code)

`run_osdr_sweep.py`, `root_replication*.py` and `power_analysis.py` query the
OSDR biodata API (`https://visualization.osdr.nasa.gov/biodata/api/v2`) and
download GeneLab-processed counts into `results/sweep/`. The v1 OSD-37 run
reads `data/osdr/GLDS-37_rna_seq_VST_Counts_GLbulkRNAseq.csv`.

Accessions the sweep enumerated (Table S10): OSD-37, 38, 120, 193, 208, 217, 218, 219, 223,
251, 281, 313, 314, 321, 346, 406, 411, 427, 437, 476, 480, 498, 502, 508, 510,
518, 519, 522, 565, 624, 658, 678, 782. Root replication (Table S27) used
OSD-120, 193, 218, 406 and 624.

Sample metadata for OSD-193, 218, 281 and 406 was recovered from NCBI
BioSample through E-utilities (`recover_sample_meta.py`). The recovered tables
are in `data/derived/`.

## Methylation — OSD-217 WGBS (fetched by code)

`osd217_methylation.py` and `osd217_promoter_methylation.py` download the
author-processed tables from `https://osdr.nasa.gov/geode-py/ws/studies/OSD-217/`
(GSE95594; Ws ecotype; root and leaf; FLT vs GC) into `data/osd217/`.

## SOG1 ChIP-seq — OSD-496 / GEO GSE112529 (SuperSeries GSE112773) (manual + HPC)

These are six SRA runs: SRR6919520–SRR6919524 and SRR6919526. Alignment and peak calling are described in
[`code/hpc/sog1_alignment.md`](code/hpc/sog1_alignment.md). The outputs go in `data/osd496/bam/` and
`data/osd496/peaks/`.
