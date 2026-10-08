# SOG1 ChIP-seq (OSD-496 / GEO GSE112529 (SuperSeries GSE112773)) — alignment and peak calling

This is the provenance for §2.7 / Table S26 / fig33–34. The alignments ran as
Biomni HPC jobs. `worker-0.ipynb` in this folder is the original job log;
`PLAN.md` is the plan from the later v10 session, which happened to be saved
alongside it. The BAMs, SAMs and peak files were **not retained** (≈21 GB of SAM).
Everything below regenerates them from public data.

## Inputs

| Library | SRA run | GEO sample | Role | Reads (SRA) |
|---|---|---|---|---|
| `input_20min` | SRR6919520 | GSM3072264 | input, 20 min after irradiation | 10.6 M |
| `input_1h` | SRR6919521 | GSM3072265 | input, 1 h | 14.9 M |
| `IP_SOG1_20min` | SRR6919522 | GSM3072266 | SOG1-3xFLAG IP, 20 min | 12.8 M |
| `IP_SOG1_1h` | SRR6919523 | GSM3072267 | SOG1-3xFLAG IP, 1 h | 14.8 M |
| `IP_wt_20min` | SRR6919524 | GSM3072268 | wild-type (no FLAG) IP, 20 min | 17.8 M |
| `IP_wt_1h` | SRR6919526 | GSM3072269 | wild-type IP, 1 h | 14.9 M |

Per GEO, all samples are 10-day-old **Col-0** seedlings (pSOG1::SOG1-3xFLAG or wild type), given **75 Gy of Co-60
gamma irradiation** at 3.75 Gy/min and sampled 20 min or 1 h later. ChIP used anti-FLAG M2. Read counts are SRA
`total_spots` (BioProject SRP136806).

In the analysis the FASTQs were the OSDR copies named
`GLDS-496_chip-seq_<SRR>.fastq.gz`. SRR6919525 was not used. The reference is TAIR10
with Ensembl-style chromosome names (`1`..`5`, `Mt`, `Pt`). `sog1_chipseq.py`
renames these to `chr1`..`chr5` when it reads the peaks.

## Commands

Alignment. These are the exact commands from `worker-0.ipynb`, run with HISAT2 2.2.1:

```bash
hisat2-build -p 8 TAIR10.fa TAIR10_index
hisat2 -p 8 --no-spliced-alignment -x TAIR10_index -U SRR6919522.fastq.gz -S IP_SOG1_20min.sam
# ...repeated for each library in the table above
```

Sorting/indexing and peak calling. These ran after the notebook, so they are reconstructed from the parameters in
the report (§16) and the manuscript Methods. Treat them as the specification, not a
verbatim log:

```bash
samtools sort -o bam/IP_SOG1_1h.bam IP_SOG1_1h.sam && samtools index bam/IP_SOG1_1h.bam
# MACS2 2.2.9.1, narrow peaks, IP vs matched-timepoint input
macs2 callpeak -t bam/IP_SOG1_20min.bam -c bam/input_20min.bam -f BAM -g 1.2e8 -q 0.01 -n SOG1_20min --outdir peaks
macs2 callpeak -t bam/IP_SOG1_1h.bam    -c bam/input_1h.bam    -f BAM -g 1.2e8 -q 0.01 -n SOG1_1h    --outdir peaks
# specificity: SOG1-IP vs wild-type IP
macs2 callpeak -t bam/IP_SOG1_20min.bam -c bam/IP_wt_20min.bam -f BAM -g 1.2e8 -q 0.01 -n SOG1vsWT_20min --outdir peaks
macs2 callpeak -t bam/IP_SOG1_1h.bam    -c bam/IP_wt_1h.bam    -f BAM -g 1.2e8 -q 0.01 -n SOG1vsWT_1h    --outdir peaks
```

`code/sog1_chipseq.py` expects these outputs under `$CHROMDECODE_ROOT/data/osd496/`:
`bam/<library>.bam` and `peaks/{SOG1_20min,SOG1_1h,SOG1vsWT_20min,SOG1vsWT_1h}_peaks.narrowPeak`.

Expected outcome (report §16): 56 peaks for SOG1 at 20 min and 100 at 1 h. SOG1-vs-WT gives "160/180 peaks", as the report words it.
