# references.bib verification report (2026-10-07)

All 24 entries in `references.bib` were built from live Crossref `/works/{DOI}` records. Each DOI was found with a Crossref bibliographic query, and its BibTeX was also pulled through doi.org content negotiation. Title, first author and year were checked against what was intended. Year is the print/issue year, not the online-first year. The only non-Crossref field is the Bourbousse 2018 issue and pages, which came from PubMed esummary (PMID 30541889). Abstracts are not included. Some author names contain UTF-8 accents (JASPAR, MACS, Villacampa), so compile with UTF-8 input (biber, or inputenc utf8).

## Requested references

| # | Item | Key | DOI | Notes |
|---|---|---|---|---|
| 1 | PCSD | `liu2018pcsd` | 10.1093/nar/gkx919 | NAR 46(D1):D1157-67, 2018 (online 2017) |
| 2 | AraENCODE | `wang2023araencode` | 10.1016/j.molp.2023.06.005 | Mol Plant 16:1113-16. A bioRxiv preprint also exists (10.1101/2023.06.10.544382), but the journal version is cited |
| 3 | PlantCADB | `ding2023plantcadb` | 10.1016/j.gpb.2022.10.005 | GPB 21:311-23, **2023** print (online 2022) |
| 4 | OSDR / GeneLab | `gebre2025osdr` (OSDR), `berrios2021genelab` (GeneLab) | 10.1093/nar/gkae1116; 10.1093/nar/gkaa887 | Gebre et al. NAR 53:D1697, 2025 is the OSDR paper. Berrios et al. NAR 49:D1515, 2021 is also included |
| 5 | Paul 2017 CARA | `paul2017genetic` | 10.1371/journal.pone.0180186 | Matches OSD-120 metadata |
| 6 | SOG1 ChIP-seq | `bourbousse2018sog1` | 10.1073/pnas.1810582115 | PNAS 115(52):E12453-62. **See discrepancy 1** |
| 7 | Yoshiyama 2009 | `yoshiyama2009sog1` | 10.1073/pnas.0810304106 | PNAS 106:12843-48 |
| 8 | Ws-0 genome | Zhou et al. 2019 BMC Genomics, 10.1186/s12864-019-5554-z | `zhou2019epigenomics` | **RESOLVED by author decision (2026-10-07)**: Ws (Wassilewskija) processing cited to the OSD-217 paper; no separate assembly citation |
| 9 | JASPAR 2024 | `rauluseviciute2024jaspar` | 10.1093/nar/gkad1059 | NAR 52:D174-82 |
| 10 | Gene Ontology | `go2023knowledgebase` | 10.1093/genetics/iyad031 | Genetics 224(1):iyad031, 2023. A "GO knowledgebase in 2026" paper also exists (10.1093/nar/gkaf1292) if a newer citation is wanted. It was not added |
| 11 | ChromHMM | `ernst2012chromhmm` | 10.1038/nmeth.1906 | Nat Methods 9:215-16 |
| 12 | MACS | `zhang2008macs` | 10.1186/gb-2008-9-9-r137 | |
| 13 | HISAT2 | `kim2019hisat2` | 10.1038/s41587-019-0201-4 | Nat Biotechnol 37:907-15 |
| 14 | Benjamini-Hochberg | `benjamini1995fdr` | 10.1111/j.2517-6161.1995.tb02031.x | JRSS-B 57:289-300 |
| 15 | OSD-217 methylome | `zhou2019epigenomics` | 10.1186/s12864-019-5554-z | Confirmed: BMC Genomics 20:205, 2019. First author Zhou M. Matches OSD-217 metadata |
| 17 | TAIR | `berardini2015tair` | 10.1002/dvg.22877 | Genesis 53:474-85. Lamesch et al. 2012 (10.1093/nar/gkr1090) was also verified but not added |
| 18 | bsmap / MOABS | `xi2009bsmap`, `sun2014moabs` | 10.1186/1471-2105-10-232; 10.1186/gb-2014-15-2-r38 | Both are named in the OSD-217 protocol text |

## OSDR datasets (from `https://osdr.nasa.gov/osdr/data/osd/meta/{n}`, publications then checked on Crossref)

| OSD | Study title (OSDR) | Publication(s) + DOI | Key |
|---|---|---|---|
| 37 | Comparison of the spaceflight transcriptome of four commonly used Arabidopsis thaliana ecotypes | Choi et al. 2019 Am J Bot 106:123, 10.1002/ajb2.1223; Barker et al. 2020 Front Plant Sci 11:147, 10.3389/fpls.2020.00147 | `choi2019variation`, `barker2020tost` |
| 120 | Genetic dissection of the Arabidopsis spaceflight transcriptome... | Paul et al. 2017 PLoS ONE, 10.1371/journal.pone.0180186 | `paul2017genetic` |
| 193 | During development, the Sku6 mutant roots engage different genes than wild type Col-0 roots... | Califar et al. 2020 Front Plant Sci 11:239, 10.3389/fpls.2020.00239 | `califar2020root` |
| 218 | Plant development on ISS differs from the development on the ground and is influenced by the genetic background | **None listed in OSDR** (publications array empty) | -- (cite dataset as OSD-218 only) |
| 314 | Adaptive response of Arabidopsis seedlings in microgravity and Mars reduced gravity... red light | Villacampa et al. 2021 IJMS 22:899, 10.3390/ijms22020899 | `villacampa2021mars` |
| 406 | Transcriptomic responses of Arabidopsis WS and sku5 to Blue Origin NS-12 and Virgin Galactic VP-03 | Califar et al. 2021 Grav Space Res 9:13-29, 10.2478/gsr-2021-0002 | `califar2021suborbital` |
| 496 | The SOG1 transcriptional activator and the MyB3R family of repressors... [SOG1 ChIP-seq] | Bourbousse et al. 2018 PNAS, 10.1073/pnas.1810582115 | `bourbousse2018sog1` |
| 624 | Arabidopsis transcriptome in Virgin Galactic human-tended suborbital spaceflight | Ferl et al. 2023 npj Microgravity 9:95, 10.1038/s41526-023-00340-w | `ferl2023virgin` |
| 217 | Characterization of Epigenetic Regulation in an Extraterrestrial Environment: The Arabidopsis Spaceflight Methylome | Zhou et al. 2019 BMC Genomics, 10.1186/s12864-019-5554-z | `zhou2019epigenomics` |

## Discrepancies with the draft

1. **GSE107980 is the wrong GEO accession.** GEO lists GSE107980 as a colorectal-cancer RNF6 RNA-seq series (PMID 29288235, Clin Cancer Res). OSD-496 metadata points to **GSE112529** ("...[SOG1 ChIP-seq]"). That series is a SubSeries of GSE112773, and GEO links GSE112773 to PMID 30541889, which is Bourbousse et al. 2018 PNAS. The paper is right, but the draft's accession must be changed to GSE112529 (or GSE112773 for the SuperSeries).
2. **The "Mott et al. 2011" Ws genome citation is not supported.** No 2011 paper with Mott as first author turned up. The likely intended paper is Gan X, ..., Mott R (2011) "Multiple reference genomes and transcriptomes for Arabidopsis thaliana", Nature, 10.1038/nature10414 (Crossref-verified). Mott is the *last* author, and the first author is Gan. However, see the unresolved item below.
3. PlantCADB print year is 2023, not 2022 (online 2022). PCSD print year is 2018, as the draft says (online 2017). JASPAR 2024 is a 2024 print (online 2023). Gebre OSDR is a 2025 print (online 2024). Berrios GeneLab is a 2021 print (online 2020).
4. The OSD-217 author protocol in OSDR says sites need coverage of 20 in at least 2 of 3 replicates. The published Zhou 2019 methods say at least 10x. Worth checking whichever the manuscript quotes.

## UNRESOLVED

- **Ws-0 reference genome ("ws_0.v7")**. What I tried:
  - The OSD-217 OSDR protocol and the full text of Zhou 2019 (Europe PMC, PMC6416986). Both name the `ws_0.v7.allPlusChlMito` cscall index and say only "(WS) reference genome", with **no citation** for the assembly.
  - Crossref searches for a Ws-0 assembly and for "Mott 2011". These found Gan et al. 2011 Nature (the 19-genomes project, which includes Ws-0) and later chromosome-scale ecotype assemblies, but nothing ties "v7" to any of them.
  - The 19-genomes download server (mus.well.ox.ac.uk/19genomes) did not respond, so I could not check its version labels.

  Gan et al. 2011 (10.1038/nature10414) is the most plausible source, but it was **not added** to the .bib because the link is unconfirmed. Options: confirm with the UF ICBR / Riva group (cscall authors) which Ws-0 build "v7" is, or describe the reference only as "the Ws-0 index distributed with cscall" without an assembly citation.
