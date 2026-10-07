"""Run the v2 sweep pipeline on Col-0 root RNA-seq candidates (v7 task).

Reuses run_osdr_sweep.process_dataset verbatim; adds Polycomb_repressed
odds/padj extraction for down genes per dataset.
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import io
import json
import sys
import time
import urllib.request

sys.path.insert(0, f"{ROOT}/code")
import pandas as pd
import chromdecode as cd
import run_osdr_sweep as sw

OUT = f"{ROOT}/results/root_replication"
API = "https://visualization.osdr.nasa.gov/biodata/api/v2"
import os
os.makedirs(OUT, exist_ok=True)

# candidates with counts (verified via files endpoint)
CANDS = {
    "OSD-120": "OSD-120_transcription-profiling_rna-sequencing-(rna-seq)_illumina",
    "OSD-208": "OSD-208_transcription-profiling_rna-sequencing-(rna-seq)_illumina",
    "OSD-218": "OSD-218_transcription-profiling_rna-sequencing-(rna-seq)_illumina",
    "OSD-281": "OSD-281_transcription-profiling_rna-sequencing-(rna-seq)_illumina",
    "OSD-406": "OSD-406_transcription-profiling_rna-sequencing-(rna-seq)_illumina",
    "OSD-411": "OSD-411_transcription-profiling_rna-sequencing-(rna-seq)_illumina",
    "OSD-624": "OSD-624_transcription-profiling_rna-sequencing-(rna-seq)_Illumina",
}

# state assignments (same as sweep)
tss_state = pd.read_csv(f"{ROOT}/results/gene_state_assignments.csv")

rows = []
for acc, assay in CANDS.items():
    # reuse sweep checkpoints where present
    sw.SWEEP = f"{ROOT}/results/sweep"
    res = sw.process_dataset(acc, assay, tss_state, None)
    res["accession"] = acc
    rows.append(res)
    print(acc, res.get("status"), res.get("contrast"),
          "up:", res.get("n_de_up"), "down:", res.get("n_de_down"))

    # Polycomb row for down genes
    edn = f"{ROOT}/results/sweep/{acc}_enrichment_down.csv"
    pc = None
    if os.path.exists(edn):
        e = pd.read_csv(edn)
        m = e[e["st"] == "Polycomb_repressed"]
        if len(m):
            pc = m.iloc[0].to_dict()
    elif res.get("status") == "ok" and res.get("n_de_down", 0) > 0:
        # process_dataset wrote enrichment only if n_up>0? check both files
        pass
    if pc is not None:
        print("  Polycomb_down odds:", round(pc["odds_ratio"], 2),
              "padj:", pc["padj"])
    rows[-1]["pc_odds_down"] = pc["odds_ratio"] if pc else None
    rows[-1]["pc_padj_down"] = pc["padj"] if pc else None

pd.DataFrame(rows).to_csv(f"{OUT}/replication_summary.csv", index=False)
print("saved", f"{OUT}/replication_summary.csv")
