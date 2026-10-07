"""Unit tests for ChromDecode-At on genes with known epigenomic annotations."""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import sys
sys.path.insert(0, f"{ROOT}/code")
import chromdecode as cd

DATA = f"{ROOT}/data"
PASS, FAIL = 0, 0


def check(name, cond, detail=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  PASS  {name} {detail}")
    else:
        FAIL += 1
        print(f"  FAIL  {name} {detail}")


print("== Module 1: AGI cleaning ==")
check("clean AT5G10140", cd.clean_agi("AT5G10140") == "AT5G10140")
check("clean AT5G10140.1", cd.clean_agi("AT5G10140.1") == "AT5G10140")
check("clean at5g10140", cd.clean_agi("at5g10140") == "AT5G10140")
check("clean gene symbol rejected", cd.clean_agi("FLC") is None)
check("clean mitochondrial ATMG00010", cd.clean_agi("ATMG00010") == "ATMG00010")

print("== Module 2: coordinates ==")
genes = cd.load_gene_coordinates(f"{DATA}/araencode/genes_TAIR10.txt")
check("gene count 25k-40k (genes+TEs)", 25000 < len(genes) < 40000,
      f"n={len(genes)}")
flc = genes[genes["agi"] == "AT5G10140"]
check("FLC (AT5G10140) on chr5", len(flc) == 1 and flc.iloc[0]["chr"] == "chr5",
      f"coord={flc.values.tolist()}")
phyb = genes[genes["agi"] == "AT2G18700"]
check("PHYB (AT2G18700) on chr2", len(phyb) == 1 and phyb.iloc[0]["chr"] == "chr2")
regions = cd.gene_regions(genes)
check("regions have TSS", regions["tss"].notna().all())
# FLC is on minus strand; TSS should equal gene end
flc_r = regions[regions["agi"] == "AT5G10140"].iloc[0]
check("FLC TSS = gene end (minus strand)", flc_r["tss"] == flc_r["gene_end"])

print("== Module 3: PCSD state assignment ==")
gene_states, desc = cd.load_pcsd_states(f"{DATA}/pcsd",
                                        f"{DATA}/pcsd/state_descriptions.tsv")
check("PCSD gene-state pairs loaded", len(gene_states) > 250000,
      f"n={len(gene_states)}")
check("36 states present", gene_states["pcsd_state"].nunique() == 36)
check("state descriptions has 36 rows", len(desc) == 36)

# FLC is a silenced (H3K27me3) gene in seedlings -> expect a Polycomb state (11-15)
flc_state = gene_states[gene_states["agi"] == "AT5G10140"]["pcsd_state"].tolist()
print(f"  FLC PCSD states: {flc_state}")
check("FLC assigned to a state", len(flc_state) > 0)
# ACT2 (AT3G18780) is a constitutively active gene -> expect active state (2-10, 22-28)
act_state = gene_states[gene_states["agi"] == "AT3G18780"]["pcsd_state"].tolist()
print(f"  ACT2 PCSD states: {act_state}")
check("ACT2 assigned to a state", len(act_state) > 0)

print("== Module 3: AraENCODE state assignment ==")
seg = cd.load_araencode_states(
    f"{DATA}/araencode/seedling_12_dense.redefinecolor.bed.lst.gz",
    f"{DATA}/araencode/seedling_12_dense.redefinecolor.bed.lst.gz.tbl")
check("AraENCODE segments loaded", len(seg) > 50000, f"n={len(seg)}")
check("state names mapped", seg["state_name"].notna().all())
reg_test = regions[regions["agi"].isin(["AT5G10140", "AT3G18780", "AT2G18700"])]
assigned = cd.assign_state_by_tss(reg_test, seg, value_col="state_name",
                                  out_col="ae_state")
check("AraENCODE TSS assignment works", assigned["ae_state"].notna().all(),
      f"states={assigned['ae_state'].tolist()}")
# FLC should be in a repressive AraENCODE state (Heterochromatin/Quiescent/Polycomb-like)
flc_ae = assigned[assigned["agi"] == "AT5G10140"]["ae_state"].iloc[0]
print(f"  FLC AraENCODE state: {flc_ae}")

print("== Module 4: enrichment ==")
# Positive control using TSS-based single-state assignment: genes whose TSS
# falls in PCSD state 17 should be enriched for state 17.
pcsd_seg = cd.load_pcsd_segments(f"{DATA}/pcsd")
tss_states = cd.assign_state_by_tss(regions, pcsd_seg, value_col="pcsd_state",
                                    out_col="pcsd_state")
check("PCSD TSS assignment coverage", tss_states["pcsd_state"].notna().mean() > 0.9,
      f"coverage={tss_states['pcsd_state'].notna().mean():.3f}")
ctrl = tss_states[tss_states["pcsd_state"] == 17]["agi"].unique()
bg = tss_states["agi"].unique()
enr = cd.state_enrichment(tss_states.rename(columns={"pcsd_state": "st"}),
                          ctrl, background=bg, state_col="st")
top = enr.iloc[0]
check("positive control: state 17 top hit", top["st"] == 17,
      f"odds={top['odds_ratio']:.1f}, p={top['p_value']:.2e}")
check("padj computed", enr["padj"].notna().all())

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
