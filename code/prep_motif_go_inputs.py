"""Download reference data for motif/GO enrichment:
TAIR10 genome FASTA, Arabidopsis GO (QuickGO GAF + OBO), JASPAR plant PFMs.
"""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import concurrent.futures as cf
import io
import json
import os
import sys
import urllib.request

DATA = f"{ROOT}/data"
MOTIF = f"{DATA}/motif_go"
os.makedirs(MOTIF, exist_ok=True)


def fetch(url, dest=None, timeout=600):
    req = urllib.request.Request(url, headers={"User-Agent": "python"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = r.read()
    if dest:
        open(dest, "wb").write(data)
    return data


# 1. Genome FASTA
fa = f"{MOTIF}/TAIR10.fa.gz"
if not os.path.exists(fa):
    print("Downloading TAIR10 genome ...")
    fetch("https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/"
          "fasta/arabidopsis_thaliana/dna/"
          "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz", fa)
print(f"genome: {os.path.getsize(fa):,} bytes")

# 2. GO annotations + ontology
gaf = f"{MOTIF}/goa_arabidopsis.gaf"
if not os.path.exists(gaf):
    print("Downloading Arabidopsis GO annotations (QuickGO) ...")
    req = urllib.request.Request(
        "https://www.ebi.ac.uk/QuickGO/services/annotation/"
        "downloadSearch?taxonId=3702&geneNodeType=protein",
        headers={"Accept": "text/gaf", "User-Agent": "python"})
    with urllib.request.urlopen(req, timeout=300) as r:
        data = r.read()
    open(gaf, "wb").write(data)
print(f"GAF: {os.path.getsize(gaf):,} bytes")

obo = f"{MOTIF}/go.obo"
if not os.path.exists(obo):
    print("Downloading GO ontology ...")
    fetch("https://purl.obolibrary.org/obo/go.obo", obo, timeout=300)
print(f"OBO: {os.path.getsize(obo):,} bytes")

# 3. JASPAR plant PFMs (latest version per base_id)
pfm_dir = f"{MOTIF}/pfms"
os.makedirs(pfm_dir, exist_ok=True)
if len(os.listdir(pfm_dir)) < 100:
    print("Fetching JASPAR plant matrix list ...")
    base = "https://jaspar.elixir.no/api/v1/matrix"
    matrices = []
    page = f"{base}/?tax_group=plants&collection=CORE&page_size=500"
    while page:
        d = json.loads(fetch(page, timeout=120))
        matrices += [(r["matrix_id"], r["base_id"], r["name"], r["version"])
                     for r in d["results"]]
        page = d.get("next")
    print(f"  {len(matrices)} matrices total")
    # latest version per base_id
    latest = {}
    for mid, bid, name, ver in matrices:
        if bid not in latest or ver > latest[bid][1]:
            latest[bid] = (mid, ver, name)
    print(f"  {len(latest)} unique TFs (latest versions)")

    def get_pfm(item):
        bid, (mid, ver, name) = item
        dest = f"{pfm_dir}/{mid}.jaspar"
        if os.path.exists(dest):
            return mid, "cached"
        try:
            data = fetch(f"{base}/{mid}/?format=jaspar", timeout=60)
            if b"No PFM data" in data or len(data) < 50:
                return mid, "no_data"
            open(dest, "wb").write(data)
            return mid, "ok"
        except Exception as e:
            return mid, f"fail:{e}"

    ok = fail = 0
    with cf.ThreadPoolExecutor(max_workers=12) as ex:
        for mid, status in ex.map(get_pfm, sorted(latest.items())):
            if status in ("ok", "cached"):
                ok += 1
            else:
                fail += 1
    print(f"  PFMs: {ok} ok, {fail} failed")
print(f"PFM dir: {len(os.listdir(pfm_dir))} files")
print("Done.")
