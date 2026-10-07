"""Recover sample metadata for OSD-193/218/281/406 from NCBI Biosample."""
import os as _os
# Project root (holds code/, data/, results/, figures/); override with CHROMDECODE_ROOT.
ROOT = _os.environ.get("CHROMDECODE_ROOT",
                       _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import pandas as pd

OUT = f"{ROOT}/results/root_replication"
SWEEP = f"{ROOT}/results/sweep"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def fetch(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2)

def biosample_title(acc):
    uid = ET.fromstring(fetch(
        f"{EUTILS}/esearch.fcgi?db=biosample&term={acc}"))
    ids = [e.text for e in uid.iter("Id")]
    if not ids:
        return None
    time.sleep(0.4)
    summ = ET.fromstring(fetch(
        f"{EUTILS}/esummary.fcgi?db=biosample&id={ids[0]}"))
    title = None
    for d in summ.iter("DocumentSummary"):
        t = d.find("Title")
        if t is not None and t.text:
            title = t.text
    time.sleep(0.4)
    return title

def parse_title(title):
    """'Spaceflight, 8 days old, Col-0 Rep4' -> structured fields."""
    if not title:
        return None
    cond = ("Spaceflight" if re.search(r"space\s*flight|spaceflight|flight|ISS",
                                       title, re.I)
            else "Ground Control" if re.search(r"ground|control", title, re.I)
            else None)
    age = re.search(r"(\d+)\s*days? old", title, re.I)
    rep = re.search(r"Rep\s*(\d+)", title, re.I)
    geno = ("Col-0" if re.search(r"col\s*[-]?0", title, re.I)
            else "sku6" if re.search(r"sku\s*[-]?6", title, re.I) else None)
    return {"title": title, "condition": cond,
            "age_days": int(age.group(1)) if age else None,
            "genotype": geno, "rep": int(rep.group(1)) if rep else None}

results = {}
for acc in ["OSD-193", "OSD-218", "OSD-281", "OSD-406"]:
    counts = pd.read_csv(f"{SWEEP}/{acc}_counts.csv", index_col=0, nrows=1)
    samples = [c for c in counts.columns
               if c.startswith(("GSM", "SAMN", "SRS"))]
    print(acc, len(samples), "samples")
    rows = []
    for s in samples:
        t = biosample_title(s)
        p = parse_title(t)
        rows.append({"sample": s, **(p or {"title": t})})
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/{acc}_sample_meta.csv", index=False)
    results[acc] = df
    print(df["condition"].value_counts(dropna=False).to_dict(),
          "| ages:", sorted(df["age_days"].dropna().unique()),
          "| genotypes:", df["genotype"].dropna().unique().tolist())
    time.sleep(1)
print("done")
