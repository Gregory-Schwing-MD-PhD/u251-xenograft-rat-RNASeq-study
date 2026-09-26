# -*- coding: utf-8 -*-
"""Annotate the differentially expressed genes with gene type and official name.

The DESeq2 tables carry only Ensembl IDs and symbols, and the GTF used for
quantification is not in this repository, so gene class ("is this locus actually
protein coding?") is looked up from MyGene.info rather than typed from memory.
Writes figures/gene_annotation.json for the deck builder to read.
"""
import io
import json
import os
import urllib.parse
import urllib.request

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
XLSX = os.path.join(ROOT, "manuscript", "Supplementary_Data.xlsx")
OUT = os.path.join(HERE, "figures", "gene_annotation.json")

# MyGene's type_of_gene vocabulary -> the phrase used on a slide
READABLE = {
    "protein-coding": "protein coding",
    "ncRNA": "long non-coding RNA",
    "snRNA": "small nuclear RNA",
    "snoRNA": "small nucleolar RNA",
    "rRNA": "ribosomal RNA",
    "pseudo": "pseudogene",
    "biological-region": "unclassified locus",
    "unknown": "unclassified locus",
    "other": "unclassified locus",
}


def query(ids, scopes):
    body = urllib.parse.urlencode({
        "q": ",".join(ids), "scopes": scopes,
        "fields": "symbol,name,type_of_gene", "species": "human"}).encode()
    req = urllib.request.Request("https://mygene.info/v3/query", data=body,
                                 headers={"Content-Type":
                                          "application/x-www-form-urlencoded"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


de = pd.read_excel(XLSX, "S2_DE_significant")
de["ens"] = de.gene_id.astype(str).str.split(".").str[0]

ann = {}
for hit in query(de.ens.tolist(), "ensembl.gene"):
    q, sym = hit.get("query"), hit.get("symbol")
    if hit.get("notfound") or not sym:
        continue
    ann[q] = {"symbol": sym, "name": hit.get("name", ""),
              "kind": hit.get("type_of_gene", "unknown")}

# anything Ensembl-unmatched: try again by symbol before giving up
missing = [r.ens for _, r in de.iterrows() if r.ens not in ann]
by_sym = {r.ens: r.symbol for _, r in de.iterrows()}
if missing:
    for hit in query([by_sym[m] for m in missing], "symbol"):
        for ens, sym in by_sym.items():
            if ens in missing and hit.get("query") == sym and hit.get("symbol"):
                ann[ens] = {"symbol": hit["symbol"], "name": hit.get("name", ""),
                            "kind": hit.get("type_of_gene", "unknown")}

out = {}
for _, r in de.iterrows():
    a = ann.get(r.ens, {"symbol": r.symbol, "name": "", "kind": "unknown"})
    out[str(r.symbol)] = {
        "lfc": round(float(r.log2FoldChange), 2),
        "padj": float(r.padj),
        "name": a["name"],
        "kind": a["kind"],
        "kind_label": READABLE.get(a["kind"], a["kind"]),
        "coding": a["kind"] == "protein-coding",
    }

json.dump(out, io.open(OUT, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
n_cod = sum(1 for v in out.values() if v["coding"])
print("gene_annotation.json: %d genes, %d protein coding, %d not"
      % (len(out), n_cod, len(out) - n_cod))
print("\ntop |log2FC| and how each locus is classified:")
for sym, v in sorted(out.items(), key=lambda kv: -abs(kv[1]["lfc"]))[:10]:
    print("  %-16s %+7.2f  %-20s %s" % (sym, v["lfc"], v["kind_label"],
                                        v["name"][:46]))
