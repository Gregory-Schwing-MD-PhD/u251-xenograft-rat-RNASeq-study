#!/usr/bin/env python
"""A tissue-complete host signature: Zhang 2014 non-myeloid brain cell types + Bowman 2016 myeloid states.

Why: the Bowman 2016 signature alone (four myeloid states) cannot describe bulk rat brain. On its 189 genes in the
host mixture every sample correlates with every myeloid column at |r| < 0.08 (the three controls negatively with all
four), because the genes that make the sorted-microglia column distinct include transcripts the sort carried from
neighbouring cells (S100B, NRGN, KIF5A, NEFL, MLC1, TPPP), and in whole tissue those come from astrocytes, neurons
and oligodendrocytes. CIBERSORT v1.04 then either fits nothing (best r 0.02-0.08 in tumours, below 0 in controls) or
stops: a permutation draw no nu-SVR can fit with a positive weight makes CoreAlg divide 0 by 0, which.min() of three
NaN RMSEs is empty, and out[[mn]] ends the run.

Here the myeloid states are fitted alongside the cells that make up the rest of the tissue. Zhang's microglia column
is left out (Bowman's normal microglia stand for it; keeping both would put two near-copies of one cell type from
two batches in one signature). Both mean matrices are mapped mouse -> human by 1:1 orthologs exactly as in
host_deconv_prep.py, joined on the genes both references measured, and reduced with host_deconv_prep.signature(), so a
myeloid marker is kept only if it is >= 2x every Zhang non-myeloid cell type too; that is what removes the carried
neuronal and glial transcripts. The two references differ in platform (Zhang: Cufflinks FPKM -> TPM; Bowman: STAR
counts -> TPM), so a residual batch difference between the two blocks of columns is part of this signature.

    python build_combined_signature.py --zhang zhang2014_brain_mean7_tpm.txt --bowman bowman2016_myeloid_mean4.txt \
        --mouse2human mouse_human_1to1.tsv --mixture mixture_host_hs.txt --out results
"""
import argparse
import json
from pathlib import Path

import pandas as pd

from host_deconv_prep import one_to_one, signature, to_human

ZHANG_KEEP = ["astrocyte", "neuron", "OPC", "newly_formed_oligodendrocyte", "myelinating_oligodendrocyte", "endothelial"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zhang", required=True)
    ap.add_argument("--bowman", required=True)
    ap.add_argument("--mouse2human", required=True)
    ap.add_argument("--mixture", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--name", default="zhang6_bowman4")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    m2h = one_to_one(a.mouse2human, "mouse")
    read = lambda f: pd.read_csv(f, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
    z = to_human(read(a.zhang), m2h)[ZHANG_KEEP]
    b = to_human(read(a.bowman), m2h)
    genes = z.index.intersection(b.index)
    mean = pd.concat([z.loc[genes], b.loc[genes]], axis=1)
    sig, info = signature(mean)
    sig.index.name = "GeneSymbol"
    p = out / f"signature_host_{a.name}.txt"
    sig.to_csv(p, sep="\t", float_format="%.4f")
    mix = pd.read_csv(a.mixture, sep="\t", index_col=0)
    info.update({"zhang_genes_human": int(len(z)), "bowman_genes_human": int(len(b)), "shared_genes": int(len(genes)),
                 "columns": list(sig.columns), "genes_in_mixture": int(sig.index.isin(mix.index).sum()),
                 "genes_per_column": sig.idxmax(axis=1).value_counts().to_dict()})
    carried = ["S100B", "NRGN", "KIF5A", "NEFL", "MLC1", "TPPP", "COL3A1"]
    info["carried_transcripts_still_in_signature"] = [g for g in carried if g in sig.index]
    (out / f"signature_host_{a.name}_log.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))
    print(p.resolve())


if __name__ == "__main__":
    main()
