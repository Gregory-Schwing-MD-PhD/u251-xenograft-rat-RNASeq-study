#!/usr/bin/env python
"""Deconvolve the bulk mixture with signature matrices already written by cibersort_neftel.py, one signature per
process, each result written as soon as it exists (the first run held all three in memory until the end, and the
third signature's permutations ran long).

    python deconv_only.py --results results --cohorts ../sample_cohorts.csv --perm 1000 neftel4_confident neftel4_allcells ...
"""
import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cibersort_neftel import deconvolve  # noqa: E402


def run(args):
    name, res, perm, coh = args
    sig = pd.read_csv(Path(res) / f"signature_{name}.txt", sep="\t", index_col=0)
    mix = pd.read_csv(Path(res) / "mixture_u251_tpm.txt", sep="\t", index_col=0)
    d = deconvolve(sig, mix, perm=perm)
    d.insert(1, "cohort", d["sample"].map(coh).fillna(""))
    d.insert(0, "signature", name)
    d.to_csv(Path(res) / f"fractions_nusvr_{name}.tsv", sep="\t", index=False, float_format="%.4f")
    return name, d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--cohorts", required=True)
    ap.add_argument("--perm", type=int, default=1000)
    ap.add_argument("names", nargs="+")
    a = ap.parse_args()
    coh = pd.read_csv(a.cohorts).set_index("sample")["cohort"].to_dict()
    with ProcessPoolExecutor(len(a.names)) as ex:
        for name, d in ex.map(run, [(n, a.results, a.perm, coh) for n in a.names]):
            print("=== " + name, flush=True)
            print(d.to_string(float_format=lambda v: f"{v:.3f}"), flush=True)


if __name__ == "__main__":
    main()
