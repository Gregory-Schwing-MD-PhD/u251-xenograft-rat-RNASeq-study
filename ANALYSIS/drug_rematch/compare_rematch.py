#!/usr/bin/env python
"""Old ChEMBL matching (first search hit) against exact-name matching, for each arm of the drug ranking.

  published   September ranking (subtypes/drug_ranking_final.csv), today's rerun with the old matching (the hold-out
              control arm, inputs/IL68B_control_ranking_final.csv), and today's rerun with exact matching
  IL68B       hold-out arm with the old matching (inputs/IL68B_holdout_ranking_final.csv) vs exact matching
  IL66B       the same for IL66B
The difference between 'old today' and 'exact' is the matching alone: same software, same day, same inputs.
Writes REMATCH.md and rematch.json next to this script.
"""
import json
from pathlib import Path

import pandas as pd

R = Path(__file__).resolve().parent
U = R.parents[1]
ARMS = {"published": {"September": U / "subtypes/drug_ranking_final.csv", "old matching": R / "inputs/IL68B_control_ranking_final.csv",
                      "exact matching": R / "runs/published/subtypes/drug_ranking_final.csv"},
        "IL68B": {"old matching": R / "inputs/IL68B_holdout_ranking_final.csv", "exact matching": R / "runs/IL68B/subtypes/drug_ranking_final.csv"},
        "IL66B": {"old matching": R / "inputs/IL66B_holdout_ranking_final.csv", "exact matching": R / "runs/IL66B/subtypes/drug_ranking_final.csv"}}
WATCH = ["ciclopirox", "deferoxamine", "pyrantel", "vandetanib", "mefloquine", "bromocriptine", "pentetrazol", "ifosfamide"]


def summary(p):
    d = pd.read_csv(p)
    d["key"] = d.Drug.str.strip().str.lower()
    nob = d.sort_values("score_nobbb", ascending=False)
    out = {"n_clinical": int(len(d)), "n_both_agree": int(d.both_agree.astype(bool).sum()),
           "top10_bbb": d.sort_values("rank_bbb").Drug.head(10).tolist(), "top10_nobbb": nob.Drug.head(10).tolist(),
           "table": d.sort_values("rank_bbb").head(15)[["Drug", "NES", "BBB_Martins", "both_agree", "phase", "rank_bbb", "rank_nobbb", "score_bbb", "score_nobbb"]]
           .assign(ChEMBL_ID=d.sort_values("rank_bbb").head(15).ChEMBL_ID).to_dict("records"),
           "watch": {}}
    for w in WATCH:
        r = d[d.key == w]
        out["watch"][w] = None if r.empty else {"rank_bbb": int(r.rank_bbb.iloc[0]), "rank_nobbb": int(r.rank_nobbb.iloc[0]),
                                                "NES": float(r.NES.iloc[0]), "BBB_Martins": float(r.BBB_Martins.iloc[0]),
                                                "both_agree": bool(r.both_agree.iloc[0]), "phase": float(r.phase.iloc[0]),
                                                "score_bbb": float(r.score_bbb.iloc[0]), "score_nobbb": float(r.score_nobbb.iloc[0])}
    out["_set"] = set(d.key)
    return out


res, md = {}, ["# Drug ranking with exact-name ChEMBL matching\n",
               "Old matching = the first hit of ChEMBL's free-text search (the published chain). Exact matching = "
               "`resolve_chembl.py`: pref_name or synonym equal to the DSigDB name, mapped to the parent molecule, audited "
               "(`chembl_mapping.csv`, `overrides.tsv`). Everything downstream is the September chain.\n"]
for arm, runs in ARMS.items():
    res[arm] = {}
    md.append(f"\n## {arm}\n")
    md.append("| run | clinical compounds | both barrier models | ciclopirox rank (weighted / unweighted) | deferoxamine rank |\n|---|---|---|---|---|\n")
    sets = {}
    for lab, p in runs.items():
        if not p.exists():
            md.append(f"| {lab} | missing: {p.name} | | | |\n")
            continue
        s = summary(p)
        sets[lab] = s.pop("_set")
        res[arm][lab] = s
        c, dfo = s["watch"]["ciclopirox"], s["watch"]["deferoxamine"]
        md.append(f"| {lab} | {s['n_clinical']} | {s['n_both_agree']} | "
                  f"{'%d / %d' % (c['rank_bbb'], c['rank_nobbb']) if c else 'absent'} | "
                  f"{'%d / %d (BBB %.2f)' % (dfo['rank_bbb'], dfo['rank_nobbb'], dfo['BBB_Martins']) if dfo else 'absent'} |\n")
    if "old matching" in sets and "exact matching" in sets:
        gained, lost = sorted(sets["exact matching"] - sets["old matching"]), sorted(sets["old matching"] - sets["exact matching"])
        res[arm]["entered"], res[arm]["left"] = gained, lost
        md.append(f"\nEntered the clinical list with exact matching ({len(gained)}): {', '.join(gained) or 'none'}\n\n"
                  f"Left it ({len(lost)}): {', '.join(lost) or 'none'}\n")
    if "exact matching" in res[arm]:
        md.append("\nTop 15 with exact matching (rank by |NES|^1.5 x ADMET-AI BBB):\n\n| rank | compound | NES | BBB | both models | phase | unweighted rank |\n|---|---|---|---|---|---|---|\n")
        for r in res[arm]["exact matching"]["table"]:
            md.append(f"| {r['rank_bbb']} | {r['Drug']} | {r['NES']:+.2f} | {r['BBB_Martins']:.2f} | {'yes' if r['both_agree'] else 'no'} | {r['phase']:.0f} | {r['rank_nobbb']} |\n")
(R / "REMATCH.md").write_text("".join(md), encoding="utf-8")
json.dump(res, open(R / "rematch.json", "w"), indent=1, default=str)
print("".join(md))
