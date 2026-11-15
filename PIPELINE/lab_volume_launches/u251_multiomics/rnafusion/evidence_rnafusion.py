#!/usr/bin/env python3
"""Pass-2 evidence for every whitelisted fusion in every library, and the recurrence group claim with its 10-split
relabelling null (standard library only).

Per fusion x library, from FusionInspector (pass 2):
  present         FusionInspector reports it (it lists a fusion from 2 fragments by default, so 'not reported' can
                  still hide 1 fragment)
  absent          not reported AND the 5'-partner expression of that library predicts >= 5 fragments at the carrier
                  arm's mean FFPM:  expected = mean_arm(FFPM) / mean_arm(CPM of the 5' partner) x partner count here
  not assessable  not reported and the prediction is < 5 (or the partner is not expressed)
Group claim: for each of the 10 balanced splits of the six tumours, count fusions present in all 3 of one arm and
assessably absent in all 3 of the other (both directions). The rank of the true split among the 10 is the only
p-value this design allows (floor 0.10).
"""
import argparse, csv, json, os, sys, time
from collections import defaultdict

LIBS = ["C2B", "IL64B", "IL66B", "IL67B", "IL68B", "IL69B", "N168B", "N269B", "NL70B", "NL71B"]
PRIMARY = ["IL67B", "IL68B", "IL69B"]
RECURRENT = ["IL66B", "NL70B", "NL71B"]
TUMOURS = PRIMARY + RECURRENT
CONTROLS = ["IL64B", "N168B"]
GRAFT = {"C2B": 84519445, "IL64B": 314037, "IL66B": 25459283, "IL67B": 31886295, "IL68B": 32225887,
         "IL69B": 49327614, "N168B": 435428, "N269B": 3823923, "NL70B": 31293679, "NL71B": 29210358}
# the ten balanced 3-v-3 splits; IL67B stays in A so no split appears with its mirror; split 1 is the true labelling
SPLITS = [(["IL67B", "IL68B", "IL69B"], ["IL66B", "NL70B", "NL71B"]),
          (["IL67B", "IL68B", "IL66B"], ["IL69B", "NL70B", "NL71B"]),
          (["IL67B", "IL68B", "NL70B"], ["IL69B", "IL66B", "NL71B"]),
          (["IL67B", "IL68B", "NL71B"], ["IL69B", "IL66B", "NL70B"]),
          (["IL67B", "IL69B", "IL66B"], ["IL68B", "NL70B", "NL71B"]),
          (["IL67B", "IL69B", "NL70B"], ["IL68B", "IL66B", "NL71B"]),
          (["IL67B", "IL69B", "NL71B"], ["IL68B", "IL66B", "NL70B"]),
          (["IL67B", "IL66B", "NL70B"], ["IL68B", "IL69B", "NL71B"]),
          (["IL67B", "IL66B", "NL71B"], ["IL68B", "IL69B", "NL70B"]),
          (["IL67B", "NL70B", "NL71B"], ["IL68B", "IL69B", "IL66B"])]
MIN_EXPECTED = 5.0


def sname(lib):
    return lib + "_humclip"


def read_table(path):
    with open(path, newline="") as f:
        header = f.readline().rstrip("\n").lstrip("#").split("\t")
        rows = [dict(zip(header, l.rstrip("\n").split("\t"))) for l in f if l.strip()]
    return header, rows


def to_f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def star_total(path):
    for line in open(path):
        if "Number of input reads" in line:
            return float(line.split("|")[1])
    return 0.0


def gene_names(gene_info):
    """STAR geneInfo.tab (first line = count; then gene_id, name, type) -> list of names in ReadsPerGene order."""
    lines = open(gene_info).read().splitlines()
    return {l.split("\t")[0]: l.split("\t")[1] for l in lines[1:] if "\t" in l}


def partner_counts(rpg, id2name):
    """Reverse-stranded library: ReadsPerGene column 4. Summed over gene ids that share a name. Also returns the
    column 3 and column 4 totals (a reverse-stranded library has column 4 far above column 3)."""
    out = defaultdict(float)
    c3 = c4 = 0.0
    for line in open(rpg):
        p = line.rstrip("\n").split("\t")
        if p[0].startswith("N_") or len(p) < 4:
            continue
        out[id2name.get(p[0], p[0])] += to_f(p[3])
        c3 += to_f(p[2]); c4 += to_f(p[3])
    return out, c3, c4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--d", required=True, help="M/rnafusion")
    ap.add_argument("--gene-info", required=True)
    a = ap.parse_args()
    D, R = a.d, f"{a.d}/results"
    E = f"{D}/evidence"
    os.makedirs(E, exist_ok=True)
    notes = []
    start = os.path.getmtime(f"{D}/.pass2_start") if os.path.isfile(f"{D}/.pass2_start") else 0
    wl = [l.strip() for l in open(f"{D}/whitelist_pass2.tsv") if l.strip()]
    id2name = gene_names(a.gene_info)
    named = sum(1 for k, v in id2name.items() if k != v)
    if named < 0.5 * max(1, len(id2name)):
        notes.append(f"geneInfo.tab names only {named} of {len(id2name)} genes: partner lookups fall back to ids, so "
                     "most fusions will read 'not assessable' rather than 'absent'")
    fi, total, pc, stale = {}, {}, {}, []
    for lib in LIBS:
        s = sname(lib)
        total[lib] = star_total(f"{R}/star/{s}.Log.final.out")
        pc[lib], c3, c4 = partner_counts(f"{R}/star/{s}.ReadsPerGene.out.tab", id2name)
        if lib not in CONTROLS and not c4 > 5 * c3:
            notes.append(f"{lib}: ReadsPerGene column 4 ({c4:.0f}) is not > 5x column 3 ({c3:.0f}); the library was "
                         "verified reverse-stranded (ISR), so the partner counts are suspect")
        p = f"{R}/fusioninspector/{s}/{s}.FusionInspector.fusions.abridged.tsv"
        d = defaultdict(lambda: [0.0, 0.0])
        if os.path.isfile(p):
            if os.path.getmtime(p) < start:
                stale.append(lib)
            for r in read_table(p)[1]:
                d[r.get("FusionName", "")][0] += to_f(r.get("JunctionReadCount")) + to_f(r.get("SpanningFragCount"))
                d[r.get("FusionName", "")][1] += to_f(r.get("FFPM"))
        else:
            notes.append(f"{lib}: no FusionInspector abridged table (FusionInspector reported nothing, or it failed)")
        fi[lib] = d
    if stale:
        print("STALE FusionInspector tables (older than the pass-2 start):", stale)
        open(f"{E}/EVIDENCE_FAILED", "w").write(f"stale FusionInspector tables: {stale}\n")
        return 1

    def cpm5(lib, f):
        g5 = f.split("--")[0]
        return pc[lib].get(g5, 0.0) / (total[lib] / 1e6) if total[lib] else 0.0

    def status(lib, f, carriers):
        """present / absent / not_assessable for lib, given the carrier libraries that set the expectation."""
        frags = fi[lib].get(f, [0.0, 0.0])[0]
        if frags > 0:
            return "present", frags, None
        ffpm = sum(fi[c].get(f, [0, 0])[1] for c in carriers) / len(carriers)
        cp = sum(cpm5(c, f) for c in carriers) / len(carriers)
        g5 = f.split("--")[0]
        if cp > 0:
            exp = ffpm / cp * pc[lib].get(g5, 0.0)
        else:
            exp = ffpm * total[lib] / 1e6 if pc[lib].get(g5, 0.0) > 0 else 0.0
        return ("absent" if exp >= MIN_EXPECTED else "not_assessable"), 0.0, exp

    # artefact flags from the pass-1 audit: (a) or (d) anywhere, (b) or (c) in any library that called it
    flagged = defaultdict(set)
    sc = f"{D}/audit/scored_calls.tsv"
    if os.path.isfile(sc):
        for r in read_table(sc)[1]:
            for k in ("a", "b", "c", "d"):
                if r.get(k) == "1":
                    flagged[r["key"]].add(k)
    for f in wl:
        for lib in CONTROLS:
            if fi[lib].get(f, [0])[0] > 0:
                flagged[f].add("a_fi")

    # ---- evidence matrix (expectation for each library from the six tumours' carriers, for information)
    with open(f"{E}/evidence_matrix.tsv", "w") as o:
        o.write("fusion\tflags\t" + "\t".join(f"{l}_frags\t{l}_ffpm" for l in LIBS) + "\n")
        for f in wl:
            o.write(f"{f}\t{','.join(sorted(flagged[f])) or '-'}\t" +
                    "\t".join(f"{fi[l].get(f, [0, 0])[0]:.0f}\t{fi[l].get(f, [0, 0])[1]:.4g}" for l in LIBS) + "\n")

    # ---- group claim over the 10 splits
    def arm_only(f, X, Y):
        if not all(fi[x].get(f, [0])[0] > 0 for x in X):
            return None
        st = [status(y, f, X) for y in Y]
        if all(s[0] == "absent" for s in st):
            return st
        return None

    results = {}
    for scope in ("all_candidates", "artefact_filtered"):
        cand = [f for f in wl if scope == "all_candidates" or not flagged[f]]
        per_split = []
        for i, (A_, B_) in enumerate(SPLITS, 1):
            a_only = [f for f in cand if arm_only(f, A_, B_)]
            b_only = [f for f in cand if arm_only(f, B_, A_)]
            per_split.append(dict(split=i, A=A_, B=B_, A_only=a_only, B_only=b_only, stat=len(a_only) + len(b_only)))
        s1 = per_split[0]["stat"]
        hi = sum(1 for p in per_split if p["stat"] > s1) + 1
        lo_rank_hi = sum(1 for p in per_split if p["stat"] >= s1)
        results[scope] = dict(n_candidates=len(cand), splits=per_split, rank=(hi, lo_rank_hi),
                              p=max(0.10, lo_rank_hi / 10))
    with open(f"{E}/group_claim_splits.tsv", "w") as o:
        o.write("scope\tsplit\tarm_A\tarm_B\tA_only\tB_only\tstatistic\n")
        for scope, r in results.items():
            for p in r["splits"]:
                o.write(f"{scope}\t{p['split']}\t{'+'.join(p['A'])}\t{'+'.join(p['B'])}\t{','.join(p['A_only']) or '-'}"
                        f"\t{','.join(p['B_only']) or '-'}\t{p['stat']}\n")

    # ---- per-library status for the true split's recurrence-only and primary-only candidates
    rec_only = results["all_candidates"]["splits"][0]["B_only"]
    pri_only = results["all_candidates"]["splits"][0]["A_only"]
    with open(f"{E}/true_split_candidates.tsv", "w") as o:
        o.write("fusion\tdirection\tflags\tin_C2B\tin_IL64B\tin_N168B\tin_N269B\t" +
                "\t".join(f"{l}_frags\t{l}_expected_if_absent" for l in TUMOURS) + "\n")
        for f, dirn, X in [(f, "recurrence_only", RECURRENT) for f in rec_only] + \
                          [(f, "primary_only", PRIMARY) for f in pri_only]:
            cells = []
            for l in TUMOURS:
                s = status(l, f, X)
                cells.append(f"{fi[l].get(f, [0])[0]:.0f}\t{'' if s[2] is None else f'{s[2]:.1f}'}")
            o.write(f"{f}\t{dirn}\t{','.join(sorted(flagged[f])) or '-'}\t" +
                    "\t".join("yes" if fi[l].get(f, [0])[0] > 0 else "no" for l in ("C2B", "IL64B", "N168B", "N269B"))
                    + "\t" + "\t".join(cells) + "\n")

    # ---- N269B at clonal level: DHX38--AMPD1 per million graft pairs, against the tumours and the controls
    k = "DHX38--AMPD1"
    rate = {l: fi[l].get(k, [0])[0] / (GRAFT[l] / 1e6) for l in LIBS}
    md = ["# rnafusion pass-2 evidence and the recurrence group claim", "",
          f"Written {time.strftime('%Y-%m-%d %H:%M')}. FusionInspector (pass 2, whitelist of {len(wl)} fusions) in all ten "
          "libraries. 'Absent' = not reported by FusionInspector AND >= 5 fragments predicted from the 5' partner's "
          "expression at the carrier arm's mean FFPM; otherwise 'not assessable'. FusionInspector does not list a "
          "fusion below its own floor (2 fragments by default), so 'not reported' can hide one fragment.", ""]
    for scope, r in results.items():
        s = r["splits"][0]
        md.append(f"- **{scope}** ({r['n_candidates']} fusions): true split has {len(s['B_only'])} recurrence-only and "
                  f"{len(s['A_only'])} primary-only fusions (statistic {s['stat']}); rank {r['rank'][0]}"
                  + (f"-{r['rank'][1]}" if r['rank'][1] != r['rank'][0] else "") + f" of 10, p = {r['p']:.2f} "
                  "(floor 0.10). Per split: " + ", ".join(str(p["stat"]) for p in r["splits"]) + ".")
    md += ["", "Recurrence-only in the true split (all fusions; flags from the pass-1 audit, a_fi = reported by "
           "FusionInspector in IL64B or N168B):"]
    md += [f"- {f}: flags {','.join(sorted(flagged[f])) or 'none'}; in C2B "
           f"{'yes (pre-existing in culture)' if fi['C2B'].get(f, [0])[0] > 0 else 'no'}" for f in rec_only] or ["- none"]
    md += ["", f"## N269B and {k}", ""]
    md += [f"- {l}: {fi[l].get(k, [0])[0]:.0f} fragments, {rate[l]:.2f} per million graft pairs" for l in
           ["C2B"] + TUMOURS + ["N269B"] + CONTROLS]
    md += ["", "Fragments in N269B show U251-lineage transcripts in the contralateral library. Resident cells and index "
           "hopping both predict them (M/reads/INDEX_HOPPING.txt: N269B shares its i7 with IL66B and IL69B and its i5 "
           "with NL71B and C2B). The controls IL64B and N168B carry the same exposure to hopping from human-rich "
           "libraries; a rate in N269B far above theirs is what hopping alone does not explain, and it is the only "
           "comparison here that bears on it."]
    if notes:
        md += ["", "## Notes", ""] + [f"- {n}" for n in notes]
    open(f"{E}/EVIDENCE.md", "w").write("\n".join(md) + "\n")
    json.dump({s: dict(n=r["n_candidates"], rank=r["rank"], p=r["p"],
                       stats=[p["stat"] for p in r["splits"]]) for s, r in results.items()},
              open(f"{E}/summary.json", "w"), indent=1)
    open(f"{E}/EVIDENCE_OK", "w").write(time.strftime("%Y-%m-%dT%H:%M:%S") + "\n")
    print(open(f"{E}/EVIDENCE.md").read())
    return 0


if __name__ == "__main__":
    sys.exit(main())
