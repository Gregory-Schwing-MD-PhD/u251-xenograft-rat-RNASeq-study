#!/usr/bin/env python3
"""Pass B gate for human_de_null (standard library only; runs on the compute node before any container starts).

Pass B reads the graft-only counts of rnaseq_graft. It may run only when
  1. rnaseq_graft's Nextflow run finished with exit 0 (NEXTFLOW_EXIT holds "0");
  2. checks.json was written after that exit (a checks.json older than NEXTFLOW_EXIT belongs to an earlier run);
  3. every required check in checks.json passed. The default set is 0 (same release and tool versions as the
     published run), 1 (clip applied exactly once), 2 (input counts carried from xengsort to STAR), 3 (reverse
     stranded), 4 (GENCODE v44, same gene set) and 10 (no empty outputs): the checks on which the counts depend.
     Check 9 (the saved human index) does not touch the counts and is reported, not required;
  4. pass A's R0 passed (pass A is the analysis of record, and pass B uses its universe rule and settings);
  5. the counts and lengths tables exist and carry every passB_column of the samplesheet.
Writes --out (gate.json) and exits 0 when the gate is open, 3 when it is closed.
"""
import argparse
import csv
import json
import os
import sys
import time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checks", required=True)
    ap.add_argument("--nextflow-exit", required=True)
    ap.add_argument("--required", default="0,1,2,3,4,10")
    ap.add_argument("--r0", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--lengths", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    req = [x.strip() for x in a.required.split(",") if x.strip()]
    g = {"written": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "checks_json": a.checks, "required": req, "reasons": []}

    def fail(msg):
        g["reasons"].append(msg)
        print("GATE: " + msg)

    # 1-2. Nextflow exit and freshness of checks.json
    nfx = None
    if os.path.isfile(a.nextflow_exit):
        with open(a.nextflow_exit) as fh:
            nfx = fh.read().strip()
    g["nextflow_exit"] = nfx
    if nfx != "0":
        fail("rnaseq_graft NEXTFLOW_EXIT is %r, not '0'" % nfx)
    if not os.path.isfile(a.checks):
        fail("rnaseq_graft checks.json does not exist")
        C = {}
    else:
        if os.path.isfile(a.nextflow_exit) and os.path.getmtime(a.checks) < os.path.getmtime(a.nextflow_exit):
            fail("checks.json is older than NEXTFLOW_EXIT (it describes an earlier run)")
        with open(a.checks) as fh:
            C = json.load(fh)
    g["all_hard_pass"] = C.get("all_hard_pass")
    g["failed_hard"] = C.get("failed_hard")
    # 3. required checks
    got = {}
    for k in req:
        c = C.get("checks", {}).get(k)
        if c is None:
            fail("required check %s is missing from checks.json" % k)
            got[k] = None
            continue
        got[k] = {"title": c.get("title"), "hard": c.get("hard"), "pass": c.get("pass")}
        if c.get("pass") is not True:
            fail("required check %s (%s) did not pass: %r" % (k, c.get("title"), c.get("pass")))
    g["required_checks"] = got
    g["other_hard_checks"] = {k: v.get("pass") for k, v in C.get("checks", {}).items() if v.get("hard") and k not in req}
    # 4. pass A R0
    r0 = None
    if os.path.isfile(a.r0):
        with open(a.r0) as fh:
            try:
                r0 = json.load(fh)
            except ValueError as e:
                fail("pass A r0_check.json is not valid JSON: %s" % e)
    g["passA_r0_pass"] = None if r0 is None else r0.get("pass")
    if r0 is None:
        fail("pass A r0_check.json does not exist (pass A did not reach R0)")
    elif r0.get("pass") is not True:
        fail("pass A R0 did not pass; nothing from this item is used")
    # 5. inputs
    with open(a.samplesheet) as fh:
        cols = [r["passB_column"] for r in csv.DictReader(fh)]
    for f in (a.counts, a.lengths):
        if not os.path.isfile(f) or os.path.getsize(f) == 0:
            fail("missing or empty: %s" % f)
            continue
        with open(f) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
        miss = [c for c in cols if c not in hdr]
        if miss:
            fail("%s lacks columns %s" % (f, ",".join(miss)))
    g["open"] = not g["reasons"]
    with open(a.out + ".part", "w") as fh:
        json.dump(g, fh, indent=1)
    os.replace(a.out + ".part", a.out)
    print("GATE %s (all_hard_pass of rnaseq_graft: %s; failed hard: %s)" % (
        "OPEN" if g["open"] else "CLOSED", g["all_hard_pass"], g["failed_hard"]))
    return 0 if g["open"] else 3


if __name__ == "__main__":
    sys.exit(main())
