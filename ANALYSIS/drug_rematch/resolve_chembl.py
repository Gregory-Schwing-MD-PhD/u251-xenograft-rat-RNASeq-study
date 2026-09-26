#!/usr/bin/env python3
"""Resolve every DSigDB compound name to a ChEMBL molecule by EXACT name, not by the first free-text search hit.

Why: the published chain (create_publication_figure_600_dpi.R, query_chembl_with_api) took molecule/search?q=<name>
and kept result [0]. For "deferoxamine" that is CHEMBL4635234, an unnamed analogue with no clinical phase, so the
phase filter dropped an approved drug; vandetanib, mefloquine and bromocriptine went the same way, and pyrantel kept a
stale phase (claim check of 2026-09-26, SLIDES/deck_v4 notes).

Rule, per distinct name (case-insensitive):
  1. molecules whose pref_name equals the name (ChEMBL filter pref_name__iexact);
  2. else molecules that list the name as a synonym (molecule_synonyms__molecule_synonym__iexact);
  3. every hit is mapped to its parent (molecule_hierarchy.parent_chembl_id), so salts and forms collapse;
  4. one distinct parent -> that parent; several -> the one with the highest max_phase (ties: the one whose own
     pref_name equals the name, then the lowest ChEMBL number), flagged 'ambiguous' for review; none -> unresolved.
Unresolved names get no ChEMBL ID and no phase: they cannot enter the clinical ranking by a guessed record.
Overrides (overrides.tsv: drug, chembl_id, reason) win over the rule; they come from the audit, never from a search.

    python3 resolve_chembl.py resolve <table.csv> [<table.csv> ...] --out mapping.csv   # needs network (ChEMBL REST)
    python3 resolve_chembl.py apply mapping.csv <in.csv> <out.csv>                        # offline

Only the standard library, so it runs in any Python 3.
"""
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "https://www.ebi.ac.uk/chembl/api/data"
UA = {"User-Agent": "u251-litt-study/1.0 (mailto:go2432@wayne.edu)"}
CACHE = os.path.join(HERE, "chembl_exact_cache.json")
_cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}


def get(url):
    if url in _cache:
        return _cache[url]
    for attempt in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                d = json.load(r)
            _cache[url] = d
            if len(_cache) % 25 == 0:
                json.dump(_cache, open(CACHE, "w"))
            time.sleep(0.1)
            return d
        except urllib.error.HTTPError as e:
            if e.code == 404:
                _cache[url] = None
                return None
            time.sleep(2 * (attempt + 1))
        except Exception:
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("ChEMBL unreachable: " + url)


def molecule(cid):
    return get("%s/molecule/%s.json" % (BASE, cid))


def search(field, name):
    d = get("%s/molecule.json?%s=%s&limit=100" % (BASE, field, urllib.parse.quote(name)))
    return (d or {}).get("molecules", []) or []


def parent_of(m):
    p = (m.get("molecule_hierarchy") or {}).get("parent_chembl_id")
    if p and p != m["molecule_chembl_id"]:
        pm = molecule(p)
        return pm or m
    return m


def phase(m):
    try:
        return float(m.get("max_phase")) if m.get("max_phase") not in (None, "") else None
    except (TypeError, ValueError):
        return None


def num(cid):
    return int(re.sub(r"\D", "", cid) or 0)


def resolve(name):
    n = name.strip()
    how, hits = "none", []
    for field, label in (("pref_name__iexact", "pref_name"), ("molecule_synonyms__molecule_synonym__iexact", "synonym")):
        hits = search(field, n)
        if hits:
            how = label
            break
    parents = {}
    for h in hits:
        p = parent_of(h)
        parents[p["molecule_chembl_id"]] = p
    if not parents:
        return {"match": "unresolved", "n_parents": 0, "candidates": ""}
    ranked = sorted(parents.values(), key=lambda m: (-(phase(m) if phase(m) is not None else -1),
                                                     0 if (m.get("pref_name") or "").lower() == n.lower() else 1,
                                                     num(m["molecule_chembl_id"])))
    best = ranked[0]
    st = best.get("molecule_structures") or {}
    return {"match": how + ("_ambiguous" if len(parents) > 1 else ""), "n_parents": len(parents),
            "candidates": ";".join("%s|%s|%s" % (m["molecule_chembl_id"], m.get("pref_name") or "", phase(m)) for m in ranked),
            "chembl_id": best["molecule_chembl_id"], "pref_name": best.get("pref_name") or "",
            "max_phase": phase(best), "first_approval": best.get("first_approval") or "",
            "molecule_type": best.get("molecule_type") or "", "smiles": st.get("canonical_smiles") or ""}


def load_overrides():
    p = os.path.join(HERE, "overrides.tsv")
    if not os.path.exists(p):
        return {}
    return {r["drug"].strip().lower(): r for r in csv.DictReader(open(p, encoding="utf-8"), delimiter="\t")}


def cmd_resolve(tables, out):
    old = {}
    for t in tables:
        for r in csv.DictReader(open(t, encoding="utf-8")):
            k = r["Drug"].strip().lower()
            old.setdefault(k, {"drug": r["Drug"].strip(), "old_chembl_id": (r.get("ChEMBL_ID") or "").strip(),
                               "old_max_phase": (r.get("Max_Phase") or "").strip()})
    ov = load_overrides()
    rows = []
    for i, (k, o) in enumerate(sorted(old.items())):
        r = dict(o)
        r.update(resolve(o["drug"]))
        if k in ov:                                        # an audited decision replaces the rule's pick
            m = molecule(ov[k]["chembl_id"]) if ov[k]["chembl_id"] else None
            if m:
                st = m.get("molecule_structures") or {}
                r.update({"chembl_id": m["molecule_chembl_id"], "pref_name": m.get("pref_name") or "", "max_phase": phase(m),
                          "first_approval": m.get("first_approval") or "", "molecule_type": m.get("molecule_type") or "",
                          "smiles": st.get("canonical_smiles") or "", "match": "override"})
            else:
                r.update({"chembl_id": "", "pref_name": "", "max_phase": None, "smiles": "", "match": "override_none"})
            r["override_reason"] = ov[k].get("reason", "")
        oid = o["old_chembl_id"]
        if oid.upper().startswith("CHEMBL"):
            om = molecule(oid) or {}
            r["old_pref_name"] = om.get("pref_name") or ""
            r["old_phase_today"] = phase(om) if om else ""
            r["old_parent"] = (om.get("molecule_hierarchy") or {}).get("parent_chembl_id") or oid if om else ""
        r["id_changed"] = bool(r.get("chembl_id")) and r.get("chembl_id") != r.get("old_parent", oid)
        try:
            op = float(o["old_max_phase"]) if o["old_max_phase"] not in ("", "NA") else None
        except ValueError:
            op = None
        r["phase_changed"] = (op if op is not None else -99) != (r.get("max_phase") if r.get("max_phase") is not None else -99)
        rows.append(r)
        if (i + 1) % 20 == 0:
            print("  %d/%d" % (i + 1, len(old)), file=sys.stderr)
    json.dump(_cache, open(CACHE, "w"))
    cols = ["drug", "match", "chembl_id", "pref_name", "max_phase", "first_approval", "molecule_type", "n_parents", "candidates",
            "old_chembl_id", "old_pref_name", "old_max_phase", "old_phase_today", "old_parent", "id_changed", "phase_changed",
            "override_reason", "smiles"]
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c, "")) for c in cols})
    ch = [r for r in rows if r["id_changed"] or r["phase_changed"]]
    print("names %d | exact pref_name %d | synonym %d | ambiguous %d | unresolved %d | override %d | id or phase changed %d"
          % (len(rows), sum(r["match"].startswith("pref_name") for r in rows), sum(r["match"].startswith("synonym") for r in rows),
             sum("ambiguous" in r["match"] for r in rows), sum(r["match"] == "unresolved" for r in rows),
             sum(r["match"].startswith("override") for r in rows), len(ch)))


def cmd_apply(mapping, src, dst):
    m = {r["drug"].strip().lower(): r for r in csv.DictReader(open(mapping, encoding="utf-8"))}
    rows = list(csv.DictReader(open(src, encoding="utf-8")))
    cols = list(rows[0].keys()) + ["ChEMBL_ID_search", "Max_Phase_search", "chembl_match", "chembl_pref_name"]
    miss = 0
    for r in rows:
        k = r["Drug"].strip().lower()
        r["ChEMBL_ID_search"], r["Max_Phase_search"] = r.get("ChEMBL_ID", ""), r.get("Max_Phase", "")
        x = m.get(k)
        if x is None:
            miss += 1
            r["chembl_match"] = "not in mapping"
            continue
        r["ChEMBL_ID"] = x["chembl_id"] or ""
        r["Max_Phase"] = x["max_phase"] if x["max_phase"] not in ("", "None") else "NA"
        r["chembl_match"], r["chembl_pref_name"] = x["match"], x["pref_name"]
    assert miss == 0, "%d names are missing from the mapping" % miss
    with open(dst, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, quoting=csv.QUOTE_MINIMAL)
        w.writeheader()
        w.writerows(rows)
    print("applied to %d rows -> %s" % (len(rows), dst))


if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "resolve":
        out = a[a.index("--out") + 1]
        cmd_resolve([t for t in a[1:] if t not in ("--out", out)], out)
    elif a and a[0] == "apply":
        cmd_apply(a[1], a[2], a[3])
    else:
        sys.exit(__doc__)
