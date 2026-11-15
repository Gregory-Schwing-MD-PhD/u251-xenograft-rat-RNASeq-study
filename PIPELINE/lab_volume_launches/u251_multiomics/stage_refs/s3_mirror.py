#!/usr/bin/env python3
"""Mirror a subset of a public S3 bucket with anonymous ListObjectsV2 + curl, then verify every object
against the ETag S3 computed at upload (MD5 for single-part uploads; MD5-of-part-MD5s for multipart).

Checks are written as TSV rows (check, expected, observed, result, outside_reference) to --checks.
Exit 0 only if every listed object is present at its listed size AND its content matches its ETag.

Idempotent: objects already at their listed size are not re-fetched; partial files resume with curl -C -.
A verification cache (key, size, mtime_ns -> computed ETag) avoids re-hashing unchanged files on a requeue.
Standard library only (runs on the grid's nextflow-env python).
"""
import argparse, concurrent.futures as cf, hashlib, os, subprocess, sys, time
import urllib.parse, urllib.request, xml.etree.ElementTree as ET

NS = {"s": "http://s3.amazonaws.com/doc/2006-03-01/"}
MiB = 1024 * 1024


def list_prefix(bucket, prefix, curl="curl"):
    objs, token = [], None
    while True:
        q = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if token:
            q["continuation-token"] = token
        url = f"https://{bucket}.s3.amazonaws.com/?" + urllib.parse.urlencode(q)
        for attempt in range(6):
            try:
                body = urllib.request.urlopen(url, timeout=120).read()
                break
            except Exception as e:  # network flake, or no CA bundle for this python: try curl, then back off
                r = subprocess.run([curl, "-s", "--fail", "--connect-timeout", "60", url], capture_output=True)
                if r.returncode == 0 and r.stdout:
                    body = r.stdout
                    break
                if attempt == 5:
                    raise
                time.sleep(10 * (attempt + 1))
        root = ET.fromstring(body)
        for c in root.findall("s:Contents", NS):
            objs.append((c.find("s:Key", NS).text, int(c.find("s:Size", NS).text),
                         c.find("s:ETag", NS).text.strip('"'), c.find("s:LastModified", NS).text))
        if root.find("s:IsTruncated", NS).text == "true":
            token = root.find("s:NextContinuationToken", NS).text
        else:
            return objs


def etag_of(path, etag):
    """Compute the ETag S3 would give this file, using the part count encoded in the listed ETag."""
    size = os.path.getsize(path)
    if "-" not in etag:
        h = hashlib.md5()
        with open(path, "rb") as f:
            for blk in iter(lambda: f.read(8 * MiB), b""):
                h.update(blk)
        return h.hexdigest()
    nparts = int(etag.split("-")[1])
    # aws cli / boto default is 8 MiB parts. Every common part size that gives the listed part count is tried, so a
    # file uploaded with another part size is not called corrupt; the first candidate's value is returned on mismatch.
    cands = [m * MiB for m in (8, 16, 5, 10, 15, 32, 50, 64, 100, 128, 256, 512) if -(-size // (m * MiB)) == nparts]
    if not cands:
        return "UNVERIFIABLE_PARTSIZE"
    first = None
    for ps in cands:
        digests = []
        with open(path, "rb") as f:
            while True:
                blk = f.read(ps)
                if not blk:
                    break
                digests.append(hashlib.md5(blk).digest())
        got = hashlib.md5(b"".join(digests)).hexdigest() + f"-{len(digests)}"
        if got == etag:
            return got
        first = first or got
    return first


def _etag_job(args):
    path, etag = args
    try:
        return etag_of(path, etag)
    except Exception as e:
        return f"ERROR:{e}"


def fetch(curl, bucket, key, size, local):
    os.makedirs(os.path.dirname(local), exist_ok=True)
    if os.path.exists(local):
        have = os.path.getsize(local)
        if have == size:
            return "present"
        if have > size:
            os.remove(local)
    url = f"https://{bucket}.s3.amazonaws.com/" + urllib.parse.quote(key, safe="/")
    cmd = [curl, "--silent", "--show-error", "--fail", "--location", "--retry", "5", "--retry-delay", "15",
           "--connect-timeout", "60", "-o", local]
    if size > 0 and os.path.exists(local) and os.path.getsize(local) > 0:
        cmd[1:1] = ["-C", "-"]
    cmd.append(url)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return f"curl_rc={r.returncode}:{r.stderr.strip()[:200]}"
    if size == 0 and not os.path.exists(local):
        open(local, "ab").close()  # curl does not always create a file for an empty body
    return "fetched"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bucket", required=True)
    ap.add_argument("--root", required=True, help="key prefix stripped when writing locally, e.g. rnafusion/references/")
    ap.add_argument("--prefix", action="append", required=True, help="key prefix below --root to mirror (repeatable)")
    ap.add_argument("--dest", required=True)
    ap.add_argument("--expected", help="snapshot listing (size etag key) taken when the plan was made")
    ap.add_argument("--checks", required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--curl", default="curl")
    a = ap.parse_args()

    os.makedirs(a.dest, exist_ok=True)
    rows = []

    def chk(name, exp, obs, res, ref):
        rows.append("\t".join(str(x) for x in (name, exp, obs, res, ref)))

    objs = []
    for p in a.prefix:
        got = list_prefix(a.bucket, a.root + p, a.curl)
        objs += got
        chk(f"s3_list:{p}", ">0 objects", f"{len(got)} objects, {sum(o[1] for o in got)} bytes",
            "PASS" if got else "FAIL", f"s3://{a.bucket}/{a.root}{p} ListObjectsV2")
    objs = sorted(set(objs))
    print(f"[s3_mirror] listed {len(objs)} objects, {sum(o[1] for o in objs) / 1e9:.2f} GB", flush=True)
    with open(os.path.join(a.dest, "LISTING.tsv"), "w") as f:
        f.write("key\tsize\tetag\tlast_modified\n")
        for k, s, e, m in objs:
            f.write(f"{k}\t{s}\t{e}\t{m}\n")
    total = sum(o[1] for o in objs)
    chk("s3_list:total", "122138849929 bytes in 322 objects (listing of 2026-09-26)",
        f"{total} bytes in {len(objs)} objects", "INFO", "planning snapshot")

    if a.expected and os.path.exists(a.expected):
        snap = {}
        wanted = tuple(a.root + p for p in a.prefix)
        for line in open(a.expected):
            s, e, k = line.split()
            if k.startswith(wanted):
                snap[k] = (int(s), e.strip('"'))
        live = {k: (s, e) for k, s, e, _ in objs}
        added = sorted(set(live) - set(snap))
        removed = sorted(set(snap) - set(live))
        changed = sorted(k for k in set(live) & set(snap) if live[k] != snap[k])
        res = "PASS" if not (added or removed or changed) else "INFO"
        chk("s3_listing_vs_planning_snapshot", "identical", f"added={len(added)} removed={len(removed)} changed={len(changed)}"
            + (f" first_changed={(changed + added + removed)[0]}" if res != "PASS" else ""), res,
            "listing captured 2026-09-26 when the pipeline was evaluated (a change is reported, not blocked)")

    local = {k: os.path.join(a.dest, k[len(a.root):]) for k, _, _, _ in objs}
    cache_path = os.path.join(a.dest, ".verify_cache.tsv")
    cache = {}
    if os.path.exists(cache_path):
        for line in open(cache_path):
            k, s, mt, et = line.rstrip("\n").split("\t")
            cache[k] = (int(s), int(mt), et)

    bad = {}
    for rnd in (1, 2):
        todo = [o for o in objs if rnd == 1 or o[0] in bad]
        on_disk = lambda: sum(os.path.getsize(local[o[0]]) for o in todo if os.path.exists(local[o[0]]))
        b0, t0 = on_disk(), time.time()
        with cf.ThreadPoolExecutor(a.jobs) as ex:
            fres = dict(zip([o[0] for o in todo],
                            ex.map(lambda o: fetch(a.curl, a.bucket, o[0], o[1], local[o[0]]), todo)))
        dt, db = time.time() - t0, on_disk() - b0
        print(f"[s3_mirror] round {rnd}: {len(todo)} objects, fetched {db / 1e9:.2f} GB in {dt:.0f} s "
              f"({db / 1e6 / max(dt, 1):.1f} MB/s with {a.jobs} parallel curl)", flush=True)
        need = []
        for k, s, e, _ in todo:
            p = local[k]
            if not os.path.exists(p) or os.path.getsize(p) != s:
                continue
            st = os.stat(p)
            c = cache.get(k)
            if c and c[0] == st.st_size and c[1] == st.st_mtime_ns and c[2] == e:
                continue
            need.append((k, p, e))
        # threads, not processes: hashlib releases the GIL on large buffers, and a process pool would need
        # multiprocessing IPC sockets in a temp directory (node-local /tmp is swept on this grid)
        t_hash = time.time()
        with cf.ThreadPoolExecutor(a.jobs) as ex:
            comp = dict(zip([n[0] for n in need], ex.map(_etag_job, [(n[1], n[2]) for n in need])))
        hb = sum(os.path.getsize(n[1]) for n in need)
        print(f"[s3_mirror] round {rnd}: hashed {len(need)} objects, {hb / 1e9:.2f} GB in {time.time() - t_hash:.0f} s",
              flush=True)
        for k, p, e in need:
            if comp[k] == e:
                st = os.stat(p)
                cache[k] = (st.st_size, st.st_mtime_ns, e)
            else:
                cache.pop(k, None)  # a stale entry must never vouch for changed content
        with open(cache_path + ".part", "w") as f:
            for k, (s, mt, et) in sorted(cache.items()):
                f.write(f"{k}\t{s}\t{mt}\t{et}\n")
        os.replace(cache_path + ".part", cache_path)
        bad = {}
        for k, s, e, _ in objs:
            p = local[k]
            st = os.stat(p) if os.path.exists(p) else None
            ls = st.st_size if st else -1
            ok = (ls == s and k in cache and cache[k][2] == e
                  and cache[k][0] == st.st_size and cache[k][1] == st.st_mtime_ns)
            if not ok:
                bad[k] = (s, ls, fres.get(k, ""), comp.get(k, ""))
        if not bad:
            break
        for k in bad:  # wrong content at the right size: delete so round 2 fetches afresh
            p = local[k]
            if os.path.exists(p) and os.path.getsize(p) == bad[k][0]:
                os.remove(p)
                cache.pop(k, None)

    with open(os.path.join(a.dest, "MANIFEST.tsv.part"), "w") as f:
        f.write("key\tlisted_size\tlocal_size\tlisted_etag\tcomputed_etag\tok\n")
        for k, s, e, _ in objs:
            p = local[k]
            ls = os.path.getsize(p) if os.path.exists(p) else -1
            ce = cache[k][2] if k not in bad and k in cache else (bad.get(k, (0, 0, "", ""))[3] or "NA")
            ok = int(ls == s and ce == e and k not in bad)
            f.write(f"{k}\t{s}\t{ls}\t{e}\t{ce}\t{ok}\n")
    os.replace(os.path.join(a.dest, "MANIFEST.tsv.part"), os.path.join(a.dest, "MANIFEST.tsv"))

    n_ok = len(objs) - len(bad)
    chk("s3_objects_present_at_listed_size_and_etag", f"{len(objs)}/{len(objs)}", f"{n_ok}/{len(objs)}",
        "PASS" if not bad else "FAIL", "S3 ETag (MD5 / multipart MD5, 8 MiB parts) computed by S3 at upload")
    for k in sorted(bad)[:20]:
        s, ls, fr, ce = bad[k]
        chk(f"s3_object:{k[len(a.root):]}", f"size={s}", f"size={ls} fetch={fr} etag={ce}", "FAIL", "S3 listing")
    n_multi = sum(1 for o in objs if "-" in o[2])
    chk("s3_etag_kinds", "", f"{len(objs) - n_multi} single-part (MD5) + {n_multi} multipart", "INFO", "")
    with open(a.checks, "a") as f:
        for r in rows:
            f.write(r + "\n")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
