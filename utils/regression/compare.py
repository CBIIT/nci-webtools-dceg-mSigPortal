#!/usr/bin/env python3
"""Compare two regression runs produced by run_cases.py.

  python3 utils/regression/compare.py data/regression/baseline data/regression/candidate [--json report.json]

Per-file verdicts
  identical   byte-identical
  tolerated   numeric table differences within tolerance, or row/column order only
  visual      image/PDF/SVG differs within the visual threshold (dimensions and page counts equal)
  log         log file; not diffed, scanned for Traceback/Error/Warning counts instead
  varies      output that differs run to run on every version (SigProfilerClusters results); line overlap reported
  FAIL        content differs beyond tolerance, or structure differs
  missing     present in baseline, absent in candidate (FAIL)
  extra       present in candidate only (informational)

Stdlib only; uses Pillow for perceptual image hashing when it is importable.
"""
import argparse
import io
import json
import math
import re
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

try:
    from PIL import Image  # optional
except Exception:  # pragma: no cover
    Image = None

LOG_PATTERNS = re.compile(r"(_log\.txt$|\.log$|\.out$|\.err$|JOB_METADATA\.txt$|_case\.json$|_env\.json$)")
TABLE_EXT = {".txt", ".tsv", ".csv", ".all", ".bed", ".exome", ".region", ".vcf"}
DATE_RE = re.compile(r"\d{4}-\d{2}-\d{2}|\d{14}")
NUM_RE = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")
# outputs embed their own run directory, e.g. /data/regression/<label>/...
LABEL_PATH_RE = re.compile(rb"/data/regression/[^/\s]+/")
# SigProfilerClusters assigns subclasses nondeterministically (worker order, unseeded random)
VARIES_RE = re.compile(r"Cluster/Result/[^/]*\.txt$")
PROVENANCE_PREFIXES = ("##bcftools",)


# --------------------------------------------------------------------------- file discovery

def normalize_rel(rel: str) -> str:
    return DATE_RE.sub("<DATE>", rel)


def index_dir(root: Path):
    files = {}
    for p in root.rglob("*"):
        if p.is_file():
            files[normalize_rel(p.relative_to(root).as_posix())] = p
    return files


# --------------------------------------------------------------------------- comparers

def is_num(s):
    return bool(NUM_RE.match(s.strip())) if s else False


def num_equal(a, b, rtol, atol):
    fa, fb = float(a), float(b)
    if math.isnan(fa) and math.isnan(fb):
        return True
    if "." not in a and "." not in b and "e" not in a.lower() and "e" not in b.lower():
        return fa == fb  # integers must match exactly
    return math.isclose(fa, fb, rel_tol=rtol, abs_tol=atol)


def parse_table(text: str):
    """Return (header, rows, delimiter) or None if the text is not a consistent delimited table."""
    lines = [l.rstrip("\r\n") for l in text.splitlines() if l.strip()]
    if len(lines) < 2:
        return None
    for delim in ("\t", ","):
        widths = {len(l.split(delim)) for l in lines}
        if len(widths) == 1 and widths.pop() >= 2:
            rows = [l.split(delim) for l in lines]
            return rows[0], rows[1:], delim
    return None


def compare_tables(ta, tb, rtol, atol):
    ha, ra, _ = ta
    hb, rb, _ = tb
    if sorted(ha) != sorted(hb):
        return "FAIL", f"header differs: {sorted(set(ha) ^ set(hb))[:5]}"
    if len(ra) != len(rb):
        return "FAIL", f"row count {len(ra)} vs {len(rb)}"
    reorder = ha != hb
    col_b = [hb.index(c) for c in ha]
    rb = [[y[j] for j in col_b] for y in rb]
    sa, sb = sorted(ra), sorted(rb)
    if sa != ra or sb != rb:
        reorder = True
    diffs = []
    for i, (x, y) in enumerate(zip(sa, sb)):
        for j, (cx, cy) in enumerate(zip(x, y)):
            if cx == cy:
                continue
            if is_num(cx) and is_num(cy):
                if num_equal(cx, cy, rtol, atol):
                    reorder = reorder or False
                    continue
            diffs.append(f"row {x[0]!r} col {ha[j]!r}: {cx!r} != {cy!r}")
    if diffs:
        return "FAIL", f"{len(diffs)} cell(s) differ; e.g. " + "; ".join(diffs[:3])
    return ("tolerated", "order/format differences only") if reorder else ("tolerated", "numeric within tolerance")


def compare_text(a: bytes, b: bytes, rtol, atol):
    try:
        sa, sb = a.decode("utf-8"), b.decode("utf-8")
    except UnicodeDecodeError:
        return "FAIL", "binary content differs"
    if sa.startswith("##") or sb.startswith("##"):
        keep = lambda s: "\n".join(l for l in s.splitlines() if not l.startswith(PROVENANCE_PREFIXES))
        sa, sb = keep(sa), keep(sb)
        if sa == sb:
            return "tolerated", "provenance header lines only"
    ta, tb = parse_table(sa), parse_table(sb)
    if ta and tb:
        return compare_tables(ta, tb, rtol, atol)
    la = [l.strip() for l in sa.splitlines() if l.strip()]
    lb = [l.strip() for l in sb.splitlines() if l.strip()]
    if la == lb:
        return "tolerated", "whitespace only"
    # token-wise numeric tolerance for free text such as Statistics.txt
    if len(la) == len(lb):
        bad = []
        for x, y in zip(la, lb):
            tx, ty = x.split(), y.split()
            if len(tx) != len(ty):
                bad.append((x, y))
                continue
            for u, v in zip(tx, ty):
                if u != v and not (is_num(u) and is_num(v) and num_equal(u, v, rtol, atol)):
                    bad.append((x, y))
                    break
        if not bad:
            return "tolerated", "numeric within tolerance"
        return "FAIL", f"{len(bad)} line(s) differ; e.g. {bad[0][0][:80]!r} vs {bad[0][1][:80]!r}"
    return "FAIL", f"line count {len(la)} vs {len(lb)}"


def png_size(data: bytes):
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return int.from_bytes(data[16:20], "big"), int.from_bytes(data[20:24], "big")


def ahash(data: bytes):
    img = Image.open(io.BytesIO(data)).convert("L").resize((16, 16))
    px = list(img.getdata())
    mean = sum(px) / len(px)
    return sum(1 << i for i, v in enumerate(px) if v > mean)


def compare_image(a, b, threshold):
    sa, sb = png_size(a), png_size(b)
    if sa != sb:
        return "FAIL", f"dimensions {sa} vs {sb}"
    if Image is None:
        return "visual", "pixels differ (install Pillow for perceptual hash)"
    try:
        d = bin(ahash(a) ^ ahash(b)).count("1")
    except Exception as e:  # pragma: no cover
        return "visual", f"pixels differ (hash failed: {e})"
    return ("visual", f"hash distance {d}") if d <= threshold else ("FAIL", f"hash distance {d} > {threshold}")


PDF_PAGE_RE = re.compile(rb"/Type\s*/Page(?!s)")


def compare_pdf(a, b):
    pa, pb = len(PDF_PAGE_RE.findall(a)), len(PDF_PAGE_RE.findall(b))
    if pa != pb:
        return "FAIL", f"page count {pa} vs {pb}"
    return "visual", f"{pa} page(s); content differs"


def compare_zip(a, b, ctx):
    za, zb = zipfile.ZipFile(io.BytesIO(a)), zipfile.ZipFile(io.BytesIO(b))
    na = {normalize_rel(n): n for n in za.namelist() if not n.endswith("/")}
    nb = {normalize_rel(n): n for n in zb.namelist() if not n.endswith("/")}
    if set(na) != set(nb):
        return "FAIL", f"member lists differ: {sorted(set(na) ^ set(nb))[:5]}"
    worst, notes = "identical", []
    rank = {"identical": 0, "tolerated": 1, "visual": 2, "log": 2, "varies": 2, "FAIL": 3}
    for key in na:
        v, note = compare_bytes(Path(na[key]), za.read(na[key]), zb.read(nb[key]), ctx)
        if rank[v] > rank[worst]:
            worst, notes = v, [f"{na[key]}: {note}"]
    return worst, "; ".join(notes) if notes else f"{len(na)} members"


def scan_log(data: bytes):
    text = data.decode("utf-8", "replace")
    c = Counter()
    for line in text.splitlines():
        l = line.lower()
        if "traceback" in l:
            c["traceback"] += 1
        elif "deprecat" in l:
            c["deprecation"] += 1
        elif "warning" in l:
            c["warning"] += 1
        elif "error" in l:
            c["error"] += 1
    return c


def compare_varies(a: bytes, b: bytes):
    key_cols = ("samples", "chr", "start", "ref", "alt", "subclass")

    def rows(data):
        lines = data.decode("utf-8", "replace").splitlines()
        head = lines[0].split("\t") if lines else []
        idx = [head.index(c) for c in key_cols] if all(c in head for c in key_cols) else None
        out = []
        for l in lines[1:]:
            f = l.split("\t")
            out.append("|".join(f[i] for i in idx) if idx and len(f) > max(idx) else l)
        return out, idx is not None

    (la, keyed), (lb, _) = rows(a), rows(b)
    common = sum((Counter(la) & Counter(lb)).values())
    pct = 100 * common / max(len(la), len(lb), 1)
    what = "same mutation+subclass" if keyed else "identical rows"
    return "varies", f"{len(la)} vs {len(lb)} rows; {pct:.1f}% {what} (order ignored)"


def compare_bytes(path: Path, a: bytes, b: bytes, ctx):
    a = LABEL_PATH_RE.sub(b"/data/regression/<LABEL>/", a)
    b = LABEL_PATH_RE.sub(b"/data/regression/<LABEL>/", b)
    if a == b:
        return "identical", ""
    name = path.name
    ext = path.suffix.lower()
    if LOG_PATTERNS.search(name):
        ca, cb = scan_log(a), scan_log(b)
        note = f"baseline {dict(ca)} candidate {dict(cb)}"
        return ("FAIL", "new traceback: " + note) if cb["traceback"] > ca["traceback"] else ("log", note)
    if VARIES_RE.search(path.as_posix()):
        return compare_varies(a, b)
    if ext == ".png":
        return compare_image(a, b, ctx["image_threshold"])
    if ext == ".pdf":
        return compare_pdf(a, b)
    if ext == ".svg":
        return "visual", "svg differs"
    if ext == ".zip":
        return compare_zip(a, b, ctx)
    if ext in TABLE_EXT or ext == "" or ext == ".json":
        return compare_text(a, b, ctx["rtol"], ctx["atol"])
    return "FAIL", "binary content differs"


# --------------------------------------------------------------------------- driver

def compare_case(base_dir: Path, cand_dir: Path, ctx):
    fa, fb = index_dir(base_dir), index_dir(cand_dir)
    out = {}
    for rel in sorted(set(fa) | set(fb)):
        if rel not in fb:
            out[rel] = ("missing", "")
            continue
        if rel not in fa:
            out[rel] = ("extra", "")
            continue
        out[rel] = compare_bytes(fa[rel], fa[rel].read_bytes(), fb[rel].read_bytes(), ctx)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("baseline", type=Path)
    p.add_argument("candidate", type=Path)
    p.add_argument("--suite", choices=["profiler", "extractor"], help="limit to one suite")
    p.add_argument("--case", nargs="*", help="limit to case names")
    p.add_argument("--rtol", type=float, default=1e-6)
    p.add_argument("--atol", type=float, default=1e-9)
    p.add_argument("--image-threshold", type=int, default=12, help="max 256-bit average-hash distance for PNGs")
    p.add_argument("--json", type=Path, help="write a machine-readable report")
    p.add_argument("--verbose", "-v", action="store_true", help="list every file, not just non-identical ones")
    args = p.parse_args()
    ctx = {"rtol": args.rtol, "atol": args.atol, "image_threshold": args.image_threshold}

    report, failures = {}, 0
    suites = [args.suite] if args.suite else ["profiler", "extractor"]
    for suite in suites:
        bs, cs = args.baseline / suite, args.candidate / suite
        if not bs.exists():
            continue
        print(f"\n=== {suite} ===")
        for env_name in ("_env.json",):
            try:
                eb = json.loads((bs / env_name).read_text())["python"]
                ec = json.loads((cs / env_name).read_text())["python"]
                print(f"baseline python : {eb[0] if isinstance(eb, list) else eb}")
                print(f"candidate python: {ec[0] if isinstance(ec, list) else ec}")
            except Exception:
                pass
        for case_dir in sorted(d for d in bs.iterdir() if d.is_dir()):
            name = case_dir.name
            if args.case and name not in args.case:
                continue
            cand = cs / name
            if not cand.exists():
                print(f"\n- {name}: MISSING in candidate")
                report[f"{suite}/{name}"] = {"status": "missing"}
                failures += 1
                continue
            res = compare_case(case_dir, cand, ctx)
            counts = Counter(v for v, _ in res.values())
            bad = counts["FAIL"] + counts["missing"]
            failures += 1 if bad else 0
            verdict = "FAIL" if bad else "pass"
            print(f"\n- {name}: {verdict}  " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
            for rel, (v, note) in res.items():
                if v == "identical" and not args.verbose:
                    continue
                if v in ("FAIL", "missing") or args.verbose or v in ("visual", "log", "tolerated", "extra", "varies"):
                    print(f"    [{v:9s}] {rel}" + (f"  -- {note}" if note else ""))
            report[f"{suite}/{name}"] = {"status": verdict, "counts": dict(counts),
                                         "files": {k: {"verdict": v, "note": n} for k, (v, n) in res.items()}}

    if args.json:
        args.json.write_text(json.dumps(report, indent=2))
        print(f"\nreport written to {args.json}")
    print(f"\n{failures} case(s) failed" if failures else "\nall cases passed")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
