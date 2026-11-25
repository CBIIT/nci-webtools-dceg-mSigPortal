#!/usr/bin/env python3
"""Run mSigPortal Python regression cases inside the backend / extraction images.

Stdlib only. Shells out to `docker`. Outputs go to data/regression/<label>/<suite>/<case>/.

Examples
  # capture baselines on the current images
  python3 utils/regression/run_cases.py profiler  --label baseline
  python3 utils/regression/run_cases.py extractor --label baseline --image msigportal-extraction:py39

  # rerun on the new images, then compare
  python3 utils/regression/run_cases.py profiler  --label candidate
  python3 utils/regression/run_cases.py extractor --label candidate --image msigportal-extraction:py313
  python3 utils/regression/compare.py data/regression/baseline data/regression/candidate

  # build the wide signature database TSV the extractor needs from a long-format psql export
  python3 utils/regression/run_cases.py pivot-signatures sig_long.tsv data/regression/inputs/signature_SBS96.tsv
"""
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
CASES_FILE = HERE / "cases.json"
DATA_DIR = REPO / "data"
REG_DIR = DATA_DIR / "regression"
INPUTS_DIR = REG_DIR / "inputs"
CONTAINER_DATA = "/data"  # both images mount ./data here


def load_cases():
    with open(CASES_FILE) as f:
        return json.load(f)


def container_path(host_path: Path) -> str:
    """Map a path under ./data to its in-container path."""
    rel = host_path.resolve().relative_to(DATA_DIR.resolve())
    return f"{CONTAINER_DATA}/{rel.as_posix()}"


def stage_inputs(suite_cfg):
    """Copy example inputs into data/regression/inputs so the container can see them."""
    INPUTS_DIR.mkdir(parents=True, exist_ok=True)
    aliases = suite_cfg.get("input_aliases", {})
    for rel in suite_cfg["inputs"]:
        src = REPO / rel
        if not src.exists():
            print(f"  [inputs] missing (optional cases using it will be skipped): {rel}")
            continue
        name = aliases.get(src.name, src.name)
        dst = INPUTS_DIR / name
        if not dst.exists() or dst.stat().st_size != src.stat().st_size:
            shutil.copy2(src, dst)
    # signature databases: use the COSMIC reference files shipped with the Extractor checkout
    # unless the user already exported a table from Postgres (see README)
    for ctx, name in suite_cfg.get("signature_db", {}).items():
        dst = INPUTS_DIR / name
        src = REPO / suite_cfg.get("signature_db_fallback", {}).get(ctx, "")
        if not dst.exists() and src.is_file():
            shutil.copy2(src, dst)
            print(f"  [inputs] {name} <- {src.relative_to(REPO)}")
    if "seed" in suite_cfg:
        (INPUTS_DIR / "Seeds.txt").write_text(f"Replicate\tSeed\n1\t{suite_cfg['seed']}\n")


DRY_RUN = False


class _FakeProc:
    returncode = 0
    stdout = "dry-run\n"
    stderr = ""


def run(cmd, log_path: Path = None, check=False):
    print("  $ " + " ".join(cmd))
    if DRY_RUN:
        return _FakeProc()
    if log_path:
        with open(log_path, "w") as log:
            proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
    else:
        proc = subprocess.run(cmd, text=True, capture_output=True)
    if check and proc.returncode != 0:
        raise SystemExit(f"command failed ({proc.returncode}): {' '.join(cmd)}")
    return proc


def docker_prefix(args, suite):
    """Return the docker command prefix that executes inside the right image."""
    data_mount = f"{DATA_DIR.resolve()}:{CONTAINER_DATA}"
    if args.image:
        workdir = "/deploy/server" if suite == "profiler" else "/deploy/app"
        cmd = ["docker", "run", "--rm", "-v", data_mount, "-w", workdir]
        if args.gpus:
            cmd += ["--gpus", args.gpus]
        if args.docker_args:
            cmd += args.docker_args.split()
        return cmd + [args.image]
    if suite != "profiler":
        raise SystemExit("--image is required for the extractor suite (no compose service exists for it)")
    # docker compose: backend service already mounts ./data:/data and loads server/.env
    return ["docker", "compose", "-f", str(REPO / "docker-compose.yml"),
            "run", "--rm", "--no-deps", "-T", "backend"]


# the py3.13 images keep system python3 and select the interpreter via ENV PYTHON_BIN, like the Node spawn sites
PYTHON = ["sh", "-c", 'exec "${PYTHON_BIN:-python3}" "$@"', "python"]


def capture_env(prefix, out_dir: Path, label: str, image_desc: str):
    """Record interpreter and package versions of the image under test."""
    py = run(prefix + PYTHON + ["-c",
                       "import sys,platform;print(sys.version.replace(chr(10),' '));print(platform.platform())"])
    freeze = run(prefix + PYTHON + ["-m", "pip", "freeze"])
    env = {
        "label": label,
        "image": image_desc,
        "captured_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "python": py.stdout.strip().splitlines() if py.returncode == 0 else py.stderr.strip(),
        "pip_freeze": sorted(freeze.stdout.strip().splitlines()) if freeze.returncode == 0 else freeze.stderr.strip(),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "_env.json", "w") as f:
        json.dump(env, f, indent=2)
    print(f"  [env] {env['python'][0] if isinstance(env['python'], list) else env['python']}")


def select_cases(cases: dict, args):
    tags_in = set(args.include_tags.split(",")) if args.include_tags else set()
    for name, case in cases.items():
        if args.case and name not in args.case:
            continue
        tags = set(case.get("tags", []))
        if "slow" in tags and "slow" not in tags_in and not args.all:
            print(f"- {name}: skipped (slow; pass --include-tags slow or --all)")
            continue
        if "gpu" in tags and "gpu" not in tags_in:
            print(f"- {name}: skipped (gpu; pass --include-tags gpu and --gpus all)")
            continue
        yield name, case


# --------------------------------------------------------------------------- profiler

def profiler_cli(case_name, case, out_dir: Path):
    """Mirror server/services/api/visualization/profilerExtraction.js + userForm.jsx arg assembly."""
    input_path = INPUTS_DIR / case["input"]
    if not input_path.exists():
        return None
    output_dir = out_dir / "profilerExtraction"
    args = {
        "Input_Format": case["Input_Format"],
        "Input_Path": container_path(input_path),
        "Project_ID": case_name,
        "Output_Dir": container_path(output_dir),
        "Genome_Building": case["Genome_Building"],
        "Data_Type": case["Data_Type"],
        "Cluster": case.get("Cluster", "False"),
        "SeqInfo": case.get("SeqInfo", "False"),
    }
    if case["Input_Format"] in ("vcf", "csv", "tsv"):
        if case.get("Filter"):
            args["Filter"] = case["Filter"]
        args["Collapse"] = case.get("Collapse", "False")
        args["vcf_split_all_filter"] = case.get("vcf_split_all_filter", "False")
        if case.get("Bed"):
            args["Bed"] = container_path(INPUTS_DIR / case["Bed"])
    cli = []
    for k, v in args.items():
        cli += [f"--{k}", str(v)]
    return cli


def cmd_profiler(args):
    cfg = load_cases()["profiler"]
    stage_inputs(cfg)
    suite_dir = REG_DIR / args.label / "profiler"
    suite_dir.mkdir(parents=True, exist_ok=True)
    prefix = docker_prefix(args, "profiler")
    capture_env(prefix, suite_dir, args.label, args.image or "docker compose service: backend")

    results = {}
    for name, case in select_cases(cfg["cases"], args):
        out_dir = suite_dir / name
        if out_dir.exists():
            shutil.rmtree(out_dir)
        out_dir.mkdir(parents=True)
        cli = profiler_cli(name, case, out_dir)
        if cli is None:
            print(f"- {name}: skipped (input {case['input']} not staged)")
            results[name] = "skipped"
            continue
        print(f"- {name} [{case.get('package')}]")
        t0 = time.time()
        proc = run(prefix + PYTHON + [cfg["script"]] + cli,
                   log_path=out_dir / "profiler_extraction_log.txt")
        dt = time.time() - t0
        status = "ok" if proc.returncode == 0 else f"exit {proc.returncode}"
        if proc.returncode != 0 and "optional" in case.get("tags", []):
            status += " (optional)"
        results[name] = status
        print(f"  -> {status} in {dt:.0f}s")
        with open(out_dir / "_case.json", "w") as f:
            json.dump({"case": name, "args": cli, "exit_code": proc.returncode,
                       "seconds": round(dt, 1), "package": case.get("package")}, f, indent=2)
    print_summary(results)


# --------------------------------------------------------------------------- extractor

def extractor_cli(case_name, case, cfg, out_dir: Path, sig_db_override: str):
    input_path = INPUTS_DIR / case["input"]
    if not input_path.exists():
        return None, f"input {case['input']} not staged"
    ctx = case["context_type"]
    sig_db = Path(sig_db_override) if sig_db_override else INPUTS_DIR / cfg["signature_db"][ctx]
    if not sig_db.exists():
        return None, (f"signature database {sig_db} missing; export it with psql and "
                      f"`run_cases.py pivot-signatures` (see README)")
    merged = dict(cfg["defaults"])
    merged.update({k: v for k, v in case.items() if k not in ("package", "tags", "_note", "input")})
    seeds = INPUTS_DIR / merged["seeds"]
    merged["seeds"] = container_path(seeds) if seeds.exists() else "random"
    merged["input_data"] = container_path(input_path)
    merged["output"] = container_path(out_dir)
    merged["signature_database"] = container_path(sig_db) if sig_db.resolve().is_relative_to(DATA_DIR.resolve()) \
        else sig_db.as_posix()
    cli = []
    for k, v in merged.items():
        cli += [f"--{k}", str(v)]
    return cli, None


def cmd_extractor(args):
    cfg = load_cases()["extractor"]
    stage_inputs(cfg)
    suite_dir = REG_DIR / args.label / "extractor"
    suite_dir.mkdir(parents=True, exist_ok=True)
    prefix = docker_prefix(args, "extractor")
    capture_env(prefix, suite_dir, args.label, args.image)

    results = {}
    for name, case in select_cases(cfg["cases"], args):
        out_dir = suite_dir / name
        if out_dir.exists():
            shutil.rmtree(out_dir)
        out_dir.mkdir(parents=True)
        cli, why = extractor_cli(name, case, cfg, out_dir, args.signature_db)
        if cli is None:
            print(f"- {name}: skipped ({why})")
            results[name] = "skipped"
            continue
        print(f"- {name} [{case.get('package')}]")
        t0 = time.time()
        proc = run(prefix + PYTHON + [cfg["script"]] + cli,
                   log_path=out_dir / "extraction_log.txt")
        dt = time.time() - t0
        status = "ok" if proc.returncode == 0 else f"exit {proc.returncode}"
        results[name] = status
        print(f"  -> {status} in {dt:.0f}s")
        with open(out_dir / "_case.json", "w") as f:
            json.dump({"case": name, "args": cli, "exit_code": proc.returncode,
                       "seconds": round(dt, 1), "package": case.get("package")}, f, indent=2)
    print_summary(results)


# --------------------------------------------------------------------------- helpers

def cmd_pivot(args):
    """Long psql export (signatureName, mutationType, contribution) -> wide TSV with a Type column.

    Mirrors extraction-service/services/extraction.js (signature.tsv written before the Python call).
    """
    wide = {}
    names = []
    with open(args.long_tsv, newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        reader = csv.DictReader(f, dialect=dialect)
        for row in reader:
            sig, mt, c = row["signatureName"], row["mutationType"], row["contribution"]
            if sig not in names:
                names.append(sig)
            wide.setdefault(mt, {})[sig] = c
    out = Path(args.wide_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["Type"] + names)
        for mt in wide:
            w.writerow([mt] + [wide[mt].get(s, "0") for s in names])
    print(f"wrote {out} ({len(wide)} mutation types x {len(names)} signatures)")


def print_summary(results):
    print("\nsummary")
    for k, v in results.items():
        print(f"  {k:32s} {v}")
    failed = [k for k, v in results.items() if v.startswith("exit") and "optional" not in v]
    if failed:
        print(f"\n{len(failed)} case(s) failed: {', '.join(failed)}")
        sys.exit(1)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--label", required=True, help="output folder under data/regression (e.g. baseline, candidate)")
        sp.add_argument("--case", nargs="*", help="run only these case names")
        sp.add_argument("--include-tags", default="", help="comma list of tags to include: slow,gpu")
        sp.add_argument("--all", action="store_true", help="include slow cases")
        sp.add_argument("--image", help="run via `docker run <image>` instead of docker compose")
        sp.add_argument("--gpus", help="value for `docker run --gpus` (e.g. all); only with --image")
        sp.add_argument("--docker-args", help="extra args for `docker run` (quoted string)")
        sp.add_argument("--dry-run", action="store_true", help="print the docker commands without running them")

    sp = sub.add_parser("profiler", help="backend image: mSigPortal_Profiler_Extraction.py cases")
    common(sp)
    sp.set_defaults(func=cmd_profiler)

    sp = sub.add_parser("extractor", help="extraction image: mSigPortal-SigProfilerExtractor.py cases")
    common(sp)
    sp.add_argument("--signature-db", help="override wide signature TSV (default per context in data/regression/inputs)")
    sp.set_defaults(func=cmd_extractor)

    sp = sub.add_parser("pivot-signatures", help="long psql export -> wide signature TSV")
    sp.add_argument("long_tsv")
    sp.add_argument("wide_tsv")
    sp.set_defaults(func=cmd_pivot)

    args = p.parse_args()
    global DRY_RUN
    DRY_RUN = getattr(args, "dry_run", False)
    os.chdir(REPO)
    args.func(args)


if __name__ == "__main__":
    main()
