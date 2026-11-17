#!/usr/bin/env python3
"""Run the R suite against a built backend image. Python standard library only."""
import argparse
from pathlib import Path
import shlex
import subprocess
import sys
import uuid

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent


def command(args, output):
    container = getattr(args, "container", None)
    suite = getattr(args, "suite", "backend")
    suite_path = f"{args.staging}/suite" if container else "/opt/msigportal-r-regression"
    report_path = f"{args.staging}/reports" if container else "/regression"
    baseline_path = f"{args.staging}/baseline" if container else "/baseline"
    suite_mount = f"{HERE / 'r'}:/opt/msigportal-r-regression:ro"
    output_mount = f"{output}:/regression"
    if container:
        cmd = ["docker", "exec", container, "Rscript"]
    elif args.image:
        cmd = ["docker", "run", "--rm", "--network", "none", "--entrypoint", "Rscript",
               "--volume", suite_mount, "--volume", output_mount]
        if args.baseline:
            cmd += ["--volume", f"{args.baseline.resolve()}:/baseline:ro"]
        cmd += [args.image]
    else:
        cmd = ["docker", "compose", "run", "--rm", "--no-deps", "--entrypoint", "Rscript",
               "--volume", suite_mount, "--volume", output_mount]
        if args.baseline:
            cmd += ["--volume", f"{args.baseline.resolve()}:/baseline:ro"]
        cmd += [args.service]
    script = "run.R" if suite == "backend" else "refitting_run.R"
    cmd += ["--vanilla", f"{suite_path}/{script}", "--server", args.server,
            "--output", report_path, "--timeout", str(args.timeout)]
    if args.baseline:
        cmd += ["--baseline", baseline_path]
    if args.case:
        cmd += ["--case", args.case]
    if args.list:
        cmd += ["--list"]
    return cmd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True, help="Fresh output label under data/regression/<label>/r")
    target = parser.add_mutually_exclusive_group()
    target.add_argument("--image", help="Built image; otherwise use docker compose run backend")
    target.add_argument("--container", help="Run fresh R processes in an existing container, preserving manually installed packages")
    parser.add_argument("--suite", choices=("backend", "refitting"), default="backend")
    parser.add_argument("--service", default="backend", help="Compose service (default: backend)")
    parser.add_argument("--server", help="Application directory (default: /deploy/server or /refitting-service)")
    parser.add_argument("--baseline", type=Path, help="Previous passing R output directory for numeric/table comparison")
    parser.add_argument("--case", help="R regular expression selecting case names")
    parser.add_argument("--timeout", default=600, type=int, help="Timeout in seconds per R case")
    parser.add_argument("--list", action="store_true", help="List cases without running them")
    parser.add_argument("--dry-run", action="store_true", help="Print the Docker command only")
    args = parser.parse_args()
    args.server = args.server or ("/deploy/server" if args.suite == "backend" else "/refitting-service")
    args.staging = f"/tmp/msigportal-r-regression-{uuid.uuid4().hex}"
    if args.suite == "refitting" and not (args.image or args.container):
        parser.error("--suite refitting requires --image or --container for the separate refitting service")
    if not args.label or Path(args.label).name != args.label or args.label in (".", ".."):
        parser.error("--label must be a single directory name")
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    output = REPO / "data" / "regression" / args.label / ("r" if args.suite == "backend" else "r-refitting")
    if args.baseline and not (args.baseline / "summary.csv").is_file():
        parser.error("--baseline must contain an R summary.csv")
    cmd = command(args, output)
    print(shlex.join(cmd), flush=True)
    if args.dry_run:
        return 0
    if output.exists() and any(output.iterdir()):
        parser.error(f"Output directory is not empty; use a new --label: {output}")
    output.mkdir(parents=True, exist_ok=True)
    with (output / "runner.log").open("w") as log:
        try:
            if args.container:
                prepare = [["docker", "exec", args.container, "mkdir", "-p", f"{args.staging}/suite"],
                           ["docker", "cp", f"{HERE / 'r'}/.", f"{args.container}:{args.staging}/suite"]]
                if args.baseline:
                    prepare += [["docker", "exec", args.container, "mkdir", "-p", f"{args.staging}/baseline"],
                                ["docker", "cp", f"{args.baseline.resolve()}/.", f"{args.container}:{args.staging}/baseline"]]
                for step in prepare:
                    subprocess.run(step, check=True, stdout=log, stderr=subprocess.STDOUT)
            with subprocess.Popen(cmd, cwd=REPO, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  text=True, bufsize=1) as process:
                for line in process.stdout:
                    print(line, end="", flush=True)
                    log.write(line)
                    log.flush()
                status = process.wait()
            if args.container:
                copied = subprocess.run(["docker", "cp", f"{args.container}:{args.staging}/reports/.", str(output)],
                                        stdout=log, stderr=subprocess.STDOUT)
                if copied.returncode:
                    status = status or copied.returncode
                print(f"Container scratch reports: {args.staging}/reports")
        except subprocess.CalledProcessError as exc:
            print(f"Container staging failed; see {output / 'runner.log'}", file=sys.stderr)
            return exc.returncode
        except FileNotFoundError:
            print("Docker is not installed or is not on PATH", file=sys.stderr)
            return 2
    print(f"Reports: {output}")
    return status


if __name__ == "__main__":
    sys.exit(main())
