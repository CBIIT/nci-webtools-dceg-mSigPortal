# Python 3.13 upgrade of the backend and extraction-batch images

Status: in progress (started 2026-10-28). Branch `python313` in this repo and in each xtmgah fork.

## Why

Both Python-bearing images run Python 3.9 (Amazon Linux 2023 system Python in the backend image,
Rocky Linux 9 system Python in the GPU extraction image). The pinned stack has no Python 3.13 wheels:

| Pin in the Dockerfiles | First version with cp313 wheels |
|---|---|
| numpy==1.26.4 | numpy 2.1.0 |
| pandas==1.3.5 | pandas 2.2.3 |
| torch==2.1.0 | torch 2.5.x |
| matplotlib<=3.4.3 (SigProfilerExtractor cap) | matplotlib 3.9.x |

Moving to 3.13 therefore also moves the SigProfiler packages onto numpy 2 and pandas 2.

## Decisions

- The xtmgah forks stay the source of truth and are patched in place on a `python313` branch each,
  using the upstream SigProfilerSuite / AlexandrovLab repositories as the reference for fixes. Custom
  fork changes are identified against the upstream base version first and preserved.
- The extraction image moves from `nvidia/cuda:12.1.1-runtime-rockylinux9` to an
  `nvidia/cuda:12.x-runtime-amzn2023` base, because Rocky Linux 9 AppStream ships Python 3.11, 3.12
  and 3.14 but not 3.13, while Amazon Linux 2023 ships `python3.13` packages.
- Regression baselines are captured on the current images before the Dockerfiles change
  (see `utils/regression/README.md`).

## Which features exercise the Python packages

Only two application features reach Python. See `utils/regression/README.md` for the full map.

- Visualization, user data upload: `mSigPortal_Profiler_Extraction.py` with SigProfilerMatrixGenerator
  (vcf, maf, csv, tsv), sigProfilerPlotting (catalog inputs), SigProfilerClusters plus
  SigProfilerSimulator (vcf with Cluster on).
- Extraction page: `mSigPortal-SigProfilerExtractor.py` with SigProfilerExtractor and, inside it,
  torch, SigProfilerAssignment, sigProfilerPlotting, PyPDF2 and reportlab.

## Per-package upgrade notes

| Package | Fork version | Document | Branch |
|---|---|---|---|
| SigProfilerPlotting | 1.1.6 | [python313-sigprofilerplotting.md](python313-sigprofilerplotting.md) | `xtmgah/SigProfilerPlotting@python313` |
| SigProfilerMatrixGenerator | 1.1.30 | [python313-sigprofilermatrixgenerator.md](python313-sigprofilermatrixgenerator.md) | `xtmgah/SigProfilerMatrixGenerator@python313` |
| SigProfilerExtractor | 1.1.13 | [python313-sigprofilerextractor.md](python313-sigprofilerextractor.md) | `xtmgah/SigProfilerExtractor@python313` |
| SigProfilerClusters | 1.0.11 | [python313-sigprofilerclusters.md](python313-sigprofilerclusters.md) | `xtmgah/SigProfilerClusters@python313` |

## Application changes in this repository

- `server/services/api/visualization/profilerExtraction.js` and
  `extraction-service/services/extraction.js` spawn the interpreter named by the `PYTHON_BIN`
  environment variable, defaulting to `python3`. The Dockerfiles set `PYTHON_BIN=python3.13` so the
  Amazon Linux system `python3` (3.9) is left untouched, as AL2023 requires.
- `docker/backend.dockerfile` and `docker/extraction-batch.dockerfile`: see the Dockerfile section
  below (filled in after the package branches are final).
- `extraction-service/services/python/mSigPortal_Profiler_Extraction.py` (an unused, stale 2023 copy
  of the backend script; gitignored) was removed from the local checkout. The private scripts repo
  copy step in `.github/workflows/deploy-extraction.yml` never copied it, so nothing changes in CI.
- `utils/regression/`: regression harness used to compare the old and new images.

## Dockerfile changes

### `docker/backend.dockerfile`

| Before | After |
|---|---|
| `dnf install python3-devel python3-pip` (system 3.9) | `dnf install python3.13 python3.13-devel python3.13-pip`; system `python3` untouched |
| `pip3 install -e git+.../SigProfilerClusters` first (pulled PyPI SPMG/SPP/Extractor/torch, then overwritten) | order from [python313-sigprofilerclusters.md](python313-sigprofilerclusters.md#install-order-for-the-dockerfile): SPP fork, SPMG fork, `--no-deps SigProfilerSimulator==1.2.2`, `fastrand>=1.2`, `--no-deps` Clusters fork, then `statsmodels scipy` |
| fork `master`/`main` branches | `@python313` branches |
| `seaborn==0.13.2`, then `--force-reinstall numpy==1.26.4 pandas==1.3.5` | `numpy>=2.1,<3 pandas>=2.2.3,<3 matplotlib>=3.9 seaborn==0.13.2` (build arg `PY_PINS`) on every install step |
| `echo ... MutationMatrixGenerator >> scripts/__init__.py` | removed; the SPMG branch ships `scripts/MutationMatrixGenerator.py` |
| (none) | `ENV PYTHON_BIN=python3.13` |
| `GENOME_PATH=/src/sigprofilermatrixgenerator/.../tsb` | unchanged: editable `git+` installs still clone into `/src/<egg name lowercased>` (pip's default `src` dir under the working directory `/`), and SPMG resolves `references/` from its source checkout |

The `cp fonts/* /usr/share/fonts` step stays. SPP now registers its bundled fonts itself, but R plots
in the same image may rely on the system copies.

### `docker/extraction-batch.dockerfile`

| Before | After |
|---|---|
| `nvidia/cuda:12.1.1-runtime-rockylinux9` | `nvidia/cuda:12.9.2-runtime-amzn2023` |
| `dnf -y update` | `dnf -y --releasever=latest update`: the nvidia amzn2023 bases are locked to an older AL2023 release with no `python3.13` or `nodejs24` packages |
| `dnf module enable nodejs:20`, `nodejs npm` | `nodejs24` (same package as the backend image; Node 20 is end of life) |
| `python3-devel python3-pip python3-setuptools` | `python3.13 python3.13-devel python3.13-pip` |
| `torch==2.1.0` (PyPI default, cu121) | `torch>=2.5,<3` from `https://download.pytorch.org/whl/cu129` |
| Extractor fork `main`; `PyPDF2==2.11.2 SigProfilerAssignment==0.0.14 sigProfilerPlotting==1.2.2`; upstream SPMG from git; forced `numpy==1.26.4 pandas==1.3.5 reportlab<5.1` | one install: Extractor fork `@python313` (editable) + `SigProfilerAssignment>=1.1.5,<1.2 sigProfilerPlotting>=1.4.3,<1.5 SigProfilerMatrixGenerator>=1.4,<1.5 pypdf>=6,<7 reportlab>=4,<6 numpy>=2.1,<3 pandas>=2.2.3,<3 matplotlib>=3.9,<4` ([python313-sigprofilerextractor.md](python313-sigprofilerextractor.md#recommended-dockerfile-pins-and-order-extraction-batch-image)) |
| `LD_LIBRARY_PATH` += `/usr/local/lib/python3.9/site-packages/nvidia/{cudnn,cuda_cupti}/lib/` | `/usr/local/lib/python3.13/site-packages/nvidia/{cudnn,cuda_cupti}/lib/` (checked in the built image) |
| (none) | `ENV PYTHON_BIN=python3.13` |
| `GENOME_PATH=/deploy/app/src/sigprofilermatrixgenerator/...` (SPMG was an editable git install) | `/usr/local/lib/python3.13/site-packages/SigProfilerMatrixGenerator/references/chromosomes/tsb`: SPMG is now a regular PyPI install (only used for `input_type=vcf`, which the UI never sends) |

CUDA choice: 12.9, matching `stage`, which moved the current image to
`nvidia/cuda:12.9.2-runtime-rockylinux9` (#125), so it relies on the same CUDA 12.9-capable driver on
the Batch GPU hosts. torch wheels bundle their own CUDA user-space libraries; the cu129 index carries
cp313 builds (resolved to torch 2.13.0+cu129 on 2026-11-03). The first local build used 12.6.3 / cu126
and passed the same checks. Merging into
`stage` will conflict in `docker/extraction-batch.dockerfile`; take this branch's version.

### Local build

The extraction image was first built before the fork branches were pushed, from a copy of the
Dockerfile with the `git+https://github.com/xtmgah/<Repo>` URLs pointed at a local `git daemon`
(`~/tmp/py313/docker/build.sh`). After the push on 2026-11-02 both images build from the unmodified
Dockerfiles.

## Running the regression suite

Full details in `utils/regression/README.md`. Short form, from the repository root (any host Python
3.9+, e.g. `python3.13`; the harness is stdlib only):

```bash
# 1. baselines: build from dev's Dockerfiles (on python313, docker/*.dockerfile are already the 3.13 ones)
git show dev:docker/backend.dockerfile | docker build -t msigportal-backend:py39 -f - .
git show dev:docker/extraction-batch.dockerfile | docker build -t msigportal-extraction:py39 -f - .
python3 utils/regression/run_cases.py profiler  --label baseline --all --image msigportal-backend:py39
python3 utils/regression/run_cases.py extractor --label baseline --image msigportal-extraction:py39

# 2. candidates on this branch
docker build -t msigportal-backend:py313 -f docker/backend.dockerfile .
docker build -t msigportal-extraction:py313 -f docker/extraction-batch.dockerfile .
python3 utils/regression/run_cases.py profiler  --label candidate --all --image msigportal-backend:py313
python3 utils/regression/run_cases.py extractor --label candidate --image msigportal-extraction:py313

# 3. compare
python3 utils/regression/compare.py data/regression/baseline data/regression/candidate --json data/regression/report.json
```

Expected: matrices and tables `identical` or `tolerated`; images `visual`; logs with deprecation counts
but no new tracebacks. Known, accepted differences from the package docs:

- Extractor de novo signatures differ by up to ~4e-6 (torch convergence check), so some de novo files
  will be `FAIL` at the 1e-6 tolerance; COSMIC decomposition tables are identical.
- Extractor: SigProfilerAssignment 1.x renames de novo labels (`96-A` → `96-SBS96A`), renames
  `Decomposed_Mutation_Probabilities.txt`, adds `JOB_METADATA_SPA.txt`, and changes the assignment log
  format. These show as `missing`/new files.
- Clusters output is nondeterministic run to run on every version.
- On an Apple-silicon host the extraction image runs under amd64 emulation and is slow.

## Regression results

### Profiler suite (backend image), 2026-11-04

Run locally (Apple silicon, Docker Desktop 8 GB). Baseline `msigportal-backend:py39` built from `dev`
(Python 3.9.25, numpy 1.26.4, pandas 1.3.5, fork default branches); candidate
`msigportal-backend:py313` (Python 3.13.16, numpy 2.5.4, pandas 2.3.3, forks at `python313`, Clusters
`4354a91`). Report: `data/regression/report-profiler.json`. **All 15 cases pass.**

| Case(s) | Result |
|---|---|
| vcf GRCh37/GRCh38/mm10, maf, tsv split-filters, bed, seqinfo, wes | every matrix, statistics and SeqInfo file `identical` |
| csv_grch37_collapse | matrices `identical`; the converted `*_mSigPortal_{SNV,INDEL}[_Collapse].txt` inputs hold the same rows in a different order (`tolerated`) |
| catalog tsv/csv, CNV48, RS32 | matrix/list files `identical`; plots (PDF/SVG) `visual`, same page counts |
| vcf_grch37_cluster, vcf_grch38_cluster | matrices `identical`; per-sample VCFs differ only in the `##bcftools` header (`tolerated`); `*_clustered_class_All.txt` `varies`: 99.8% / 99.7% of (sample, chr, start, ref, alt, subclass) rows agree, within Clusters' run-to-run variation plus the stdev and write-race fixes |
| all logs | no tracebacks; error/warning counts equal baseline vs candidate |

`compare.py` changes made for this run: paths containing the run label (`/data/regression/<label>/`)
are normalized; table rows are compared order-independently; `.vcf` files are compared as text with
`##bcftools` provenance lines ignored; Cluster result tables get the `varies` verdict.

### Extractor suite (extraction image)

Not run locally: the old image's torch 2.1.0 fails in `nmf_cpu.py` (`could not create a primitive
descriptor for a matmul primitive`) under amd64 emulation on Apple silicon. To be run on x86_64
(see Push checklist).

## Push checklist

1. Done 2026-11-02: `python313` pushed to xtmgah/SigProfilerPlotting (`7c3b1f2`),
   SigProfilerMatrixGenerator (`caace75`), SigProfilerClusters (`6b8a7ff`) and SigProfilerExtractor
   (`e97d4db`). The Dockerfiles install from these branches, so any later fork commit changes the next
   image build.
2. `git push -u origin python313` in this repository.
3. Build both images, run the regression suite, deploy to the dev tier, run the UI checks listed in
   `utils/regression/README.md`.
