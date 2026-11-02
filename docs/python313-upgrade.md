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

To be completed once the sub-agent reports provide the working dependency versions.

## Regression results

To be filled in after `utils/regression/compare.py` is run against the baseline and candidate images.

## Push checklist

1. `gh auth login -h github.com` (the stored token was invalid on 2026-10-28).
2. In each of `src/sigprofilerplotting`, `src/sigprofilermatrixgenerator`, `src/sigprofilerclusters`,
   `src/sigprofilerextractor`: `git push -u origin python313`.
3. `git push -u origin python313` in this repository.
4. Build both images, run the regression suite, deploy to the dev tier, run the UI checks listed in
   `utils/regression/README.md`.
