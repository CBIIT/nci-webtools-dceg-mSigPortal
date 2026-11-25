# Python 3.13 upgrade: per-package sub-agent plan (revision 2)

## Context

Revision of the approved plan at `~/.claude/plans/i-need-to-upgrade-staged-snail.md`. The package-to-feature map and the regression harness (`utils/regression/`) are done. Baselines have not been captured yet. This revision adds the user's new requirements for the package upgrade phase:

1. Each xtmgah fork is upgraded by its own sub-agent, using the upstream repo as the reference for how numpy 2 / pandas 2 / matplotlib 3.9+ / Python 3.13 issues were fixed.
2. Every custom change in the fork (the delta against the upstream version it was cut from) is identified first and preserved.
3. Work is committed to a new branch named `python313` in each fork repo.
4. Sub-agent model is chosen per package complexity: Fable 5.1 1M for the two hardest, Opus 5.5 1M for the other two.
5. Every change is documented in markdown under `docs/` in this project (the folder exists and is empty).

Decisions already made and unchanged: patch forks in place (not move to upstream), extraction image moves to `nvidia/cuda:12.x-runtime-amzn2023`, baselines before Dockerfile changes.

## Environment facts that constrain execution

- Local fork clones: `src/sigprofilerplotting` (branch master), `src/sigprofilermatrixgenerator` (master), `src/sigprofilerclusters` (main). All have origin `https://github.com/xtmgah/<repo>` and a few dirty build artifacts (egg-info, install.log, version.py) that must not be committed. `packages/SigProfilerExtractor/SigProfilerExtractor` is a plain vendored copy, not a git repo; the Extractor sub-agent clones `https://github.com/xtmgah/SigProfilerExtractor` into `src/sigprofilerextractor` (gitignored via `/src`).
- Host Pythons: 3.9, 3.10, 3.11, 3.12, 3.14 and `uv`. No 3.13 installed. `uv` needs `UV_CACHE_DIR=$TMPDIR/uv-cache` and `UV_PYTHON_INSTALL_DIR=$TMPDIR/uv-python` in this sandbox; `uv python install 3.13` downloads from GitHub. Fallback if the download is blocked: test on python3.12 with numpy>=2.1 / pandas>=2.2.3 (same API surface) and syntax-check on 3.14.
- Network: `api.github.com` was denied by the user in the last command batch; `github.com` and `pypi.org` were reachable earlier. Sub-agents need `git clone`/`git fetch` from github.com and `pip`/`uv` from pypi.org and files.pythonhosted.org. `gh` token is invalid, so pushing needs `gh auth login` or the osxkeychain git credential first.
- `packages/SigProfilerExtractor/Reference_Signatures/GRCh37/COSMIC_v3.3_{SBS,DBS,ID}_GRCh37.txt` already have the `Type` + signature-column format the extractor needs. The regression README's psql export becomes optional; copy these as `data/regression/inputs/signature_{SBS96,DBS78,ID83}.tsv`.
- Docker is not reachable from this sandbox. Image builds and regression runs are commands handed to the user.

## Sub-agent design

### Common brief (identical for all four, parameterized)

Each sub-agent gets this contract in its prompt:

1. **Workspace**: `cd <clone>`; `git fetch origin`; `git stash --include-untracked` or discard build artifacts only (never commit egg-info, dist, install.log, version.py bumps); `git checkout -b python313 origin/<default>`.
2. **Upstream reference**: clone upstream into `$TMPDIR/upstream-<pkg>` (SigProfilerSuite/SigProfilerPlotting, SigProfilerSuite/SigProfilerMatrixGenerator, AlexandrovLab/SigProfilerClusters, AlexandrovLab/SigProfilerExtractor; if an AlexandrovLab repo redirects to SigProfilerSuite use that). Find the upstream tag matching the fork's `setup.py` version (SPP 1.1.6, SPMG 1.1.30, Clusters 1.0.11, Extractor 1.1.13) and produce `fork_delta.diff = diff -ru upstream@<tag> fork@HEAD` excluding build artifacts. That diff is the list of custom changes to preserve; summarize each hunk in the doc.
3. **Compatibility pass**: grep the fork for the known breakers (list below) and for every hit look at how upstream HEAD rewrote the same code; apply the minimal equivalent change without reverting any hunk in `fork_delta.diff`. Relax `install_requires` caps; add `python_requires=">=3.9"`; do not change the package version number unless the Dockerfile pin needs it (it does not; Dockerfiles will reference the branch).
4. **Verify**: create a 3.13 venv (`uv venv --python 3.13 $TMPDIR/venv-<pkg>`), `pip install -e .` with numpy>=2.1 pandas>=2.2.3 matplotlib>=3.9, then run the package-specific smoke tests below with `-W error::DeprecationWarning` where feasible. Record versions used.
5. **Commit** on `python313` with a message per logical change. Do not push; report the branch SHA. Push happens in the finalization step once auth is refreshed.
6. **Document**: write `docs/python313-<package>.md` with sections: fork origin and upstream base tag; custom changes preserved (from `fork_delta.diff`); compatibility changes made with file:line and the upstream commit or line used as reference; dependency pin changes; how it was tested and results; open risks. Keep tables short; link upstream commits by URL.
7. **Report back**: branch name, SHA, files changed, test results, anything that could not be verified.

Known breakers to grep for (from the earlier fork analysis): `font_manager._rebuild`, `np.mat(` / `np.matrix`, `DataFrame.append(`, `.iteritems(`, `np.int`/`np.float`/`np.bool`/`np.object`/`np.NaN`/`np.product`/`np.in1d`, `delim_whitespace=`, `error_bad_lines=`, `line_terminator=`, `squeeze=`, `applymap`, `matplotlib.cm.get_cmap`, `register_cmap`, `distutils`, `pkg_resources`, `imp`, `collections.Mapping`, `scipy.misc`, `sklearn.externals`, `torch.FloatTensor`/`torch.cuda.set_device` legacy use, `setup.py` `install_requires` caps.

### Per-package assignments

| Package | Clone | Model | Why | Package-specific verification |
|---|---|---|---|---|
| SigProfilerPlotting (xtmgah, 1.1.6, 68 ahead / 278 behind SigProfilerSuite) | `src/sigprofilerplotting` | **fable** | 5k-line monolith, largest custom delta (SBS288, 384_extended, zero-mutation handling, bundled fonts), matplotlib API drift is the main 3.13 risk, `_rebuild()` in `setup.py:59` | run `plot_example.py`; call `plotSBS/plotDBS/plotID/plotCNV/plotSV` with `percentage=True` on `client/public/assets/exampleInput/demo_input_catalog.tsv` and the CNV48/RS32 `.all` files; confirm bundled font is registered via `font_manager.fontManager.addfont` and appears in one rendered PNG |
| SigProfilerExtractor (xtmgah, 1.1.13, standalone repo) | `src/sigprofilerextractor` (new clone) | **fable** | NMF math on `np.mat` in `nmf_cpu.py`/`nmf_gpu.py` (matrix `*` semantics), torch >=2.5 migration, matplotlib `<=3.4.3` cap, PyPDF2/reportlab/SigProfilerAssignment coupling, GPU path cannot be run locally | CPU run of `sigProfilerExtractor('matrix', ..., 'packages/SigProfilerExtractor/Demo_input_Test/vcftest.SBS96.all', minimum_signatures=1, maximum_signatures=2, nmf_replicates=3, ...)` with `Reference_Signatures/GRCh37/COSMIC_v3.3_SBS_GRCh37.txt`; unit-compare `nmf_cpu` before/after on a fixed seed (store old outputs first); import `nmf_gpu` and exercise its code paths with `torch` CPU tensors where the class allows; pick SigProfilerAssignment version that imports and runs (try newest 0.0.x, then 1.1.x) |
| SigProfilerMatrixGenerator (xtmgah, 1.1.30, 7 ahead / 247 behind) | `src/sigprofilermatrixgenerator` | **opus** | moderate: pandas/numpy drift in `SigProfilerMatrixGeneratorFunc.py` and CNV generator; add `scripts/MutationMatrixGenerator.py` alias module so the Dockerfile echo patch can be removed | `SigProfilerMatrixGeneratorFunc('t', 'GRCh37', <dir with demo_input_multi.vcf>, plot=False, seqInfo=True)` using `data/genomes/GRCh37` symlinked into the package `references/chromosomes/tsb`; byte-compare the SBS96/DBS78/ID83 matrices against a pre-change run in the same venv on python3.12 |
| SigProfilerClusters (xtmgah, 1.0.11, standalone, 3 commits: VAF calculation change) | `src/sigprofilerclusters` | **opus** | small code base; main work is `matplotlib>=3.3,<3.5` cap, pandas drift in `classifyFunctions.py`/`hotspot.py`, and making it resolve against the upgraded fork SPP/SPMG plus PyPI SigProfilerSimulator; depends on wave 1 | install fork SPP and SPMG `python313` branches first; `hp.analysis(...)` on a single-sample VCF from `demo_input_multi.vcf` with `simulations=5` to keep runtime minutes; confirm `clustered_class_All.txt` columns match a python3.12/old-pins run |

Waves: wave 1 runs SigProfilerPlotting, SigProfilerExtractor, SigProfilerMatrixGenerator in parallel. Wave 2 runs SigProfilerClusters after SPP and SPMG report, because its tests import them. The orchestrator (this session) reviews each report, reads the docs file, and spot-checks `git diff origin/<default>..python313` for accidental reverts of fork hunks before starting wave 2.

Each sub-agent is launched with `Agent(subagent_type="general-purpose", model="fable"|"opus", run_in_background=true)` and a prompt containing the common brief, its row of the table, the breaker list, and the environment facts (uv cache dirs, no Docker, network hosts to allow).

## After the sub-agents: integration (this session)

1. **Dockerfiles** on a new msigportal branch `python313` (cut from `dev`):
   - `docker/backend.dockerfile`: `dnf install python3.13 python3.13-pip python3.13-devel`; use `python3.13 -m pip`; install forks with `-e 'git+https://github.com/xtmgah/<repo>@python313#egg=<name>'`; install Clusters with `--no-deps` then explicit deps so PyPI SPMG/SPP do not shadow the forks; pins `numpy>=2.1,<3 pandas>=2.2.3,<3 matplotlib>=3.9 seaborn==0.13.2`; drop the line-34 alias echo if SPMG ships the alias module; recheck `GENOME_PATH`.
   - `docker/extraction-batch.dockerfile`: base `nvidia/cuda:12.6.3-runtime-amzn2023` (or 12.8.2/12.9.2 depending on the Batch host driver); python3.13 packages; `torch>=2.5` from the matching cu12x index; Extractor from `@python313`; `SigProfilerAssignment` and `sigProfilerPlotting` versions chosen by the Extractor sub-agent; fix `LD_LIBRARY_PATH` to `python3.13/site-packages/nvidia/...`; recheck `GENOME_PATH`.
   - Node spawns: add `PYTHON_BIN` env (default `python3`) read at [profilerExtraction.js:111](server/services/api/visualization/profilerExtraction.js:111) and [extraction.js:233](extraction-service/services/extraction.js:233); Dockerfiles set `PYTHON_BIN=python3.13`.
2. **Scripts** (local gitignored copies; user syncs to the private repo): no changes expected; any fix found by regression lands in `server/services/python/` and the matching `extraction-service/` copy. Delete the stale unused `extraction-service/services/python/mSigPortal_Profiler_Extraction.py`.
3. **Regression harness**: update `utils/regression/README.md` and `cases.json` to source signature tables from `packages/SigProfilerExtractor/Reference_Signatures/GRCh37/` (copy step in `run_cases.py`), keeping psql as the alternative.
4. **Docs**: `docs/python313-upgrade.md` as the index: goals, decisions, per-package links, Dockerfile diffs, how to run the regression suite, results table filled in after the user runs it, and a push checklist.
5. **Push**: after the user runs `gh auth login` (or confirms the keychain credential), push `python313` in the four fork repos, then push the msigportal branch. Dockerfile builds depend on the fork branches being on GitHub.

## Verification

1. Per package: sub-agent smoke tests above in a 3.13 venv (or 3.12 + numpy 2 fallback, stated in the doc).
2. User captures baselines on current images (commands in `utils/regression/README.md`), then builds the `python313` images and runs `run_cases.py ... --label candidate` and `compare.py`. Expected: matrices/tables `identical` or `tolerated`; images `visual`; logs with deprecation counts but no new tracebacks.
3. Deploy to dev tier; one UI run each of Visualization upload (vcf with Cluster, catalog_tsv) and Extraction SBS96 end to end, including the auto-created Exploration sessions.
4. Results recorded in `docs/python313-upgrade.md` before promotion.
