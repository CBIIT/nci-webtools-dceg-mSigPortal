# Python 3.13 upgrade: SigProfilerMatrixGenerator (xtmgah fork)

## Summary

The xtmgah fork of SigProfilerMatrixGenerator was updated so it installs and
runs on Python 3.13 with numpy >= 2.1, pandas >= 2.2.3 and matplotlib >= 3.9,
while keeping every mSigPortal-specific customization. The work is on branch
`python313` of the fork clone. All changes are small, targeted compatibility
fixes plus a new alias module; no behavior relevant to mSigPortal changed, and
the package version was not bumped.

Code areas touched:

| Area | Why |
|------|-----|
| `DataFrame.append` -> `pd.concat` | pandas 2.0 removed `DataFrame.append` |
| `scipy.stats.binom_test` -> shim over `binomtest` | removed in SciPy 1.12 |
| invalid regex escape `'\/'` -> raw strings | SyntaxWarning on Python 3.12+ |
| new `MutationMatrixGenerator.py` alias | lets Simulator import the module and removes a Dockerfile patch |
| `setup.py` `python_requires` | resolve modern wheels on 3.13 |
| `DataFrame.from_dict(DataFrame)` -> `DataFrame(...)` | Pandas4Warning on pandas 3 (added in the re-verification) |

Verified on the real target stack (Python 3.13.16, numpy 2.1.3, pandas 2.2.3,
matplotlib 3.10.9) and on the latest releases (numpy 2.5.4, pandas 3.1.0,
matplotlib 3.11.2). Every matrix, `vcf_files` seqInfo file and `Statistics.txt`
is byte-identical to a Python 3.9 / numpy 1.26.4 / pandas 1.3.5 run of
`origin/master` in all four scenarios. Matrices stay integer (no `.0`).

The first pass of this work ran without network access and was only checked on
a substitute Python 3.9 / pandas 1.5.3 stack; that is superseded by the
Verification section below.

## Fork origin and upstream base

| Item | Value |
|------|-------|
| Fork clone | `/Users/chop2/Projects/msigportal/src/sigprofilermatrixgenerator` |
| Fork origin | `https://github.com/xtmgah/SigProfilerMatrixGenerator` (default branch `master`) |
| Fork HEAD at start | `0e5cd42` "updated to v2.0" (Tongwu Zhang, 2021-11-05) |
| setup.py VERSION | 1.1.30 (unchanged) |
| Upstream parent | `https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator` |
| Upstream base | [`585cc62`](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator/commit/585cc62b0cecd1a97c8cd0c86f0de6410031e14d) "Merge branch 'genome_updates'" (2021-05-28), `setup.py` VERSION 1.1.30 |

Upstream has no 1.1.30 tag (tags jump from `v1.2` to `v1.2.18`). The base is the
`git merge-base` of upstream `master` and fork `master`; its `setup.py` says
1.1.30, matching the fork. The fork branches were fetched into
`~/tmp/py313/upstream-spmg` from the local clone.

Delta files (excluding `*.egg-info`, `dist/`, `build/`, `z-bak/`,
`references/`, `install.log`, `__pycache__`):

| File | Range | Content |
|------|-------|---------|
| `~/tmp/py313/spmg-fork_delta.diff` | `585cc62..origin/master` | 1 file, +18/-1: the `Statistics.txt` block |
| `~/tmp/py313/spmg-fork313_delta.diff` | `585cc62..python313` | the above plus this branch's 5 files |

`git diff origin/master python313` touches only the 5 files listed under Branch
and commits; `SigProfilerMatrixGeneratorFunc.py` is not among them, so the fork
hunk is not reverted. The fork also tracks
`SigProfilerMatrixGenerator/z-bak/SigProfilerMatrixGeneratorFunc.py` (added in
fork commit `0e5cd42`, not by this branch).

## Custom changes preserved

The upstream delta has exactly one fork customization:

| Customization | Location | Notes |
|---------------|----------|-------|
| mSigPortal `Statistics.txt` output | `SigProfilerMatrixGeneratorFunc.py` end of `SigProfilerMatrixGeneratorFunc` (the "Modified By Jian Sang" block) | writes `Statistics.txt` with the "Total of N SNVs, N DINUCs, N INDELs" summary into the input directory. Fork commit `0e5cd42`. mSigPortal reads this file. Unchanged on `python313`; byte-identical output verified. |

Correction to the first pass: the SigProfilerPlotting version log line
(`SigProfilerMatrixGeneratorFunc.py:317` upstream), the `df2csv` matrix writer
(`SigProfilerMatrixGenerator.py:29` upstream) and `scripts/CNVMatrixGenerator.py`
are all present in upstream `585cc62`. They are upstream 1.1.30 code, not fork
customizations. They are still unchanged on `python313`.

No customization conflicted with an upstream fix.

## Compatibility changes made

Line numbers are in the final committed `SigProfilerMatrixGenerator.py`.
Upstream references (all in SigProfilerSuite/SigProfilerMatrixGenerator):

| Change | Upstream commit |
|--------|-----------------|
| `append` -> `pd.concat` | [`0983bad`](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator/commit/0983badb3add7a18441211d5c9e4b72193ee8c9d) "Upgrade to pandas and numpy >= 2.0.0" (2025-02-03): same `trans = pd.concat([trans, a, b])` rewrite of the strand-bias collapse code in `MutationMatrixGenerator.py` |
| `binom_test` -> `binomtest` | [`f6cddc9`](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator/commit/f6cddc98dc2570d8c2f0e98472dd1159fddad831) "Update scipy dependency and binom_test changes" (2024-04-25): `stats.binomtest(k=..., n=...).pvalue` |
| `'\/scripts$'` regex | [`5b7105f`](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator/commit/5b7105f6f19041ca87929219b2d717e1b0481a87) "Change how fasta and tsb folders are referred to" (2023-06-09) removed the regex in favor of `ref_install` path helpers; a raw string is the minimal equivalent here |
| `MutationMatrixGenerator` module name | [`507483d`](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator/commit/507483d50563b7bbbb4ae0bf6bbcc624366cf312) "CLI Updates ... Rename SigProfilerMatrixGenerator.py to MutationMatrixGenerator.py" (2023-09-02) |
| `python_requires=">=3.9"` | [`66c9859`](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator/commit/66c98594d31b5e30a9ae2066b02ef4b37e1ca282) "Upgrade to v1.2.14" (2023-03-07); still `>=3.9` at HEAD `e67c5c4` (v1.4.0) |
| `DataFrame.from_dict(DataFrame)` | no upstream fix: HEAD `e67c5c4` (v1.4.0) still has `pd.DataFrame.from_dict(mut_4992)` at `MutationMatrixGenerator.py:3863` |

| File:line | Change | Reason / upstream equivalent |
|-----------|--------|------------------------------|
| `scripts/SigProfilerMatrixGenerator.py` (18 sites: 1878, 1879, 1882, 1898, 1899, 1902, 1915, 1916, 1919, 2313, 2314, 2317, 2487, 2488, 2491, 2504, 2505, 2508) | `X = X.append(a)` -> `X = pd.concat([X, a])` | pandas 2.0 removed `DataFrame.append`. Upstream replaced the same strand-bias collapse code with `pd.concat`. These sites are inside `matrix_generator` (1758), `matrix_generator_INDEL` (2186) and `matrix_generator_DINUC` (2380) and run on the standard mSigPortal path (SBS 384->288, 24->18, 6144->4608; DBS 186->150, 2976->2400; INDEL TSB collapse). |
| `scripts/SigProfilerMatrixGenerator.py:29` (new `binom_test` helper); call sites 596, 606, 1285, 1295, 2041, 2090, 2092, 2094 | add module-level `binom_test(counts, p=0.5)` wrapping `stats.binomtest(...).pvalue`; route `stats.binom_test([x, y])` call sites through it | `scipy.stats.binom_test` was deprecated in SciPy 1.10 and removed in 1.12. Upstream switched to `scipy.stats.binomtest`. The shim reproduces the old `binom_test([successes, failures])` semantics (and returns 1.0 for the empty `[0, 0]` case). These call sites only execute with `tsb_stat=True` or `gs=True`, both of which mSigPortal does not set. |
| `scripts/save_chrom_tsb_separate.py:93` | `re.sub('\/scripts$', ...)` -> `re.sub(r'/scripts$', ...)` | Python 3.12 raises `SyntaxWarning: invalid escape sequence '\/'`. Install-time helper only. Regex is unchanged. |
| `scripts/save_context_distribution.py:549` | same raw-string fix | same reason. |
| `scripts/SigProfilerMatrixGenerator.py:2416` (`matrix_generator_DINUC`) | `pd.DataFrame.from_dict(mut_4992)` -> `pd.DataFrame(mut_4992)` | pandas 3 emits `Pandas4Warning: Passing a DataFrame to DataFrame.from_dict is deprecated` on every run with DBS output. For the default `orient='columns'`, `from_dict` calls the constructor, so behavior is the same. Not hit on pandas 2.2.3. Commit `caace75`. |

Scanned and found clean (no change needed): `np.int` / `np.float` / `np.bool` /
`np.object` / `np.NaN` / `np.product` / `np.in1d` / `np.mat`, `DataFrame.iteritems`,
`applymap`, `.ix`, `delim_whitespace=` / `error_bad_lines=` / `line_terminator=`
/ `squeeze=` in `read_csv`, `matplotlib.cm.get_cmap` / `register_cmap`,
`distutils` / `pkg_resources` / `import imp`, `collections.Mapping`,
`scipy.misc`, `sklearn.externals`. Chained assignments found
(`gene_counts[gene][sample] = ...`, etc.) are all on plain Python dicts, not
pandas objects, so pandas Copy-on-Write does not affect them.

Out of scope: `install.py` and `save_context_distribution.py` use `os.system`
with `wget` / `rsync` / `gzip` / `tar` / `sort` to download and build genome
references. This is the genome-install path, which mSigPortal does not use at
runtime (genomes are mounted from the host into
`references/chromosomes/tsb`). These calls are not Python-version dependent and
were left unchanged.

## MutationMatrixGenerator alias and the Dockerfile patch

Upstream renamed `scripts/SigProfilerMatrixGenerator.py` to
`MutationMatrixGenerator.py`. Newer SigProfilerSimulator and SigProfilerClusters
import the module as
`from SigProfilerMatrixGenerator.scripts import MutationMatrixGenerator as matGen`.
This fork keeps the original file name, so the backend Dockerfile currently
patches around it:

```
# docker/backend.dockerfile (current, lines 33-34)
RUN echo 'from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGenerator as MutationMatrixGenerator' \
    >> /src/sigprofilermatrixgenerator/SigProfilerMatrixGenerator/scripts/__init__.py
```

A new real submodule `scripts/MutationMatrixGenerator.py` now re-exports the
full public namespace from `SigProfilerMatrixGenerator`:

```python
from SigProfilerMatrixGenerator.scripts.SigProfilerMatrixGenerator import *
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGenerator as _spmg
globals().update({k: v for k, v in vars(_spmg).items() if not k.startswith("__")})
```

With this file present,
`from SigProfilerMatrixGenerator.scripts import MutationMatrixGenerator`
imports a genuine module that exposes every function Simulator/Clusters call by
attribute (`catalogue_generator_single`, `catalogue_generator_INDEL_single`,
`matrix_generator`, `matrix_generator_INDEL`, `matrix_generator_DINUC`,
`exome_check`, `panel_check`, `reference_paths`, `gene_range`, `BED_filtering`,
`df2csv`, `perm`, ...). Verified on Python 3.13 with PyPI SigProfilerSimulator
1.2.2 (see Verification).

Result: the Dockerfile `echo ... >> scripts/__init__.py` patch (lines 33-34) is
no longer needed and can be deleted.

SigProfilerSimulator 1.2.2 declares `SigProfilerMatrixGenerator>=1.3.2` and
`sigProfilerPlotting>=1.4.0`. A plain `pip install` of Simulator, or of
Clusters (which pulls Simulator), would replace both forks with PyPI releases.
The Dockerfile must install Clusters/Simulator with `--no-deps` (as the plan
already says) and then install `fastrand` explicitly.

## Dependency changes (setup.py)

Only `python_requires` was added. No caps existed to relax, and the floors are
left permissive so the resolver picks modern wheels on 3.13. The
`sigProfilerPlotting` floor is deliberately kept at `>=1.0.1` so the xtmgah
SigProfilerPlotting fork (1.1.6) that the backend image installs on top still
satisfies it.

| Field | Before | After |
|-------|--------|-------|
| `python_requires` | (absent) | `>=3.9` |
| `matplotlib` | `>=2.2.2` | `>=2.2.2` |
| `sigProfilerPlotting` | `>=1.0.1` | `>=1.0.1` |
| `statsmodels` | `>=0.9.0` | `>=0.9.0` |
| `scipy` | `>=1.1.0` | `>=1.1.0` |
| `numpy` | `>=1.18.5` | `>=1.18.5` |
| `pandas` | `>=0.23.4` | `>=0.23.4` |

Tested resolutions on Python 3.13 (`uv pip install` of the SPP fork
`python313`, then this branch):

| Package | Pinned stack (recommended) | Latest stack |
|---------|----------------------------|--------------|
| numpy | 2.1.3 | 2.5.4 |
| pandas | 2.2.3 | 3.1.0 |
| matplotlib | 3.10.9 | 3.11.2 |
| scipy | 1.18.1 | 1.18.1 |
| statsmodels | 0.15.0 | 0.15.0 |
| pyparsing | 3.3.3 | 3.3.3 |
| seaborn | 0.13.2 | 0.13.2 |
| sigProfilerPlotting | 1.1.6 (xtmgah fork, `python313`) | same |

Recommended Dockerfile pins for this package: same as SigProfilerPlotting
(`numpy==2.1.3 pandas==2.2.3 matplotlib==3.10.9 seaborn==0.13.2`). The plan's
`numpy<3 pandas<3` caps are not needed by this package (pandas 3.1.0 also
passes), but are harmless. Full freezes: `~/tmp/py313/venv-spmg.freeze.txt`,
`venv-spmg-latest.freeze.txt`, `venv-spmg39.freeze.txt`.

## Verification

### Environments

All venvs are under `~/tmp/py313/`, installed non-editable from `git archive`
exports (`spmg-src-313`, `spmg-src-39`, `spp-src-313`, `spp-src-39`), so the
clones were not touched. In each venv,
`site-packages/SigProfilerMatrixGenerator/references/chromosomes/tsb` is a
symlink to `data/genomes`. Setup scripts: `spmg-setup.sh`, `spmg-setup313.sh`.

| Env | Interpreter | numpy | pandas | scipy | matplotlib | statsmodels | SPMG source | SPP source |
|-----|-------------|-------|--------|-------|------------|-------------|-------------|------------|
| `venv-spmg39` (baseline) | CPython 3.9.25 x86_64 (Rosetta) | 1.26.4 | 1.3.5 | 1.13.1 | 3.9.4 | 0.14.6 | fork `origin/master` (`0e5cd42`) | fork `origin/master` (`8fbef65`) |
| `venv-spmg` (pinned) | CPython 3.13.16 arm64 | 2.1.3 | 2.2.3 | 1.18.1 | 3.10.9 | 0.15.0 | `python313` (`caace75`) | fork `python313` (`7c3b1f2`) |
| `venv-spmg-latest` | CPython 3.13.16 arm64 | 2.5.4 | 3.1.0 | 1.18.1 | 3.11.2 | 0.15.0 | same | same |

Baseline order mirrors `docker/backend.dockerfile`: SPP fork, SPMG fork,
`seaborn==0.13.2`, then reinstall `numpy==1.26.4 pandas==1.3.5`. Deviation:
the prompt asked for `/opt/homebrew/bin/python3.9 -m venv`, but pandas 1.3.5 has
no macOS arm64 wheel and building it from source fails (no `clang`; Xcode tools
not installed). The baseline uses `uv`'s x86_64 CPython 3.9.25 under Rosetta,
which has a pandas 1.3.5 wheel. Both are CPython 3.9 with the same pinned
numpy/pandas as the image.

### Inputs and scenarios

mSigPortal never passes a raw VCF to SigProfilerMatrixGenerator.
`mSigPortal_Profiler_Extraction.py` first converts it with
`vcf_Multiple_Convert` to `<project>_mSigPortal_SNV.txt` and
`<project>_mSigPortal_INDEL.txt`, then calls `SigProfilerMatrixGeneratorFunc`
on that directory. The comparison uses that path:
`~/tmp/py313/spmg-convert.py` runs the real `vcf_Multiple_Convert` (loaded from
the server script via `ast`, not copied) on
`client/public/assets/exampleInput/demo_input_multi.vcf` and writes the two txt
files into `<scenario>/input/`. Result: 3 samples (SB749362, SC349747, SC420396),
4086 SNV lines, 740 INDEL lines.

One line of the demo file (line 118, `"##contig=<ID=GL000241.1,...>"`) is
wrapped in quotes. mSigPortal's own `Normalize_Delimited_Header` rejects the file
because of it (`Error 2727: The column of "#CHROM" can not be found`), on any
Python version. The converter drops lines starting with `"` before converting.
See Open risks.

The prompt's literal setup (raw VCF copied into `input/`) was also run on all
three envs. SPMG's VCF reader yields 0 samples / 0 mutations for this file on
every stack, so that run only shows `Statistics.txt` is identical
(`Total of 0 SNVs...`); no matrices are written.

Scenarios (`~/tmp/py313/spmg-run.zsh`), each in a fresh project directory,
project `test`, `exome=False, chrom_based=False, plot=False, tsb_stat=False`:

| Scenario | genome | seqInfo | bed_file | SPMG summary (all three envs) |
|----------|--------|---------|----------|-------------------------------|
| grch37 | GRCh37 | False | None | 3 samples, 0 errors, 4085 SNVs, 919 DINUCs, 739 INDELs |
| grch37_seq | GRCh37 | True | None | 3 samples, 0 errors, 4085 SNVs, 774 DINUCs, 739 INDELs |
| grch37_bed | GRCh37 | False | `demo_input_bed.bed` | 3 samples, 0 errors, 4085 SNVs, 919 DINUCs, 739 INDELs |
| grch38 | GRCh38 | False | None | 3 samples, 3433 errors, 1245 SNVs, 856 DINUCs, 209 INDELs |

The GRCh38 errors are expected (GRCh37 coordinates checked against the GRCh38
reference). The lower DINUC count with `seqInfo=True` is existing SPMG behavior,
identical on 3.9.

### Byte comparison

`~/tmp/py313/spmg-compare.zsh` runs `diff -rq` per directory and `cmp` on
`Statistics.txt` (log: `spmg-compare2.log`, after the last commit).

| Scenario | SBS | DBS | ID | vcf_files | Statistics.txt |
|----------|-----|-----|----|-----------|----------------|
| grch37 | identical (9) | identical (6) | identical (6) | n/a | identical |
| grch37_seq | identical (9) | identical (6) | identical (6) | identical (100) | identical |
| grch37_bed | identical (9) | identical (6) | identical (6) | identical (5) | identical |
| grch38 | identical (9) | identical (6) | identical (6) | n/a | identical |

Same result for `venv-spmg` and `venv-spmg-latest` against `venv-spmg39`. Zero
diffs. The 9 SBS files include SBS18/SBS288/SBS4608, the 6 DBS files include
DBS150/DBS2400, and the 6 ID files include ID332: these are the outputs built by
the `pd.concat` collapse code.

Integer-dtype risk (from the first pass) is closed:

| Check | pandas 2.2.3 | pandas 3.1.0 |
|-------|--------------|--------------|
| `grep -E '\.0(\s|$)'` over every matrix file | no matches | no matches |
| dtype of every DataFrame returned by `SigProfilerMatrixGeneratorFunc` (11 per run) | all integer | all integer |

### Warnings

| Check | `venv-spmg` | `venv-spmg-latest` |
|-------|-------------|--------------------|
| `PYTHONWARNINGS=error::DeprecationWarning python -c "import SigProfilerMatrixGenerator.scripts.SigProfilerMatrixGeneratorFunc"` | clean | clean |
| Strict import of `SigProfilerMatrixGeneratorFunc`, `MutationMatrixGenerator`, `CNVMatrixGenerator`, `save_chrom_tsb_separate`, `save_context_distribution`, `install` with `-W error::DeprecationWarning -W error::FutureWarning -W error::PendingDeprecationWarning` and `warnings.filterwarnings` stubbed out | clean (a probe `DeprecationWarning` raised, proving the filter was live) | clean |
| `python -W error::SyntaxWarning -m compileall -f` on the installed package | clean | n/a |
| Runtime, all warnings recorded (`spmg-warnrun.py`; grch37_seq, grch37_bed, grch38) before `caace75` | no Deprecation/Future warnings | 1 `Pandas4Warning` at `SigProfilerMatrixGenerator.py:2416` (`from_dict`) |
| Same, after `caace75` | none | none |

The `filterwarnings` stub matters: the SigProfilerPlotting fork calls
`warnings.filterwarnings("ignore")` at import, and SPMG imports it, so
`PYTHONWARNINGS` alone stops seeing anything after that import.

Remaining warnings, all `ResourceWarning` (unclosed files, existing upstream
code, not version-related, left as is): `SigProfilerMatrixGeneratorFunc.py:444`
(seqInfo and BED temp files) and `:532` (log file handle). 25 to 50 per run.

### Simulator and Clusters imports

PyPI `SigProfilerSimulator` 1.2.2 installed with `--no-deps` into `venv-spmg`,
plus `fastrand` 3.1.1. `~/tmp/py313/spmg-simcheck.py` resolves every
`from SigProfilerMatrixGenerator.scripts import X as Y` /
`import SigProfilerMatrixGenerator as Y` in Simulator and in the Clusters fork
(`src/sigprofilerclusters`, `main` `ffa8298`) and checks each `Y.attr` they use.

| Consumer | SPMG module (alias) | Names used | Missing |
|----------|---------------------|------------|---------|
| Simulator `SigProfilerSimulator.py`, `mutational_simulator.py` | `MutationMatrixGenerator` (`matRef`) | `reference_paths` | none |
| Simulator | `SigProfilerMatrixGeneratorFunc` (`matGen`) | `SigProfilerMatrixGeneratorFunc` | none |
| Simulator | `save_context_distribution` (`context_dist`) | `context_distribution`, `context_distribution_BED` | none |
| Simulator | package (`sig`) | `__version__` | none |
| Clusters fork `SigProfilerClusters.py` | `SigProfilerMatrixGenerator` (`matRef`, old name) | `reference_paths` | none |
| Clusters fork | `convert_input_to_simple_files` (`convertIn`) | `convertVCF`, `convertTxt`, `convertMAF`, `convertICGC` | none |
| Clusters fork | `SigProfilerMatrixGeneratorFunc` (`matGen`) | `SigProfilerMatrixGeneratorFunc` | none |

`import SigProfilerSimulator.SigProfilerSimulator` and
`SigProfilerSimulator.mutational_simulator` succeed, and
`MutationMatrixGenerator.matrix_generator is SigProfilerMatrixGenerator.matrix_generator`
is `True` (same function objects). `scripts/__init__.py` in the installed package
is empty, so this works without the Dockerfile line-34 patch. Lines 33-34 of
`docker/backend.dockerfile` can be deleted. A Simulator run was not done (not
part of the mSigPortal path).

## Open risks and follow-ups

- Install order in the image: Simulator 1.2.2 requires
  `SigProfilerMatrixGenerator>=1.3.2` and `sigProfilerPlotting>=1.4.0`. Without
  `--no-deps` for Clusters/Simulator, pip replaces both forks with PyPI
  releases (this is how the current image ends up fragile).
- `binom_test` shim still not exercised at runtime: `tsb_stat=True` and
  `gs=True` are not used by mSigPortal and were not run.
- Demo input: `client/public/assets/exampleInput/demo_input_multi.vcf` line 118
  is a quoted meta line that makes mSigPortal's VCF header check fail with
  Error 2727. This is independent of the upgrade and was not changed (project
  files are out of scope for this task); the comparison removed that line.
- Baseline env ran under Rosetta (x86_64) while the 3.13 envs are arm64. Outputs
  are text counts, so this does not affect the byte comparison, but it is not
  the Linux image. A Docker run of the backend image is still needed after the
  Dockerfile change.
- `install.py` genome download path (`wget`/`rsync`) was left as is; mSigPortal
  mounts genomes and does not use it.

## Branch and commits

Branch `python313` (from `origin/master`), not pushed. HEAD `caace75`
(`caace757a853f33e956073cb234e3a8a752cc01a`).

| Commit | Summary |
|--------|---------|
| `117c35c` | Replace removed pandas DataFrame.append with pd.concat |
| `d8ea031` | Add scipy binomtest shim for removed stats.binom_test |
| `127f877` | Use raw strings for regex escapes to silence Python 3.12+ SyntaxWarning |
| `7ee9dc8` | Add MutationMatrixGenerator alias module |
| `e6905f2` | setup.py: declare python_requires>=3.9 |
| `caace75` | Use DataFrame constructor instead of DataFrame.from_dict(DataFrame) (added during re-verification) |

Files changed: `SigProfilerMatrixGenerator/scripts/SigProfilerMatrixGenerator.py`,
`SigProfilerMatrixGenerator/scripts/MutationMatrixGenerator.py` (new),
`SigProfilerMatrixGenerator/scripts/save_chrom_tsb_separate.py`,
`SigProfilerMatrixGenerator/scripts/save_context_distribution.py`, `setup.py`.
No artifacts (install.log, version.py, egg-info, build, `__pycache__`, `z-bak/`,
genome data) were committed.
