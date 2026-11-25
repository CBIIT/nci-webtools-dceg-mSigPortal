# Python 3.13 upgrade: SigProfilerPlotting (xtmgah fork)

## Summary

The xtmgah fork of SigProfilerPlotting (1.1.6) now installs and runs on Python 3.13 with numpy 2,
pandas 2.2 / 3.x and matplotlib 3.9 to 3.11. Every fork customization is kept. Three small commits
on branch `python313`; the package version is unchanged.

| Area | Why |
|---|---|
| `Series.append` -> `pd.concat` in `plotCNV` | pandas 2.0 removed `Series.append` |
| Bundled fonts registered at import time | the `setup.py` install hook used `font_manager._rebuild()` (gone in matplotlib 3.6) and never ran for wheel installs |
| `setup.py` metadata | declare `python_requires`, real dependencies, stop deleting the tracked `dist/` on every build |

All 12 plot cases render on Python 3.13 and produce the same PDF page counts, SVG counts and
pixel-identical first pages as the unmodified fork on Python 3.9 / matplotlib 3.4.3.

## Fork origin and upstream base

| Item | Value |
|---|---|
| Fork | `https://github.com/xtmgah/SigProfilerPlotting`, branch `master`, HEAD `8fbef65` "updated to 2.0" (Tongwu Zhang, 2021-11-05) |
| Upstream | `https://github.com/SigProfilerSuite/SigProfilerPlotting` |
| Upstream base | commit `10c478b` "Fixed tick_param conflict with matplotlib>=3.1.0" (2020-04-27), the only upstream commit with `VERSION = '1.1.6'` (no tag exists). Its `sigProfilerPlotting.py` is byte-identical to fork commit `b94023e`. |
| Fork-only commits | `aa84804` "updated sigProfilerPlotting" (2020-06-18), `8fbef65` "updated to 2.0" (2021-11-05) |
| Delta file | `~/tmp/py313/spp-fork_delta.diff` (`diff -ru` upstream `10c478b` vs fork `origin/master`, excluding egg-info, dist, build, z-bak, `__pycache__`, fonts, .DS_Store) |

Upstream history was rewritten after the fork: the fork's 2019-2020 commits by Erik Bergstrom
(e.g. `74b6304`, `e92670d`) exist upstream under different SHAs (`b285e19`, `16d0bfb`), with
identical content.

Correction to the earlier session's notes: SBS288, SBS288_Normalized, 384_extended, the SBS1536
label fix, zero-mutation sample handling and the upper-right legend placement are upstream
1.1.6 features (upstream commits `f6e2581`, `b285e19`, `16d0bfb`, `fb29da5`, `dad61a0`,
`f54f483`), not fork customizations. They are still present and verified below.

## Custom changes preserved

The real fork delta (from `spp-fork_delta.diff`):

| Customization | Location (python313 branch) | Status |
|---|---|---|
| `plotSV()` for RS32 rearrangement matrices | `sigProfilerPlotting.py:62` | unchanged, verified (RS32 case) |
| `plotCNV()` for CNV48 matrices, with `read_from_file` / `write_to_file` (PNG `BytesIO` buffers) | `sigProfilerPlotting.py:240` | one line changed (`Series.append`), verified (CNV48 case) |
| SBS4608 plot type (1536 heatmap plus TSB panel) | `sigProfilerPlotting.py:2312` | unchanged, verified (SBS4608 case) |
| Per-sample SVG copy of every plot in `<output>/svg/<TYPE>_plots_<project>_<sample>.svg` | every plot branch | unchanged, verified (SVG count per case) |
| ID83 label spelling ("Microhomology", ">1bp Insertion at Repeats") | `plotID` | unchanged |
| `import pandas`, global `warnings.filterwarnings("ignore")` | module top | unchanged |
| README typo fix, `.gitignore` | repo root | unchanged |
| `z-bak/sigProfilerPlotting.py` backup copy | `sigProfilerPlotting/z-bak/` | unchanged, excluded from the delta |

`git diff origin/master..python313` touches only `setup.py` and adds lines to
`sigProfilerPlotting.py`; no fork hunk was reverted.

Pre-existing fork quirk kept as is: the SBS288 branch calls `pp.savefig(plot1)` twice per sample
(`sigProfilerPlotting.py:3377` and `:3379`), so `SBS_288_plots_*.pdf` has two pages per sample.
Same before and after.

## Compatibility changes made

| File:line | Change | Upstream reference |
|---|---|---|
| `sigProfilerPlotting/sigProfilerPlotting.py:278` | `lab_homdel.append(lab_LOH.append(lab_het))` -> `pd.concat([lab_homdel, lab_LOH, lab_het])` | same line in upstream `caf1ac3` "Support pandas and numpy >= 2.0.0": https://github.com/SigProfilerSuite/SigProfilerPlotting/commit/caf1ac3672414a8eebf1c8fb3acbfd1a5242d1ad |
| `sigProfilerPlotting/sigProfilerPlotting.py:32-60` | new `_register_bundled_fonts()`, called at import: `fontManager.addfont()` for each bundled TTF, looking in the installed `sigProfilerPlotting/fonts/` first, then the repo-root `fonts/` (editable installs) | upstream `load_custom_fonts()` uses `fontManager.addfont`, introduced in `82ca6e6` "Modify font installation": https://github.com/SigProfilerSuite/SigProfilerPlotting/commit/82ca6e63ff64e6c61cd3b8a0211f142bdbdcf639 |
| `setup.py` | removed `move_ttf` install command (`cmdclass`), which copied fonts into matplotlib and called `font_manager._rebuild()` | upstream removed it in the same commit `82ca6e6` |
| `setup.py` | `packages` adds `sigProfilerPlotting.fonts`, `package_dir={'sigProfilerPlotting.fonts': 'fonts'}`, `package_data` `*.ttf` | fonts stay at the repo root, where `docker/backend.dockerfile` copies them from |
| `setup.py` | removed the `shutil.rmtree("dist")` that ran on every import of `setup.py` | it deleted the tracked `dist/sigProfilerPlotting-1.1.6.tar.gz` on each build |
| `setup.py` | `python_requires=">=3.9"` | upstream `ec680a2` "v1.4.0: Pandas and Numpy >=": https://github.com/SigProfilerSuite/SigProfilerPlotting/commit/ec680a2b25f7e989a0375e8d1ea5aa2d1f0eaaa3 |

Scanned and clean: `np.mat`/`np.matrix`, `DataFrame.append`, `.iteritems`, `np.int`/`np.float`/
`np.bool`/`np.object`/`np.NaN`/`np.product`/`np.in1d`, `delim_whitespace`, `error_bad_lines`,
`line_terminator`, `squeeze=`, `applymap`, `cm.get_cmap`, `register_cmap`, `distutils`,
`pkg_resources`, `imp`, `collections.Mapping`, `scipy.misc`, PyPDF2. The remaining `.append(`
hits are on Python lists. `_rebuild` survives only in comments.

## Dependency changes

| Field | Before | After |
|---|---|---|
| `python_requires` | (absent) | `>=3.9` |
| `install_requires` | `matplotlib` | `matplotlib>=3.2`, `numpy`, `pandas` (no upper caps) |

The `>=3.2` floor is for `fontManager.addfont`. numpy and pandas were imported but undeclared.

## Verification

### Environments

| Env | Where | Versions | Package source |
|---|---|---|---|
| before | `python:3.9-slim` container (Linux arm64) | Python 3.9.25, numpy 1.24.4, pandas 1.5.3, matplotlib 3.4.3 | fork `origin/master`, fonts copied to `/usr/share/fonts` as in `docker/backend.dockerfile` |
| after, minimum | macOS uv venv | Python 3.13.16, numpy 2.1.3, pandas 2.2.3, matplotlib 3.9.4, pyparsing 3.3.3, seaborn 0.13.2 | `python313` (git archive, non-editable) |
| after, recommended | macOS uv venv | Python 3.13.16, numpy 2.1.3, pandas 2.2.3, matplotlib 3.10.9, pyparsing 3.3.3 | same |
| after, latest | macOS uv venv | Python 3.13.16, numpy 2.5.4, pandas 3.1.0, matplotlib 3.11.2, seaborn 0.13.2 | same |
| font check | `python:3.13-slim` container | Python 3.13.18, numpy 2.1.3, pandas 2.2.3, matplotlib 3.9.4 | same |

The before env could not be built on macOS (no arm64 wheel for matplotlib 3.4.3 and the Xcode
license blocks compiling it), so it ran in Docker.

### Plot cases

Harness: `~/tmp/py313/spp-test/run_plots.py`. Every case is called the way
`mSigPortal_Profiler_Extraction.py` calls it (`percentage=True`; CNV/SV with `aggregate=False`).
A case passes when exactly one PDF is written, its page count matches the sample count, one SVG
per sample is written, and nothing is printed to stdout (the package swallows exceptions and
prints "There may be an issue with the formatting..." instead).

Inputs: `extraction_sample_SBS96.all` (21 samples), `Example_DBS78_SigProfileExtractor/DBS78/Samples.txt`
(251), `Example_ID83_SigProfileExtractor/ID83/Samples.txt` (260), `breast_cancer_samples_example.CNV48.all`
(4), `breast_cancer_samples_example.RS32.all` (4). SBS288 (T/U/N), SBS384 (T/U/B/N), SBS1536 and
SBS4608 matrices were synthesized from the first 3 SBS96 samples with a seeded random split of
each count. The zero case is SBS96 with an added all-zero sample.

| Case | before (3.9) | 3.13 min | 3.13 mpl 3.10 | 3.13 latest |
|---|---|---|---|---|
| SBS96 | 21 pages, 21 svg | same | same | same |
| DBS78 | 251 / 251 | same | same | same |
| ID83 | 260 / 260 | same | same | same |
| CNV48 | 4 / 4 | same | same | same |
| RS32 | 4 / 4 | same | same | same |
| SBS288 | 6 / 3 (double page, see above) | same | same | same |
| SBS288_Normalized | 3 / 3 | same | same | same |
| SBS384 | 3 / 3 | same | same | same |
| SBS384_extended | 3 / 3 | same | same | same |
| SBS1536 | 3 / 3 | same | same | same |
| SBS4608 (fork) | 3 / 3 | same | same | same |
| SBS96 with zero sample | 3 / 3 | same | same | same |

First pages of SBS96, CNV48, SBS4608 and ID83, rasterized at 1400 px with `sips`, are
pixel-identical between before and after (matplotlib 3.10.9) (PIL `ImageChops.difference` bbox
is empty).

### Fonts on a Linux host

In a clean `python:3.13-slim` container (0 system TTFs), with the branch installed non-editable:

| Query | Before `import sigProfilerPlotting` | After import |
|---|---|---|
| Arial | `ValueError` (not found) | `site-packages/sigProfilerPlotting/fonts/Arial.ttf` |
| Arial bold | n/a | `.../fonts/Arial Bold.ttf` |
| Times New Roman, regular and bold | n/a | `.../fonts/Times New Roman.ttf`, `Times New Roman Bold.ttf` |
| Courier New, regular and bold | n/a | `.../fonts/Courier New.ttf`, `Courier New Bold.ttf` |

Editable install (`uv pip install -e`, as the backend Dockerfile does with `pip install -e git+...`):
`_register_bundled_fonts()` returns the repo-root `fonts/` and all six TTFs are in
`fontManager.ttflist`.

For comparison, the unmodified fork on `python:3.9-slim` without the Dockerfile's
`cp fonts/* /usr/share/fonts` step also fails `findfont('Arial')`: its `move_ttf` hook does not
run for wheel builds. The new code does not need that Dockerfile step.

### Deprecation warnings

| Run | Result |
|---|---|
| `PYTHONWARNINGS=error::DeprecationWarning`, full harness, latest stack | all 12 cases pass |
| same, import only, matplotlib 3.9.4 + pyparsing 3.3.3 | fails inside matplotlib: `PyparsingDeprecationWarning: 'oneOf' deprecated` (`matplotlib/_fontconfig_pattern.py`) |
| same, import only, matplotlib 3.9.4 + pyparsing 3.2.5 | clean |
| same, import only, matplotlib 3.10.9 + pyparsing 3.3.3 | clean |

The package calls `warnings.filterwarnings("ignore")` at import (a fork customization), which
hides everything after import, so `PYTHONWARNINGS` alone only tests import. A second mode in the
harness (`SPP_WARN=record`) resets the filters after import and records every warning by origin:

| Stack | Deprecation/Future/PendingDeprecation warnings |
|---|---|
| latest (numpy 2.5.4, pandas 3.1.0, mpl 3.11.2) | none |
| mpl 3.10.9 | none |
| mpl 3.9.4 | `PyparsingDeprecationWarning` `parseString`, `resetCache` from matplotlib 3.9.4 (third-party) |

Remaining non-deprecation warnings, all pre-existing and cosmetic, left unchanged:
`UserWarning: set_ticklabels() should only be used with a fixed number of ticks` at
`sigProfilerPlotting.py:185`, `:409`, `:3010`, `:3373`, `:3695`, and
`Attempting to set identical low and high ylims` at `:660` for the all-zero sample. Nothing
originating in this package needed a fix.

## Recommended pins (Dockerfile)

| Package | Pin |
|---|---|
| Python | 3.13 |
| numpy | 2.1.3 |
| pandas | 2.2.3 |
| matplotlib | 3.10.9 (or 3.9.4 with `pyparsing<3.3`) |
| seaborn | 0.13.2 |

numpy 2.5.4 / pandas 3.1.0 / matplotlib 3.11.2 also pass every case; the other SigProfiler forks
set the final upper limits.

## Open risks

- GPU and the extraction image are not involved: `docker/extraction-batch.dockerfile` installs
  `sigProfilerPlotting==1.2.2` from PyPI, not this fork.
- The after runs were on macOS, where system Arial exists. Linux font resolution was checked
  with `findfont` in a container; full plot rendering on Linux was only done for the before env.
- SBS288/384/1536/4608 inputs were synthesized; real SigProfilerMatrixGenerator output of those
  types was not used.
- `os.system("mkdir %s")` for the `svg/` folder (fork code) is unchanged; output paths come from
  the server, not users.

## Branch and commits

Branch `python313` in `src/sigprofilerplotting`, off `origin/master` (`8fbef65`). Not pushed.
HEAD `7c3b1f207c69b604db24268321e9162a932a5ceb`.

| Commit | Summary |
|---|---|
| `ad4358f` | Replace removed pandas Series.append in plotCNV |
| `11a0b29` | Register bundled fonts at import time instead of a setup.py install hook |
| `7c3b1f2` | setup.py: declare python_requires and the real runtime dependencies (also drops the `dist/` deletion) |

Files changed: `setup.py`, `sigProfilerPlotting/sigProfilerPlotting.py`.
