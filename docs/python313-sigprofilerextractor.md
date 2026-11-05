# Python 3.13 upgrade: SigProfilerExtractor (xtmgah fork)

## Summary

The xtmgah SigProfilerExtractor 1.1.13 fork now installs and runs on Python 3.13 with numpy 2,
pandas 2, matplotlib 3.11, torch 2.15 and SigProfilerAssignment 1.1.5. The work is on branch
`python313` (8 commits, version still 1.1.13). All mSigPortal customizations are kept.

The branch code was checked in two ways:

- **Same code on the old stack.** The `python313` code, run in the old Python 3.9 stack, gives
  byte-identical output to the unpatched fork for every non-PDF file.
- **Old vs new stack.** Remaining differences come from library versions only:
  - SigProfilerAssignment 1.x changes some labels and output file names.
  - Float32 NMF convergence lands one check earlier in some replicates, which shifts de novo
    signatures by up to 3.8e-6.
  - COSMIC signatures, COSMIC activities and the decomposition mapping are numerically identical.

The GPU path was exercised on CPU tensors only. There is no GPU on this host.

## Fork origin and upstream base

| Item | Value |
|---|---|
| Fork | `https://github.com/xtmgah/SigProfilerExtractor`, default branch `main`, HEAD `2425fdd` |
| Clone | `src/sigprofilerextractor` (gitignored) |
| Upstream | `https://github.com/AlexandrovLab/SigProfilerExtractor` (redirects to SigProfilerSuite), tag `v1.1.13`; HEAD `v1.5.0` used as the fix reference |
| Delta | `~/tmp/py313/spe-fork_delta.diff` (`diff -ru` upstream v1.1.13 vs fork HEAD, excluding egg-info, `__pycache__`, `.DS_Store`) |

## Custom changes preserved

| Hunk (fork commit) | Location | Status |
|---|---|---|
| `Signature_Database` parameter replaces `cosmic_version` (03d41c3, "Modify By Jian 2022-10-27") | `sigpro.py` signature, `execution_parameters`, JOB_METADATA writer | unchanged |
| `decomp.spa_analyze(..., signature_database=Signature_Database, ...)` | `sigpro.py:886` | unchanged; SPA 1.1.5 still accepts `signature_database` |
| `"192"` added to the SBS context `elif` (16c40d9) | `sigpro.py:517` | unchanged (the condition is still always true, as before) |
| `version.py`, `.gitattributes` | repo root / package | unchanged |
| Data files under `data/` and `src/FormatFiles/` | | differ from upstream only by CRLF; unchanged |

`sigpro.py` still differs from upstream 1.1.13 by exactly the fork hunks, plus the two lines this
branch adds. No customization conflicted with an upstream fix.

## Compatibility changes made

Upstream links use `https://github.com/AlexandrovLab/SigProfilerExtractor/commit/<sha>`.

| File:line | Change | Why | Upstream reference |
|---|---|---|---|
| new `nndsvd.py`; `nmf_cpu.py:14`, `nmf_gpu.py:7`, `_initialise_wh` | Use a vendored `Nndsvd` instead of `nimfa`; pass ndarrays instead of `np.mat(...)`; remove `nimfa` from `setup.py` | nimfa 1.4.0 calls `np.mat` internally and fails under numpy 2 even with `np.asmatrix` input | [89b5518](https://github.com/AlexandrovLab/SigProfilerExtractor/commit/89b551886abce876beb553a9c24d7559d89c9105), [68eac2c](https://github.com/AlexandrovLab/SigProfilerExtractor/commit/68eac2c1cca331556f4f195e7d1aa09227a1a41b) |
| `nndsvd.py` body | Line-by-line port of nimfa's dense NNDSVD, not upstream's rewrite. Reproduces three nimfa behaviours: `float*np.matrix` routes through `np.dot` in float64; numpy-1 value-based casting of the float64 scalar; NNDSVDar always draws `np.random.uniform(1, 1)` because `len()` of a matrix selection is 1 | The UI exposes nndsvd/nndsvda/nndsvdar/nndsvd_min. Upstream's rewrite changes their results (`eps` term, generator RNG) | fork behaviour kept over upstream |
| `subroutines.py:25`, `:1413-1418`; `PlotDecomposition_*.py:14` and `crop_margins` | PyPDF2 1.x/2.x API → `pypdf` (`PdfReader/PdfWriter`, `.pages`, `add_page`, `mediabox.lower_left` etc.); unused `PdfFileMerger` import dropped; `setup.py` `PyPDF2` → `pypdf>=3.0.0` | PyPDF2 3.x raises on the old names; SPA ≥0.2 already depends on pypdf | upstream deleted these modules instead ([a185f6f](https://github.com/AlexandrovLab/SigProfilerExtractor/commit/a185f6f)). Rename-only edit, checked by reverse-mapping to HEAD |
| `subroutines.py:743` | `clusters.append(...)` → collect frames, then `pd.concat(..., ignore_index=True)` | `DataFrame.append` removed in pandas 2.0 | [68eac2c](https://github.com/AlexandrovLab/SigProfilerExtractor/commit/68eac2c1cca331556f4f195e7d1aa09227a1a41b) |
| `subroutines.py:199,206,219` | `np.array(totalMutations)[idx][:, np.newaxis]` | `totalMutations` is a Series; `Series[:, None]` raises in pandas 2.0 | [68eac2c](https://github.com/AlexandrovLab/SigProfilerExtractor/commit/68eac2c1cca331556f4f195e7d1aa09227a1a41b) |
| `sigpro.py:867` | Add `devopts['make_decomposition_plots']` | SPA ≥0.0.33 raises KeyError without it | [9b6d241](https://github.com/AlexandrovLab/SigProfilerExtractor/commit/9b6d241) |
| `sigpro.py:460` | `np.array(seeds["Seed"]).item()` | `int()` of a 1-element 1-D array raises in numpy 2.5 (deprecated since 1.25); seeds-file runs failed | minimal form of [6ef047c](https://github.com/AlexandrovLab/SigProfilerExtractor/commit/6ef047c) |
| `subroutines.py:1293`; `setup.py` | `set_xticklabels(list(...), rotation=30)`; matplotlib cap `<=3.4.3` removed | `fontdict` was passed positionally (TypeError on current matplotlib) | [f500883](https://github.com/AlexandrovLab/SigProfilerExtractor/commit/f5008831b7ef3b274804517c3fac2fc3e4c5647a) |
| `subroutines.py:1033-1088` (6 calls) | `to_csv(path, "\t", ...)` → `sep="\t"` | Positional arguments are deprecated in pandas 2 and keyword-only in pandas 3. Output verified unchanged | upstream HEAD uses keywords |
| `setup.py` | Add `python_requires=">=3.9"`; drop the dead win32 torch 1.5.1 block and the now-unused `sys`/`subprocess` imports | The block would raise ValueError and install a torch with no wheels | — |

Checked and left as is:

- **torch legacy API.** `torch.FloatTensor`/`DoubleTensor` as `.type()` targets, `.cuda()`,
  `torch.cuda.set_device` and `generator.random(..., dtype=np.float64)` all still work under
  torch 2.15 with every warning raised as an error (see test 3). Upstream HEAD still uses them.
- **`np.matrix`.** `subroutines.nnmf_cpu/nnmf_gpu` still use it. It exists in numpy 2 and upstream
  HEAD keeps it.
- **Other breakers.** No hits for `iteritems`, `np.int/float/bool/object/NaN/product/in1d`, the
  removed `read_csv` arguments, `applymap`, `cm.get_cmap`, `register_cmap`, `distutils`,
  `pkg_resources`, `imp`, `collections.Mapping`, `scipy.misc`, `binom_test`,
  `sklearn.externals`. The two `font_manager._rebuild` hits are comments.

## Dependency changes (setup.py)

| Requirement | Before | After |
|---|---|---|
| `python_requires` | none | `>=3.9` |
| matplotlib | `>=3.4.2,<=3.4.3` | `>=3.4.2` |
| nimfa | `>=1.1.0` | removed |
| PyPDF2 | `>=1.26.0` | replaced by `pypdf>=3.0.0` |
| others (scipy, torch, numpy, pandas, SPMG, sigProfilerPlotting, SPA, pillow, statsmodels, scikit-learn, psutil, reportlab) | | unchanged floors |

Version selection, tested empirically:

| Package | Result |
|---|---|
| SigProfilerAssignment 0.0.x (0.0.14–0.0.33) | not usable: every 0.0.x pins `pandas<2.0.0`, uses PyPDF2 2.x names, and 0.0.33 adds PyMuPDF |
| SigProfilerAssignment 1.1.5 | **works**: `spa_analyze(signature_database=..., devopts=...)` API compatible; uses pypdf |
| SigProfilerAssignment 0.2.6 | same API and requirements as 1.1.5 (pypdf, numpy/pandas ≥2); not run |
| sigProfilerPlotting 1.4.3 | works. The current image already gets 1.4.3, not 1.2.2: installing SPMG from git upgrades it |
| SigProfilerMatrixGenerator 1.4.0 | works (only imported by Extractor here) |
| pypdf 6.21.4, reportlab 5.1.0 | work on 3.13. The `reportlab<5.1` pin was a Python 3.9 problem only |

## Recommended Dockerfile pins and order (extraction-batch image)

```dockerfile
# 1. torch first, from the CUDA index matching the base image (cu126 for CUDA 12.6, etc.)
RUN python3.13 -m pip install "torch>=2.5,<3" --index-url https://download.pytorch.org/whl/cu126
# 2. the fork, editable (see Open risks: a non-editable install picks up the stale committed build/lib)
RUN python3.13 -m pip install -e 'git+https://github.com/xtmgah/SigProfilerExtractor@python313#egg=SigProfilerExtractor' \
    "SigProfilerAssignment>=1.1.5,<1.2" "sigProfilerPlotting>=1.4.3,<1.5" "SigProfilerMatrixGenerator>=1.4,<1.5" \
    "pypdf>=6,<7" "reportlab>=4,<6" "numpy>=2.1,<3" "pandas>=2.2.3,<3" "matplotlib>=3.9,<4"
```

- Drop `PyPDF2==2.11.2`, `SigProfilerAssignment==0.0.14`, `sigProfilerPlotting==1.2.2`, the
  forced numpy/pandas reinstall, and `reportlab<5.1`.
- The separate SPMG git install is not needed; Extractor only imports it.
- Update `LD_LIBRARY_PATH` to `.../python3.13/site-packages/nvidia/...` for the new interpreter.

This exact command, minus the CUDA index, was run in a fresh Python 3.13 venv and the SBS96 run
matched the after run (`venv-spe-fresh`, below).

## Verification

Environments (all in `~/tmp/py313/`):

| Venv | Python | Key versions |
|---|---|---|
| `venv-spe-before` (Dockerfile order; freeze in `spe-before-freeze.txt`) | 3.9.22 | torch 2.1.0, numpy 1.26.4, pandas 1.3.5 (built from source), matplotlib 3.4.3, SPA 0.0.14, sigProfilerPlotting 1.4.3, SPMG 1.4.0 (git e67c5c4), PyPDF2 2.11.2, nimfa 1.4.0, reportlab 5.0.1, scipy 1.13.1, scikit-learn 1.6.1 |
| `venv-spe-after` (freeze in `spe-after-freeze.txt`) | 3.13.16 | torch 2.15.0 (CPU), numpy 2.5.4, pandas 2.3.3, matplotlib 3.11.2, SPA 1.1.5, sigProfilerPlotting 1.4.3, SPMG 1.4.0, pypdf 6.21.4, reportlab 5.1.0, scipy 1.18.1, scikit-learn 1.9.1 |
| `venv-spe-fresh` | 3.13.16 | Recommended pins, fork installed with `pip install -e git+file://...@python313` (resolved to the same versions) |

Setup notes for the before venv:

- pandas 1.3.5 has no macOS arm64 cp39 wheel. It was built with
  `CC=/Library/Developer/CommandLineTools/usr/bin/clang pip install --no-build-isolation`
  because the Xcode clang shim is broken on this host.
- The fork was installed editable (`-e`), as in the Dockerfile.

The example `data/examples/extraction/*/Seeds.txt` files use the column header `Seeds` and have 500
rows, but 1.1.13 requires one `Seed` row. The runs used the first seed with the header renamed
(`spe-seeds-{SBS96,DBS78}.txt`: 5177184 and 3688686).

| # | Test | Result |
|---|---|---|
| 1 | `spe_nmf_test.py`: `nmf_cpu.NMF` on a fixed 96×12 matrix, all 5 inits × single/double, plus `subroutines.nnmf_cpu`; before saved to `spe_nmf_before.npz` | 52 of 53 arrays (all W, H, initial W/H, convergence) **bit-identical**. `nnmf_cpu` similarity summary differs by ≤3.3e-6 absolute (≤1e-7 relative): NEP 50 keeps `round(np.float32*100, 3)` in float32 instead of float64. Same 3-decimal values |
| 1b | `spe_nndsvd_vs_nimfa.py`: vendored `nndsvd.py` vs nimfa in the numpy-1 venv, 96 cases (4 shapes × float32/64 × 4 ranks × 3 flags, incl. rank > min dim) | all bit-identical; out-of-range rank raises IndexError in both |
| 2a | Branch code run on the **old** stack, SBS96 + DBS78 | every non-PDF output byte-identical to the unpatched baseline |
| 2b | Old stack run twice | reproducible (byte-identical) |
| 2c | Before vs after, SBS96 (collapse_to_SBS96=True) and DBS78 (False), `spe_compare.py` | no tracebacks; PDF page counts equal (14/14 each); see the next table |
| 3 | `spe_gpu_harness.py`: `nmf_gpu.NMF` with `torch.cuda.set_device`/`Tensor.cuda` patched in the harness, all inits × precisions, plus a batch of 2; `warnings.simplefilter("error")` | GPU class == CPU class exactly (max diff 0.0); batched shapes correct; no warnings |
| 4 | `PYTHONWARNINGS=error::DeprecationWarning` import of `sigpro`, `nmf_gpu`, `PlotDecomposition`; also `-W error` on `sigpro` | clean |

Before vs after extraction details (2c):

| Output | SBS96 | DBS78 |
|---|---|---|
| `Suggested_Solution/*De-Novo_Solution/Signatures/*_De-Novo_Signatures.txt` | max abs diff 3.8e-6 (exceeds the 1e-6 target) | identical |
| De novo activities (`*_De-Novo_Activities_refit.txt`), All_Solutions NMF activities | identical | identical |
| `COSMIC_*_Signatures.txt`, `COSMIC_*_Activities.txt`, `COSMIC_*_Samples_Stats.txt`, decomposition log | identical | identical |
| `De_Novo_map_to_COSMIC_*.csv` | same signatures, percentages and errors; de novo label `Signature 96-A` → `Signature 96-SBS96A` | same, `78-A` → `78-DBS78A` |
| All_Solutions signatures, Wall, SEM, sample stats, `All_solutions_stat.csv` | ≤8e-6 on signatures; convergence iterations 400→300 and 300→200 in two replicates; one `Maximum Sample L1%` 40.57→40.56 | identical |
| Signature assignment logs | different format (SPA 1.x rounding of initial composition, fewer add/remove lines); written in nondeterministic order across runs | same |
| Probability files | `Decomposed_Mutation_Probabilities.txt` → `De_Novo_MutationType_Probabilities[_refit].txt`; new `JOB_METADATA_SPA.txt` | same |

Why the de novo values move:

- Test 2a shows the branch code reproduces the old outputs exactly on the old stack.
- Test 1 shows NMF is bit-identical when it runs to `max_iterations`.
- The remaining difference is the float32 KL-loss convergence check under torch 2.15 vs 2.1. It
  stops two replicates one `nmf_test_conv` (100 iterations) earlier, and the averaged signatures
  shift in the 6th decimal.
- The old stack can't run on Python 3.13, so this can't be removed without pinning torch 2.1.

mSigPortal impact checked:

- The client derives the de novo name with `denovoSigString.slice(-1)` and
  `endsWith(split('-')[1])` (`client/src/components/pages/extraction/apiSlice.js:113,146`). Both
  work with the new label.
- Nothing in `server/`, `extraction-service/` or `client/src` reads
  `Decomposed_Mutation_Probabilities.txt` or the assignment logs. The private scripts repo was not
  checked.

Remaining warnings:

- None at import with `-W error`.
- At runtime, `subroutines.py` calls `warnings.filterwarnings("ignore")` at import (unchanged), so
  library warnings are muted inside the package.

## Open risks

- **GPU not run.** The GPU path was not run on a GPU. Still to verify on a real CUDA host:
  - `torch.cuda.set_device`, `.cuda(gpu_id)` and `.type(torch.FloatTensor)` moving tensors back
    to CPU in `nmf_gpu.py:39-120`.
  - `subroutines.nnmf_gpu` worker GPU assignment (`current_process()._identity % device_count()`).
  - The `LD_LIBRARY_PATH` for cuDNN/cupti under python3.13.
  - Float32 GPU results compared with the current image.
- **Stale `build/lib` in git.** The fork commits a `build/lib/` copy of the package whose
  `sigpro.py` predates the mSigPortal changes. A non-editable `pip install` copies that stale file
  (it lacks `Signature_Database`). Keep `pip install -e` in the Dockerfile. `build/` and egg-info
  were not touched, per the repo rules.
- **New SPA output format.** SigProfilerAssignment 1.x changes labels, probability file names and
  log format (see above). Check the private scripts repo and any downstream parsing.
- **De novo drift.** De novo signatures differ from the current image by ~4e-6, which is above the
  1e-6 target. The cause is the torch upgrade, not the patches.
- **Seed example files.** The `Seeds.txt` files under `data/examples/extraction/` are not valid
  1.1.13 seed files (`Seeds` header, 500 rows). Uploading one raises "Please set valid seeds"
  before and after this change.
- **Untested inputs.** ID83 and the non-default normalizations were not run end to end. ID83
  shares all code paths with SBS96/DBS78 except plotting by sigProfilerPlotting.

## Branch and commits

Branch `python313` in `src/sigprofilerextractor`, based on `origin/main` 2425fdd, not pushed.
HEAD `e97d4db84eb1fca4b3d2b8ba35eafc90918de0fb`.

| Commit | Subject |
|---|---|
| da68da9 | Replace nimfa NNDSVD with a vendored NumPy 2 port |
| 7cfc81d | Move from PyPDF2 1.x/2.x API to pypdf |
| e8b63b1 | pandas 2: replace DataFrame.append and Series multi-dim indexing |
| 2b535fd | Pass make_decomposition_plots to SigProfilerAssignment via devopts |
| 83259af | Read the seeds-file seed as a Python int |
| bffba4c | matplotlib >=3.5: fix set_xticklabels call and lift the <=3.4.3 cap |
| be5a3be | setup.py: add python_requires and drop the dead Windows torch 1.5.1 hack |
| e97d4db | Pass to_csv separator by keyword |

Files changed: `setup.py`, `SigProfilerExtractor/{nndsvd (new), nmf_cpu, nmf_gpu, sigpro, subroutines}.py`,
`SigProfilerExtractor/PlotDecomposition_{SBS96,SBS288,SBS1536,DBS78,ID83,CNV48}.py`.

Scratch artifacts in `~/tmp/py313/`:

- Scripts: `spe_nmf_test.py`, `spe_nndsvd_vs_nimfa.py`, `spe_run.py`, `spe_compare.py`,
  `spe_gpu_harness.py`.
- Outputs: `spe-out-before{,2}/`, `spe-out-branch-oldstack/`, `spe-out-after{,2}/`,
  `spe-out-fresh/`, `spe-compare-{SBS96,DBS78}.txt`.
