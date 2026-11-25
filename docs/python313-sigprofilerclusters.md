# Python 3.13 upgrade: SigProfilerClusters (xtmgah fork)

## Summary

The xtmgah fork of SigProfilerClusters (1.0.11) now installs and runs on Python 3.13 with
numpy 2.1, pandas 2.2 and matplotlib 3.10, on top of the `python313` branches of the
SigProfilerPlotting and SigProfilerMatrixGenerator forks and PyPI SigProfilerSimulator 1.2.2.
Three commits on branch `python313` (two for 3.13, one bug fix for parallel file writes); the fork's
VAF code is untouched and the version is unchanged.

| Area | Why |
|---|---|
| `statistics.stdev()` gets Python ints in `hotspot.first_run` | on Python 3.12+ it raises `AttributeError: 'numpy.int64' object has no attribute 'bit_length'`, so every clusters run failed |
| `setup.py` | `matplotlib<3.5` cap blocked any Python > 3.9; SPP/SPMG floors were above the fork versions; seaborn was imported but not declared; no `python_requires` |

mSigPortal's real `Run_SigProfilerClusters` (bcftools split, Simulator, `hp.analysis`, result
assembly) runs end to end on 3.13. With the run made deterministic (see Verification),
`test_clustered_class_All.txt`, every subclass table, every matrix and 1927 text files are
byte-identical to the Python 3.9 / numpy 1.26.4 / pandas 1.3.5 / matplotlib 3.4.3 baseline.
Two internal pickles differ because the old stdev was wrong on Python 3.9 (details below).

## Fork origin and upstream base

| Item | Value |
|---|---|
| Fork | `https://github.com/xtmgah/SigProfilerClusters`, branch `main`, HEAD `ffa8298` "Create .gitignore" |
| Fork history | `12b478e` "Initial commit" (`.gitattributes`), `faa2387` "Updated VAF calculation" (Tongwu Zhang, 2022-08-09; adds the whole package), `ffa8298` |
| Upstream | `https://github.com/AlexandrovLab/SigProfilerClusters` (redirects to SigProfilerSuite) |
| Upstream base | [`804a633`](https://github.com/SigProfilerSuite/SigProfilerClusters/commit/804a633) "Allow single line MNVs with CCFs" (2022-05-03), the commit that set `VERSION = '1.0.11'`. No 1.0.11 tag exists (tags start at v1.0.0, then v1.1.2). |
| Base check | PyPI sdist `SigProfilerClusters==1.0.11`: package files and `setup.py` identical to `804a633` |
| Delta file | `~/tmp/py313/spc-fork_delta.diff` (`diff -ru` upstream `804a633` vs fork `origin/main`, excluding `__pycache__`, egg-info, `examples/`, dist, build) |

## Custom changes preserved

The fork delta is one source file plus `.gitignore` / `.gitattributes`:

| Customization | Location | Status |
|---|---|---|
| `is_number()` helper | `classifyFunctions.py:134` | unchanged |
| `pullVaf` standardVC branch: detects input type from column count / FORMAT (`Sanger` 8 columns, `non-Sanger-VAF`, `non-Sanger-AF`, `non-Sanger-AD`); VAF from INFO `VAF=`/`AF=`, FORMAT `AF`, or `AD` ratio rounded to 3 places; `-2.5` for multi-allelic, `-1.5` for missing | `classifyFunctions.py:238-349` | unchanged; exercised by the test (demo VCF FORMAT has `AF` -> `non-Sanger-AF`) |
| `pullVaf` writes `<project>_nonClustered_Vaf.txt` (non-clustered SNVs with VAF), which mSigPortal reads as `cN_Path` | `classifyFunctions.py:360-398` | unchanged; file produced and byte-identical before/after |
| Debug prints (`!!!!` banners, file paths) | `pullVaf` | unchanged |
| `.gitignore` reduced to `.DS_Store`; tracked `egg-info/` and `__pycache__/*.cpython-37.pyc` | repo root | unchanged (not touched by this branch) |

`git diff origin/main..python313` touches only `setup.py` and two lines of `hotspot.py`;
`classifyFunctions.py` is not changed, so no fork hunk is reverted.

Pre-existing fork quirk kept: the `non-Sanger-VAF` branch looks up `FORMAT_arr.index("AF")`, so a
FORMAT with `VAF` but no `AF` raises `ValueError`. Not version related; left as is.

## Compatibility changes made

| File:line | Change | Upstream reference |
|---|---|---|
| `SigProfilerClusters/hotspot.py:332` | `stdev(bin_std_dev)` -> `stdev([int(x) for x in bin_std_dev])` | identical change in [`57fe3ca`](https://github.com/SigProfilerSuite/SigProfilerClusters/commit/57fe3ca05de17867402bbc605bf54011ae1a112a) "Updated to support python 3.12" |
| `SigProfilerClusters/hotspot.py:339` | `stdev(bin_std_dev_cum)` -> `stdev([int(x) for x in bin_std_dev_cum])` | same commit |
| `setup.py` | `python_requires=">=3.9"` | [`f610ba6`](https://github.com/SigProfilerSuite/SigProfilerClusters/commit/f610ba6) "v1.2.0: Require Pandas and Numpy >= 2.0.0 and Python >= 3.9" |
| `setup.py` | `matplotlib>=3.3.0,<3.5` -> `matplotlib>=3.3.0` | upstream dropped the cap in [`2f19cfd`](https://github.com/SigProfilerSuite/SigProfilerClusters/commit/2f19cfd4321c84824b2a8f5f713b7890da96d4c5) (`>=3.5`); the floor is left as is |
| `setup.py` | `sigProfilerPlotting>=1.1.8` -> `>=1.1.6`; `SigProfilerMatrixGenerator>=1.2.4` -> `>=1.1.30` | fork versions; upstream raised these instead, which would exclude the forks |
| `setup.py` | add `seaborn>=0.13.0` | [`74977e9`](https://github.com/SigProfilerSuite/SigProfilerClusters/commit/74977e99658e358788fb1a64895cca2e6cd4f252) "Add seaborn dependency" |

## Bug fix: parallel appends to the clustered files (added 2026-11-04)

Not a 3.13 issue; found by the image regression run. `hotspot.first_run` runs in `mp.Pool` workers,
and every worker appends to the same `<project>_clustered.txt` and `<project>_nonClustered.txt`.
Their buffered writes interleaved: rows from different samples were spliced together and cut
mid-field, and a split on a line boundary left a blank line that crashed `classifyFunctions.pullVaf`
(`IndexError: list index out of range` at `sample = lines[1]`). Observed on `vcf_grch37_cluster`
(3 samples, `max_cpu=6`). Upstream HEAD has the same code and the same race.

| File | Change |
|---|---|
| `SigProfilerClusters/hotspot.py` | `hotSpotAnalysis` creates the pool with `initializer=_init_write_lock, initargs=(mp.Lock(),)`; `first_run` holds that lock (`_append_lock()`) for each of the two `with open(..., 'a')` blocks, so a sample's rows are written as one block. Outside a pool `_append_lock()` is a no-op. |

Effect: identical rows when no interleaving happened; corrupted rows from the old code no longer
occur. Sample order in the files still follows worker completion order (unchanged). Verified by
running both regression Cluster cases (`vcf_grch37_cluster`, `vcf_grch38_cluster`) on the
`msigportal-backend:py313` image with the patched file: both exit 0, result tables have no blank
or truncated rows. The 17-column `ClassIII` rows in `*_clustered_class_All.txt` (extra `vaf`
field) are pre-existing and also appear in the Python 3.9 baseline.

Behavior note on the stdev fix: on Python 3.9, `statistics.stdev()` of `numpy.int64` values
casts the variance back to the input type (`statistics._convert`), so it is truncated to an
integer. Example on the baseline env: `[0]*9 + [1]` gives `0.0` (true 0.316); `[0]*7 + [2,3,3]`
gives `1.0` (true 1.317). The fixed code computes the exact value, as upstream does. In this test
it changed the corrected IMD cutoff of 4 density-correction regions that each contain one
original mutation (`imds_corrected` 1 -> 10000, 10000, 2048, 32); no mutation changed class. Feeding
the dumped 3.9 inputs through the 3.9 code with only this change reproduces the 3.13 values
exactly (`~/tmp/py313/spc_dbg.py`), so this is the only source of the difference.

Scanned and clean: `np.mat`/`np.matrix`, `DataFrame.append`/`Series.append` (all `.append` hits
are Python lists), `.iteritems`, `np.int`/`np.float`/`np.bool`/`np.object`/`np.NaN`/
`np.product`/`np.in1d`, `delim_whitespace`, `error_bad_lines`, `line_terminator`, `squeeze=`,
`applymap`, `cm.get_cmap`, `register_cmap`, `distutils`, `pkg_resources`, `imp`,
`collections.Mapping`, `scipy.misc`, `binom_test`, `sklearn.externals`, PyPDF2, chained pandas
assignment. `font_manager._rebuild()` appears only in comments. `plottingFunctions.rainfall` and
`eventProbability` use pandas `read_csv`/`.loc`/`iterrows` in ways that work on pandas 2.2.

Expected pandas drift in `classifyFunctions.py` did not exist: that module does not use pandas.

## MutationMatrixGenerator import and the Dockerfile patch

The Clusters fork imports the old module name
(`from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGenerator as matRef`), which the
SPMG fork still has, so Clusters itself never needed the alias. SigProfilerSimulator 1.2.2 imports
`MutationMatrixGenerator`; on the 3.13 stack that resolves to the real alias module on the SPMG
`python313` branch (`scripts/__init__.py` is empty), and the Simulator ran. The
`echo ... >> scripts/__init__.py` patch in `docker/backend.dockerfile` (lines 33-34) is not needed
with the SPMG `python313` branch. The 3.9 baseline does need it: without it the Simulator import
fails with `ImportError: cannot import name 'MutationMatrixGenerator'`.

## Dependency changes

| Field | Before | After |
|---|---|---|
| `python_requires` | (absent) | `>=3.9` |
| `matplotlib` | `>=3.3.0,<3.5` | `>=3.3.0` |
| `sigProfilerPlotting` | `>=1.1.8` | `>=1.1.6` |
| `SigProfilerMatrixGenerator` | `>=1.2.4` | `>=1.1.30` |
| `SigProfilerSimulator` | `>=1.1.3` | unchanged |
| `seaborn` | (absent, but imported) | `>=0.13.0` |
| `statsmodels`, `scipy`, `pandas`, `numpy` | floors | unchanged |
| `SigProfilerExtractor` | `>=1.0.17` | unchanged (its only import is commented out; not installed in the test envs) |

### Install order for the Dockerfile

SigProfilerSimulator 1.2.2 requires `SigProfilerMatrixGenerator>=1.3.2` and
`sigProfilerPlotting>=1.4.0`, and Clusters requires `SigProfilerExtractor`. Any install with
dependency resolution of Simulator or Clusters replaces the forks with PyPI releases (and pulls
torch). This order keeps the forks:

```sh
PINS="numpy==2.1.3 pandas==2.2.3 matplotlib==3.10.9 seaborn==0.13.2"
pip install <sigprofilerplotting@python313> $PINS
pip install <sigprofilermatrixgenerator@python313> $PINS
pip install --no-deps SigProfilerSimulator==1.2.2
pip install "fastrand>=1.2"
pip install --no-deps <sigprofilerclusters@python313>
pip install statsmodels scipy $PINS
```

Then link `SigProfilerMatrixGenerator/references/chromosomes/tsb` to the genome data. Check with
`pip show sigProfilerPlotting SigProfilerMatrixGenerator` (1.1.6 / 1.1.30). `pip check` will
report the Simulator and Extractor requirements above; that is expected.

## Verification

### Environments

All under `~/tmp/py313/`, built by `spc-setup.sh` from `git archive` exports (clones untouched),
installed in the order above. `references/chromosomes/tsb` is a symlink to `data/genomes`.

| Env | Interpreter | numpy | pandas | matplotlib | scipy | statsmodels | seaborn | SPP / SPMG / SPC source | Simulator |
|---|---|---|---|---|---|---|---|---|---|
| `venv-spc39` (before) | CPython 3.9.25 x86_64 (Rosetta) | 1.26.4 | 1.3.5 | 3.4.3 | 1.13.1 | 0.14.1 | 0.13.2 | `origin/master` / `origin/master` / `origin/main`, plus the Dockerfile alias patch | 1.2.2 |
| `venv-spc` (after) | CPython 3.13.16 arm64 | 2.1.3 | 2.2.3 | 3.10.9 | 1.18.1 | 0.15.0 | 0.13.2 | `python313` (`7c3b1f2` / `caace75` / `6b8a7ff`) | 1.2.2 |

`pip show` in `venv-spc`: sigProfilerPlotting 1.1.6, SigProfilerMatrixGenerator 1.1.30,
SigProfilerClusters 1.0.11, all from the fork exports. fastrand 3.1.1 in both. Freezes:
`venv-spc.freeze.txt`, `venv-spc39.freeze.txt`.

Deviation: the prompt asked for `/opt/homebrew/bin/python3.9`; pandas 1.3.5 has no arm64 wheel
and cannot be built here (no clang), so the baseline uses uv's x86_64 CPython 3.9 (same as the
SPMG doc). SigProfilerSimulator 1.2.2 is the latest release and supports Python 3.9, so it is
what the current Dockerfile's unpinned install resolves to; both envs use it.

### Input

`demo_input_multi.vcf` cannot be read by bcftools as shipped: line 118 is a quoted
`"##contig..."` line, and the header `Description=` values lost their quotes (CSV export).
`spc-inputs/demo_clean.vcf` strips the quotes and re-quotes `Description` values. Sample
SB749362 was extracted with `bcftools view -c1 -s SB749362` (1662 records) into
`spc-inputs/SB749362.vcf`.

### Harness

`~/tmp/py313/spc_run.py` loads mSigPortal's real `Run_SigProfilerClusters` from
`server/services/python/mSigPortal_Profiler_Extraction.py` (via `ast`, not copied):

1. `sim`: same bcftools split and `sigSim.SigProfilerSimulator("test", <dir>, "GRCh37",
   contexts=['96'], chrom_based=True, simulations=100)` as mSigPortal, plus a fixed `seed_file`.
   Run once on 3.13 (`spc-sim`, 9 s): 1396 SNVs, 268 DINUCs, 249 INDELs, 100 simulation files.
2. `clusters`: copies of `spc-sim` per env; runs `Run_SigProfilerClusters` unchanged except the
   simulator call is skipped (precomputed) and temp-file cleanup is disabled. So `hp.analysis` is
   called with mSigPortal's exact arguments (`analysis="all", sortSims=True, subClassify=True,
   correction=True, calculateIMD=True, max_cpu=6, standardVC=True, TCGA=False, sanger=False`) and
   the result assembly that writes `Result/test_clustered_class_All.txt` runs.

`simulations` cannot be small: `analysis()` exits with "Please simulate a minimum of 100
simulations" below 100. Both envs analyze the same simulation files, so differences come from
Clusters and its stack only. Multiprocessing start method is `fork` (the Linux default in the image).

Clusters is not deterministic run to run, on either Python version: worker processes append to
the same per-sample IMD file in arbitrary order, and `hotspot.densityCorrection` calls
`random.sample(..., 10)` in a worker, where `random` is reseeded from urandom after fork. Two
3.13 runs with `max_cpu=6` differed in the class of a few mutations; two 3.9 runs did too, and so did two `max_cpu=1` runs without the re-seed. For the
comparison the harness sets `max_cpu=1` (`SPC_MAX_CPU=1`) and re-seeds `random` in each forked
child (`os.register_at_fork`). With that, two runs on each env are byte-identical.

### Results

Unmodified mSigPortal call (`max_cpu=6`), before the stdev fix: 3.13 fails in
`hotspot.first_run` (`'numpy.int64' object has no attribute 'bit_length'`). After the fix: 3.13
and 3.9 both exit 0 and produce `Result/test_clustered_class_All.txt` (1394 rows),
`rainfallPlots_clustered_test_corrected.pdf` and `test_intradistance_plots_96_corrected.pdf`.

Deterministic comparison (`spc-compare.zsh`, `spc_cmp_pickles.py`; `spc-39-s1` vs `spc-313-f1`):

| Output | Result |
|---|---|
| `Result/test_clustered_class_All.txt` | byte-identical (1394 rows, same 16 columns) |
| `output/vcf_files_corrected/test_clustered/subclasses/class{1,1a,1b,1c,2,3}/*.txt` | byte-identical |
| `output/vcf_files_corrected/test_nonClustered/SNV/test_nonClustered_Vaf.txt` (fork file) | byte-identical |
| `output/clustered/{DBS,MBS,omikli,kataegis,other}/*.vcf`, `output/nonClustered/SBS/*.vcf` | byte-identical |
| all subclass matrices (`SBS6`..`SBS6144`, DBS) and `Statistics.txt` | byte-identical |
| `output/simulations/` IMD files (original and simulated) | byte-identical |
| total text files compared | 1927 identical; differing: 13 SPMG/Clusters log files (versions, paths, timestamps), `SB749362.vcf` (bcftools header records the output path), 19 pickles |
| 19 pickles in `output/simulations/data/` | 17 equal by value (bytes differ only by numpy pickle format); `imds_corrected` and `interval_lines_corrected` differ for 4 regions, explained under the stdev fix |
| PDFs | same page count (1 each); not compared visually |

### Warnings

| Check | Result |
|---|---|
| `PYTHONWARNINGS=error::DeprecationWarning python -c "import SigProfilerClusters.SigProfilerClusters, .hotspot, .classifyFunctions, .plottingFunctions, .convertToVCF, SigProfilerSimulator.SigProfilerSimulator"` | clean (there is no `SigProfilerClusters.hp` module; `hp` is mSigPortal's alias for `SigProfilerClusters.SigProfilerClusters`) |
| `python -W error::SyntaxWarning -m compileall -f` on the installed package | clean |
| Full run with `warnings.simplefilter("always")` and `filterwarnings` stubbed (SPP and `hotspot.py` install ignore filters at import) | no Deprecation/Future/PendingDeprecation warnings |

Remaining warnings, pre-existing and left unchanged: `RuntimeWarning: invalid value encountered
in (scalar) divide` at `hotspot.py:50` (`z_test` with sigma 0), `:397`, `:400`, `:413`
(`hotspot.py:34` filters the numpy 1.x message "long_scalars"; numpy 2 words it differently, so
this filter no longer matches, but SPP's global ignore still hides it in normal runs);
`UserWarning: Glyph 65534 missing from font(s)` at `hotspot.py:1245` and
`plottingFunctions.py:1535`; `set_ticklabels()`/`FixedFormatter` at `plottingFunctions.py:249`, `:358`.

## Recommended pins (Dockerfile)

| Package | Pin |
|---|---|
| Python | 3.13 |
| numpy | 2.1.3 |
| pandas | 2.2.3 |
| matplotlib | 3.10.9 |
| seaborn | 0.13.2 |
| SigProfilerSimulator | 1.2.2 (`--no-deps`) |
| fastrand | 3.1.1 (`>=1.2`) |

## Open risks

- Not verified: newer stacks (numpy 2.5 / pandas 3) for Clusters; `probability=True`,
  `includedCCFs`, `TCGA`/`sanger` VAF modes, ID contexts, GRCh38/mm10. mSigPortal uses none of these.
- Run-to-run nondeterminism (worker completion order and unseeded `random` in workers) is upstream
  behavior on all versions. mSigPortal results can vary slightly between runs of the same input.
  The interleaved-write corruption is fixed on this branch (`4354a91`) but still present in the
  current production image, so baseline cluster files can contain corrupted rows.
- Behavior change from the stdev fix is real but small: density-corrected regions with very few
  mutations can get a different IMD cutoff than on the old image.
- The demo input `demo_input_multi.vcf` is not valid for bcftools, so the clustering path in
  mSigPortal fails on it at the `bcftools query -l` step on any Python. Project files were out
  of scope and not changed.
- The run was on macOS (arm64 after, x86_64 Rosetta before), not the Linux image. A run in the
  rebuilt backend image is still needed.
- `SigProfilerExtractor` stays in `install_requires` though unused; installs must use `--no-deps`.

## Branch and commits

Branch `python313` in `src/sigprofilerclusters`, off `origin/main` (`ffa8298`). Pushed to
`xtmgah/SigProfilerClusters` 2026-11-02; `4354a91` pushed 2026-11-04.
HEAD `4354a91`.

| Commit | Summary |
|---|---|
| `a9db281` | setup.py: declare python_requires, drop matplotlib<3.5 cap, relax fork floors (and declare seaborn) |
| `6b8a7ff` | Pass Python ints to statistics.stdev in first_run |
| `4354a91` | Serialize worker appends to the shared clustered/nonClustered files |

Files changed: `setup.py`, `SigProfilerClusters/hotspot.py`.

Scratch files: `~/tmp/py313/spc-setup.sh`, `spc_run.py`, `spc-runall.zsh`, `spc-compare.zsh`,
`spc_cmp_table.py`, `spc_cmp_pickles.py`, `spc_dbg.py`, `spc-inputs/`, run dirs `spc-sim`,
`spc-39-s1`, `spc-313-f1`.
