# Python package regression harness

Supports the Python 3.13 upgrade of `docker/backend.dockerfile` and `docker/extraction-batch.dockerfile`.
The idea: run every Python code path the application can reach on the **current** images, keep the outputs
as a baseline, run the same cases on the **new** images, and diff.

Outputs land in `data/regression/<label>/<suite>/<case>/` (the `data/` folder is gitignored).

## 1. Which features use which Python package

Only two application features invoke Python. Everything else (Exploration, Association, Signature,
Etiology, Pattern, Publications, Refitting) runs through R or Node only.

### Backend image: Visualization → user data upload

`POST /submitVisualization/:id` → `server/services/api/visualization/profilerExtraction.js` →
`python3 services/python/mSigPortal_Profiler_Extraction.py`. Runs on the web task (local) or the ECS
Fargate worker (`server/worker.js`). Form options that select code paths: input format
(vcf, maf, csv, tsv, catalog_tsv, catalog_csv), genome (GRCh37, GRCh38, mm10), WGS/WES,
Cluster, SeqInfo, and for vcf/csv/tsv: Filter, Collapse, split-by-filter, Bed file.

| Package (xtmgah fork) | Reached when | What the user sees |
|---|---|---|
| SigProfilerMatrixGenerator 1.1.30 | input format vcf, maf, csv, tsv | SBS/DBS/ID matrices in all sizes, `matrix_files_list.txt`, Statistics, Bed filtering, SeqInfo files |
| sigProfilerPlotting 1.1.6 | input format catalog_tsv or catalog_csv; plot function picked by matrix row count (SBS, DBS, ID, CNV48, RS32) | profile plot images zipped per type under `output/` |
| SigProfilerClusters 1.0.11 + SigProfilerSimulator (PyPI, transitive) | vcf **and** Cluster on | `Cluster/Result/<id>_clustered_class_All.txt`, cluster class files and plots |
| seaborn, numpy, pandas, matplotlib | transitive only | n/a |

### Extraction image: Extraction page

`POST /submitExtraction/:id` → `server/services/workers.js` submits an AWS Batch job →
`extraction-service/worker.js` → `extraction-service/services/extraction.js` →
`python3 services/python/mSigPortal-SigProfilerExtractor.py`. GPU is on for every Batch run.
The script makes one call, `sigpro.sigProfilerExtractor(...)`; all other packages are used inside it.

| Package | Role | Form options that matter |
|---|---|---|
| SigProfilerExtractor 1.1.13 (xtmgah) | NMF (`nmf_cpu.py` / `nmf_gpu.py`), solution selection, De-Novo and Decomposed outputs | context type SBS96/DBS78/ID83, genome, exome, min/max signatures, replicates, iterations |
| torch 2.1.0 | GPU NMF | always in Batch |
| SigProfilerAssignment 0.0.14 | COSMIC decomposition, `De_Novo_map_to_COSMIC_*.csv` | collapse_to_SBS96 |
| sigProfilerPlotting 1.2.2 (PyPI) + PyPDF2 + reportlab | decomposition plots and PDFs | make_decomposition_plots |
| SigProfilerMatrixGenerator (upstream) | only for `input_type=vcf`, which the UI never sends; provides the genome path | none |

After the Python run, Node imports the results into SQLite and submits De-Novo and Decomposed
Exploration sessions through the API, so one end-to-end UI run per tier is still required.

## 2. Cases

Defined in `cases.json`. Profiler cases cover every input format, all three genomes, WES, Collapse,
split-by-filter, Bed, SeqInfo, both catalog formats, CNV48 and RS32 plotting, and two Cluster runs
(tagged `slow`). Extractor cases cover SBS96/DBS78/ID83 on CPU, SBS96 on GPU (tagged `gpu`), the
no-plot/no-collapse path, and exome normalization. Extractor runs use a fixed one-row seeds file written
by `run_cases.py` (seed 5177184, the first entry of
`data/examples/extraction/Example_SBS96_SigProfileExtractor/Seeds.txt`, whose 500-row `Seeds` format
the Extractor rejects) and small NMF settings so a
case finishes in minutes.

## 3. Running

Prerequisites: Docker, the genome data in `data/genomes/` (GRCh37, GRCh38, mm10), and the gitignored
Python scripts present under `server/services/python/` and `extraction-service/services/python/`.

### Backend (profiler) suite

Build the current image and capture the baseline before touching the Dockerfile (on the `python313`
branch, use `git show dev:docker/backend.dockerfile | docker build -t msigportal-backend:py39 -f - .`
and `--image msigportal-backend:py39` instead):

```bash
docker compose build backend
```

```bash
python3 utils/regression/run_cases.py profiler --label baseline --all
```

Drop `--all` to skip the two Cluster cases while iterating. After the Dockerfile upgrade, rebuild and run
the same command with `--label candidate`. To keep both images side by side, tag the builds and pass
`--image` instead of using compose (the profiler script does not need `server/.env`):

```bash
docker build -t msigportal-backend:py313 -f docker/backend.dockerfile .
python3 utils/regression/run_cases.py profiler --label candidate --all --image msigportal-backend:py313
```

`run_cases.py` runs `${PYTHON_BIN:-python3}` inside the container, so the py3.13 images (which set
`PYTHON_BIN=python3.13` and leave the system `python3` alone) and the old images both work unchanged.

### Extraction (extractor) suite

There is no compose service for this image, so tag the builds explicitly. On the `python313` branch the
files in `docker/` are already the 3.13 versions, so build the baseline from `dev`'s Dockerfile:

```bash
git show dev:docker/extraction-batch.dockerfile | docker build -t msigportal-extraction:py39 -f - .
```

The extractor needs a wide signature table (`Type` column plus one column per signature), the same file
`extraction.js` writes as `signature.tsv`. By default `run_cases.py` copies the COSMIC v3.3 GRCh37
reference tables from `packages/SigProfilerExtractor/Reference_Signatures/GRCh37/` into
`data/regression/inputs/signature_{SBS96,DBS78,ID83}.tsv`. To test against the exact table the
application builds from Postgres instead, export it in long format and pivot (this overrides the copy):

```bash
psql "$DATABASE_URL" -c "\copy (select \"signatureName\",\"mutationType\",contribution from signature where \"signatureSetName\"='COSMIC_v3.3_Signatures_GRCh37_SBS96' and study='Reference') to 'sig_SBS96_long.tsv' with (format csv, delimiter E'\t', header)"
```

```bash
python3 utils/regression/run_cases.py pivot-signatures sig_SBS96_long.tsv data/regression/inputs/signature_SBS96.tsv
```

Repeat for `COSMIC_v3.3_Signatures_GRCh37_DBS78` → `signature_DBS78.tsv` and
`COSMIC_v3.3_Signatures_GRCh37_ID83` → `signature_ID83.tsv`. Then:

```bash
python3 utils/regression/run_cases.py extractor --label baseline --image msigportal-extraction:py39
```

On a GPU host add `--include-tags gpu --gpus all`.

### Compare

```bash
python3 utils/regression/compare.py data/regression/baseline data/regression/candidate --json data/regression/report.json
```

Verdicts per file: `identical`, `tolerated` (numeric within 1e-6 relative, or ordering only), `visual`
(PNG/PDF/SVG differ but dimensions and page counts match), `log` (scanned for new tracebacks, counts of
warnings and deprecations reported), `varies` (SigProfilerClusters result tables, which differ run to run on
every version; agreement on sample/position/subclass is reported), `FAIL`, `missing`. Exit code is non-zero if any case fails.
Install Pillow on the host (`pip install pillow`) to get a perceptual-hash threshold on PNGs instead of
a plain "pixels differ".

## 4. What to expect to change on Python 3.13

- Matrices, statistics, cluster tables and extractor signature/activity tables should be `identical` or
  `tolerated`. Anything else is a defect in a fork, a pin, or a script.
- Plot images will shift pixels (matplotlib 3.4 → 3.9+, font registration). Review one SBS96 image per
  suite by eye and confirm the bundled font is still used.
- Log files will show new deprecation warnings from pandas 2 / numpy 2; a new `Traceback` is a failure.
- Each suite writes `_env.json` with the interpreter version and `pip freeze` so the two runs are
  traceable.
