# Backend R regression suite

Run the R code shipped in a built backend image, independently of the web server,
database, and S3. The harness is mounted at run time, so **the image does not need
to be rebuilt to pick up these tests**. It never installs packages or edits the
application source. The host launcher needs only Python 3's standard library and
Docker.

From the repository root, after the backend finishes building:

```sh
python3 utils/regression/run_r.py --label r-candidate
```

This uses `docker compose run --rm --no-deps backend` and the existing Compose
configuration. To select a tagged image, avoid Compose's environment-file
requirement, and disable container networking:

```sh
python3 utils/regression/run_r.py \
  --label r-candidate --image msigportal-backend:latest
```

Each invocation needs a fresh label. Output is saved under
`data/regression/<label>/r/`, which is outside the source tree tracked by Git.
The container is disposable. This does not execute tests inside the running web
container, share its R process, or start its normal entrypoint.

To use packages installed manually in an existing container, select it explicitly:

```sh
python3 utils/regression/run_r.py --label r-live --container msigportal-backend-1
```

This copies the harness into a unique `/tmp` directory, executes fresh R
processes in that container, and copies reports back to the host. Application
source and installed packages are read without modification. The scratch
directory remains in the container for inspection. The container's existing
network configuration applies; the backend suite still substitutes local S3
fixtures. A new container from the image does not include packages installed
manually in an existing container.

## Scope and assertions

The backend suite discovers top-level function definitions in the five active
backend R scripts. Currently there are **113 definitions**, with **147 cases**.
`server/services/R/refittingWrapper.R` is excluded: refitting executes through
`refitting-service/refitting.R` in the separate service.

| File | Functions | Exercised behavior |
| --- | ---: | --- |
| `Sigvisualfunc.R` | 52 | Numerical reconstruction, probabilities, profile formatting and rendering, clustering, prevalence, burden, motif analysis, genomic contexts, association and regression |
| `explorationWrapper.R` | 21 | Landscape JSON, single-signature and too-few-sample branches, decomposition, cosine similarity, public/user exposure flows, legacy plot entry points, dispatcher serialization |
| `visualizeWrapper.R` | 28 | Public/upload profile comparisons, cosine heatmaps, PCA, motif analysis, summaries, downloads, tree/hierarchy data, dispatcher serialization |
| `associationWrapper.R` | 7 | Variable selection, collapse options, univariable and multivariable analysis, dispatcher serialization |
| `utils.R` | 5 | Known-error conditions, empty-reference guards, RData loading and reference selection |

Every discovered function must have a registered case. A new function without a
case, or a stale target after a rename, fails discovery. Workers trace function
entry, including helpers redefined by nested `source()` calls. A case cannot pass
without invoking its declared target. This is **function coverage**, not branch
or line coverage; it does not establish that every input or scientific result is
correct.

Cases use seeded synthetic profiles and exposures with known reconstruction,
orthogonal/identical-vector cosine checks, normalized probability checks, expected
sample counts, validation conditions, and FDR calculations. Helpers render real
SVG/PDF outputs so lazy ggplot errors are evaluated. Wrapper responses containing
`error`, `uncaughtError`, or a named plot error fail even when the R function
catches the exception and returns normally. Expected validation errors have
explicit assertions.

Each case gets its own R process, working directory, deterministic random seed,
and timeout (default 600 seconds). One failure does not stop subsequent cases.
The suite itself uses base R; production package requirements are exercised by
the application scripts. A full run can take several minutes; motif enumeration,
association plots, and genome annotation are the more expensive cases.

### Explicit boundaries

- Public-data tests replace S3 loading with local fixture tables. The qualified
  `aws.s3::object_exists()` probe returns false in disposable workers. They test
  downstream processing, not AWS authentication, real S3 objects, or study-specific
  reference overrides. No real S3 bucket is used.
- SATS and refitting are outside the backend suite. The separate service
  preflight is described below.
- Genomic-context cases use installed human genome packages. The rainfall case
  supplies a small local chromosome-layout RData fixture.
- Installer scripts and offline `utils/*.R` data-maintenance jobs are outside
  this backend runtime suite.

## Baseline comparison

A single run checks contracts and successful execution. To detect numerical and
table changes across dependency upgrades, run the same harness on a trusted old
image, then the new image:

```sh
python3 utils/regression/run_r.py \
  --label r-baseline --image msigportal-backend:before-upgrade

python3 utils/regression/run_r.py \
  --label r-candidate --image msigportal-backend:after-upgrade \
  --baseline data/regression/r-baseline/r
```

Each candidate case requires a passing baseline case and compares its return
snapshot, plot-layer data where available, and generated text tables with
`all.equal(tolerance = 1e-6)`. Row order, column names, factor levels, and output
structure matter. Temporary paths, response timestamps, file sizes, and random
Plotly trace/widget IDs are excluded. Package versions, source hashes, and R
session details are retained for auditing.

SVG/PDF/PNG byte differences are listed in `visual-changes.txt` for visual review;
they do not by themselves fail a case. Font and graphics-device changes often
alter those bytes. Composite plots have a render check but no pixel-diff oracle.
Warnings are captured for review rather than automatically failing all cases.

Do not treat an upgraded-image snapshot as proof of unchanged behavior without
comparing it with a trusted baseline. Existing application failures remain
failures; the harness does not silently mark them expected or skip them.

## Reports and focused runs

- `summary.csv`: case status, target, scope, duration, warning count, failure text.
- `coverage.csv`: every function, invocation status, and the cases that called it.
  A function called only by failed cases is `INVOKED_IN_FAILED_CASE`, not
  successfully exercised.
- `packages.csv`, `sources.csv`, `session-info.txt`: runtime and source provenance.
- `warnings.txt`: warnings grouped by case.
- `<case>/run.log`: source, package, plot, and error diagnostics.
- `<case>/snapshot.rds`: normalized return value and generated table data.
- `<case>/artifacts/`: rendered plots and generated output files.

Exit status is nonzero if any selected case fails, times out, crashes, or differs
from its baseline. A focused run reports unselected functions as `NOT_EXERCISED`;
it is not a full-suite pass.

```sh
# List registered cases.
python3 utils/regression/run_r.py --label r-list --image msigportal-backend:latest --list

# Reproduce the original Landscape path, including its edge cases.
python3 utils/regression/run_r.py --label r-landscape --image msigportal-backend:latest \
  --case 'explorationWrapper[.]msigportal[.]msLandscape'

# Association cases, including grouped statistics and regression.
python3 utils/regression/run_r.py --label r-association --image msigportal-backend:latest \
  --case 'associationWrapper|mSigPortal_associaiton'
```

For a shell already inside a backend container, mount this directory at
`/opt/msigportal-r-regression` and a writable output directory, then run:

```sh
Rscript --vanilla /opt/msigportal-r-regression/run.R \
  --server /deploy/server --output /tmp/r-regression
```

Use `--server` if the image stores its server at another path. Use `--timeout`
to increase the per-case timeout. Do not modify case fixtures while a run is in
progress; rerun after fixture changes with a fresh label.

## Extending the suite

Add cases in `r/cases.R`. A case names its source file and function and contains
an explicit invocation plus assertions. Add shared inputs in `r/helpers.R`.
`r/dbs78.txt` contains the mutation-type column from the repository's example
`client/public/assets/exampleInput/extraction_sample_DBS78.all`; counts are
generated deterministically, not copied from real samples.

Use the `scope` field when a test exercises validation or a mocked boundary.
Keep deployment-load checks separate from stubbed wrapper contracts, and keep
real fitting/reference-data integration tests in the service that owns them.

`server/services/R/Sigvisualfunc.R` is excluded by the repository's existing
`.gitignore`. The two comparison-helper fixes are also saved in
[`patches/profile-comparison-fixes.patch`](patches/profile-comparison-fixes.patch)
so they can be reviewed and applied to the separately managed R script.

## Separate refitting-service preflight

```sh
python3 utils/regression/run_r.py --suite refitting \
  --label refitting-candidate --image msigportal-refitting:latest
```

Or use `--container <running-refitting-container>` to check a running service.
The default application directory for this suite is `/refitting-service`.
Reports go to `data/regression/<label>/r-refitting/`.

This suite sources the real `refitting.R`, loads its actual SATS dependencies,
and checks missing-file and invalid-genome validation for both
`run_sbs_refitting()` and `run_dbs_refitting()` (five checks total). It does not
use model stubs or the unused backend wrapper's `/app/...` source path. These
are service load/validation checks; full scientific fitting regression still
requires representative MAF, genomic, clinical, and shared reference fixtures.
Numerical `--baseline` comparison is supported by the backend suite only.

The harness's own assertion, snapshot, compressed-table, and command-construction
checks run without Docker or the plotting packages:

```sh
Rscript --vanilla utils/regression/r/self_test.R
python3 -m unittest discover -s utils/regression -p test_run_r.py
```
