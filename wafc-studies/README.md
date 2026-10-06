# wafc-studies

Reproducibility compendium for the simulation study of the article on
regression with additive functional coefficients estimated by a block
lasso on wavelet bases (WAFC):

$$Y = \sum_{\ell=1}^{p} \beta_\ell(U) X_\ell + \varepsilon, \qquad
\beta_\ell(u) = c_\ell + \sum_{m=1}^{q} g_{\ell m}(u_m).$$

One command, run from the root of this folder in the locked environment,
fits every replicate of every design and writes every table the article
and its supplementary material report.

## 1. Requirements

- R 4.6.1 (the version recorded in `renv.lock`) and a C/C++ toolchain:
  `WaveBased` compiles from source, and so may some dependencies when no
  binary is available.
- [`renv`](https://rstudio.github.io/renv/). `renv.lock` pins every
  package, `WaveBased` at a commit of `michelcias/WaveBased` on GitHub.
- The WAFC code (the estimator and the competitors), which this
  compendium loads and does not reimplement: see section 6 and
  `PROVENANCE.md`.

## 2. Setup

From the root of the compendium (the folder holding this file), in a fresh
R session, which activates the project library through `.Rprofile`:

```r
renv::restore()
```

This installs the pinned versions, `WaveBased` from GitHub at the commit
recorded in `renv.lock`. Do not run a bare `renv::install()`, which would
look for `WaveBased` on CRAN.

## 3. Reproduce everything

```bash
Rscript run_all.R
```

For each design listed in `designs` of `config/study.yaml`, this fits every
unit and writes the tables under `outputs/<design>/tables/`. A unit is one
method on one replicate; each is cached under `outputs/<design>/units/` as
soon as it finishes, so an interrupted run is resumed by running the same
command again. The number of parallel workers is `workers` of
`config/study.yaml` (forked processes, Unix only); results do not depend
on it.

The full study takes several hundred processor-hours, most of it in
`gam.gcv` and `bsgl` at the largest sample size. The two steps can be run
separately, design by design and restricted to part of it, with the
scripts:

```bash
Rscript scripts/01_simulate.R --config=config/core.yaml --workers=8
Rscript scripts/02_aggregate.R --config=config/core.yaml
```

| Option | Meaning |
|---|---|
| `--config=PATH` | design file (default `config/core.yaml`) |
| `--study=PATH` | study file (default `config/study.yaml`) |
| `--out=DIR` | output directory (default `outputs/<design>`) |
| `--workers=N` | forked workers |
| `--cells=a,b`, `--sizes=250,500`, `--methods=a,b` | only these cells, sample sizes, methods |
| `--reps=R` or `--reps=FROM:TO` | only these replicates |
| `--refresh` | refit cached units whose settings no longer match the configuration |
| `--list` | print the units and stop |

A restriction selects units and never changes what a unit draws: a
restricted run computes exactly the units the full run computes.

## 4. The study

### Designs

| File | Cells | Sample sizes | Replicates |
|---|---|---|---|
| `config/core.yaml` | `smooth`, `inhomogeneous`, `mixed`, `null`, `uneven` | 250, 500, 1000, 2000 | 100 |
| `config/arms.yaml` | `inhomogeneous` and `mixed` with, one at a time, signal to noise ratio 1, correlated modulators, linear covariates dependent on the modulators | 1000 | 100 |
| `config/scale.yaml` | `(p, q) = (6, 4)`, 6 of 24 blocks active | 1000 | 100 |

Each design file describes its cells; `config/study.yaml` holds what the
designs share: the seeds, the test sample, the methods and their
options, the tuning grids by sample size, and the labels of the tables.

### Methods

| Row | What it is |
|---|---|
| `wafc.cv1se` | the WAFC estimator: `cv.wafc()` with the block lasso, `(J, lambda)` by 10-fold cross-validation over `J` in 2, ..., 8, then the threshold on the norms of the blocks at one standard error |
| `wafc.cv` | the same fit with the threshold of the smallest cross-validated error |
| `lasso`, `lasso.cv1se` | the coordinatewise lasso on the same design, without and with the threshold |
| `gam.reml`, `gam.gcv` | `mgcv` with one smooth `s(u_m, by = x_l)` per block, the basis dimension `k` in 5, 10, 20, 40, 80 chosen by REML or by GCV, with `mgcv::k.check()` on the final fit |
| `bsgl` | cubic B-splines with a group lasso by block, `2^J` basis functions |
| `aspline` | a spline with knots chosen adaptively per block |
| `klopp` | the block lasso of Klopp and Pensky (2015) on the wavelet design |
| `oracle` | the lasso on the wavelet design restricted to the truly active blocks (a reference) |
| `linear` | constant coefficients, least squares on `X` (a reference) |
| `vcbart` | VCBART (Deshpande et al.), Bayesian trees for varying coefficients |

Every method is tuned in its own terms, and every method of a replicate
sees the same training sample, the same folds and the same test sample.

### What is measured

On a test sample of 2000 observations from the same law: the prediction
error against the true regression function (`rmse_f`), against the
responses (`rmse_y`), and the error of the functional coefficients
(`mse_beta`). On a grid of 256 points over the observed range of each
modulator: the integrated squared error of the centred components, in
total and split into active and zero blocks. The blocks kept, from which
`P(S-hat = S)`, the false positives and the false negatives. The
resolution level or basis dimension chosen, and whether it is the top of
its grid, read on the levels the design is built at: the WAFC code caps
the level of a modulator at `floor(log2(d))`, `d` its distinct points on
the circle, so at n = 250 the level 8 is built at 7 and the grid of `J` is
2, ..., 7 in effect. The elapsed time of the fit, with its search over the
grids.
The definitions are in `R/metrics.R`.

## 5. Outputs

`outputs/<design>/` (not versioned):

- `units/<cell>/n<n>/<method>/repNNN.rds`: one file per unit, with its
  rows, its side tables, the settings that determine it and the
  provenance of the code that fitted it;
- `sessionInfo.txt`: the R session, the files of the configuration, the
  commit of the WAFC code and the `WaveBased` build of the last run;
- `tables/`: `replicates.csv` (one row per fit), `summary.csv` (by cell,
  sample size and method), `paired.csv` (ratios to `wafc.cv1se` within
  replicate), `slopes.csv` (log error against log n), `components.csv`
  (error of every block), and the side tables `gam_k.csv`, `kcheck.csv`,
  `convergence.csv`, `threshold.csv`, `curves.csv`; `missing.csv` lists
  the units not yet fitted.

## 6. Layout

```
wafc-studies/
├── README.md          this file
├── INSTRUCTIONS.md    how to work on the compendium
├── CLAUDE.md          rules for an assistant working here
├── PROVENANCE.md      where the WAFC code and WaveBased come from
├── DESCRIPTION        the dependencies, for renv
├── renv.lock          the locked environment
├── run_all.R          the whole study
├── config/
│   ├── study.yaml     seeds, methods, tuning grids by sample size, labels
│   ├── core.yaml      the main design
│   ├── arms.yaml      one factor at a time at n = 1000
│   └── scale.yaml     the larger model
├── R/
│   ├── config.R       reading and checking the configuration
│   ├── seeds.R        one seed per cell, sample size, replicate and stream
│   ├── data.R         training sample, test sample, grid and folds
│   ├── methods.R      the methods and how each one is read
│   ├── metrics.R      what is measured on every fit
│   ├── run.R          units, cache and the parallel driver
│   ├── aggregate.R    from the units to the tables
│   └── cli.R          the command line
└── scripts/
    ├── 00_setup.R     checks the environment and sources R/
    ├── 01_simulate.R  fits the units of a design
    └── 02_aggregate.R writes the tables of a design
```

The estimator and its competitors are loaded from the WAFC code named by
`code` in `config/study.yaml` (the path of its loader, or a package that
exports it); `PROVENANCE.md` records the version.

## 7. Protocol

- **Randomness.** A single master seed (`seeds$master` of
  `config/study.yaml`). The seeds of the training sample, the test sample
  and the folds of a replicate are an injective function of the master
  seed, the cell's seed key, the sample size and the replicate
  (`R/seeds.R`), so no two units share a seed by accident and a unit
  draws the same data whatever else is run. The methods of a replicate
  start from the random stream the folds leave behind.
- **Pairing.** Every method of a replicate is fitted on the same data and
  folds and read on the same test sample and grid. An arm takes the seed
  key of the core cell it varies.
- **Tuning.** Declared per sample size in `config/study.yaml`; a grid can
  change at one sample size, for one method, without touching the rest.
- **The cache guards the protocol.** A unit records the settings that
  determine it; a run whose configuration no longer matches a cached unit
  stops instead of mixing the two, unless told to refit (`--refresh`).

## 8. License

The code of the compendium is released under the GNU General Public
License, version 3 or later (see `DESCRIPTION`).
