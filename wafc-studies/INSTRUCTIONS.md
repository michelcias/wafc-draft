# Working on the compendium

What someone who changes this compendium needs to know. `README.md` says
what it reproduces and how to run it; this file says how it is kept
reproducible.

## Purpose

The compendium reproduces the simulation study and the data applications
reported in the article
and its supplementary material, and nothing else. A configuration, a
method or a table that does not appear there does not belong here. The
technical reason behind a constant (a cost that bounds a grid, a failure a
guard prevents) is part of the reproduction and stays in the comment next
to it.

## Environment

- Start R from the root of the compendium: `.Rprofile` activates `renv`,
  and every script checks that it runs from there.
- `renv::restore()` installs what `renv.lock` pins. The dependencies are
  declared in `DESCRIPTION` (`renv` snapshots explicitly), not detected
  from the code, because the code that uses most of them is the WAFC code
  loaded from outside `R/`.
- A system `Rprofile.site` that loads `bspm` makes `renv` warn that a
  package was loaded before it activated; the warning does not affect what
  `renv` installs.
- To move to another `WaveBased`, install it by commit and snapshot:

  ```r
  renv::install("michelcias/WaveBased@<sha>")
  renv::snapshot()
  ```

  Only a commit that is on GitHub: a local commit cannot be downloaded,
  and `renv::restore()` fails on the next machine. A new `WaveBased` that
  changes any number invalidates every cached unit (see "The cache").

## The WAFC code

The estimator and the competitors are loaded, never copied into `R/`:
`code` in `config/study.yaml` is either the path of the loader of the WAFC
code, which `scripts/01_simulate.R` sources before fitting, or
`package:<name>` for a package that exports the same functions (the
standalone release takes them from the pinned `WaveBased`). Every run
records in `outputs/<design>/sessionInfo.txt`, and every unit in its file,
what was loaded: the commit and whether it had uncommitted changes, or the
version and commit of the package.
`PROVENANCE.md` records the version the reported numbers come from;
numbers produced from a working tree with uncommitted changes are not
reported.

## Rules of the code

- **One place for each setting.** Cells, sample sizes and replicates in
  the design files; seeds, methods, method options, tuning grids and
  labels in `config/study.yaml`. The code reads them and holds no study
  constant of its own.
- **Seeds only through `unit_seeds()`** (`R/seeds.R`). Nothing else calls
  `set.seed()`. `seeds$keys` and `seeds$sample_sizes` are append-only:
  their positions are part of every seed.
- **The data of a replicate depend on the cell, the sample size and the
  replicate only,** never on the method or on what else is run. The
  methods start from the random stream the folds leave behind
  (`restore_stream()`).
- **Every method is read the same way** (`R/metrics.R`): prediction on the
  test sample, functional coefficients there, components on the grid,
  blocks kept. A method that cannot answer one of these leaves the
  column `NA`; it is not given a substitute.
- **A failed fit is a row,** with its message in `error` and no numbers,
  never a missing row: a column that disappears from a table reads as a
  method that was never run.
- **A new method** is an entry of `study_methods` (`R/methods.R`), with its
  row labels, and a name in `methods` of `config/study.yaml`. Adding it
  touches no cached unit of the other methods.

## The cache

A unit is one method on one replicate,
`outputs/<design>/units/<cell>/n<n>/<method>/repNNN.rds`, written
atomically when it finishes. It carries the settings that determine it
(the cell, the tuning of its method at its sample size, the sizes of the
test sample and of the grid, the seeds). When the configuration changes a
setting of a cached unit, the run stops and names the first stale unit;
`--refresh` refits the stale units, and removing their files does the
same. The settings do not include the code: after a change of the WAFC
code or of `WaveBased` that moves any number, remove the units by hand
(the whole `outputs/<design>/units/`, or the methods affected).

## Running long

- A full design runs for many hours. Launch it detached and log it, for
  instance

  ```bash
  nohup Rscript scripts/01_simulate.R --config=config/core.yaml \
    --workers=8 > outputs/core.log 2>&1 &
  ```

  and read the log; each finished unit prints one line.
- The peak memory is that of the largest unit times the number of
  workers; the largest units are `bsgl` and `gam.reml` in the cells with
  the most blocks (`mixed`, `scale`) at the largest sample size.
- A run with a configuration that is not part of the study (to time a
  grid, or to test a change) is given its own design file and its own
  `--out` outside `outputs/`, so it never mixes with the units of the
  study.

## Producing the reported numbers

1. `renv::restore()`; the WAFC code at the commit of `PROVENANCE.md`,
   without uncommitted changes.
2. `Rscript run_all.R` (or the two scripts, design by design).
3. The tables in `outputs/<design>/tables/` are the ones the article
   reports; `sessionInfo.txt` beside them records the session.
