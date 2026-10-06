# Provenance

Where the code that produces the numbers comes from.

## The WAFC code

The estimator (`cv.wafc()`, `wafc_threshold()`), the data generating
process (`simulate_wafc()`) and the competitors (`wafc_competitor()`) are
the WAFC code, loaded as `code` in `config/study.yaml` says. While the
compendium lives beside it, the code is the folder `wafc/` of the
repository `michelcias/wafc-draft`, sourced through its loader
(`code: ../wafc/R/load.R`):

| | |
|---|---|
| repository | `michelcias/wafc-draft` |
| commit (last change to `wafc/R`) | `f14e6f21702d842dc98ed2f307cf97f36abc051d` |
| tree of `wafc/R` at that commit | `e262faac78b2af83011408fb0dae98a62977c67d` |

In the standalone release of the compendium the same functions come from
the version of `WaveBased` that contains them, pinned in `renv.lock`, and
`code` is `package:WaveBased`; this table then names that version and the
commit of `wafc/R` it was built from.

Every run writes what it actually loaded to
`outputs/<design>/sessionInfo.txt` (for a path in a git working tree: the
commit, the tree of the directory and whether it had uncommitted changes;
for a package: its version and commit), and every unit keeps the same line
in its file. The numbers reported in the article come from runs whose line
names the code above with no uncommitted changes.

## WaveBased

The wavelet bases are evaluated by `WaveBased` (`wbasis()`, `wtable()`),
pinned in `renv.lock`:

| | |
|---|---|
| repository | `michelcias/WaveBased` |
| commit | `e494b0eb64a605c9ed8d0a50f2552aea650321eb` |
| version | 2.6-0 |

## Data

The data of the applications are in `data/`, built by `data-raw/` from
the files listed in `data-raw/sources.yaml`, whose SHA-256 digests are
those of the files the data were built from; `data/README.md` gives the
sources and their licences.

## R and the other packages

R 4.6.1 and the versions recorded in `renv.lock` (`glmnet` 5.1, `grpreg`
3.6.0, `mgcv` 1.9-4, `Matrix` 1.7-6, `VCBART` 1.2.5, `yaml` 2.3.12 and
their dependencies). The session of each run is in
`outputs/<design>/sessionInfo.txt`.
