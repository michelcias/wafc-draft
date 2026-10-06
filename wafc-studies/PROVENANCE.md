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
| commit (last change to `wafc/R`) | `9b11273bd56130a692af03e794084d7b73fc8df5` |
| tree of `wafc/R` at that commit | `0818de514a2412a832f082fbe98eed3cada3adae` |

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

## R and the other packages

R 4.6.1 and the versions recorded in `renv.lock` (`glmnet` 5.1, `grpreg`
3.6.0, `mgcv` 1.9-4, `Matrix` 1.7-6, `VCBART` 1.2.5, `yaml` 2.3.12 and
their dependencies). The session of each run is in
`outputs/<design>/sessionInfo.txt`.
