## ============================================================
## wafc/scripts/14-junction.R -- the junction across machines (E4.4)
## ============================================================
## Before the production runs on another machine, a small set of units is
## fitted there and here with the same command, and the two caches must
## agree in every field but the time. The units come from the pilot's
## design file (its own master seed, so no unit is a unit of the
## production): three cells, n = 250 and 1000, replicate 1, every method.
##
## Fit, from the root of wafc-studies/ (on each machine, with its own DIR):
##   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript scripts/01_simulate.R \
##     --config=../wafc/scripts/11-pilot.yaml --out=DIR --workers=8 \
##     --cells=smooth,mixed,inhomogeneous.xu --sizes=250,1000 --reps=1
##
## Compare, from the root of wafc-draft:
##   Rscript wafc/scripts/14-junction.R REFERENCE_DIR OTHER_DIR
##
## Prints one line per unit and "JUNCTION OK" when every unit present in
## both agrees in its rows and side tables, the columns of time excepted.
## ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: Rscript wafc/scripts/14-junction.R REF OTHER")
ref <- args[[1L]]
oth <- args[[2L]]

units <- function(dir) {
  f <- list.files(file.path(dir, "units"), pattern = "\\.rds$",
                  recursive = TRUE)
  sort(f)
}
drop_time <- function(d) {
  if (!is.data.frame(d)) return(d)
  d[, setdiff(names(d), c("time", "elapsed", "design")), drop = FALSE]
}
same <- function(a, b) {
  ra <- drop_time(a[["rows"]])
  rb <- drop_time(b[["rows"]])
  ok <- isTRUE(all.equal(ra, rb, tolerance = 0, check.attributes = FALSE))
  sa <- a[["side"]]
  sb <- b[["side"]]
  if (!identical(sort(names(sa)), sort(names(sb)))) return(FALSE)
  for (nm in names(sa)) {
    ok <- ok && isTRUE(all.equal(drop_time(sa[[nm]]), drop_time(sb[[nm]]),
                                 tolerance = 0, check.attributes = FALSE))
  }
  ok
}

fa <- units(ref)
fb <- units(oth)
common <- intersect(fa, fb)
cat(sprintf("units: %d in the reference, %d in the other, %d in both\n",
            length(fa), length(fb), length(common)))
if (length(common) == 0L) stop("no unit in common")
res <- vapply(common, function(f) {
  same(readRDS(file.path(ref, "units", f)), readRDS(file.path(oth, "units", f)))
}, logical(1))
for (f in common) cat(sprintf("%-60s %s\n", f, if (res[[f]]) "identical" else "DIFFERS"))
cat(sprintf("\n%d of %d identical (time excepted)\n", sum(res), length(res)))
cat(if (all(res) && length(fa) == length(fb)) "JUNCTION OK\n" else "JUNCTION FAILED\n")
