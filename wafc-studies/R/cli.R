## R/cli.R -- the command line of the scripts.
##
## Arguments are --name=value; lists are comma separated. Restrictions
## (--cells, --sizes, --methods, --reps) select units and never change what
## a unit draws, so a restricted run computes exactly the units a full run
## would.
##
##   --config=PATH    design file (default config/core.yaml)
##   --study=PATH     study file (default config/study.yaml)
##   --out=DIR        output directory (default outputs/<design>)
##   --workers=N      forked workers (default `workers` of the study file)
##   --cells=a,b      only these cells
##   --sizes=250,500  only these sample sizes
##   --methods=a,b    only these methods
##   --reps=R         replicates 1, ..., R; or --reps=FROM:TO
##   --refresh        refit cached units whose settings changed
##   --list           print the units and stop

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (a in args) {
    if (!startsWith(a, "--")) stop("unrecognized argument '", a, "'; use ",
                                   "--name=value.", call. = FALSE)
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1L]]
    out[[kv[1L]]] <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=")
                     else TRUE
  }
  known <- c("config", "study", "out", "workers", "cells", "sizes",
             "methods", "reps", "refresh", "list")
  bad <- setdiff(names(out), known)
  if (length(bad) > 0L) stop("unknown option(s): ",
                             paste0("--", bad, collapse = ", "),
                             call. = FALSE)
  split <- function(v) if (is.null(v)) NULL else
    strsplit(v, ",", fixed = TRUE)[[1L]]
  reps <- out[["reps"]]
  if (!is.null(reps)) {
    reps <- if (grepl(":", reps, fixed = TRUE)) {
      r <- as.integer(strsplit(reps, ":", fixed = TRUE)[[1L]])
      seq(r[1L], r[2L])
    } else seq_len(as.integer(reps))
  }
  list(config = if (is.null(out[["config"]])) file.path("config", "core.yaml")
                else out[["config"]],
       study = if (is.null(out[["study"]])) file.path("config", "study.yaml")
               else out[["study"]],
       out = out[["out"]],
       workers = if (is.null(out[["workers"]])) NULL else
         as.integer(out[["workers"]]),
       cells = split(out[["cells"]]),
       sizes = if (is.null(out[["sizes"]])) NULL else
         as.integer(split(out[["sizes"]])),
       methods = split(out[["methods"]]),
       reps = reps,
       refresh = isTRUE(out[["refresh"]]),
       list = isTRUE(out[["list"]]))
}

## The output directory of a run: the one given, or outputs/<design>.
run_dir <- function(cli, cfg) {
  if (!is.null(cli[["out"]])) cli[["out"]] else
    file.path("outputs", cfg[["design"]])
}

run_workers <- function(cli, cfg) {
  w <- if (!is.null(cli[["workers"]])) cli[["workers"]] else cfg[["workers"]]
  if (is.null(w)) w <- 1L
  max(1L, as.integer(w))
}
