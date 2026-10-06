## R/config.R -- reading and checking the configuration.
##
## A run reads two files: the study file (config/study.yaml), with what
## every design shares, and one design file (config/core.yaml,
## config/arms.yaml, config/scale.yaml, or one kept elsewhere), with the
## cells, the sample sizes and the replicates. The design file is merged
## over the study file with utils::modifyList(), so it may also override a
## tuning grid or the list of methods. Both paths are arguments of the
## scripts, which is how a run with another configuration and another
## output directory reuses this pipeline unchanged.

## The entries of a cell that are not arguments of simulate_wafc().
cell_reserved <- c("seed_key", "exclude")

## The methods the code knows; see R/methods.R.
known_methods <- function() names(study_methods)

#' Read the configuration of a run
#'
#' @param design Path of the design file.
#' @param study Path of the study file.
#' @return The merged configuration, a list, with every cell completed
#'   (its `name`, its `seed_key`, its `exclude`) and checked.
study_config <- function(design, study = file.path("config", "study.yaml")) {
  for (f in c(study, design)) {
    if (!file.exists(f)) stop("configuration file not found: ", f,
                              call. = FALSE)
  }
  st <- yaml::read_yaml(study, merge.precedence = "override")
  de <- yaml::read_yaml(design, merge.precedence = "override")
  cfg <- utils::modifyList(st, de)
  cfg[["files"]] <- list(study = normalizePath(study),
                         design = normalizePath(design))
  for (key in c("design", "replicates", "sample_sizes", "cells", "methods",
                "tuning", "seeds", "test_size", "grid_size", "code")) {
    if (is.null(cfg[[key]])) {
      stop("the configuration has no '", key, "' (files ", study, " and ",
           design, ").", call. = FALSE)
    }
  }
  cfg[["methods"]] <- unlist(cfg[["methods"]])
  cfg[["sample_sizes"]] <- as.integer(unlist(cfg[["sample_sizes"]]))
  cfg[["replicates"]] <- as.integer(cfg[["replicates"]])
  unknown <- setdiff(cfg[["methods"]], known_methods())
  if (length(unknown) > 0L) {
    stop("unknown method(s) in the configuration: ",
         paste(unknown, collapse = ", "), ". Known: ",
         paste(known_methods(), collapse = ", "), ".", call. = FALSE)
  }
  for (nm in names(cfg[["cells"]])) {
    cfg[["cells"]][[nm]] <- study_cell(cfg[["cells"]][[nm]], nm, cfg)
  }
  for (n in cfg[["sample_sizes"]]) {
    if (is.null(cfg[["tuning"]][[as.character(n)]])) {
      stop("no tuning entry for n = ", n, " (tuning: \"", n, "\").",
           call. = FALSE)
    }
  }
  cfg
}

## Completes and checks one cell.
study_cell <- function(cell, name, cfg) {
  cell[["name"]] <- name
  if (is.null(cell[["seed_key"]])) cell[["seed_key"]] <- name
  cell[["exclude"]] <- unlist(cell[["exclude"]])
  for (key in c("scenario", "p", "q")) {
    if (is.null(cell[[key]])) {
      stop("cell '", name, "' has no '", key, "'.", call. = FALSE)
    }
  }
  if (is.null(cell[["snr"]]) == is.null(cell[["sigma"]])) {
    stop("cell '", name, "' must give exactly one of 'snr' and 'sigma'.",
         call. = FALSE)
  }
  bad <- setdiff(cell[["exclude"]], known_methods())
  if (length(bad) > 0L) {
    stop("cell '", name, "' excludes unknown method(s): ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.null(cell[["structure"]])) {
    s <- do.call(rbind, lapply(cell[["structure"]],
                               function(r) as.character(unlist(r))))
    if (!identical(dim(s), as.integer(c(cell[["p"]], cell[["q"]])))) {
      stop("cell '", name, "': 'structure' must have p = ", cell[["p"]],
           " rows of q = ", cell[["q"]], " entries.", call. = FALSE)
    }
    cell[["structure"]] <- s
  }
  cell
}

## Checks, once the WAFC code is loaded, that simulate_wafc() accepts every
## entry of every cell. A cell that asks for an argument the loaded code
## does not have is reported by name before anything is fitted.
check_cells <- function(cfg) {
  args <- names(formals(simulate_wafc))
  bad <- character(0)
  for (cell in cfg[["cells"]]) {
    miss <- setdiff(setdiff(names(cell), c(cell_reserved, "name")), args)
    if (length(miss) > 0L) {
      bad <- c(bad, sprintf("%s (%s)", cell[["name"]],
                            paste(miss, collapse = ", ")))
    }
  }
  if (length(bad) > 0L) {
    stop("simulate_wafc() of the loaded code has no argument for: ",
         paste(bad, collapse = "; "), ". Leave these cells out (--cells) ",
         "or load a version of the code that has them.", call. = FALSE)
  }
  invisible(TRUE)
}

#' The tuning of one method at one sample size
#'
#' The entry of `tuning` for the sample size, with the method's own entry
#' (`tuning$<n>$methods$<method>`) merged over it. The folds and the basis
#' are shared by every method of a replicate and are read from the entry
#' of the sample size only.
#'
#' @return A list with `J`, `k`, `nfolds` and `basis`.
method_tuning <- function(cfg, n, method) {
  tn <- cfg[["tuning"]][[as.character(n)]]
  own <- if (method %in% names(tn[["methods"]])) tn[["methods"]][[method]]
  tn[["methods"]] <- NULL
  if (!is.null(own)) {
    own[["nfolds"]] <- NULL
    own[["basis"]] <- NULL
    tn <- utils::modifyList(tn, own)
  }
  tn[["J"]] <- as.integer(unlist(tn[["J"]]))
  tn[["k"]] <- as.integer(unlist(tn[["k"]]))
  tn[["nfolds"]] <- as.integer(tn[["nfolds"]])
  tn
}

## The tuning of a sample size without any method's own entry: the folds
## and the basis every method of a replicate shares.
size_tuning <- function(cfg, n) method_tuning(cfg, n, "")

## The methods a cell runs, in the order of the configuration.
cell_methods <- function(cfg, cell) {
  setdiff(cfg[["methods"]], cell[["exclude"]])
}

#' Load the WAFC code named in the configuration
#'
#' `code` is either the path of the loader of the WAFC code, which is
#' sourced, or "package:<name>" for a package that exports the same
#' functions, which is attached.
#'
#' @return What was loaded, invisibly.
load_wafc_code <- function(cfg) {
  path <- cfg[["code"]]
  if (startsWith(path, "package:")) {
    pkg <- sub("^package:", "", path)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    return(invisible(path))
  }
  if (!file.exists(path)) {
    stop("the WAFC code was not found at '", path, "' (entry 'code' of ",
         "the study configuration).", call. = FALSE)
  }
  source(path)
  invisible(normalizePath(path))
}
