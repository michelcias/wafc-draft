## R/run.R -- units, cache and the parallel driver.
##
## The unit of work is one method on one replicate: (cell, n, method,
## replicate). Each unit is cached in its own file,
##
##   <out>/units/<cell>/n<n>/<method>/rep<rrr>.rds,
##
## written atomically once the unit has finished, so a run that stops is
## resumed by running it again: finished units are read, not refitted. A
## unit also records the settings that determine it (the cell, the tuning
## of its method at its n, the sizes of the test sample and the grid, the
## seeds); a cached unit whose settings differ from the current
## configuration is stale, and the run stops before fitting anything unless
## it is told to refit the stale units. A method that fails inside a unit
## leaves a row with its message and is cached like any other row; a unit
## whose worker dies leaves no file and is refitted on the next run.

unit_path <- function(out, cell, n, method, rep) {
  file.path(out, "units", cell, paste0("n", n), method,
            sprintf("rep%03d.rds", rep))
}

#' The units of a run
#'
#' Every (cell, n, method, replicate) of the configuration, restricted by
#' the optional arguments, which never change what a unit draws.
#'
#' @return A data frame with one row per unit.
study_units <- function(cfg, cells = NULL, sample_sizes = NULL,
                        methods = NULL, reps = NULL) {
  cl <- names(cfg[["cells"]])
  if (!is.null(cells)) {
    bad <- setdiff(cells, cl)
    if (length(bad) > 0L) stop("unknown cell(s): ", paste(bad, collapse = ", "),
                               call. = FALSE)
    cl <- intersect(cl, cells)
  }
  ns <- cfg[["sample_sizes"]]
  if (!is.null(sample_sizes)) ns <- intersect(ns, sample_sizes)
  if (is.null(reps)) reps <- seq_len(cfg[["replicates"]])
  if (!is.null(methods)) {
    bad <- setdiff(methods, cfg[["methods"]])
    if (length(bad) > 0L) {
      stop("method(s) not in the configuration: ",
           paste(bad, collapse = ", "), call. = FALSE)
    }
  }
  out <- list()
  for (c in cl) {
    ms <- cell_methods(cfg, cfg[["cells"]][[c]])
    if (!is.null(methods)) ms <- intersect(ms, methods)
    if (length(ms) == 0L || length(ns) == 0L) next
    out[[length(out) + 1L]] <- expand.grid(
      cell = c, n = ns, method = ms, rep = reps, stringsAsFactors = FALSE,
      KEEP.OUT.ATTRS = FALSE)
  }
  units <- do.call(rbind, out)
  if (is.null(units)) stop("no unit to run with these restrictions.",
                           call. = FALSE)
  units
}

#' What determines a unit
unit_settings <- function(cfg, cell, n, method, rep) {
  cell[["exclude"]] <- NULL
  list(cell = cell, n = n, method = method, rep = rep,
       tuning = method_tuning(cfg, n, method),
       options = cfg[["method_options"]][[method]],
       test_size = cfg[["test_size"]], grid_size = cfg[["grid_size"]],
       seeds = unit_seeds(cfg, cell, n, rep))
}

## The wavelet tables of a run, one per basis, built once in the main
## process and shared by the workers.
basis_key <- function(b) {
  paste(b[["family"]], b[["filter_size"]], b[["prec_wavelet"]], sep = "/")
}

basis_tables <- function(cfg, ns) {
  tabs <- list()
  for (n in ns) {
    b <- cfg[["tuning"]][[as.character(n)]][["basis"]]
    key <- basis_key(b)
    if (is.null(tabs[[key]])) {
      tabs[[key]] <- WaveBased::wtable(
        family = b[["family"]], filter.size = as.integer(b[["filter_size"]]),
        prec.wavelet = as.integer(b[["prec_wavelet"]]), check = FALSE)
    }
  }
  tabs
}

#' Fit one unit
#'
#' @param code The provenance of the WAFC code (code_provenance()), kept
#'   with the unit.
#' @return The list cached for the unit: `rows`, `side`, `settings`, the
#'   seconds `elapsed` and `code`.
run_unit <- function(cfg, cell, n, method, rep, tables, code = NULL) {
  t0 <- elapsed()
  cl <- cfg[["cells"]][[cell]]
  st <- unit_settings(cfg, cl, n, method, rep)
  tn_n <- size_tuning(cfg, n)
  ctx <- replicate_context(cfg, cl, n, rep, st[["seeds"]], tn_n,
                           tables[[basis_key(tn_n[["basis"]])]])
  opts <- if (is.null(st[["options"]])) list() else st[["options"]]
  res <- study_methods[[method]][["fit"]](ctx, st[["tuning"]], opts)
  rows <- do.call(rbind, res[["rows"]])
  rows <- cbind(design = cfg[["design"]], rows, stringsAsFactors = FALSE)
  side <- res[["side"]]
  cv <- cfg[["curves"]]
  if (n %in% unlist(cv[["sample_sizes"]])) {
    for (lab in intersect(names(res[["ghat"]]), unlist(cv[["methods"]]))) {
      side[["curves"]] <- rbind(side[["curves"]],
                                curves_table(ctx, lab, res[["ghat"]][[lab]]))
    }
  }
  list(rows = rows, side = side, settings = st,
       elapsed = elapsed() - t0, code = code)
}

## Writes a unit atomically: a file that exists is a finished unit.
save_unit <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(obj, tmp)
  file.rename(tmp, path)
}

## The order units are dispatched in: the costliest first, so that the
## last units to finish are short ones.
unit_order <- function(units) {
  cost <- c(gam.gcv = 1, bsgl = 2, wafc = 3, klopp = 4, lasso = 5,
            oracle = 6, gam.reml = 7, vcbart = 8, aspline = 9, linear = 10)
  units[order(-units[["n"]], cost[units[["method"]]], units[["cell"]],
              units[["rep"]]), , drop = FALSE]
}

#' Run (or resume) the units of a configuration
#'
#' @param out The output directory of the run.
#' @param workers Number of forked workers.
#' @param refresh If TRUE, cached units whose settings differ from the
#'   configuration are refitted; if FALSE, their presence stops the run.
#' @param ... Restrictions passed to study_units().
#' @return The data frame of units, with their status, invisibly.
run_study <- function(cfg, out, workers = 1L, refresh = FALSE, ...) {
  check_cells(cfg)
  units <- study_units(cfg, ...)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  code <- code_provenance(cfg)
  record_session(cfg, out, code)
  units[["path"]] <- unit_path(out, units[["cell"]], units[["n"]],
                               units[["method"]], units[["rep"]])
  units[["status"]] <- "todo"
  have <- file.exists(units[["path"]])
  for (i in which(have)) {
    st <- unit_settings(cfg, cfg[["cells"]][[units[["cell"]][i]]],
                        units[["n"]][i], units[["method"]][i],
                        units[["rep"]][i])
    old <- readRDS(units[["path"]][i])[["settings"]]
    units[["status"]][i] <- if (identical(old, st)) "cached" else "stale"
  }
  nst <- sum(units[["status"]] == "stale")
  if (nst > 0L && !refresh) {
    stop(nst, " cached unit(s) were run with other settings (for instance ",
         units[["path"]][units[["status"]] == "stale"][1L], "). Rerun with ",
         "--refresh to refit them, or remove them.", call. = FALSE)
  }
  todo <- unit_order(units[units[["status"]] != "cached", , drop = FALSE])
  message(sprintf("[run] design '%s': %d unit(s), %d cached, %d to fit, %d worker(s)",
                  cfg[["design"]], nrow(units), sum(units[["status"]] == "cached"),
                  nrow(todo), workers))
  if (nrow(todo) == 0L) return(invisible(units))
  tables <- basis_tables(cfg, unique(todo[["n"]]))
  t0 <- elapsed()
  one <- function(i) {
    u <- todo[i, ]
    res <- try(run_unit(cfg, u[["cell"]], u[["n"]], u[["method"]],
                        u[["rep"]], tables, code), silent = TRUE)
    if (inherits(res, "try-error")) {
      message(sprintf("[run] %s n=%d %s rep %d: UNIT FAILED: %s", u[["cell"]],
                      u[["n"]], u[["method"]], u[["rep"]],
                      conditionMessage(attr(res, "condition"))))
      return("failed")
    }
    save_unit(res, u[["path"]])
    bad <- res[["rows"]][["error"]]
    message(sprintf("[run] %s n=%d %s rep %d: %.1f s%s", u[["cell"]],
                    u[["n"]], u[["method"]], u[["rep"]], res[["elapsed"]],
                    if (any(!is.na(bad))) paste0(" (row error: ",
                                                  bad[!is.na(bad)][1L], ")")
                    else ""))
    "done"
  }
  st <- if (workers > 1L) {
    parallel::mclapply(seq_len(nrow(todo)), one, mc.cores = workers,
                       mc.preschedule = FALSE)
  } else lapply(seq_len(nrow(todo)), one)
  st <- vapply(st, function(s) if (is.character(s)) s else "failed", "")
  units[["status"]][match(todo[["path"]], units[["path"]])] <- st
  message(sprintf("[run] %d unit(s) fitted, %d failed, %.1f min",
                  sum(st == "done"), sum(st == "failed"),
                  (elapsed() - t0) / 60))
  invisible(units)
}

#' Where the WAFC code of a run comes from
#'
#' For a package, its version and the commit it was installed from. For a
#' path in a git working tree: the commit, the tree of its directory at
#' that commit, and whether the directory has uncommitted changes.
#' Otherwise only the path.
#'
#' @return A one-line description.
code_provenance <- function(cfg) {
  if (startsWith(cfg[["code"]], "package:")) {
    pkg <- sub("^package:", "", cfg[["code"]])
    d <- utils::packageDescription(pkg)
    return(sprintf("package %s %s%s", pkg, d[["Version"]],
                   if (is.null(d[["RemoteSha"]])) "" else
                     paste0(" @ ", d[["RemoteSha"]])))
  }
  dir <- dirname(normalizePath(cfg[["code"]]))
  git <- function(...) {
    out <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(dir), ...),
                                             stdout = TRUE, stderr = FALSE)),
                    error = function(e) character(0))
    if (!is.null(attr(out, "status"))) character(0) else out
  }
  commit <- git("rev-parse", "HEAD")
  if (length(commit) == 0L) return(paste0(dir, " (not in a git working tree)"))
  tree <- git("rev-parse", "HEAD:./")
  dirty <- length(git("status", "--porcelain", "--", ".")) > 0L
  sprintf("%s at commit %s (tree %s)%s", dir, commit, tree,
          if (dirty) ", with uncommitted changes" else "")
}

## The session that produced a run, rewritten at every run.
record_session <- function(cfg, out, code = code_provenance(cfg)) {
  wb <- utils::packageDescription("WaveBased")
  writeLines(c(
    sprintf("# session of a run of design '%s', %s", cfg[["design"]],
            format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("# study file: %s", cfg[["files"]][["study"]]),
    sprintf("# design file: %s", cfg[["files"]][["design"]]),
    sprintf("# WAFC code: %s", code),
    sprintf("# WaveBased: %s%s", wb[["Version"]],
            if (is.null(wb[["RemoteSha"]])) "" else
              paste0(" @ ", wb[["RemoteSha"]])),
    sprintf("# master seed: %s", cfg[["seeds"]][["master"]]),
    "",
    utils::capture.output(utils::sessionInfo())
  ), file.path(out, "sessionInfo.txt"))
}
