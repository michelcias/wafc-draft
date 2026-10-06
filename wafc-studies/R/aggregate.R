## R/aggregate.R -- from the cached units to the tables.
##
## Reads every unit of a configuration from the output directory of its
## run and writes, under <out>/tables/:
##
##   replicates.csv   one row per fit (R/metrics.R lists the columns);
##   components.csv   the integrated squared error of every block of every
##                    fit, with whether the block is active and was kept;
##   summary.csv      by cell, n and row label: the replicates, the failed
##                    fits, mean, standard error and median of the errors,
##                    P(S-hat = S) with its standard error, the mean false
##                    positive and false negative blocks, the median time,
##                    the mean and the distribution of the J (or k) chosen,
##                    and the fraction of replicates at the top of the grid;
##   paired.csv       by cell, n and row label: the median, within
##                    replicate, of the ratio of each error to the one of
##                    the reference row (`reference` of the configuration),
##                    and the fraction of replicates in which the row is
##                    below the reference;
##   slopes.csv       by cell and row label, the least squares slope of the
##                    log of the mean squared prediction error and of the
##                    mean integrated squared error against log n, where
##                    the design has more than one sample size;
##   gam_k.csv, kcheck.csv, convergence.csv, threshold.csv, curves.csv:
##                    the side tables of the fits, bound as they are.
##
## Units missing from the cache are counted and listed in missing.csv, and
## the tables are written from the units present.

#' Read the cached units of a configuration
#'
#' @return A list with `rows`, the side tables by name, and `missing`, the
#'   units without a file.
collect_units <- function(cfg, out, ...) {
  units <- study_units(cfg, ...)
  units[["path"]] <- unit_path(out, units[["cell"]], units[["n"]],
                               units[["method"]], units[["rep"]])
  have <- file.exists(units[["path"]])
  rows <- list()
  side <- list()
  stale <- 0L
  for (i in which(have)) {
    u <- readRDS(units[["path"]][i])
    st <- unit_settings(cfg, cfg[["cells"]][[units[["cell"]][i]]],
                        units[["n"]][i], units[["method"]][i],
                        units[["rep"]][i])
    if (!identical(u[["settings"]], st)) stale <- stale + 1L
    rows[[length(rows) + 1L]] <- u[["rows"]]
    for (nm in names(u[["side"]])) {
      side[[nm]][[length(side[[nm]]) + 1L]] <- u[["side"]][[nm]]
    }
  }
  if (stale > 0L) {
    warning(stale, " unit(s) in the cache were run with settings other than ",
            "the configuration's; they are read as they are.", call. = FALSE)
  }
  list(rows = bind_rows(rows),
       side = lapply(side, bind_rows),
       missing = units[!have, c("cell", "n", "method", "rep")])
}

## rbind over data frames whose columns may differ, filling with NA.
bind_rows <- function(dfs) {
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0L) return(NULL)
  all <- unique(unlist(lapply(dfs, names)))
  do.call(rbind, lapply(dfs, function(d) {
    for (v in setdiff(all, names(d))) d[[v]] <- NA
    d[all]
  }))
}

se <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) < 2L) NA_real_ else stats::sd(v) / sqrt(length(v))
}

## Summary of one group of rows (one cell, n and row label).
summarise_group <- function(d) {
  ok <- is.na(d[["error"]])
  e <- d[ok, , drop = FALSE]
  exact <- e[["n_true"]] == e[["n_active"]] & e[["n_false"]] == 0L
  out <- data.frame(design = d[["design"]][1L], cell = d[["cell"]][1L],
                    n = d[["n"]][1L], method = d[["method"]][1L],
                    replicates = nrow(d), failed = sum(!ok),
                    stringsAsFactors = FALSE)
  for (v in c("rmse_f", "rmse_y", "mse_beta", "ise", "ise_active",
              "ise_null")) {
    x <- e[[v]]
    out[[paste0(v, "_mean")]] <- if (all(is.na(x))) NA_real_ else
      mean(x, na.rm = TRUE)
    out[[paste0(v, "_se")]] <- se(x)
    out[[paste0(v, "_median")]] <- if (all(is.na(x))) NA_real_ else
      stats::median(x, na.rm = TRUE)
  }
  out[["exact_support"]] <- if (all(is.na(exact))) NA_real_ else
    mean(exact, na.rm = TRUE)
  out[["exact_support_se"]] <- se(as.numeric(exact))
  out[["false_pos_mean"]] <- if (all(is.na(e[["n_false"]]))) NA_real_ else
    mean(e[["n_false"]], na.rm = TRUE)
  out[["false_neg_mean"]] <- if (all(is.na(e[["n_true"]]))) NA_real_ else
    mean(e[["n_active"]] - e[["n_true"]], na.rm = TRUE)
  out[["time_median"]] <- if (all(is.na(e[["time"]]))) NA_real_ else
    stats::median(e[["time"]], na.rm = TRUE)
  dim <- if (all(is.na(e[["J"]]))) e[["k"]] else e[["J"]]
  out[["chosen_mean"]] <- if (all(is.na(dim))) NA_real_ else
    mean(dim, na.rm = TRUE)
  out[["chosen_table"]] <- if (all(is.na(dim))) NA_character_ else {
    tb <- table(dim)
    paste(names(tb), tb, sep = ":", collapse = " ")
  }
  out[["top_fraction"]] <- if (all(is.na(e[["top"]]))) NA_real_ else
    mean(e[["top"]], na.rm = TRUE)
  out[["kcheck_flagged"]] <- if (all(is.na(e[["kcheck_n_low"]]))) NA_real_
    else mean(e[["kcheck_n_low"]] > 0, na.rm = TRUE)
  out
}

summarise_rows <- function(rows, labels) {
  key <- paste(rows[["cell"]], rows[["n"]], rows[["method"]], sep = "\r")
  out <- do.call(rbind, lapply(split(rows, key), summarise_group))
  out[order(out[["cell"]], out[["n"]],
            match(out[["method"]], labels)), , drop = FALSE]
}

#' Paired comparison with the reference row, replicate by replicate
paired_rows <- function(rows, reference, labels) {
  vars <- c("rmse_f", "mse_beta", "ise")
  ok <- rows[is.na(rows[["error"]]), , drop = FALSE]
  base <- ok[ok[["method"]] == reference, c("cell", "n", "rep", vars)]
  if (nrow(base) == 0L) return(NULL)
  names(base)[-(1:3)] <- paste0(vars, "_ref")
  m <- merge(ok, base, by = c("cell", "n", "rep"))
  key <- paste(m[["cell"]], m[["n"]], m[["method"]], sep = "\r")
  out <- do.call(rbind, lapply(split(m, key), function(d) {
    r <- data.frame(cell = d[["cell"]][1L], n = d[["n"]][1L],
                    method = d[["method"]][1L], pairs = nrow(d),
                    stringsAsFactors = FALSE)
    for (v in vars) {
      ref <- d[[paste0(v, "_ref")]]
      ## a reference error of zero (every block zeroed in the null cell)
      ## leaves the ratio undefined, and the pair out of the median
      ratio <- ifelse(ref > 0, d[[v]] / ref, NA_real_)
      r[[paste0(v, "_ratio_median")]] <- if (all(is.na(ratio))) NA_real_
        else stats::median(ratio, na.rm = TRUE)
      r[[paste0(v, "_below_ref")]] <- if (all(is.na(d[[v]]))) NA_real_
        else mean(d[[v]] < ref, na.rm = TRUE)
    }
    r
  }))
  out[order(out[["cell"]], out[["n"]],
            match(out[["method"]], labels)), , drop = FALSE]
}

#' Slope of the log error against log n
slope_rows <- function(summary) {
  key <- paste(summary[["cell"]], summary[["method"]], sep = "\r")
  out <- lapply(split(summary, key), function(d) {
    d <- d[order(d[["n"]]), , drop = FALSE]
    if (nrow(d) < 2L) return(NULL)
    fit_slope <- function(v) {
      ok <- is.finite(v) & v > 0
      if (sum(ok) < 2L) return(NA_real_)
      unname(stats::coef(stats::lm(log(v[ok]) ~ log(d[["n"]][ok])))[2L])
    }
    data.frame(cell = d[["cell"]][1L], method = d[["method"]][1L],
               sizes = paste(d[["n"]], collapse = ","),
               slope_mse_f = fit_slope(d[["rmse_f_mean"]]^2),
               slope_ise = fit_slope(d[["ise_mean"]]),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Write the tables of a run
#'
#' @return The path of the tables directory, invisibly.
aggregate_study <- function(cfg, out, ...) {
  res <- collect_units(cfg, out, ...)
  dir <- file.path(out, "tables")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  write <- function(x, name) {
    if (is.null(x) || nrow(x) == 0L) return(invisible(NULL))
    utils::write.csv(x, file.path(dir, name), row.names = FALSE)
    message("[aggregate] wrote ", file.path(dir, name), " (", nrow(x),
            " rows)")
  }
  labels <- row_labels(cfg[["methods"]])
  rows <- res[["rows"]]
  if (is.null(rows)) stop("no unit of this configuration in ", out, ".",
                          call. = FALSE)
  rows <- rows[order(rows[["cell"]], rows[["n"]], rows[["rep"]],
                     match(rows[["method"]], labels)), , drop = FALSE]
  write(rows, "replicates.csv")
  summ <- summarise_rows(rows, labels)
  summ[["label"]] <- unlist(cfg[["labels"]])[summ[["method"]]]
  write(summ, "summary.csv")
  write(paired_rows(rows, cfg[["reference"]], labels), "paired.csv")
  write(slope_rows(summ), "slopes.csv")
  for (nm in names(res[["side"]])) {
    write(res[["side"]][[nm]], paste0(nm, ".csv"))
  }
  write(res[["missing"]], "missing.csv")
  if (nrow(res[["missing"]]) > 0L) {
    message("[aggregate] ", nrow(res[["missing"]]), " unit(s) missing; ",
            "see missing.csv")
  }
  invisible(dir)
}
