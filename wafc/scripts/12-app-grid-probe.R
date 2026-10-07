## wafc/scripts/12-app-grid-probe.R -- step E6.2b (decision D75(a)): do the
## grids of the data applications cap the methods? From wafc-studies/ (so
## that renv is active):
##
##     Rscript ../wafc/scripts/12-app-grid-probe.R [parts] [workers]
##
## parts: "fit", "read" or "fit,read" (default); workers: default 2.
## Outputs, not versioned, in wafc/cache/e62b/: units/<job>.rds (the run
## resumes from them), probe.csv, readings.csv, and the log on stdout.
##
## The fits of the chain (wafc-studies/outputs/application/units/, grids of
## config/application.yaml: k in 5, 10, 20, 40, 80 for the splines, J in
## 2:8 for the WAFC) chose the top of their grids in marylebone.ukair: k = 80
## in the 21 fits of both splines, J = 8 in the WAFC on the whole sample.
## This probe refits on the same data, partitions, folds and random streams
## (application_data(), application_context() and the seeds of the
## compendium), with the grids extended:
##
## - the two splines, gam.reml and gam.gcv, with k in 80, 120, 160, on the
##   whole sample (split 0) and partitions 1 to 3 of both applications;
## - the WAFC with J in 2:9 on the whole sample and partitions 1 and 2 of
##   marylebone.ukair. Every candidate J is scored on the same folds
##   independently of the others, so 2:9 chooses 9 exactly when the
##   cross-validated error at J = 9 is below the one at the J the chain
##   chose; the probe fits that J and J = 9 only, and when J = 9 does not
##   win the fit is the chain's own and is not repeated. In beijing.heat the
##   WAFC chooses J = 3 or 4 in every fit of the chain, far from the top,
##   and is not probed.
##
## The criterion, declared before the run (docs/TAREFA.md, E6.2b): the grid
## of the applications grows if, in some application and method, the ratio
## of the test RMSE to the current grid falls below 0.99 on average over
## the partitions probed, or if the reading of the step of marylebone moves
## by more than the range of the partitions; otherwise it stays.
##
## The jobs run in the order below, the costliest of marylebone first, so
## that a run stopped early has marylebone.

args <- commandArgs(trailingOnly = TRUE)
parts <- strsplit(if (length(args) >= 1L) args[1L] else "fit,read", ",")[[1L]]
workers <- if (length(args) >= 2L) as.integer(args[2L]) else 2L
source(file.path("scripts", "00_setup.R"))
cfg <- application_config()
load_wafc_code(cfg)
chain <- file.path("outputs", "application")
out <- file.path("..", "wafc", "cache", "e62b")
dir.create(file.path(out, "units"), recursive = TRUE, showWarnings = FALSE)
cat("WAFC code:", code_provenance(cfg), "\n")

k_ext <- c(80L, 120L, 160L)
chain_unit <- function(base, method, split) {
  readRDS(application_unit_path(chain, base, method, split))
}

jobs <- list()
add <- function(base, method, split) {
  jobs[[length(jobs) + 1L]] <<- list(base = base, method = method,
                                     split = split)
}
add("marylebone.ukair", "wafc", 0L)
for (s in 0:3) add("marylebone.ukair", "gam.gcv", s)
for (s in 0:3) add("marylebone.ukair", "gam.reml", s)
add("marylebone.ukair", "wafc", 1L)
add("marylebone.ukair", "wafc", 2L)
for (s in 0:3) add("beijing.heat", "gam.reml", s)
for (s in 0:3) add("beijing.heat", "gam.gcv", s)
job_path <- function(j) {
  file.path(out, "units", sprintf("%s-%s-split%02d.rds", j[["base"]],
                                  j[["method"]], j[["split"]]))
}

## The WAFC on the grid (J of the chain, 9): the rows of application.R when
## J = 9 wins, the cross-validation table otherwise.
probe_wafc <- function(ctx, J_chain) {
  restore_stream(ctx)
  t0 <- elapsed()
  cv <- cv.wafc(ctx[["train"]][["x"]], ctx[["train"]][["u"]],
                ctx[["train"]][["y"]], J = c(J_chain, 9L), penalty = "block",
                foldid = ctx[["foldid"]], wavelet.table = ctx[["table"]],
                threshold = "none")
  res <- list(cvtab = cv[["cvtab"]], J = cv[["J.min"]],
              search = elapsed() - t0)
  if (cv[["J.min"]] != 9L) return(res)
  ff <- wafc_threshold_folds(cv)
  rows <- list()
  for (lab in c("wafc.cv1se", "wafc.cv")) {
    th <- wafc_threshold(cv, rule = sub("^wafc\\.", "", lab), fold.fits = ff)
    ex <- th[["extra"]]
    rows[[lab]] <- application_read(ctx, th, elapsed() - t0,
                                    list(J = ex[["J"]], lambda = ex[["lambda"]],
                                         t = ex[["t"]], t_frac = ex[["c"]],
                                         norm = ex[["norm"]]))
  }
  c(res, list(rows = rows))
}

run_job <- function(j) {
  path <- job_path(j)
  if (file.exists(path)) return("cached")
  t0 <- elapsed()
  data <- datas[[j[["base"]]]]
  ctx <- application_context(cfg, data, j[["split"]], table)
  if (j[["method"]] == "wafc") {
    J_chain <- chain_unit(j[["base"]], "wafc", j[["split"]])[["rows"]][[
      "wafc.cv1se"]][["J"]]
    res <- probe_wafc(ctx, J_chain)
    res[["J_chain"]] <- J_chain
  } else {
    tn <- cfg[["tuning"]]
    tn[["k"]] <- k_ext
    res <- list(rows = application_methods[[j[["method"]]]][["fit"]](
      ctx, tn, cfg[["method_options"]][[j[["method"]]]]))
  }
  res[["job"]] <- j
  res[["sizes"]] <- ctx[["sizes"]]
  res[["elapsed"]] <- elapsed() - t0
  save_unit(res, path)
  message(sprintf("[probe] %s %s split %d: %.0f s", j[["base"]], j[["method"]],
                  j[["split"]], res[["elapsed"]]))
  "done"
}

if ("fit" %in% parts) {
  todo <- Filter(function(j) !file.exists(job_path(j)), jobs)
  message(sprintf("[probe] %d job(s) to fit, %d worker(s)", length(todo),
                  workers))
  datas <- lapply(stats::setNames(nm = unique(vapply(todo, `[[`, "", "base"))),
                  function(b) application_data(cfg, b))
  b <- cfg[["tuning"]][["basis"]]
  table <- WaveBased::wtable(family = b[["family"]],
                             filter.size = as.integer(b[["filter_size"]]),
                             prec.wavelet = as.integer(b[["prec_wavelet"]]),
                             check = FALSE)
  st <- parallel::mclapply(todo, function(j) {
    tryCatch(run_job(j), error = function(e) {
      message(sprintf("[probe] %s %s split %d FAILED: %s", j[["base"]],
                      j[["method"]], j[["split"]], conditionMessage(e)))
      "failed"
    })
  }, mc.cores = workers, mc.preschedule = FALSE)
  message("[probe] fit exit")
}

if ("read" %in% parts) {
  rows <- list()
  cols <- list()
  for (j in jobs) {
    path <- job_path(j)
    if (!file.exists(path)) next
    p <- readRDS(path)
    if (isTRUE(p[["skipped"]])) next  # not fitted; the reason is in the file
    ch <- chain_unit(j[["base"]], j[["method"]], j[["split"]])
    labs <- names(ch[["rows"]])
    for (lab in labs) {
      old <- ch[["rows"]][[lab]]
      unchanged <- j[["method"]] == "wafc" && is.null(p[["rows"]])
      new <- if (unchanged) old else p[["rows"]][[lab]]
      cvt <- p[["cvtab"]]
      rows[[length(rows) + 1L]] <- data.frame(
        base = j[["base"]], label = lab, split = j[["split"]],
        chosen_chain = if (j[["method"]] == "wafc") old[["J"]] else
          paste(old[["k"]], collapse = ","),
        chosen_probe = if (j[["method"]] == "wafc") p[["J"]] else
          paste(new[["k"]], collapse = ","),
        top_probe = if (j[["method"]] == "wafc") p[["J"]] == 9L else
          isTRUE(new[["k_top"]]),
        cv_gain_J9 = if (is.null(cvt)) NA_real_ else
          cvt[[3L]][cvt[["J"]] == 9L] - cvt[[3L]][cvt[["J"]] == p[["J_chain"]]],
        kcheck_min_p = if (is.null(new[["kcheck_min_p"]])) NA_real_ else
          new[["kcheck_min_p"]],
        kcheck_n_low = if (is.null(new[["kcheck_n_low"]])) NA_integer_ else
          new[["kcheck_n_low"]],
        rmse_chain = old[["rmse"]], rmse_probe = new[["rmse"]],
        ratio = new[["rmse"]] / old[["rmse"]],
        unchanged = unchanged, seconds = p[["elapsed"]],
        stringsAsFactors = FALSE)
      key <- paste(j[["base"]], lab)
      if (is.null(cols[[key]])) {
        cols[[key]] <- list(base = j[["base"]], label = lab,
                            col = list(rows = list(), grid = ch[["grid"]]))
      }
      cols[[key]][["col"]][["rows"]][[length(cols[[key]][["col"]][["rows"]]) +
                                        1L]] <-
        list(split = j[["split"]], label = lab, row = new, sizes = p[["sizes"]])
    }
  }
  tab <- do.call(rbind, rows)
  utils::write.csv(tab, file.path(out, "probe.csv"), row.names = FALSE)
  print(format(tab[setdiff(names(tab), "base")], digits = 4), row.names = FALSE)

  ## The verdict on the error: the mean over the partitions probed (split
  ## 0 has no test sample) of the ratio to the current grid.
  pt <- tab[tab[["split"]] > 0, ]
  if (nrow(pt) > 0L) {
    v <- stats::aggregate(ratio ~ base + label, pt, mean)
    v[["partitions"]] <- stats::aggregate(ratio ~ base + label, pt,
                                          length)[["ratio"]]
    cat("\nmean ratio of the test RMSE to the current grid:\n")
    print(v, row.names = FALSE)
  }

  ## The readings of marylebone, with the definitions of the report
  ## (R/application_report.R), for the probed fits and for the same splits
  ## of the chain, and the range of the chain over its 20 partitions.
  ## The points read are those the figure of the chain draws (drawn_grid()
  ## on the chain's units), for the probed fits as well.
  data <- application_data(cfg, "marylebone.ukair")
  spec <- cfg[["report"]][["readings"]][["marylebone.ukair"]]
  full <- collect_application(cfg, chain, "marylebone.ukair")
  drawn <- drawn_grid(full, data, cfg[["report"]][["figure"]])
  rd <- list()
  for (key in names(cols)) {
    if (cols[[key]][["base"]] != "marylebone.ukair") next
    lab <- cols[[key]][["label"]]
    r <- application_readings(cols[[key]][["col"]], lab, data, spec, drawn)
    rd[[key]] <- cbind(source = "probe", r)
  }
  for (lab in unique(vapply(cols, `[[`, "", "label"))) {
    if (!any(vapply(cols, function(z) z[["base"]] == "marylebone.ukair" &&
                      z[["label"]] == lab, NA))) next
    r <- application_readings(full, lab, data, spec, drawn)
    rd[[paste("chain", lab)]] <- cbind(source = "chain", r)
  }
  rd <- do.call(rbind, rd)
  utils::write.csv(rd, file.path(out, "readings.csv"), row.names = FALSE)
  keep <- rd[["quantity"]] %in% c("rise", "fraction_at_2003-01-01",
                                  "fraction_at_2003-07-01",
                                  "peak_above_before", "wind_range")
  cat("\nreadings of marylebone (probe: split 0 and partitions probed; chain: split 0 and its 20 partitions):\n")
  print(format(rd[keep, ], digits = 4), row.names = FALSE)
}
