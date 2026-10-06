## R/application.R -- the data applications: partitions, fits and cache.
##
## The applications (R/application_data.R) are run the way the simulation
## study is (R/run.R). The unit of work is one method on one split of one
## application, cached in its own file, written atomically,
##
##   <out>/units/<application>/<method>/split<ss>.rds,
##
## with the settings that determine it; a cached unit whose settings
## differ from the configuration stops the run unless it is told to refit.
##
## A split is either the whole sample (split 0, the fit the figures show)
## or one of `partitions` random partitions by block (splits 1, 2, ...): a
## fraction `test_fraction` of the blocks, drawn at random, is the test
## sample. In both, the blocks of the training sample are spread at random
## over `nfolds` folds, which every method tuned by cross-validation uses.
## The partitions give the prediction error on held-out weeks, with a
## standard error between partitions, and the frequency with which each
## block (l, m) is kept; the training samples of two partitions share most
## of their observations, which the standard error of the report corrects
## for (R/application_report.R).
##
## Every method is tuned in its own terms. The WAFC chooses (J, lambda) by
## cross-validation on the blocked folds and the threshold on the norms of
## the blocks from the same folds; mgcv chooses the smoothing parameters
## and the basis dimension by REML or by GCV, the criteria it is used with;
## the linear model has nothing to tune. None of them models the dependence
## of the errors.

#' Read the configuration of the applications
#'
#' The application file (config/application.yaml) with three entries taken
#' from the study file unless it sets them: `code` (the WAFC code to load),
#' `workers`, and the master seed and radix of `seeds`, so that the study
#' and the applications share one source of randomness and one version of
#' the code.
#'
#' @return The configuration, a list, checked.
application_config <- function(path = file.path("config", "application.yaml"),
                               study = file.path("config", "study.yaml")) {
  for (f in c(study, path)) {
    if (!file.exists(f)) stop("configuration file not found: ", f,
                              call. = FALSE)
  }
  st <- yaml::read_yaml(study, merge.precedence = "override")
  cfg <- yaml::read_yaml(path, merge.precedence = "override")
  for (key in c("code", "workers")) {
    if (is.null(cfg[[key]])) cfg[[key]] <- st[[key]]
  }
  for (key in c("master", "rep_radix")) {
    if (is.null(cfg[["seeds"]][[key]])) {
      cfg[["seeds"]][[key]] <- st[["seeds"]][[key]]
    }
  }
  cfg[["files"]] <- list(study = normalizePath(study),
                         design = normalizePath(path))
  for (key in c("design", "bases", "data_dir", "partitions", "test_fraction",
                "methods", "tuning", "grid_size", "grid_quantiles", "code")) {
    if (is.null(cfg[[key]])) {
      stop("the configuration has no '", key, "' (files ", study, " and ",
           path, ").", call. = FALSE)
    }
  }
  cfg[["bases"]] <- unlist(cfg[["bases"]])
  cfg[["methods"]] <- unlist(cfg[["methods"]])
  cfg[["partitions"]] <- as.integer(cfg[["partitions"]])
  cfg[["tuning"]][["J"]] <- as.integer(unlist(cfg[["tuning"]][["J"]]))
  cfg[["tuning"]][["k"]] <- as.integer(unlist(cfg[["tuning"]][["k"]]))
  cfg[["tuning"]][["nfolds"]] <- as.integer(cfg[["tuning"]][["nfolds"]])
  bad <- setdiff(cfg[["bases"]], names(application_bases))
  if (length(bad) > 0L) {
    stop("unknown application(s): ", paste(bad, collapse = ", "), ". Known: ",
         paste(names(application_bases), collapse = ", "), ".", call. = FALSE)
  }
  bad <- setdiff(cfg[["bases"]], unlist(cfg[["seeds"]][["keys"]]))
  if (length(bad) > 0L) {
    stop("application(s) not in seeds$keys: ", paste(bad, collapse = ", "),
         "; append them there.", call. = FALSE)
  }
  bad <- setdiff(cfg[["methods"]], names(application_methods))
  if (length(bad) > 0L) {
    stop("unknown method(s): ", paste(bad, collapse = ", "), ". Known: ",
         paste(names(application_methods), collapse = ", "), ".",
         call. = FALSE)
  }
  if (cfg[["partitions"]] + 1L >= as.integer(cfg[["seeds"]][["rep_radix"]])) {
    stop("more partitions than seeds$rep_radix allows.", call. = FALSE)
  }
  cfg
}

#' The seed of one split
#'
#' Through unit_seeds() of the simulation study: the application is the
#' seed key (its position in `seeds$keys`, which is append-only), there is
#' one sample per application (entered as the single sample size 0), and
#' the split is the replicate, the whole sample being replicate 1 and
#' partition s replicate s + 1. The seed is that of the test stream: it
#' draws the blocks of the test sample, and then the folds.
application_seed <- function(cfg, base, split) {
  sc <- list(seeds = list(master = cfg[["seeds"]][["master"]],
                          keys = cfg[["seeds"]][["keys"]],
                          sample_sizes = list(0L),
                          rep_radix = cfg[["seeds"]][["rep_radix"]]))
  unit_seeds(sc, list(seed_key = base), 0L, split + 1L)[["test"]]
}

#' One partition by block
#'
#' From the seed: the blocks of the test sample, a fraction
#' `test_fraction` of all blocks drawn at random (none for the whole
#' sample), and then the fold of each block of the training sample, drawn
#' at random so that the folds have as equal numbers of blocks as possible.
#'
#' @return A list with the rows of the training and test samples, the fold
#'   of each training row, and the numbers of blocks.
application_split <- function(block, seed, test_fraction, nfolds) {
  set.seed(seed)
  ub <- sort(unique(block))
  te_b <- if (test_fraction > 0) {
    sample(ub, round(test_fraction * length(ub)))
  } else ub[0L]
  te <- which(block %in% te_b)
  tr <- which(!(block %in% te_b))
  tb <- sort(unique(block[tr]))
  fold_of <- stats::setNames(sample(rep_len(seq_len(nfolds), length(tb))), tb)
  list(train = tr, test = te,
       foldid = unname(fold_of[as.character(block[tr])]),
       test_blocks = length(te_b), blocks = length(ub))
}

#' Everything a method needs from one split
#'
#' Leaves the random stream where the partition leaves it: every method
#' starts from there (restore_stream()).
application_context <- function(cfg, data, split, table) {
  sp <- application_split(data[["block"]],
                          application_seed(cfg, data[["base"]], split),
                          if (split == 0L) 0 else cfg[["test_fraction"]],
                          cfg[["tuning"]][["nfolds"]])
  tr <- sp[["train"]]
  te <- sp[["test"]]
  part <- function(i) list(x = data[["x"]][i, , drop = FALSE],
                           u = data[["u"]][i, , drop = FALSE],
                           y = data[["y"]][i])
  list(base = data[["base"]], split = split, train = part(tr),
       test = if (length(te) > 0L) part(te) else NULL,
       grid = data[["grid"]], foldid = sp[["foldid"]], table = table,
       sizes = list(ntrain = length(tr), ntest = length(te),
                    test_blocks = sp[["test_blocks"]],
                    blocks = sp[["blocks"]],
                    sd_test = if (length(te) > 0L)
                      stats::sd(data[["y"]][te]) else NA_real_),
       rng = get(".Random.seed", envir = globalenv()))
}

## Norm of each resolution level of each block of a fit on the wavelet
## design: where in scale a component carries its energy (a step puts it
## in the fine levels, a smooth cycle in the coarse ones).
application_level_norms <- function(b, design) {
  out <- list()
  for (l in seq_len(design[["p"]])) {
    for (m in seq_len(design[["q"]])) {
      nm <- wafc_block_name(design, l, m)
      idx <- design[["blocks"]][[nm]]
      Jm <- design[["J"]][m]
      lv <- rep(0:(Jm - 1L), times = 2^(0:(Jm - 1L)))
      out[[nm]] <- as.numeric(tapply(b[idx], lv, function(z) sqrt(sum(z^2))))
    }
  }
  out
}

application_rmse <- function(a, b) sqrt(mean((a - b)^2))

#' How every fit of the applications is read
#'
#' The root mean squared error on the test sample (none for the whole
#' sample) and on the training sample; the blocks kept; the components
#' g_lm on the grid, centred there as wafc_grid_components() returns them
#' (every method fixes the level of a component by its own convention, so
#' the shape is what compares); and, to put beta_l back together, the
#' levels c_l, the mean on the grid that the centring removed (`offset`)
#' and the mean of each component over the training sample (`gmean`), all
#' in the convention of the fit, in which beta_l = c_l + sum_m g_lm.
#'
#' @return A list, the row of the fit.
application_read <- function(ctx, f, time, extra = list()) {
  g <- wafc_grid_components(f, ctx[["grid"]])
  means <- function(h) matrix(vapply(h, mean, 0), nrow(h), ncol(h),
                              dimnames = dimnames(h))
  offset <- if (is.null(g)) NULL else means(f[["g"]](ctx[["grid"]]))
  gmean <- if (is.null(g)) NULL else means(f[["g"]](ctx[["train"]][["u"]]))
  c(list(rmse = if (is.null(ctx[["test"]])) NA_real_ else
           application_rmse(ctx[["test"]][["y"]],
                predict(f, ctx[["test"]][["x"]], ctx[["test"]][["u"]])),
         rmse_train = application_rmse(ctx[["train"]][["y"]], predict(f)),
         time = time, blocks = f[["blocks"]], cc = f[["cc"]],
         components = g, offset = offset, gmean = gmean),
    extra)
}

application_fail <- function(err) {
  list(error = gsub("\\s+", " ", if (inherits(err, "try-error"))
    conditionMessage(attr(err, "condition")) else as.character(err)))
}

#' The WAFC: one search and its two thresholds
#'
#' cv.wafc() with the block lasso, (J, lambda) by cross-validation on the
#' blocked folds over the grid of J, then the threshold on the norms of the
#' blocks from the same folds: row "wafc.cv1se" with the threshold at one
#' standard error (the estimator; the structure the article reads) and row
#' "wafc.cv" with the threshold of the smallest cross-validated error (the
#' predictor; it keeps the small blocks the first rule drops). The folds
#' are refitted once for the two rules.
application_fit_wafc <- function(ctx, tuning, options) {
  restore_stream(ctx)
  labels <- c(wafc.cv1se = "cv1se", wafc.cv = "cv")
  rows <- list()
  t0 <- elapsed()
  cv <- try(do.call(cv.wafc, c(list(ctx[["train"]][["x"]],
                                     ctx[["train"]][["u"]],
                                     ctx[["train"]][["y"]],
                                     J = tuning[["J"]], penalty = "block",
                                     foldid = ctx[["foldid"]],
                                     wavelet.table = ctx[["table"]],
                                     threshold = "none"), options)),
            silent = TRUE)
  if (inherits(cv, "try-error")) {
    for (lab in names(labels)) rows[[lab]] <- application_fail(cv)
    return(rows)
  }
  t_search <- elapsed() - t0
  top <- at_top(ctx, cv[["J.min"]], tuning[["J"]])
  t0 <- elapsed()
  ff <- try(wafc_threshold_folds(cv), silent = TRUE)
  t_folds <- elapsed() - t0
  for (lab in names(labels)) {
    th <- if (inherits(ff, "try-error")) ff else
      try(wafc_threshold(cv, rule = labels[[lab]], fold.fits = ff),
          silent = TRUE)
    rd <- if (inherits(th, "try-error")) th else {
      ex <- th[["extra"]]
      try(application_read(ctx, th, t_search + t_folds + th[["time"]],
                           list(J = ex[["J"]], top = top,
                                lambda = ex[["lambda"]], t = ex[["t"]],
                                t_frac = ex[["c"]], nzero = ex[["nzero"]],
                                norm = ex[["norm"]],
                                candidates = ex[["candidates"]],
                                levels = application_level_norms(th[["coef"]][-1L],
                                                     th[["design"]]))),
          silent = TRUE)
    }
    rows[[lab]] <- if (inherits(rd, "try-error")) application_fail(rd) else rd
  }
  rows
}

#' A competitor through wafc_competitor()
application_fit_competitor <- function(ctx, tuning, options, label, engine,
                                       own = list()) {
  restore_stream(ctx)
  f <- try(do.call(wafc_competitor,
                   c(list(engine, ctx[["train"]][["x"]],
                          ctx[["train"]][["u"]], ctx[["train"]][["y"]],
                          active = NULL, foldid = ctx[["foldid"]],
                          wavelet.table = ctx[["table"]]),
                     own, options)),
           silent = TRUE)
  if (inherits(f, "try-error")) {
    return(stats::setNames(list(application_fail(f)), label))
  }
  ex <- f[["extra"]]
  extra <- list()
  if (engine == "gam") {
    extra <- list(k = ex[["k"]], k_top = ex[["k.top"]], edf = ex[["edf"]],
                  k_table = ex[["k.table"]])
    kc <- try(gam_kcheck(f[["fit"]]), silent = TRUE)
    if (!inherits(kc, "try-error")) {
      extra[["kcheck"]] <- kc[["table"]]
      extra[["kcheck_min_p"]] <- kc[["kcheck_min_p"]]
      extra[["kcheck_n_low"]] <- kc[["kcheck_n_low"]]
    }
  }
  rd <- try(application_read(ctx, f, f[["time"]], extra), silent = TRUE)
  stats::setNames(list(if (inherits(rd, "try-error")) application_fail(rd)
                       else rd), label)
}

## The methods of the applications: for each one, the rows it produces and
## its fit, function(ctx, tuning, options) returning a named list of rows.
application_methods <- list(
  wafc = list(rows = c("wafc.cv1se", "wafc.cv"),
              fit = application_fit_wafc),
  gam.reml = list(rows = "gam.reml", fit = function(ctx, tuning, options)
    application_fit_competitor(ctx, tuning, options, "gam.reml", "gam",
                               list(k = tuning[["k"]]))),
  gam.gcv = list(rows = "gam.gcv", fit = function(ctx, tuning, options)
    application_fit_competitor(ctx, tuning, options, "gam.gcv", "gam",
                               list(k = tuning[["k"]]))),
  linear = list(rows = "linear", fit = function(ctx, tuning, options)
    application_fit_competitor(ctx, tuning, options, "linear", "linear"))
)

application_unit_path <- function(out, base, method, split) {
  file.path(out, "units", base, method, sprintf("split%02d.rds", split))
}

#' What determines a unit
application_settings <- function(cfg, base, split, method) {
  list(base = base, split = split, method = method,
       tuning = cfg[["tuning"]], options = cfg[["method_options"]][[method]],
       test_fraction = if (split == 0L) 0 else cfg[["test_fraction"]],
       grid_size = cfg[["grid_size"]],
       grid_quantiles = cfg[["grid_quantiles"]],
       seed = application_seed(cfg, base, split),
       data = unname(tools::md5sum(file.path(
         cfg[["data_dir"]], application_bases[[base]][["file"]]))))
}

#' The units of a run
#'
#' Every (application, split, method) of the configuration, restricted by
#' the optional arguments, which never change what a unit draws.
application_units <- function(cfg, bases = NULL, splits = NULL,
                              methods = NULL) {
  b <- cfg[["bases"]]
  if (!is.null(bases)) {
    bad <- setdiff(bases, b)
    if (length(bad) > 0L) stop("application(s) not in the configuration: ",
                               paste(bad, collapse = ", "), call. = FALSE)
    b <- intersect(b, bases)
  }
  s <- 0:cfg[["partitions"]]
  if (!is.null(splits)) {
    bad <- setdiff(splits, s)
    if (length(bad) > 0L) stop("split(s) outside 0, ..., ",
                               cfg[["partitions"]], ": ",
                               paste(bad, collapse = ", "), call. = FALSE)
    s <- intersect(s, splits)
  }
  m <- cfg[["methods"]]
  if (!is.null(methods)) {
    bad <- setdiff(methods, m)
    if (length(bad) > 0L) stop("method(s) not in the configuration: ",
                               paste(bad, collapse = ", "), call. = FALSE)
    m <- intersect(m, methods)
  }
  if (length(b) == 0L || length(s) == 0L || length(m) == 0L) {
    stop("no unit to run with these restrictions.", call. = FALSE)
  }
  expand.grid(base = b, split = s, method = m, stringsAsFactors = FALSE,
              KEEP.OUT.ATTRS = FALSE)
}

#' Fit one unit
application_unit <- function(cfg, data, split, method, table, code = NULL) {
  t0 <- elapsed()
  ctx <- application_context(cfg, data, split, table)
  opts <- cfg[["method_options"]][[method]]
  if (is.null(opts)) opts <- list()
  rows <- application_methods[[method]][["fit"]](ctx, cfg[["tuning"]], opts)
  list(rows = rows, sizes = ctx[["sizes"]], grid = ctx[["grid"]],
       settings = application_settings(cfg, data[["base"]], split, method),
       elapsed = elapsed() - t0, code = code)
}

#' Run (or resume) the units of the applications
#'
#' @param out The output directory.
#' @param workers Number of forked workers.
#' @param refresh If TRUE, cached units whose settings differ from the
#'   configuration are refitted; if FALSE, their presence stops the run.
#' @param ... Restrictions passed to application_units().
#' @return The data frame of units, with their status, invisibly.
run_application <- function(cfg, out, workers = 1L, refresh = FALSE, ...) {
  units <- application_units(cfg, ...)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  code <- code_provenance(cfg)
  record_session(cfg, out, code)
  units[["path"]] <- application_unit_path(out, units[["base"]],
                                           units[["method"]], units[["split"]])
  units[["status"]] <- "todo"
  for (i in which(file.exists(units[["path"]]))) {
    st <- application_settings(cfg, units[["base"]][i], units[["split"]][i],
                               units[["method"]][i])
    old <- readRDS(units[["path"]][i])[["settings"]]
    units[["status"]][i] <- if (identical(old, st)) "cached" else "stale"
  }
  nst <- sum(units[["status"]] == "stale")
  if (nst > 0L && !refresh) {
    stop(nst, " cached unit(s) were run with other settings (for instance ",
         units[["path"]][units[["status"]] == "stale"][1L], "). Rerun with ",
         "--refresh to refit them, or remove them.", call. = FALSE)
  }
  todo <- units[units[["status"]] != "cached", , drop = FALSE]
  ## the costliest first: the WAFC search, then GCV, which refits mgcv
  ## without discretization; the larger model first; the whole sample first
  cost <- c(wafc = 1, gam.gcv = 2, gam.reml = 3, linear = 4)
  size <- vapply(todo[["base"]], function(b) {
    length(application_bases[[b]][["xlab"]]) *
      length(application_bases[[b]][["ulab"]])
  }, 0)
  todo <- todo[order(cost[todo[["method"]]], -size, todo[["split"]]), ,
               drop = FALSE]
  message(sprintf("[application] %d unit(s), %d cached, %d to fit, %d worker(s)",
                  nrow(units), sum(units[["status"]] == "cached"), nrow(todo),
                  workers))
  if (nrow(todo) == 0L) return(invisible(units))
  ## the data and the wavelet table are prepared here, once, and inherited
  ## by the forked workers
  data <- lapply(stats::setNames(nm = unique(todo[["base"]])),
                 function(b) application_data(cfg, b))
  b <- cfg[["tuning"]][["basis"]]
  table <- WaveBased::wtable(family = b[["family"]],
                             filter.size = as.integer(b[["filter_size"]]),
                             prec.wavelet = as.integer(b[["prec_wavelet"]]),
                             check = FALSE)
  t0 <- elapsed()
  one <- function(i) {
    u <- todo[i, ]
    res <- try(application_unit(cfg, data[[u[["base"]]]], u[["split"]],
                                u[["method"]], table, code), silent = TRUE)
    if (inherits(res, "try-error")) {
      message(sprintf("[application] %s split %d %s: UNIT FAILED: %s",
                      u[["base"]], u[["split"]], u[["method"]],
                      conditionMessage(attr(res, "condition"))))
      return("failed")
    }
    save_unit(res, u[["path"]])
    errs <- unlist(lapply(res[["rows"]], `[[`, "error"))
    message(sprintf("[application] %s split %d %s: %.1f s%s", u[["base"]],
                    u[["split"]], u[["method"]], res[["elapsed"]],
                    if (length(errs) > 0L) paste0(" (row error: ", errs[1L],
                                                   ")") else ""))
    "done"
  }
  st <- if (workers > 1L) {
    parallel::mclapply(seq_len(nrow(todo)), one, mc.cores = workers,
                       mc.preschedule = FALSE)
  } else lapply(seq_len(nrow(todo)), one)
  st <- vapply(st, function(s) if (is.character(s)) s else "failed", "")
  units[["status"]][match(todo[["path"]], units[["path"]])] <- st
  message(sprintf("[application] %d unit(s) fitted, %d failed, %.1f min",
                  sum(st == "done"), sum(st == "failed"),
                  (elapsed() - t0) / 60))
  invisible(units)
}

#' The command line of scripts/03_application.R
#'
#' The options of the simulation scripts (R/cli.R) that apply here
#' (--config, --study, --out, --workers, --methods, --refresh, --list), with
#' the configuration config/application.yaml by default, and three of its
#' own:
#'
#'   --bases=a,b      only these applications
#'   --splits=0:20    only these splits (0 the whole sample, 1, 2, ... the
#'                    partitions); a range FROM:TO or a list
#'   --parts=fit,report  fit the units, write the tables and figures, or
#'                    both (the default)
application_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  key <- sub("=.*$", "", sub("^--", "", args))
  own <- key %in% c("bases", "splits", "parts")
  cli <- parse_cli(args[!own])
  for (k in c("cells", "sizes", "reps")) {
    if (!is.null(cli[[k]])) {
      stop("--", k, " is an option of the simulation study; the ",
           "applications take --bases and --splits.", call. = FALSE)
    }
  }
  if (!any(key == "config")) {
    cli[["config"]] <- file.path("config", "application.yaml")
  }
  val <- function(k) {
    a <- args[own & key == k]
    if (length(a) == 0L) NULL else sub("^[^=]*=", "", a[length(a)])
  }
  sp <- val("splits")
  if (!is.null(sp)) {
    sp <- if (grepl(":", sp, fixed = TRUE)) {
      r <- as.integer(strsplit(sp, ":", fixed = TRUE)[[1L]])
      seq(r[1L], r[2L])
    } else as.integer(strsplit(sp, ",", fixed = TRUE)[[1L]])
  }
  parts <- val("parts")
  parts <- if (is.null(parts)) c("fit", "report") else
    strsplit(parts, ",", fixed = TRUE)[[1L]]
  bad <- setdiff(parts, c("fit", "report"))
  if (length(bad) > 0L) stop("unknown part(s): ", paste(bad, collapse = ", "),
                             "; use fit and report.", call. = FALSE)
  c(cli, list(bases = if (is.null(val("bases"))) NULL else
                strsplit(val("bases"), ",", fixed = TRUE)[[1L]],
              splits = sp, parts = parts))
}
