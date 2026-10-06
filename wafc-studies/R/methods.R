## R/methods.R -- the methods of the study and how each one is read.
##
## Every method is fitted on the training sample of a replicate, with the
## folds of the replicate and the tuning grids of its sample size, and read
## on the test sample and on the grid (R/metrics.R). A method is a fit; two
## of them give two rows each, because two estimators come out of one
## search:
##
##   wafc      cv.wafc() with the block lasso, (J, lambda) chosen by
##             cross-validation over the grid of J, followed by the
##             threshold on the norms of the blocks: row "wafc.cv1se" with
##             the threshold at one standard error (the estimator, the
##             default of cv.wafc()), row "wafc.cv" with the threshold of
##             the smallest cross-validated error. The folds are refitted
##             once for the two thresholds.
##   lasso     the same search with the coordinatewise lasso: row "lasso",
##             the fit at the selected pair, and row "lasso.cv1se", the
##             same threshold as the estimator.
##   gam.reml  mgcv with one smooth s(u_m, by = x_l) per block, the
##   gam.gcv   smoothing parameters and the basis dimension, common to the
##             smooths, chosen from the grid k by REML or by GCV.
##   bsgl      cubic B-splines with a group lasso by block, the dimension
##             2^J chosen by cross-validation over the grid of J.
##   aspline   a spline with knots chosen adaptively per block.
##   klopp     the block lasso of Klopp and Pensky (2015) on the wavelet
##             design: chunks of about log n coefficients, the level terms
##             penalized as blocks of their own.
##   oracle    the block lasso of wafc on the wavelet design restricted to
##             the blocks that are truly active, (J, lambda) chosen by
##             cross-validation over the grid of J, without threshold; a
##             reference, not a competitor: the distance from wafc to it is
##             the price of not knowing the structure.
##   linear    least squares on the linear covariates alone (constant
##             coefficients); the oracle of the null cell.
##   vcbart    Bayesian additive regression trees for varying coefficients
##             (Deshpande et al., 2026), with no additive restriction.
##
## The time of a row is the time of its fit, the search over the grids
## included; a thresholded row adds the fold fits and the threshold.

elapsed <- function() proc.time()[["elapsed"]]

## The levels, one per modulator, a wavelet design asked for level J is
## built at: wafc_design() caps the level of a modulator at
## floor(log2(d)), d its distinct points on the circle after the rescaling
## (at n = 250 with continuous modulators, d = 249 and J = 8 is built at 7).
## Two levels with the same capped design are the same candidate, so "the
## top of the grid" is read on the capped levels.
capped_levels <- function(ctx, J) {
  u <- wafc_rescale(ctx[["train"]][["u"]], eps = 0)[["u"]] %% 1
  d <- apply(u, 2L, function(v) length(unique(v)))
  pmin(as.integer(J), as.integer(floor(log2(pmax(d, 1L)))))
}

at_top <- function(ctx, J, grid) {
  if (is.null(J) || length(J) != 1L || is.na(J)) return(NA)
  all(capped_levels(ctx, J) == capped_levels(ctx, max(grid)))
}

## The reading of a fit that has the interface of wafc_competitor() (every
## competitor, and the thresholded WAFC fits).
read_competitor <- function(ctx, f) {
  list(f_test = as.numeric(predict(f, ctx[["test"]][["x"]],
                                   ctx[["test"]][["u"]])),
       beta_test = f[["beta"]](ctx[["test"]][["u"]]),
       ghat = wafc_grid_components(f, ctx[["grid"]]),
       blocks = f[["blocks"]])
}

## The reading of a "cv.wafc" object without threshold, at lambda.min of
## the selected J.
read_cv <- function(ctx, cv) {
  list(f_test = as.numeric(predict(cv, ctx[["test"]][["x"]],
                                   ctx[["test"]][["u"]])),
       beta_test = predict(cv, newu = ctx[["test"]][["u"]], type = "beta"),
       ghat = wafc_grid_components(cv, ctx[["grid"]]),
       blocks = wafc_blocks(cv)[["nonzero"]] > 0L)
}

## Appends one row (and its components and curves) to a list of outputs.
add_row <- function(out, ctx, method, fit, reading, time, extra = list()) {
  r <- make_row(ctx, method, fit, reading, time, extra)
  out[["rows"]][[length(out[["rows"]]) + 1L]] <- r[["row"]]
  out[["side"]][["components"]] <- rbind(out[["side"]][["components"]],
                                         r[["components"]])
  out[["ghat"]][[method]] <- reading[["ghat"]]
  out
}

add_fail <- function(out, ctx, method, fit, err) {
  r <- fail_row(ctx, method, fit, err)
  out[["rows"]][[length(out[["rows"]]) + 1L]] <- r[["row"]]
  out
}

add_side <- function(out, name, tab) {
  if (!is.null(tab)) out[["side"]][[name]] <- rbind(out[["side"]][[name]],
                                                    tab)
  out
}

## The identification columns of a side table.
side_id <- function(ctx, method) {
  data.frame(cell = ctx[["cell"]][["name"]], n = ctx[["n"]],
             rep = ctx[["rep"]], method = method, stringsAsFactors = FALSE)
}

## The WAFC search and its thresholds, for the block lasso ("wafc") and the
## coordinatewise lasso ("lasso"). `rules` maps a row label to a threshold
## rule, with NA for the fit without threshold.
fit_wavelet <- function(ctx, tuning, options, fit, penalty, rules) {
  out <- list(rows = list(), side = list(), ghat = list())
  t0 <- elapsed()
  cv <- try(do.call(cv.wafc, c(list(ctx[["train"]][["x"]],
                                     ctx[["train"]][["u"]],
                                     ctx[["train"]][["y"]],
                                     J = tuning[["J"]], penalty = penalty,
                                     foldid = ctx[["foldid"]],
                                     wavelet.table = ctx[["table"]],
                                     threshold = "none"), options)),
            silent = TRUE)
  if (inherits(cv, "try-error")) {
    for (lab in names(rules)) out <- add_fail(out, ctx, lab, fit, cv)
    return(out)
  }
  t_search <- elapsed() - t0
  top <- at_top(ctx, cv[["J.min"]], tuning[["J"]])
  conv <- try(wafc_cv_convergence(cv), silent = TRUE)
  if (is.data.frame(conv)) {
    out <- add_side(out, "convergence", cbind(side_id(ctx, fit), conv))
  }
  ff <- NULL
  t_folds <- 0
  if (any(!is.na(rules))) {
    t0 <- elapsed()
    ff <- try(wafc_threshold_folds(cv), silent = TRUE)
    t_folds <- elapsed() - t0
  }
  for (lab in names(rules)) {
    rule <- rules[[lab]]
    if (is.na(rule)) {
      rd <- try(read_cv(ctx, cv), silent = TRUE)
      if (inherits(rd, "try-error")) {
        out <- add_fail(out, ctx, lab, fit, rd)
        next
      }
      out <- add_row(out, ctx, lab, fit, rd, t_search,
                     list(J = cv[["J.min"]], top = top,
                          lambda = cv[["lambda.min"]],
                          nzero = sum(wafc_blocks(cv)[["nonzero"]])))
      next
    }
    th <- if (inherits(ff, "try-error")) ff else
      try(wafc_threshold(cv, rule = rule, fold.fits = ff), silent = TRUE)
    if (inherits(th, "try-error")) {
      out <- add_fail(out, ctx, lab, fit, th)
      next
    }
    ex <- th[["extra"]]
    out <- add_row(out, ctx, lab, fit, read_competitor(ctx, th),
                   t_search + t_folds + th[["time"]],
                   list(J = ex[["J"]], top = top, lambda = ex[["lambda"]],
                        nzero = ex[["nzero"]]))
    out <- add_side(out, "threshold",
                    cbind(side_id(ctx, lab),
                          data.frame(t = ex[["t"]], c = ex[["c"]],
                                     blocks_kept = sum(ex[["kept"]]))))
  }
  out
}

## A competitor through wafc_competitor(), with the arguments every
## competitor receives (the ones a fitter does not take are dropped there).
fit_competitor <- function(ctx, tuning, options, fit, engine, own = list()) {
  out <- list(rows = list(), side = list(), ghat = list())
  restore_stream(ctx)
  f <- try(do.call(wafc_competitor,
                   c(list(engine, ctx[["train"]][["x"]],
                          ctx[["train"]][["u"]], ctx[["train"]][["y"]],
                          active = ctx[["active"]],
                          foldid = ctx[["foldid"]],
                          wavelet.table = ctx[["table"]]),
                     own, options)),
           silent = TRUE)
  if (inherits(f, "try-error")) return(add_fail(out, ctx, fit, fit, f))
  ex <- f[["extra"]]
  extra <- list(J = ex[["J"]], lambda = ex[["lambda"]], nzero = ex[["nzero"]])
  if (engine %in% c("klopp", "oracle")) {
    extra[["top"]] <- at_top(ctx, ex[["J"]], tuning[["J"]])
  }
  if (engine == "bsgl" && !is.null(ex[["df"]])) {
    ## the dimensions bsgl keeps: at most the distinct values minus one
    d <- min(apply(ctx[["train"]][["u"]], 2L, function(v) length(unique(v))))
    df <- 2L^tuning[["J"]]
    adm <- df[df <= d - 1L]
    extra[["top"]] <- ex[["df"]] == if (length(adm)) max(adm) else min(df)
  }
  if (engine == "gam") {
    kt <- ex[["k.table"]]
    if (!is.null(kt)) {
      kt[["chosen"]] <- kt[["k.used"]] == paste(ex[["k"]], collapse = ",")
      kt[["top"]] <- kt[["chosen"]] & ex[["k.top"]]
      extra[["k"]] <- kt[["k"]][kt[["chosen"]]][1L]
      extra[["top"]] <- ex[["k.top"]]
      out <- add_side(out, "gam_k", cbind(side_id(ctx, fit), kt))
    }
    kc <- try(gam_kcheck(f[["fit"]]), silent = TRUE)
    if (!inherits(kc, "try-error")) {
      extra[["kcheck_min_p"]] <- kc[["kcheck_min_p"]]
      extra[["kcheck_n_low"]] <- kc[["kcheck_n_low"]]
      out <- add_side(out, "kcheck", cbind(side_id(ctx, fit), kc[["table"]]))
    }
  }
  if (is.data.frame(ex[["conv"]])) {
    out <- add_side(out, "convergence", cbind(side_id(ctx, fit), ex[["conv"]]))
  }
  rd <- try(read_competitor(ctx, f), silent = TRUE)
  if (inherits(rd, "try-error")) return(add_fail(out, ctx, fit, fit, rd))
  add_row(out, ctx, fit, fit, rd, f[["time"]], extra)
}

## The registry: for each method, the row labels it produces and its fit.
## Every fit is function(ctx, tuning, options) and returns a list with
## `rows` (a list of one-row data frames), `side` (named side tables) and
## `ghat` (the components on the grid of each row, for the curves).
study_methods <- list(
  wafc = list(
    rows = c("wafc.cv1se", "wafc.cv"),
    fit = function(ctx, tuning, options) {
      restore_stream(ctx)
      fit_wavelet(ctx, tuning, options, "wafc", "block",
                  c(wafc.cv1se = "cv1se", wafc.cv = "cv"))
    }),
  lasso = list(
    rows = c("lasso", "lasso.cv1se"),
    fit = function(ctx, tuning, options) {
      restore_stream(ctx)
      fit_wavelet(ctx, tuning, options, "lasso", "lasso",
                  c(lasso = NA_character_, lasso.cv1se = "cv1se"))
    }),
  gam.reml = list(
    rows = "gam.reml",
    fit = function(ctx, tuning, options) {
      fit_competitor(ctx, tuning, options, "gam.reml", "gam",
                     list(k = tuning[["k"]]))
    }),
  gam.gcv = list(
    rows = "gam.gcv",
    fit = function(ctx, tuning, options) {
      fit_competitor(ctx, tuning, options, "gam.gcv", "gam",
                     list(k = tuning[["k"]]))
    }),
  bsgl = list(
    rows = "bsgl",
    fit = function(ctx, tuning, options) {
      fit_competitor(ctx, tuning, options, "bsgl", "bsgl",
                     list(df = 2L^tuning[["J"]]))
    }),
  aspline = list(
    rows = "aspline",
    fit = function(ctx, tuning, options) {
      fit_competitor(ctx, tuning, options, "aspline", "aspline")
    }),
  klopp = list(
    rows = "klopp",
    fit = function(ctx, tuning, options) {
      fit_competitor(ctx, tuning, options, "klopp", "klopp",
                     list(J = tuning[["J"]]))
    }),
  oracle = list(
    rows = "oracle",
    fit = function(ctx, tuning, options) {
      fit_competitor(ctx, tuning, options, "oracle", "oracle",
                     list(J = tuning[["J"]]))
    }),
  linear = list(
    rows = "linear",
    fit = function(ctx, tuning, options) {
      fit_competitor(ctx, tuning, options, "linear", "linear")
    }),
  vcbart = list(
    rows = "vcbart",
    fit = function(ctx, tuning, options) {
      fit_competitor(ctx, tuning, options, "vcbart", "vcbart")
    })
)

## Every row label of the study, in the order of the tables.
row_labels <- function(methods = names(study_methods)) {
  unlist(lapply(methods, function(m) study_methods[[m]][["rows"]]),
         use.names = FALSE)
}
