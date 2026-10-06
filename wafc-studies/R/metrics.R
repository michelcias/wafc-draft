## R/metrics.R -- what is measured on every fit.
##
## Every method is asked the same questions on the test sample and on the
## grid: its prediction of the regression function, its functional
## coefficients beta_l(u), its additive components g_lm, and which blocks
## (l, m) it kept. One row per fit:
##
##   rmse_f      root mean squared error against the true regression
##               function on the test sample (the prediction error);
##   rmse_y      the same against the test responses;
##   mse_beta    mean over the test sample of the squared error of the p
##               functional coefficients, summed over l;
##   ise         integrated squared error of the components on the grid,
##               summed over the blocks, and split into ise_active and
##               ise_null (the blocks whose true component is zero);
##   n_true      active blocks kept, n_false zero blocks kept, out of
##               n_active and n_block; S-hat = S when n_true = n_active and
##               n_false = 0;
##   J, k        the resolution level, or the basis dimension of the gam
##               methods, the method chose; top, whether that choice is
##               the largest value of its grid;
##   lambda      the penalty level, nzero the nonzero penalized
##               coefficients after any threshold;
##   time        elapsed seconds of the fit, the search over the tuning
##               grids and any threshold included;
##   kcheck_*    for the gam methods, the basis dimension check of
##               mgcv::k.check() on the final fit: the smallest p-value
##               over the smooths, and the number of smooths with
##               k-index below 1 and p-value below 0.05;
##   error       the message of a fit that failed, whose row has no numbers.

#' Integrated squared error of each component on the grid
#'
#' Both the estimate and the truth are centred on the grid: every method
#' fixes the level of a component by its own convention, so only the shape
#' is comparable.
#'
#' @param ghat A p by q list of vectors on the grid, or NULL for a method
#'   without an additive decomposition.
#' @return A p by q matrix.
ise_components <- function(ghat, train, grid) {
  p <- train[["p"]]
  q <- train[["q"]]
  ise <- matrix(NA_real_, p, q)
  if (is.null(ghat)) return(ise)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      gm <- train[["g"]][[l, m]]
      gt <- if (is.null(gm)) numeric(nrow(grid)) else gm(grid[, m])
      gh <- ghat[[l, m]] - mean(ghat[[l, m]])
      gt <- gt - mean(gt)
      ise[l, m] <- mean((gh - gt)^2) * diff(range(grid[, m]))
    }
  }
  ise
}

## The empty row: the identification columns and NA everywhere else.
row_template <- function(ctx, method, fit) {
  data.frame(
    cell = ctx[["cell"]][["name"]], scenario = ctx[["cell"]][["scenario"]],
    n = ctx[["n"]], rep = ctx[["rep"]], method = method, fit = fit,
    rmse_f = NA_real_, rmse_y = NA_real_, mse_beta = NA_real_,
    ise = NA_real_, ise_active = NA_real_, ise_null = NA_real_,
    n_true = NA_integer_, n_false = NA_integer_,
    n_active = sum(ctx[["active"]]), n_block = length(ctx[["active"]]),
    J = NA_integer_, k = NA_integer_, top = NA, lambda = NA_real_,
    nzero = NA_integer_, sigma = ctx[["train"]][["sigma"]], time = NA_real_,
    kcheck_min_p = NA_real_, kcheck_n_low = NA_integer_,
    error = NA_character_, stringsAsFactors = FALSE)
}

#' One row of the table, and the components behind it
#'
#' @param reading A list with `f_test` (prediction on the test sample),
#'   `beta_test` (the n_test by p functional coefficients there), `ghat`
#'   (the components on the grid, or NULL) and `blocks` (the p by q
#'   logical matrix of blocks kept, or NULL).
#' @param extra A list with any of `J`, `k`, `top`, `lambda`, `nzero`,
#'   `kcheck_min_p`, `kcheck_n_low`.
#' @return A list with the `row` and the per-block table `components`.
make_row <- function(ctx, method, fit, reading, time, extra = list()) {
  test <- ctx[["test"]]
  active <- ctx[["active"]]
  row <- row_template(ctx, method, fit)
  f <- reading[["f_test"]]
  row[["rmse_f"]] <- sqrt(mean((f - test[["f"]])^2))
  row[["rmse_y"]] <- sqrt(mean((f - test[["y"]])^2))
  row[["mse_beta"]] <- mean(rowSums((reading[["beta_test"]] -
                                       test[["beta"]])^2))
  ise <- ise_components(reading[["ghat"]], ctx[["train"]], ctx[["grid"]])
  row[["ise"]] <- sum(ise)
  row[["ise_active"]] <- sum(ise[active])
  row[["ise_null"]] <- sum(ise[!active])
  sel <- reading[["blocks"]]
  if (!is.null(sel)) {
    row[["n_true"]] <- sum(sel[active])
    row[["n_false"]] <- sum(sel[!active])
  }
  for (v in names(extra)) {
    if (!is.null(extra[[v]]) && length(extra[[v]]) == 1L) {
      row[[v]] <- switch(typeof(row[[v]]),
                         integer = as.integer(extra[[v]]),
                         logical = as.logical(extra[[v]]),
                         as.double(extra[[v]]))
    }
  }
  row[["time"]] <- time
  comp <- data.frame(
    cell = row[["cell"]], n = row[["n"]], rep = row[["rep"]],
    method = method, l = as.vector(row(ise)), m = as.vector(col(ise)),
    active = as.vector(active), ise = as.vector(ise),
    kept = if (is.null(sel)) NA else as.vector(sel),
    stringsAsFactors = FALSE)
  list(row = row, components = comp)
}

## The row of a fit that failed: its message, and no numbers.
fail_row <- function(ctx, method, fit, err) {
  row <- row_template(ctx, method, fit)
  msg <- if (inherits(err, "try-error")) {
    conditionMessage(attr(err, "condition"))
  } else as.character(err)
  row[["error"]] <- gsub("\\s+", " ", msg)
  list(row = row, components = NULL)
}

#' The basis dimension check of a gam fit
#'
#' mgcv::k.check() on the final fit: for every smooth, its basis dimension
#' k', its effective degrees of freedom, the k-index and the p-value of the
#' randomization test of Wood (2017, section 5.9). A k-index below 1 with a
#' small p-value says the residuals keep pattern along the covariate of the
#' smooth, the sign of a basis too small.
#'
#' @return A list with the per-smooth `table` and the two summaries of the
#'   row.
gam_kcheck <- function(fit) {
  kc <- mgcv::k.check(fit)
  tab <- data.frame(smooth = rownames(kc), k_prime = kc[, 1L],
                    edf = kc[, 2L], k_index = kc[, 3L], p_value = kc[, 4L],
                    stringsAsFactors = FALSE, row.names = NULL)
  list(table = tab,
       kcheck_min_p = min(tab[["p_value"]]),
       kcheck_n_low = sum(tab[["k_index"]] < 1 & tab[["p_value"]] < 0.05))
}

#' The fitted components on the grid, with the truth, for the figures
curves_table <- function(ctx, method, ghat) {
  if (is.null(ghat)) return(NULL)
  train <- ctx[["train"]]
  grid <- ctx[["grid"]]
  out <- list()
  for (l in seq_len(train[["p"]])) {
    for (m in seq_len(train[["q"]])) {
      gm <- train[["g"]][[l, m]]
      gt <- if (is.null(gm)) numeric(nrow(grid)) else gm(grid[, m])
      out[[length(out) + 1L]] <- data.frame(
        cell = ctx[["cell"]][["name"]], n = ctx[["n"]], rep = ctx[["rep"]],
        method = method, l = l, m = m, u = grid[, m],
        ghat = ghat[[l, m]] - mean(ghat[[l, m]]), truth = gt - mean(gt),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}
