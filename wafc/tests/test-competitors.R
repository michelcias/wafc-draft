## wafc/tests/test-competitors.R -- tests of the competitors of step E2.4
## and of the penalty rule the pilot adds. Run from the root of the
## repository:
##
##     Rscript -e 'testthat::test_dir("wafc/tests")'
##
## Four tests carry the step. The first is the one docs/instrucoes.md,
## section 6, makes the first of every estimator: with the truth inside the
## model and no noise, the method recovers it. The second checks that every
## competitor answers the three questions the pilot asks in the same
## coordinates as wafc(), which is what makes the table a comparison and
## not a list. The third checks the block LASSO of Klopp and Pensky against
## a group vector built by hand, since the whole point of that column is to
## be their penalty and not another one. The fourth checks the quantile
## universal threshold against its own definition: it is the quantile of
## the smallest penalty level that kills the penalized part under the null,
## so simulating from the null and solving for that level has to reproduce
## it.

library(testthat)

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

set.seed(20260919)
n <- 250L
p <- 3L
q <- 2L
dgp <- simulate_wafc(n, p = p, q = q, scenario = "smooth", seed = 20260919L,
                     snr = 6)
x0 <- dgp[["x"]]
u0 <- dgp[["u"]]
y0 <- dgp[["y"]]
active0 <- nzchar(dgp[["structure"]])
folds <- rep_len(1:5, n)
grid0 <- cbind(seq(min(u0[, 1L]), max(u0[, 1L]), length.out = 64),
               seq(min(u0[, 2L]), max(u0[, 2L]), length.out = 64))

has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

## ---------------------------------------------------------------------------
## Exact recovery
## ---------------------------------------------------------------------------

test_that("every competitor recovers a model of its own form without noise", {
  ## A model that is inside the span of all the sieves at once: the
  ## coefficients vary with u through a single smooth function, so the
  ## spline, the wavelet and the tree method can all reproduce it, and
  ## without noise the recovery is a question about the plumbing and not
  ## about the basis.
  gg <- function(v) sin(2 * pi * v)
  beta <- cbind(1 + gg(u0[, 1L]), 2 + numeric(n), -1.5 + numeric(n))
  f <- as.numeric(rowSums(x0 * beta))
  for (mth in c("gam", "bsgl", "klopp", "aspline", "linear")) {
    if (mth %in% c("bsgl", "klopp") && !has("grpreg")) next
    if (mth == "gam" && !has("mgcv")) next
    fit <- wafc_competitor(mth, x0, u0, f, foldid = folds, J = 4L)
    fh <- predict(fit, x0, u0)
    ## the linear fit cannot reproduce the varying part and is the control
    ## that the tolerance is not vacuous
    tol <- if (mth == "linear") 10 else 0.05 * stats::sd(f)
    if (mth == "linear") {
      expect_gt(sqrt(mean((fh - f)^2)), 0.1 * stats::sd(f))
    } else {
      expect_lt(sqrt(mean((fh - f)^2)), tol)
    }
  }
})

test_that("the oracle of the structure recovers a sparse truth in the basis", {
  d0 <- wafc_design(x0, u0, J = 3L, rescale = FALSE)
  theta <- numeric(d0[["nvars"]])
  theta[d0[["unpenalized"]]] <- c(1, 2, -1.5)
  theta[c(d0[["blocks"]][["x1:u1"]][c(1L, 4L)],
          d0[["blocks"]][["x2:u1"]][2L])] <- c(1.5, -1, 0.8)
  y_exact <- as.numeric(d0[["Z"]] %*% theta)
  act <- matrix(FALSE, p, q)
  act[1L, 1L] <- TRUE
  act[2L, 1L] <- TRUE
  fit <- wafc_competitor("oracle", x0, u0, y_exact, active = act,
                         foldid = folds, J = 3L, rescale = FALSE,
                         lambda = c(1e-1, 1e-3, 1e-6, 1e-9))
  expect_lt(sqrt(mean((predict(fit, x0, u0) - y_exact)^2)), 1e-4)
  ## it cannot have used a block it was told is zero
  expect_false(any(fit[["blocks"]][!act]))
})

## ---------------------------------------------------------------------------
## The common interface
## ---------------------------------------------------------------------------

test_that("every competitor answers in the coordinates of the model", {
  mths <- c("gam", "bsgl", "klopp", "aspline", "linear", "oracle")
  if (has("VCBART")) mths <- c(mths, "vcbart")
  for (mth in mths) {
    if (mth %in% c("bsgl", "klopp") && !has("grpreg")) next
    if (mth == "gam" && !has("mgcv")) next
    fit <- wafc_competitor(mth, x0, u0, y0, active = active0, foldid = folds,
                           J = 3L, burn = 100L, nd = 100L)
    expect_s3_class(fit, "wafc_competitor")
    ## beta(u) has one column per linear covariate, in the order of the design
    b <- fit[["beta"]](u0)
    expect_equal(dim(b), c(n, p))
    expect_equal(colnames(b), colnames(x0))
    ## the prediction is the model's own: rowSums(x * beta(u))
    expect_equal(predict(fit, x0, u0), as.numeric(rowSums(x0 * b)),
                 tolerance = 1e-8)
    ## predict on the training data agrees with the stored fitted values
    expect_equal(predict(fit, x0, u0), fit[["fitted"]], tolerance = 1e-6)
    ## the components, when the method has them, are centred on the grid and
    ## add up to beta minus its level
    g <- wafc_grid_components(fit, grid0)
    if (is.null(g)) {
      expect_true(mth == "vcbart")
      next
    }
    expect_equal(dim(g), c(p, q))
    expect_true(all(vapply(g, function(v) abs(mean(v)) < 1e-8, TRUE)))
    gb <- fit[["beta"]](grid0)
    for (l in seq_len(p)) {
      s <- Reduce(`+`, g[l, ])
      expect_equal(gb[, l] - mean(gb[, l]), s - mean(s), tolerance = 1e-6)
    }
    if (!is.null(fit[["blocks"]])) {
      expect_equal(dim(fit[["blocks"]]), c(p, q))
      expect_type(fit[["blocks"]], "logical")
    }
  }
})

test_that("a competitor ignores an argument meant for another one", {
  ## the pilot calls every method with one list of arguments
  expect_s3_class(wafc_competitor("linear", x0, u0, y0, foldid = folds,
                                  J = 3L, nfolds = 5L),
                  "wafc_competitor")
  expect_error(wafc_competitor("oracle", x0, u0, y0), "needs 'active'")
})

## ---------------------------------------------------------------------------
## The block LASSO of Klopp and Pensky
## ---------------------------------------------------------------------------

test_that("the chunks of the block LASSO are the ones of Klopp and Pensky", {
  skip_if_not(has("grpreg"))
  des <- wafc_design(x0, u0, J = 4L)
  bs <- 3L
  grp <- wafc_kp_groups(des, bs, penalize.levels = TRUE)
  ## every column belongs to exactly one group, and the p level terms are a
  ## group of their own, which is their norm (3.1) and not decision D3
  expect_length(grp, des[["nvars"]])
  expect_true(all(grp > 0L))
  expect_equal(length(unique(grp[des[["unpenalized"]]])), 1L)
  ## a chunk never straddles two resolution levels: inside a block the
  ## columns are ordered by increasing j (D12), so the sizes are
  ## min(2^j, bs) repeated as the level is cut
  idx <- des[["blocks"]][["x1:u1"]]
  sizes <- as.integer(table(grp[idx])[as.character(unique(grp[idx]))])
  want <- unlist(lapply(0:3, function(j) {
    nj <- 2^j
    diff(c(seq(0L, nj - 1L, by = bs), nj))
  }))
  expect_equal(sizes, as.integer(want))
  ## and the unpenalized version leaves the level terms out
  g0 <- wafc_kp_groups(des, bs, penalize.levels = FALSE)
  expect_true(all(g0[des[["unpenalized"]]] == 0L))
  expect_equal(max(g0), max(grp) - 1L)
})

test_that("the block LASSO zeroes a whole chunk at a time", {
  skip_if_not(has("grpreg"))
  fit <- wafc_competitor("klopp", x0, u0, y0, foldid = folds, J = 3L)
  des <- fit[["design"]]
  b <- as.numeric(stats::coef(fit[["fit"]]))[-1L]
  grp <- wafc_kp_groups(des, fit[["extra"]][["block.size"]],
                        fit[["extra"]][["penalize.levels"]])
  ## Only the wavelet chunks: the group of the level terms contains the
  ## constant covariate, which grpreg drops as a column of zero variance and
  ## leaves at zero whatever the rest of the group does. Columns that are
  ## identically zero in the design (they appear from J = 8 with the default
  ## margin) are out for the same reason.
  alive <- as.numeric(Matrix::colSums(abs(des[["Z"]]))) > 0
  pen <- setdiff(seq_len(des[["nvars"]]), des[["unpenalized"]])
  for (g in unique(grp[pen])) {
    v <- b[grp == g & alive]
    if (length(v) == 0L) next
    expect_true(all(v == 0) || all(v != 0))
  }
})

## ---------------------------------------------------------------------------
## The quantile universal threshold
## ---------------------------------------------------------------------------

test_that("the QUT is the quantile of the penalty level that kills the fit", {
  des <- wafc_design(x0, u0, J = 3L)
  lam <- wafc_lambda_qut(des, y0, alpha = 0.05, nsim = 400L, seed = 11L)
  expect_true(is.finite(lam) && lam > 0)
  expect_equal(lam, wafc_lambda_qut(des, y0, alpha = 0.05, nsim = 400L,
                                    seed = 11L))
  ## the rule is free of sigma: multiplying the response by a constant (on
  ## the residual of the unpenalized part, which is what carries the scale)
  ## multiplies the penalty level by the same constant
  unp <- des[["unpenalized"]]
  W <- as.matrix(des[["Z"]][, unp, drop = FALSE])
  r0 <- as.numeric(y0 - W %*% qr.coef(qr(W), y0))
  y2 <- y0 + 4 * r0
  expect_equal(wafc_lambda_qut(des, y2, nsim = 400L, seed = 11L) / lam, 5,
               tolerance = 1e-8)
  ## and it is the quantile of a level that really kills the penalized part:
  ## fitting a null sample at the level the rule returns for that sample
  ## leaves no non-zero wavelet coefficient in at least 1 - alpha of the
  ## draws, by construction of the statistic
  kill <- logical(40L)
  for (i in seq_along(kill)) {
    yy <- as.numeric(W %*% c(1, 2, -1.5)) + stats::rnorm(n, sd = 0.5)
    li <- wafc_lambda_qut(des, yy, alpha = 0.05, nsim = 200L, seed = 100L + i)
    fi <- wafc(design = des, y = yy,
               lambda = wafc_path_to(li, des, yy))
    kill[i] <- fi[["nzero"]][length(fi[["lambda"]])] == 0L
  }
  expect_gt(mean(kill), 0.85)
})

test_that("the QUT is more conservative than the cross-validated penalty", {
  ## the property the pilot reports: the rule is calibrated on the null, so
  ## it penalizes more than the rule calibrated on prediction
  cv <- cv.wafc(x0, u0, y0, J = 3L, foldid = folds)
  lam <- wafc_lambda_qut(cv[["wafc.fit"]][["design"]], y0, nsim = 400L,
                         seed = 11L)
  expect_gt(lam, cv[["lambda.min"]])
})

## ---------------------------------------------------------------------------
## The metric the pilot reads
## ---------------------------------------------------------------------------

test_that("wafc_grid_components gives the same components as wafc_functions", {
  fit <- wafc(x0, u0, y0, J = 3L)
  s <- fit[["lambda"]][40L]
  fn <- wafc_functions(fit, s = s, grid = grid0)
  g <- wafc_grid_components(fit, grid0, s = s)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      a <- fn[["g"]][[l, m]]
      expect_equal(g[[l, m]], a - mean(a), tolerance = 1e-10)
    }
  }
})

test_that("the components of a fit reproduce a truth that is in the basis", {
  ## the integrated squared error the pilot reports is zero when the fit is
  ## the truth, which is what says the grid, the centring and the range are
  ## the same on both sides
  d0 <- wafc_design(x0, u0, J = 3L, rescale = FALSE)
  theta <- numeric(d0[["nvars"]])
  theta[d0[["unpenalized"]]] <- c(1, 2, -1.5)
  theta[d0[["blocks"]][["x1:u1"]][2L]] <- 1.5
  y_exact <- as.numeric(d0[["Z"]] %*% theta)
  fit <- wafc(design = d0, y = y_exact, lambda = c(1e-3, 1e-8))
  g <- wafc_grid_components(fit, grid0, s = 1e-8)
  d_grid <- wafc_design(matrix(1, nrow(grid0), p), grid0, spec = d0)
  truth <- as.numeric(d_grid[["Z"]][, d0[["blocks"]][["x1:u1"]][2L]]) * 1.5
  expect_equal(g[[1L, 1L]], truth - mean(truth), tolerance = 1e-4)
  expect_lt(max(abs(unlist(g[-1L]))), 1e-4)
})
