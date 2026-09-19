## wafc/tests/test-fit.R -- tests of the estimator and of the
## reconstruction of step E2.2. Run from the root of the repository:
##
##     Rscript -e 'testthat::test_dir("wafc/tests")'
##
## The first test is again the exact recovery of the plan, now through
## wafc(): with the components in the basis (theta* known) and no noise, the
## fit with lambda -> 0 returns theta*. The other two tests that carry the
## step are the KKT one, which is what establishes the scale of lambda
## instead of trusting penalty.factor, and the one on the level of the
## constant covariate, which the engine reads into its own intercept.

library(testthat)

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

has_sparsegl <- requireNamespace("sparsegl", quietly = TRUE)
skip_no_sparsegl <- function() {
  testthat::skip_if_not(has_sparsegl, "sparsegl is not installed")
}

set.seed(20260919)
n <- 300L
p <- 3L
q <- 2L
J <- 3L
NJ <- 2L^J - 1L
dgp <- simulate_wafc(n, p = p, q = q, scenario = "smooth", seed = 20260919L,
                     snr = 4)
x0 <- dgp[["x"]]
u0 <- dgp[["u"]]
y0 <- dgp[["y"]]

## A design and a sparse truth in the basis, used by the recovery tests.
d0 <- wafc_design(x0, u0, J = J, rescale = FALSE)
theta <- numeric(d0[["nvars"]])
theta[d0[["unpenalized"]]] <- c(1, 2, -1.5)
active <- c(d0[["blocks"]][["x1:u1"]][c(1L, 4L)],
            d0[["blocks"]][["x2:u2"]][c(2L, 7L)],
            d0[["blocks"]][["x3:u1"]][3L])
theta[active] <- c(1.5, -1, 0.8, 1.2, -0.7)
y_exact <- as.numeric(d0[["Z"]] %*% theta)

## ---------------------------------------------------------------------------
## The estimator
## ---------------------------------------------------------------------------

test_that("with theta* in the basis and no noise, wafc with lambda -> 0 recovers theta*", {
  fit <- wafc(design = d0, y = y_exact,
              lambda = c(0.1, 0.01, 1e-4, 1e-8))
  expect_s3_class(fit, "wafc")
  cf <- coef(fit, s = 1e-8)
  expect_equal(unname(cf["(Intercept)", 1L]), 0)
  expect_equal(max(abs(unname(cf[-1L, 1L]) - theta)), 0, tolerance = 1e-4)
  ## the level of the constant covariate comes from the intercept of the
  ## engine, which leaves its own coefficient at zero
  expect_equal(fit[["carrier"]][["index"]], 1L)
  expect_equal(fit[["carrier"]][["value"]], 1)
  expect_equal(unname(as.matrix(fit[["beta"]])[1L, ]),
               rep(0, length(fit[["lambda"]])))
  expect_equal(unname(fit[["cc"]][, ncol(fit[["cc"]])]), theta[1:3],
               tolerance = 1e-4)
  ## at a penalty of order one the wavelet part is shrunk and the levels
  ## are not
  expect_lt(fit[["nzero"]][1L], length(active) + 5L)
  expect_true(all(fit[["cc"]][, 1L] != 0))
})

test_that("the sparse group LASSO also recovers theta* as lambda -> 0", {
  skip_no_sparsegl()
  fit <- wafc(design = d0, y = y_exact, penalty = "sglasso",
              lambda = c(0.1, 0.01, 1e-4, 1e-8), thresh = 1e-14)
  cf <- coef(fit, s = 1e-8)
  expect_equal(max(abs(unname(cf[-1L, 1L]) - theta)), 0, tolerance = 1e-4)
  expect_equal(fit[["asparse"]], 0.05)
  ## the group structure is the one of D12: one group per block, plus the
  ## group of the unpenalized level terms, which carries penalty factor zero
  g <- fit[["group"]]
  expect_equal(g[["group"]][d0[["unpenalized"]]], rep(1L, p))
  expect_equal(unique(g[["group"]][d0[["blocks"]][["x2:u1"]]]), 4L)
  expect_equal(g[["pf_group"]], c(0, rep(sqrt(NJ), p * q)))
  expect_equal(g[["pf_sparse"]][d0[["unpenalized"]]], rep(0, p))
  expect_equal(g[["names"]], c("(levels)", names(d0[["blocks"]])))
})

test_that("the KKT conditions hold along the whole path, and fix the scale of lambda", {
  fit <- wafc(x0, u0, y0, J = J)
  k <- wafc_kkt(fit)
  expect_equal(nrow(k), length(fit[["lambda"]]))
  expect_true(all(k[["ok"]]))
  expect_lte(max(k[["subgradient"]]), 1 + 1e-6)

  ## the scale itself: glmnet rescales penalty.factor to sum to nvars, so
  ## the penalty it applies to a wavelet coefficient is its own lambda
  ## times nvars/npen. On the active set the gradient equals the effective
  ## lambda, and not the one glmnet reports.
  expect_equal(fit[["lambda.factor"]], fit[["nvars"]] / fit[["npen"]])
  kk <- 40L
  Z <- as.matrix(fit[["design"]][["Z"]])
  b <- as.numeric(fit[["beta"]][, kk])
  g <- as.numeric(crossprod(Z, y0 - fit[["a0"]][kk] - Z %*% b)) / n
  act <- which(b != 0 & fit[["design"]][["penalty.factor"]] > 0)
  expect_gt(length(act), 5L)
  expect_equal(mean(abs(g[act])) / fit[["lambda"]][kk], 1, tolerance = 1e-4)
  expect_equal(mean(abs(g[act])) / fit[["fit"]][["lambda"]][kk],
               fit[["nvars"]] / fit[["npen"]], tolerance = 1e-4)
  expect_false(isTRUE(all.equal(fit[["lambda"]][kk],
                                fit[["fit"]][["lambda"]][kk])))
})

test_that("the KKT conditions hold along the whole path of the group variant", {
  skip_no_sparsegl()
  fit <- wafc(x0, u0, y0, J = J, penalty = "sglasso")
  expect_equal(fit[["lambda.factor"]], 1)
  expect_equal(fit[["lambda"]], fit[["fit"]][["lambda"]])
  k <- wafc_kkt(fit)
  expect_true(all(k[["ok"]]))
})

test_that("a penalty level that kills the wavelet part leaves least squares on the levels", {
  ## sparsegl returns the zero vector at the first point of its own path,
  ## the level terms included, which is not the solution of the objective;
  ## wafc() writes the least squares fit on the level terms there, and the
  ## KKT check is what catches the difference.
  for (pen in c("lasso", if (has_sparsegl) "sglasso")) {
    fit <- wafc(x0, u0, y0, J = J, penalty = pen)
    expect_equal(fit[["nzero"]][1L], 0L)
    ls <- stats::lm.fit(cbind(1, x0[, -fit[["carrier"]][["index"]]]), y0)
    expect_equal(unname(fit[["a0"]][1L]), unname(ls[["coefficients"]][1L]),
                 tolerance = 1e-8)
    expect_equal(unname(as.matrix(fit[["beta"]])[2:3, 1L]),
                 unname(ls[["coefficients"]][2:3]), tolerance = 1e-8)
    expect_lt(wafc_kkt(fit, s = fit[["lambda"]][1L])[["unpenalized"]], 1e-8)
  }
})

test_that("the design of the fit can be reused, and sparse equals dense", {
  dd <- wafc_design(x0, u0, J = 4L, sparse = "never")
  ds <- wafc_design(x0, u0, J = 4L, sparse = "always")
  f1 <- wafc(design = dd, y = y0)
  f2 <- wafc(design = ds, y = y0)
  expect_equal(f1[["lambda"]], f2[["lambda"]], tolerance = 1e-8)
  ## glmnet drops the constant column only when the design is dense: with a
  ## sparse one it splits the level of x1 between the column and its own
  ## intercept. The split is not identified, the sum is, and it is the sum
  ## that coef.wafc reports, so everything the model sees agrees to machine
  ## precision either way.
  expect_equal(max(abs(as.matrix(f1[["beta"]])[1L, ])), 0)
  expect_gt(max(abs(as.matrix(f2[["beta"]])[1L, ])), 0.1)
  expect_equal(f1[["cc"]], f2[["cc"]], tolerance = 1e-8)
  expect_equal(as.matrix(f1[["beta"]])[-1L, ], as.matrix(f2[["beta"]])[-1L, ],
               tolerance = 1e-8)
  expect_equal(predict(f1), predict(f2), tolerance = 1e-8)
  ## building the design inside wafc() gives the same thing
  f3 <- wafc(x0, u0, y0, J = 4L, sparse = "never")
  expect_equal(as.matrix(f3[["beta"]]), as.matrix(f1[["beta"]]),
               tolerance = 1e-10)
  expect_equal(f3[["design"]][["nvars"]], dd[["nvars"]])
})

test_that("wafc validates its arguments and refuses the ambiguous designs", {
  expect_error(wafc(y = y0), "Supply 'x' and 'u'")
  expect_error(wafc(design = d0, y = y0[-1L]), "length 300")
  expect_error(wafc(x0, u0, y0, J = J, intercept = FALSE),
               "intercept = FALSE with the constant")
  expect_error(wafc(design = list(), y = y0), "wafc_design")
  expect_error(wafc(design = d0, y = y0, lambda = -1), "non-negative")
  ## two constant covariates: the level terms are not identified
  x2 <- cbind(x0, 1)
  colnames(x2) <- c(colnames(x0), "x4")
  expect_error(wafc(x2, u0, y0, J = J), "not identified")
  ## a covariate identically zero cannot carry the intercept either
  x3 <- x0
  x3[, 1L] <- 0
  expect_error(wafc(x3, u0, y0, J = J), "identically zero")
  ## no constant covariate: no intercept by default, and it may be asked for
  xv <- cbind(stats::rnorm(n), x0[, 2:3])
  fv <- wafc(xv, u0, y0, J = J)
  expect_false(fv[["intercept"]])
  expect_true(is.na(fv[["carrier"]][["index"]]))
  expect_equal(unname(coef(fv, s = fv[["lambda"]][30])["(Intercept)", 1L]), 0)
  fi <- wafc(xv, u0, y0, J = J, intercept = TRUE)
  expect_true(fi[["intercept"]])
  expect_true(all(wafc_kkt(fi)[["ok"]]))
})

test_that("a supplied path of penalty levels is the one of the objective", {
  lam <- c(0.5, 0.2, 0.05)
  fit <- wafc(x0, u0, y0, J = J, lambda = lam)
  expect_equal(fit[["lambda"]], lam)
  expect_equal(fit[["fit"]][["lambda"]], lam / fit[["lambda.factor"]])
  expect_true(all(wafc_kkt(fit)[["ok"]]))
  ## the same penalty level means the same penalty in the two engines:
  ## with asparse = 1 the sparse group LASSO is the LASSO, so the two fits
  ## reach the same value of the objective. sparsegl does not converge on
  ## that corner, which is why wafc() warns and why the agreement is only
  ## to half a percent; the LASSO is the other engine.
  skip_no_sparsegl()
  expect_warning(fg <- wafc(x0, u0, y0, J = J, penalty = "sglasso",
                            lambda = lam, asparse = 1),
                 "does not converge reliably")
  obj <- function(f, k) {
    cf <- coef(f, s = f[["lambda"]][k])[, 1L]
    r <- y0 - cbind(1, as.matrix(f[["design"]][["Z"]])) %*% cf
    pen <- seq_len(f[["nvars"]])[-f[["design"]][["unpenalized"]]]
    mean(r^2) / 2 + f[["lambda"]][k] * sum(abs(cf[-1L][pen]))
  }
  for (k in seq_along(lam)) {
    expect_equal(obj(fg, k), obj(fit, k), tolerance = 5e-3)
  }
  ## and the optimality conditions are the ones that show it does not
  expect_false(all(wafc_kkt(fg)[["ok"]]))
  expect_lt(max(wafc_kkt(fg)[["subgradient.abs"]] / fg[["lambda"]]), 0.1)
})

test_that("print.wafc reports the fit", {
  fit <- wafc(x0, u0, y0, J = J)
  out <- utils::capture.output(print(fit))
  expect_match(out[1L], "LASSO")
  expect_true(any(grepl("carried by the intercept", out)))
  expect_true(any(grepl("penalty level", out)))
})

## ---------------------------------------------------------------------------
## coef, predict and the reconstruction
## ---------------------------------------------------------------------------

test_that("coef folds the intercept into the level of the constant covariate", {
  fit <- wafc(x0, u0, y0, J = J)
  s <- fit[["lambda"]][30]
  cf <- coef(fit, s = s)
  raw <- stats::coef(fit[["fit"]], s = s / fit[["lambda.factor"]])
  expect_equal(unname(cf["x1", 1L]), unname(as.numeric(raw)[1L]),
               tolerance = 1e-10)
  expect_equal(unname(cf["(Intercept)", 1L]), 0)
  ## the rest of the vector is untouched
  expect_equal(unname(cf[-(1:2), 1L]), unname(as.numeric(raw)[-(1:2)]),
               tolerance = 1e-10)
  ## the whole path, and the interpolation at its own points
  expect_equal(dim(coef(fit)), c(1L + fit[["nvars"]], length(fit[["lambda"]])))
  expect_equal(unname(coef(fit, s = fit[["lambda"]][c(5L, 30L)])),
               unname(coef(fit)[, c(5L, 30L)]))
  ## and the interpolation between two points of the path is the one of
  ## glmnet, on the wavelet part, which the folding does not touch
  smid <- sqrt(fit[["lambda"]][30] * fit[["lambda"]][31])
  ours <- coef(fit, s = smid)[, 1L]
  theirs <- as.numeric(stats::coef(fit[["fit"]],
                                   s = smid / fit[["lambda.factor"]]))
  expect_equal(unname(ours[-(1:2)]), theirs[-(1:2)], tolerance = 1e-10)
  expect_equal(unname(ours[2L]), theirs[1L] + theirs[2L], tolerance = 1e-10)
  ## outside the path the endpoints are used
  expect_equal(unname(coef(fit, s = 1e3)), unname(coef(fit)[, 1L, drop = FALSE]))
})

test_that("predict reproduces the fitted values and extrapolates by the spec", {
  fit <- wafc(x0, u0, y0, J = J)
  s <- fit[["lambda"]][30]
  yh <- predict(fit, s = s)
  expect_equal(dim(yh), c(n, 1L))
  manual <- fit[["a0"]][30] +
    as.numeric(as.matrix(fit[["design"]][["Z"]]) %*% fit[["beta"]][, 30])
  expect_equal(as.numeric(yh), manual, tolerance = 1e-10)
  ## the folded coefficients give the same prediction, because the constant
  ## column is one in the new data as well
  expect_equal(as.numeric(yh),
               as.numeric(cbind(1, as.matrix(fit[["design"]][["Z"]])) %*%
                            coef(fit, s = s)),
               tolerance = 1e-10)
  ## new data: same columns, modulating covariates truncated to the range
  ## used for the rescaling
  te <- simulate_wafc(200L, p = p, q = q, scenario = "smooth", seed = 77L,
                      sigma = dgp[["sigma"]])
  yn <- predict(fit, te[["x"]], te[["u"]], s = s)
  expect_equal(dim(yn), c(200L, 1L))
  expect_true(all(is.finite(yn)))
  expect_lt(sqrt(mean((yn - te[["f"]])^2)), 2 * dgp[["sigma"]])
  expect_error(predict(fit, te[["x"]], s = s), "both 'newx' and 'newu'")
  ## the whole path at once
  expect_equal(dim(predict(fit)), c(n, length(fit[["lambda"]])))
  ## coefficients and non-zero indices
  expect_equal(predict(fit, s = s, type = "coefficients"), coef(fit, s = s))
  nz <- predict(fit, s = s, type = "nonzero")[[1L]]
  expect_equal(length(nz), fit[["nzero"]][30])
  expect_true(all(nz > p))
})

test_that("the reconstruction recovers the components that generated the data", {
  ## no noise and theta* in the basis: g_hat on a grid is the true
  ## expansion of the block evaluated at the same points
  fit <- wafc(design = d0, y = y_exact, lambda = c(0.1, 0.01, 1e-4, 1e-8))
  fn <- wafc_functions(fit, s = 1e-8, n_grid = 64L)
  expect_s3_class(fn, "wafc_functions")
  expect_equal(dim(fn[["grid"]]), c(64L, q))
  expect_equal(unname(fn[["cc"]]), theta[1:3], tolerance = 1e-4)
  d_grid <- wafc_design(matrix(1, 64L, p), fn[["grid"]], spec = d0)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      idx <- d0[["blocks"]][[paste0("x", l, ":u", m)]]
      truth <- as.numeric(d_grid[["Z"]][, idx, drop = FALSE] %*% theta[idx])
      expect_equal(fn[["g"]][[l, m]], truth, tolerance = 1e-4)
    }
  }
  ## the blocks that are zero in theta* are zero in the reconstruction. At
  ## lambda = 1e-8 the fit is least squares, which leaves rounding and not
  ## exact zeros, so what is checked is the size and not the count.
  expect_equal(max(abs(fn[["g"]][[2L, 1L]])), 0, tolerance = 1e-4)
  expect_equal(unname(fn[["norm"]]["x2", "u1"]), 0, tolerance = 1e-4)
  expect_equal(unname(fn[["norm"]]["x1", "u1"]), sqrt(1.5^2 + 1),
               tolerance = 1e-4)
  ## and at a penalty of order one the zero blocks are exactly zero
  fn2 <- wafc_functions(fit, s = 0.1, n_grid = 64L)
  expect_equal(unname(fn2[["nonzero"]]["x2", "u1"]), 0L)
  expect_equal(fn2[["g"]][[2L, 1L]], numeric(64L))
  out <- utils::capture.output(print(fn))
  expect_match(out[1L], "lambda")
})

test_that("the default grid is the range of the training sample", {
  fit <- wafc(x0, u0, y0, J = J, rescale = TRUE)
  fn <- wafc_functions(fit, s = fit[["lambda"]][30], n_grid = 32L)
  expect_equal(apply(fn[["grid"]], 2L, min), apply(u0, 2L, min),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(apply(fn[["grid"]], 2L, max), apply(u0, 2L, max),
               tolerance = 1e-10, ignore_attr = TRUE)
  ## a grid given as a vector is used for every modulating covariate
  gv <- seq(0.2, 0.8, length.out = 16L)
  fn2 <- wafc_functions(fit, s = fit[["lambda"]][30], grid = gv)
  expect_equal(fn2[["grid"]][, 1L], gv, ignore_attr = TRUE)
  expect_equal(fn2[["grid"]][, 2L], gv, ignore_attr = TRUE)
  expect_error(wafc_functions(fit, s = c(0.1, 0.2)), "single non-negative")
  expect_error(wafc_functions(fit, grid = matrix(0.5, 4L, 3L)), "1 or 2 column")
})

test_that("beta_hat is the level plus the components, and matches the truth", {
  fit <- wafc(x0, u0, y0, J = 4L)
  s <- fit[["lambda"]][40]
  bh <- predict(fit, newu = u0, s = s, type = "beta")
  expect_equal(dim(bh), c(n, p))
  expect_equal(colnames(bh), colnames(x0))
  ## beta_l(u) = c_l + sum_m g_lm(u_m), read off the reconstruction
  fn <- wafc_functions(fit, s = s, grid = u0)
  manual <- matrix(rep(fn[["cc"]], each = n), n, p)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) manual[, l] <- manual[, l] + fn[["g"]][[l, m]]
  }
  expect_equal(unname(bh), manual, tolerance = 1e-10)
  ## the constant coefficient beta_3 is recovered well, and the two
  ## functional ones are closer to the truth than a constant fit would be
  expect_equal(mean(bh[, 3L]), dgp[["cc"]][3L], tolerance = 0.1)
  for (l in 1:2) {
    expect_lt(mean((bh[, l] - dgp[["beta"]][, l])^2),
              0.5 * mean((dgp[["beta"]][, l] - mean(dgp[["beta"]][, l]))^2))
  }
  expect_error(predict(fit, s = s, type = "beta"), "needs 'newu'")
})

test_that("wafc_blocks reads the structure, and the group variant zeroes blocks", {
  fit <- wafc(x0, u0, y0, J = 4L)
  b <- wafc_blocks(fit, s = fit[["lambda"]][40])
  expect_equal(dim(b[["nonzero"]]), c(p, q))
  expect_equal(dimnames(b[["nonzero"]]), list(colnames(x0), colnames(u0)))
  expect_equal(sum(b[["nonzero"]]), fit[["nzero"]][40])
  ## several penalty levels give an array
  ba <- wafc_blocks(fit, s = fit[["lambda"]][c(20L, 40L)])
  expect_equal(dim(ba[["nonzero"]]), c(p, q, 2L))

  ## the numerical argument for the variant of D3: at a comparable number
  ## of non-zero coefficients the group penalty empties whole blocks, which
  ## is the selection of structure, and the LASSO does not
  skip_no_sparsegl()
  fg <- wafc(x0, u0, y0, J = 4L, penalty = "sglasso")
  target <- fit[["nzero"]][40]
  kg <- which.min(abs(fg[["nzero"]] - target))
  zl <- sum(wafc_blocks(fit, s = fit[["lambda"]][40])[["nonzero"]] == 0)
  zg <- sum(wafc_blocks(fg, s = fg[["lambda"]][kg])[["nonzero"]] == 0)
  expect_gt(zg, zl)
  expect_gt(zg, 0L)
})

test_that("the null scenario is shrunk to a linear model by both penalties", {
  nd <- simulate_wafc(400L, p = 3L, q = 2L, scenario = "null", seed = 5L,
                      sigma = 0.5)
  for (pen in c("lasso", if (has_sparsegl) "sglasso")) {
    fit <- wafc(nd[["x"]], nd[["u"]], nd[["y"]], J = 3L, penalty = pen)
    ## somewhere on the path the wavelet part is empty and the levels are
    ## the least squares ones
    k <- which(fit[["nzero"]] == 0L)
    expect_gt(length(k), 0L)
    expect_equal(unname(fit[["cc"]][, max(k)]), nd[["cc"]], tolerance = 0.15)
    fn <- wafc_functions(fit, s = fit[["lambda"]][max(k)], n_grid = 32L)
    expect_equal(sum(fn[["nonzero"]]), 0L)
    expect_equal(max(abs(unlist(fn[["g"]]))), 0)
  }
})
