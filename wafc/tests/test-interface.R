## wafc/tests/test-interface.R -- tests of the interface of step E3.1: the
## block LASSO of decision D44 as penalty = "block" of wafc() and
## cv.wafc(), the default of both, and the threshold of decision D45 inside
## cv.wafc(), read by coef(), predict(), wafc_functions(), wafc_blocks()
## and wafc_grid_components(). Run from the root of the repository:
##
##     Rscript -e 'testthat::test_dir("wafc/tests")'
##
## Three tests carry the step. The first is again the exact recovery of the
## plan, now for the block LASSO: with the components in the basis and no
## noise, the fit with lambda -> 0 returns theta*. The second is the proof
## that nothing moved: cv.wafc(penalty = "block") selects the pair of
## wafc_fit_klopp(penalize.levels = FALSE, balanced = TRUE) on the same
## folds and returns its fit, coefficient for coefficient, and its threshold
## by the rule "cv1se" is the one of wafc_threshold() on that fit, so the
## WAFC of the interface is the 'klopp.balanced+cv1se' of the tables of
## step E2.5j. The third is the KKT check of the block LASSO, which shows
## that the lambda of grpreg is the lambda of the objective written in
## ?wafc, with the norm of a chunk on the orthonormalized scale.

library(testthat)

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

has_grpreg <- requireNamespace("grpreg", quietly = TRUE)
skip_no_grpreg <- function() {
  testthat::skip_if_not(has_grpreg, "grpreg is not installed")
}

n <- 250L
p <- 3L
q <- 2L
dgp <- simulate_wafc(n, p = p, q = q, scenario = "smooth", seed = 20261003L,
                     snr = 4)
x0 <- dgp[["x"]]
u0 <- dgp[["u"]]
y0 <- dgp[["y"]]
folds <- rep_len(1:5, n)
test0 <- simulate_wafc(400L, p = p, q = q, scenario = "smooth",
                       seed = 20261004L, snr = 4)
grid0 <- cbind(seq(min(u0[, 1L]), max(u0[, 1L]), length.out = 64),
               seq(min(u0[, 2L]), max(u0[, 2L]), length.out = 64))

## One cross-validation of each kind, shared by the tests below: the WAFC of
## D44 and D45 with every argument at its default, the block LASSO of
## wafc_fit_klopp() it has to reproduce, and the LASSO of before.
if (has_grpreg) {
  cvb <- cv.wafc(x0, u0, y0, J = 2:4, foldid = folds)
  kp <- wafc_competitor("klopp", x0, u0, y0, J = 2:4, foldid = folds,
                        penalize.levels = FALSE, balanced = TRUE)
}
cvl <- cv.wafc(x0, u0, y0, J = 3:4, foldid = folds, penalty = "lasso",
               threshold = "none")

## ---------------------------------------------------------------------------
## The defaults
## ---------------------------------------------------------------------------

test_that("the defaults are the block LASSO and the threshold cv1se", {
  expect_identical(eval(formals(wafc)[["penalty"]])[1L], "block")
  expect_identical(eval(formals(cv.wafc)[["penalty"]])[1L], "block")
  expect_identical(eval(formals(cv.wafc)[["threshold"]])[1L], "cv1se")
  expect_identical(eval(formals(wafc_threshold)[["rule"]])[1L], "cv1se")
  ## the rules that were defaults elsewhere are left as they were
  expect_identical(eval(formals(wafc_tune)[["penalty"]])[1L], "lasso")
  expect_identical(eval(formals(wafc_sigma)[["penalty"]])[1L], "lasso")
  expect_identical(eval(formals(wafc_fit_oracle)[["penalty"]])[1L], "lasso")
  ## the tolerance of the block LASSO is the default of grpreg, at which
  ## steps E2.5e to E2.5j measured it, and the other two keep 1e-9
  expect_identical(eval(formals(wafc)[["thresh"]], list(penalty = "block")),
                   1e-4)
  expect_identical(eval(formals(wafc)[["thresh"]], list(penalty = "lasso")),
                   1e-9)
  ## 'threshold' comes after the dots, so 'thresh' still reaches wafc()
  expect_true(match("threshold", names(formals(cv.wafc))) >
                match("...", names(formals(cv.wafc))))
  skip_no_grpreg()
  d <- wafc_design(x0, u0, J = 3L)
  a <- wafc(design = d, y = y0)
  b <- wafc(design = d, y = y0, penalty = "block", thresh = 1e-4)
  expect_identical(a[["penalty"]], "block")
  expect_identical(a[["beta"]], b[["beta"]])
  expect_identical(a[["lambda"]], b[["lambda"]])
  expect_identical(cvb[["penalty"]], "block")
  expect_identical(cvb[["threshold"]][["rule"]], "cv1se")
})

## ---------------------------------------------------------------------------
## The block LASSO in wafc()
## ---------------------------------------------------------------------------

test_that("wafc(penalty = \"block\") recovers a truth in the basis without noise", {
  skip_no_grpreg()
  d0 <- wafc_design(x0, u0, J = 3L, rescale = FALSE)
  theta <- numeric(d0[["nvars"]])
  theta[d0[["unpenalized"]]] <- c(1, 2, -1.5)
  theta[d0[["blocks"]][["x1:u1"]][c(2L, 5L)]] <- c(1.5, -0.8)
  theta[d0[["blocks"]][["x2:u2"]][3L]] <- 0.7
  y_exact <- as.numeric(d0[["Z"]] %*% theta)
  fit <- wafc(design = d0, y = y_exact, penalty = "block",
              lambda = c(1, 0.1, 1e-3, 1e-6, 1e-9), thresh = 1e-12)
  cf <- coef(fit, s = 1e-9)[, 1L]
  ## the level of the constant covariate is in its row, not the intercept
  expect_equal(unname(cf[1L]), 0)
  expect_equal(unname(cf[-1L]), theta, tolerance = 1e-7)
  expect_equal(unname(fit[["cc"]][, 5L]), c(1, 2, -1.5), tolerance = 1e-7)
  ## and at the largest level every chunk is zero, with the levels fitted
  expect_identical(fit[["nzero"]][1L], 0L)
})

test_that("the block LASSO is grpreg on the balanced chunks of wafc_fit_klopp", {
  skip_no_grpreg()
  d <- wafc_design(x0, u0, J = 4L)
  bs <- as.integer(ceiling(log(n)))
  grp <- wafc_kp_groups(d, bs, penalize.levels = FALSE, balanced = TRUE)
  fit <- wafc(design = d, y = y0)
  expect_identical(fit[["group"]][["group"]], grp)
  expect_identical(fit[["group"]][["block.size"]], bs)
  ref <- grpreg::grpreg(as.matrix(d[["Z"]]), y0, group = grp,
                        penalty = "grLasso")
  expect_identical(fit[["lambda"]], ref[["lambda"]])
  expect_identical(unname(as.matrix(fit[["beta"]])),
                   unname(ref[["beta"]][-1L, ]))
  expect_identical(fit[["a0"]], unname(ref[["beta"]][1L, ]))
  expect_identical(fit[["group"]][["multiplier"]],
                   as.numeric(ref[["group.multiplier"]]))
  ## every chunk of a finer level has between b_n and 2 b_n - 1 columns, and
  ## the weights are the square roots of the sizes
  expect_true(all(fit[["group"]][["sizes"]] <= 2L * bs - 1L))
  expect_equal(fit[["group"]][["multiplier"]],
               sqrt(fit[["group"]][["sizes"]]))
  ## block.size moves the chunks, and the default is ceiling(log n)
  f4 <- wafc(design = d, y = y0, block.size = 4L)
  expect_identical(f4[["group"]][["group"]],
                   wafc_kp_groups(d, 4L, penalize.levels = FALSE,
                                  balanced = TRUE))
  expect_false(identical(f4[["group"]][["group"]], grp))
  ## the levels carry the intercept as with the other two engines
  expect_identical(fit[["carrier"]][["index"]], 1L)
  expect_equal(unname(fit[["cc"]][1L, ]),
               unname(as.matrix(fit[["beta"]])[1L, ] + fit[["a0"]]))
  ## what grpreg returned of the path
  cv <- fit[["conv"]]
  expect_identical(cv[["nlambda"]], 100L)
  expect_identical(cv[["nreturned"]], length(fit[["lambda"]]))
  expect_true(is.na(cv[["jerr"]]))
  expect_identical(cv[["max.iter"]], 10000L)
  expect_identical(cv[["iter.total"]], as.integer(sum(ref[["iter"]])))
})

test_that("the lambda of the block LASSO is the lambda of its objective", {
  skip_no_grpreg()
  ## at a tight tolerance every point of the path satisfies the conditions
  ## of the objective of ?wafc, the norm of a chunk taken on the
  ## orthonormalized scale and the weight the square root of its size
  d <- wafc_design(x0, u0, J = 3L)
  fit <- wafc(design = d, y = y0, thresh = 1e-8)
  k <- wafc_kkt(fit)
  expect_true(all(k[["ok"]]))
  expect_lt(max(k[["subgradient"]]), 1 + 1e-6)
  expect_identical(k[["lambda"]], k[["lambda.engine"]])
  ## and a lambda 5% off breaks them wherever a chunk is not zero and 5% of
  ## lambda is above the absolute tolerance of the check
  off <- fit
  off[["lambda"]] <- 1.05 * fit[["lambda"]]
  ko <- wafc_kkt(off)
  act <- fit[["nzero"]] > 0L & fit[["lambda"]] > 0.01
  expect_gt(sum(act), 20L)
  expect_true(all(!ko[["ok"]][act]))
  ## at the default tolerance of grpreg, the one of the measurements, the
  ## conditions hold at the top of the path and not at its small end: the
  ## reason ?wafc gives for reading the argument 'thresh'
  k4 <- wafc_kkt(wafc(design = d, y = y0))
  expect_true(all(k4[["ok"]][1:10]))
  expect_false(all(k4[["ok"]]))
})

test_that("the block LASSO refuses what grpreg cannot honour", {
  skip_no_grpreg()
  d <- wafc_design(x0, u0, J = 3L)
  expect_error(wafc(design = d, y = y0, intercept = FALSE),
               "always has an intercept")
  expect_error(wafc(design = d, y = y0, standardize = TRUE),
               "standardize = FALSE")
  expect_error(wafc(design = d, y = y0, block.size = 0L), "positive integer")
  expect_error(cv.wafc(x0, u0, y0, J = 2L, foldid = folds,
                       intercept = FALSE), "always has an intercept")
  expect_error(cv.wafc(x0, u0, y0, J = 2L, foldid = folds,
                       standardize = TRUE), "standardize = FALSE")
  ## without a constant covariate grpreg still fits an intercept, which
  ## stays in the intercept row and in the predictions
  xv <- x0[, -1L, drop = FALSE]
  fv <- wafc(xv, u0, y0, J = 3L)
  expect_true(fv[["intercept"]])
  expect_true(is.na(fv[["carrier"]][["index"]]))
  cf <- coef(fv, s = fv[["lambda"]][30L])
  expect_equal(unname(cf[1L, 1L]), fv[["a0"]][30L])
  expect_equal(as.numeric(predict(fv, xv, u0, s = fv[["lambda"]][30L])),
               as.numeric(fv[["a0"]][30L] +
                            fv[["design"]][["Z"]] %*% cf[-1L, 1L]))
})

## ---------------------------------------------------------------------------
## cv.wafc(penalty = "block") is wafc_fit_klopp(), balanced, levels free
## ---------------------------------------------------------------------------

test_that("cv.wafc(penalty = \"block\") reproduces wafc_fit_klopp on the same folds", {
  skip_no_grpreg()
  expect_identical(cvb[["J.min"]], kp[["extra"]][["J"]])
  expect_identical(cvb[["lambda.min"]], kp[["extra"]][["lambda"]])
  expect_identical(cvb[["cvm.min"]], kp[["extra"]][["cve"]])
  ## the coefficients of the engine, intercept first, at 1e-10 (they are
  ## the same numbers: the fit is the one cv.grpreg() returned)
  raw <- wafc_raw_coef(cvb[["wafc.fit"]], s = cvb[["lambda.min"]])[, 1L]
  expect_equal(unname(raw), as.numeric(coef(kp[["fit"]])), tolerance = 1e-10)
  ## the predictions on a test sample, the levels and the functional
  ## coefficients, without the threshold
  expect_equal(as.numeric(predict(cvb, test0[["x"]], test0[["u"]],
                                  thresholded = FALSE)),
               predict(kp, test0[["x"]], test0[["u"]]), tolerance = 1e-10)
  expect_equal(unname(coef(cvb, thresholded = FALSE)[1L + 1:p, 1L]),
               unname(kp[["cc"]]), tolerance = 1e-10)
  expect_equal(unname(predict(cvb, newu = test0[["u"]], type = "beta",
                              thresholded = FALSE)),
               unname(kp[["beta"]](test0[["u"]])), tolerance = 1e-10)
  ## the cross-validated curve of every J is the one of cv.grpreg()
  for (i in seq_along(cvb[["J"]])) {
    d <- wafc_design(x0, u0, J = cvb[["J"]][i])
    ref <- grpreg::cv.grpreg(as.matrix(d[["Z"]]), y0,
                             group = wafc_kp_groups(d, 6L, FALSE,
                                                    balanced = TRUE),
                             penalty = "grLasso", fold = folds)
    z <- cvb[["cv"]][[i]]
    expect_identical(z[["lambda"]], as.numeric(ref[["lambda"]]))
    expect_identical(z[["cvm"]], as.numeric(ref[["cve"]]))
    expect_identical(z[["cvsd"]], as.numeric(ref[["cvse"]]))
  }
  ## the fit on the whole sample is the one wafc() returns at that J
  f <- wafc(x0, u0, y0, J = cvb[["J.min"]])
  expect_identical(f[["beta"]], cvb[["wafc.fit"]][["beta"]])
  expect_identical(f[["lambda"]], cvb[["wafc.fit"]][["lambda"]])
  ## and the convergence table is the one wafc_fit_klopp() keeps
  expect_equal(wafc_cv_convergence(cvb), kp[["extra"]][["conv"]])
})

test_that("the reproduction does not depend on the folds having equal sizes", {
  skip_no_grpreg()
  ## 250 observations in 7 folds: the mean over the observations, which is
  ## what cv.grpreg() reports, is not the mean of the fold means there
  f7 <- rep_len(1:7, n)
  a <- cv.wafc(x0, u0, y0, J = 3:4, foldid = f7, threshold = "none")
  b <- wafc_fit_klopp(x0, u0, y0, J = 3:4, foldid = f7,
                      penalize.levels = FALSE, balanced = TRUE)
  expect_identical(a[["J.min"]], b[["extra"]][["J"]])
  expect_identical(a[["lambda.min"]], b[["extra"]][["lambda"]])
  expect_equal(unname(wafc_raw_coef(a[["wafc.fit"]], s = a[["lambda.min"]])[, 1L]),
               as.numeric(coef(b[["fit"]])), tolerance = 1e-10)
})

test_that("type.measure = \"mae\" of the block LASSO is read from the held-out predictions", {
  skip_no_grpreg()
  a <- cv.wafc(x0, u0, y0, J = 3L, foldid = folds, type.measure = "mae",
               threshold = "none")
  d <- wafc_design(x0, u0, J = 3L)
  ref <- grpreg::cv.grpreg(as.matrix(d[["Z"]]), y0,
                           group = wafc_kp_groups(d, 6L, FALSE,
                                                  balanced = TRUE),
                           penalty = "grLasso", fold = folds, returnY = TRUE)
  E <- abs(y0 - ref[["Y"]])
  z <- a[["cv"]][[1L]]
  expect_equal(z[["cvm"]], colMeans(E), tolerance = 1e-12)
  expect_equal(z[["cvsd"]], apply(E, 2L, sd) / sqrt(n), tolerance = 1e-12)
  expect_identical(z[["lambda.min"]], z[["lambda"]][which.min(z[["cvm"]])])
  i1 <- which.min(z[["cvm"]])
  expect_identical(z[["lambda.1se"]],
                   z[["lambda"]][min(which(z[["cvm"]] <=
                                             z[["cvm"]][i1] + z[["cvsd"]][i1]))])
})

## ---------------------------------------------------------------------------
## The threshold inside cv.wafc()
## ---------------------------------------------------------------------------

test_that("the WAFC is the block LASSO followed by the threshold cv1se", {
  skip_no_grpreg()
  th <- wafc_threshold(cvb, rule = "cv1se")
  ex <- th[["extra"]]
  tt <- cvb[["threshold"]]
  expect_identical(tt[["kept"]], ex[["kept"]])
  expect_identical(tt[["t"]], ex[["t"]])
  expect_identical(tt[["candidates"]], ex[["candidates"]])
  expect_identical(tt[["nzero"]], ex[["nzero"]])
  ## the thresholded fit, read by every accessor of cv.wafc
  expect_equal(as.numeric(predict(cvb, test0[["x"]], test0[["u"]])),
               predict(th, test0[["x"]], test0[["u"]]), tolerance = 1e-10)
  expect_equal(as.numeric(predict(cvb)), th[["fitted"]], tolerance = 1e-10)
  expect_equal(unname(predict(cvb, newu = test0[["u"]], type = "beta")),
               unname(th[["beta"]](test0[["u"]])), tolerance = 1e-10)
  fn <- wafc_functions(cvb, grid = grid0)
  gt <- th[["g"]](grid0)
  for (i in seq_along(gt)) {
    expect_equal(fn[["g"]][[i]], gt[[i]], tolerance = 1e-10)
  }
  expect_equal(unname(fn[["cc"]]), unname(th[["cc"]]), tolerance = 1e-10)
  expect_identical(wafc_blocks(cvb)[["nonzero"]] > 0L, th[["blocks"]])
  gc <- wafc_grid_components(cvb, grid0)
  gr <- wafc_grid_components(th, grid0)
  for (i in seq_along(gr)) expect_equal(gc[[i]], gr[[i]], tolerance = 1e-10)
  ## coef() is the fit at lambda.min with the blocks not kept set to zero
  cf <- coef(cvb)[, 1L]
  c0 <- coef(cvb, thresholded = FALSE)[, 1L]
  des <- cvb[["wafc.fit"]][["design"]]
  for (l in seq_len(p)) for (m in seq_len(q)) {
    idx <- 1L + des[["blocks"]][[wafc_block_name(des, l, m)]]
    if (tt[["kept"]][l, m]) {
      expect_identical(cf[idx], c0[idx])
    } else {
      expect_true(all(cf[idx] == 0))
    }
  }
  expect_identical(cf[1L + des[["unpenalized"]]], c0[1L + des[["unpenalized"]]])
  ## and it is the 'klopp.balanced+cv1se' of the tables of step E2.5j
  tk <- wafc_threshold(kp, rule = "cv1se", y = y0)
  expect_identical(tk[["extra"]][["kept"]], tt[["kept"]])
  expect_equal(tk[["extra"]][["t"]], tt[["t"]], tolerance = 1e-10)
  expect_equal(as.numeric(predict(cvb, test0[["x"]], test0[["u"]])),
               predict(tk, test0[["x"]], test0[["u"]]), tolerance = 1e-10)
  ## the fold fits of the threshold keep the chunks of the whole sample, as
  ## the ones of wafc_fit_klopp() do
  fa <- wafc_threshold_folds(cvb)
  fb <- wafc_threshold_folds(kp, y = y0)
  for (k in seq_along(fa[["fits"]])) {
    expect_equal(fa[["fits"]][[k]][["b"]], fb[["fits"]][[k]][["b"]],
                 tolerance = 1e-12)
    expect_equal(fa[["fits"]][[k]][["a0"]], fb[["fits"]][[k]][["a0"]],
                 tolerance = 1e-12)
  }
})

test_that("the other rules of the threshold are options, and none is the fit", {
  skip_no_grpreg()
  for (rule in c("cv", "max")) {
    a <- cv.wafc(x0, u0, y0, J = 3L, foldid = folds, threshold = rule)
    th <- wafc_threshold(a, rule = rule)
    expect_identical(a[["threshold"]][["rule"]], rule)
    expect_identical(a[["threshold"]][["kept"]], th[["extra"]][["kept"]])
    expect_identical(a[["threshold"]][["t"]], th[["extra"]][["t"]])
    expect_equal(as.numeric(predict(a, test0[["x"]], test0[["u"]])),
                 predict(th, test0[["x"]], test0[["u"]]), tolerance = 1e-10)
  }
  ## "max" is the rule of step E2.5g at c = 0.15
  expect_identical(a[["threshold"]][["c"]], 0.15)
  ## without threshold coef() and predict() are the path at lambda.min, and
  ## the pair is the one of the thresholded call
  a1 <- a
  a <- cv.wafc(x0, u0, y0, J = 3L, foldid = folds, threshold = "none")
  expect_null(a[["threshold"]])
  expect_identical(a[["J.min"]], a1[["J.min"]])
  expect_identical(a[["lambda.min"]], a1[["lambda.min"]])
  expect_identical(coef(a), coef(a[["wafc.fit"]], s = a[["lambda.min"]]))
  expect_identical(predict(a, test0[["x"]], test0[["u"]]),
                   predict(a[["wafc.fit"]], test0[["x"]], test0[["u"]],
                           s = a[["lambda.min"]]))
  expect_error(cv.wafc(x0, u0, y0, J = 2L, foldid = folds,
                       threshold = "oracle"), "cv1se")
})

test_that("away from lambda.min the accessors read the path without threshold", {
  skip_no_grpreg()
  s1 <- cvb[["lambda.1se"]]
  expect_identical(coef(cvb, s = "lambda.1se"),
                   coef(cvb[["wafc.fit"]], s = s1))
  expect_identical(predict(cvb, test0[["x"]], test0[["u"]], s = "lambda.1se"),
                   predict(cvb[["wafc.fit"]], test0[["x"]], test0[["u"]],
                           s = s1))
  expect_identical(wafc_blocks(cvb, s = "lambda.1se"),
                   wafc_blocks(cvb[["wafc.fit"]], s = s1))
  expect_identical(coef(cvb, thresholded = FALSE),
                   coef(cvb[["wafc.fit"]], s = cvb[["lambda.min"]]))
  ## a numeric s equal to lambda.min is lambda.min
  expect_identical(coef(cvb, s = cvb[["lambda.min"]]), coef(cvb))
  expect_error(coef(cvb, s = "lambda.1se", thresholded = TRUE),
               "chosen at lambda.min")
  expect_error(coef(cvl, thresholded = TRUE), "has no threshold")
  expect_error(coef(cvb, thresholded = "yes"), "NULL, TRUE or FALSE")
  ## a cv.wafc without threshold reads the path, as before step E3.1
  expect_identical(coef(cvl), coef(cvl[["wafc.fit"]], s = cvl[["lambda.min"]]))
})

test_that("the LASSO goes through the threshold of cv.wafc as well", {
  a <- cv.wafc(x0, u0, y0, J = 3:4, foldid = folds, penalty = "lasso")
  expect_identical(a[["J.min"]], cvl[["J.min"]])
  expect_identical(a[["lambda.min"]], cvl[["lambda.min"]])
  expect_identical(a[["cvtab"]], cvl[["cvtab"]])
  th <- wafc_threshold(cvl, rule = "cv1se")
  expect_identical(a[["threshold"]][["kept"]], th[["extra"]][["kept"]])
  expect_identical(a[["threshold"]][["t"]], th[["extra"]][["t"]])
  expect_equal(as.numeric(predict(a, test0[["x"]], test0[["u"]])),
               predict(th, test0[["x"]], test0[["u"]]), tolerance = 1e-10)
  ## the engine arguments of the call reach the fold fits of the threshold
  b <- cv.wafc(x0, u0, y0, J = 3L, foldid = folds, penalty = "lasso",
               thresh = 1e-7)
  thb <- wafc_threshold(cv.wafc(x0, u0, y0, J = 3L, foldid = folds,
                                penalty = "lasso", thresh = 1e-7,
                                threshold = "none"), thresh = 1e-7)
  expect_identical(b[["threshold"]][["candidates"]],
                   thb[["extra"]][["candidates"]])
})

test_that("the cv.wafc object keeps the threshold without a second design", {
  skip_no_grpreg()
  tt <- cvb[["threshold"]]
  expect_false(any(vapply(tt, is.function, TRUE)))
  expect_named(tt, c("rule", "t", "c", "norm", "kept", "a0", "b", "nzero",
                     "candidates", "time"))
  expect_length(tt[["b"]], cvb[["wafc.fit"]][["nvars"]])
  ## a saved object grows by the coefficients and the candidates only
  bare <- cvb
  bare[["threshold"]] <- NULL
  expect_lt(length(serialize(cvb, NULL)) - length(serialize(bare, NULL)),
            1e5)
})

test_that("the print methods say which estimator they hold", {
  skip_no_grpreg()
  expect_output(print(cvb[["wafc.fit"]]), "block LASSO \\(.* balanced chunks, b_n = 6\\)")
  expect_output(print(cvb), "Threshold \"cv1se\" at lambda.min: t = ")
  expect_output(print(cvl), "No threshold on the blocks")
  expect_output(print(wafc(x0, u0, y0, J = 3L, penalty = "lasso")),
                "WAFC fit: LASSO")
})
