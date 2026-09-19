## wafc/tests/test-tune.R -- tests of the selection of (J, lambda) of step
## E2.3. Run from the root of the repository:
##
##     Rscript -e 'testthat::test_dir("wafc/tests")'
##
## The three tests that carry the step are the one that checks that the
## penalty levels the rules return are on the scale of the objective, by
## wafc_kkt() rather than by inspection of the engine (decision D17); the
## one that reproduces the resolution of equation (Jn) of E1.6 by brute
## force, including the value J = 7 at n = 12800 that the numerical check of
## that step reports; and the one that recomputes the criteria of
## wafc_bic() and wafc_ebic() from the residual sum of squares by hand.

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
                     snr = 4)
x0 <- dgp[["x"]]
u0 <- dgp[["u"]]
y0 <- dgp[["y"]]
folds <- rep_len(1:5, n)

## A design and a sparse truth in the basis, for the exact recovery test.
d0 <- wafc_design(x0, u0, J = 3L, rescale = FALSE)
theta <- numeric(d0[["nvars"]])
theta[d0[["unpenalized"]]] <- c(1, 2, -1.5)
theta[c(d0[["blocks"]][["x1:u1"]][c(1L, 4L)],
        d0[["blocks"]][["x2:u2"]][2L])] <- c(1.5, -1, 0.8)
y_exact <- as.numeric(d0[["Z"]] %*% theta)

## ---------------------------------------------------------------------------
## Cross-validation
## ---------------------------------------------------------------------------

test_that("cv.wafc selects a pair (J, lambda) whose optimality conditions hold", {
  cv <- cv.wafc(x0, u0, y0, J = 2:4, foldid = folds)
  expect_s3_class(cv, "cv.wafc")
  expect_true(cv[["J.min"]] %in% 2:4)
  expect_equal(nrow(cv[["cvtab"]]), 3L)
  expect_equal(cv[["cvtab"]][["J"]], 2:4)
  ## lambda.1se is at least lambda.min, and both belong to the path of the
  ## fit the object carries
  expect_gte(cv[["lambda.1se"]], cv[["lambda.min"]])
  lam <- cv[["wafc.fit"]][["lambda"]]
  expect_true(cv[["lambda.min"]] %in% lam)
  expect_true(cv[["lambda.1se"]] %in% lam)
  expect_equal(cv[["wafc.fit"]][["design"]][["J"]], rep(cv[["J.min"]], q))
  ## the scale of lambda is the one of the objective, checked on the
  ## stationarity conditions and not on the argument of the engine
  kkt <- wafc_kkt(cv[["wafc.fit"]],
                  s = c(cv[["lambda.min"]], cv[["lambda.1se"]]))
  expect_true(all(kkt[["ok"]]))
  ## the minimum is interior to the path, not at an end of it
  z <- cv[["cv"]][[match(cv[["J.min"]], cv[["J"]])]]
  imin <- match(cv[["lambda.min"]], z[["lambda"]])
  expect_gt(imin, 1L)
  expect_lt(imin, length(z[["lambda"]]))
  expect_lt(z[["cvm.min"]], min(z[["cvm"]][1L], z[["cvm"]][length(z[["cvm"]])]))
})

test_that("the folds are fixed, so two calls with the same foldid agree", {
  a <- cv.wafc(x0, u0, y0, J = 2:3, foldid = folds)
  b <- cv.wafc(x0, u0, y0, J = 2:3, foldid = folds)
  expect_equal(a[["cvtab"]], b[["cvtab"]])
  expect_equal(a[["cv"]][[1L]][["cvm"]], b[["cv"]][[1L]][["cvm"]])
  expect_equal(a[["lambda.min"]], b[["lambda.min"]])
  ## and the same partition serves every candidate J
  expect_equal(a[["foldid"]], folds)
  expect_equal(a[["nfolds"]], 5L)
})

test_that("the fold designs are the columns of the whole sample, restricted", {
  d <- wafc_design(x0, u0, J = 3L)
  rows <- which(folds != 1L)
  ds <- wafc_subset_design(d, rows)
  expect_equal(ds[["n"]], length(rows))
  expect_equal(ds[["nvars"]], d[["nvars"]])
  expect_equal(colnames(ds[["Z"]]), colnames(d[["Z"]]))
  expect_equal(ds[["blocks"]], d[["blocks"]])
  expect_equal(ds[["constant"]], d[["constant"]])
  expect_equal(as.matrix(ds[["Z"]]), as.matrix(d[["Z"]][rows, , drop = FALSE]))
  ## location, scale and eps are the ones of the whole sample, which is what
  ## makes the fold fits comparable
  expect_equal(ds[["location"]], d[["location"]])
  expect_equal(ds[["scale"]], d[["scale"]])
})

test_that("with theta* in the basis and no noise, the cross-validation goes to the small end", {
  ## an explicit path, because with a noiseless response the engine stops
  ## its own path early on the deviance ratio and never reaches a small
  ## enough penalty level
  cv <- cv.wafc(x0, u0, y_exact, J = 3L, foldid = folds, rescale = FALSE,
                lambda = c(0.1, 0.01, 1e-4, 1e-8))
  expect_equal(cv[["J.min"]], 3L)
  lam <- cv[["wafc.fit"]][["lambda"]]
  expect_equal(cv[["lambda.min"]], min(lam))
  expect_lt(cv[["cvm.min"]], 1e-6)
  cf <- coef(cv, s = "lambda.min")
  expect_equal(max(abs(unname(cf[-1L, 1L]) - theta)), 0, tolerance = 1e-3)
})

test_that("coef and predict of a cv.wafc read the selected pair", {
  cv <- cv.wafc(x0, u0, y0, J = 2:3, foldid = folds)
  fit <- cv[["wafc.fit"]]
  expect_equal(coef(cv, s = "lambda.min"), coef(fit, s = cv[["lambda.min"]]))
  expect_equal(coef(cv, s = "lambda.1se"), coef(fit, s = cv[["lambda.1se"]]))
  expect_equal(coef(cv, s = cv[["lambda.min"]]), coef(cv, s = "lambda.min"))
  expect_equal(predict(cv, x0, u0, s = "lambda.1se"),
               predict(fit, x0, u0, s = cv[["lambda.1se"]]))
  expect_equal(predict(cv, newu = u0, s = "lambda.min", type = "beta"),
               predict(fit, newu = u0, s = cv[["lambda.min"]], type = "beta"))
  expect_output(print(cv), "Cross-validated WAFC fit")
  expect_output(print(cv), "Selected: J =")
})

test_that("cv.wafc validates its arguments", {
  expect_error(cv.wafc(x0, u0, y0, J = 2.5), "integer")
  expect_error(cv.wafc(x0, u0, y0, J = 0L), "larger than j0")
  expect_error(cv.wafc(x0, u0, y0, J = 2L, nfolds = 2L), "at least 3")
  expect_error(cv.wafc(x0, u0, y0, J = 2L, foldid = folds[-1L]),
               "one entry per observation")
  expect_error(cv.wafc(x0, u0, y0, J = 2L, foldid = rep(1L, n)),
               "at least 3 folds")
  expect_error(cv.wafc(x0, u0, y0[-1L], J = 2L), "length")
})

## ---------------------------------------------------------------------------
## The information criteria
## ---------------------------------------------------------------------------

test_that("wafc_bic and wafc_ebic are the criteria they claim to be", {
  fit <- wafc(x0, u0, y0, J = 3L)
  ic <- wafc_bic(fit)
  expect_equal(nrow(ic), length(fit[["lambda"]]))
  expect_equal(ic[["lambda"]], fit[["lambda"]])
  ## degrees of freedom: the non-zero wavelet coefficients plus the p levels
  expect_equal(ic[["nzero"]], fit[["nzero"]])
  expect_equal(ic[["df"]], fit[["nzero"]] + p)
  ## residual sum of squares, recomputed from the fitted values
  eta <- predict(fit, x0, u0)
  expect_equal(unname(ic[["rss"]]), unname(colSums((y0 - eta)^2)))
  ## the criterion itself
  expect_equal(ic[["bic"]], n * log(ic[["rss"]] / n) + ic[["df"]] * log(n))
  ## and the extended one, which adds the log of the number of models of
  ## that size and reduces to the BIC at gamma = 0
  eb <- wafc_ebic(fit, gamma = 0.5)
  expect_equal(eb[["ebic"]],
               ic[["bic"]] + 2 * 0.5 * lchoose(fit[["npen"]], ic[["nzero"]]))
  expect_equal(wafc_ebic(fit, gamma = 0)[["ebic"]], ic[["bic"]])
  ## the minimiser is reported as an attribute and is interior
  k <- attr(ic, "which.min")
  expect_equal(attr(ic, "lambda.min"), ic[["lambda"]][k])
  expect_equal(ic[["bic"]][k], min(ic[["bic"]]))
  expect_gt(k, 1L)
  expect_lt(k, nrow(ic))
  expect_error(wafc_ebic(fit, gamma = 2), "in \\[0, 1\\]")
  expect_error(wafc_bic(fit[["design"]]), "returned by wafc")
})

test_that("the criteria can be read at a chosen penalty level", {
  fit <- wafc(x0, u0, y0, J = 3L)
  s <- fit[["lambda"]][c(20L, 40L)]
  ic <- wafc_bic(fit, s = s)
  expect_equal(ic[["lambda"]], s)
  expect_equal(ic[["bic"]], wafc_bic(fit)[["bic"]][c(20L, 40L)])
})

## ---------------------------------------------------------------------------
## The rule of the theory (E1.6)
## ---------------------------------------------------------------------------

test_that("wafc_J_theory is the resolution of equation (Jn) of E1.6", {
  brute <- function(n, s) {
    target <- (n / log(n))^(1 / (2 * s + 1))
    min(which(2^(1:40) >= target))
  }
  for (s in c(0.25, 0.5, 1.5, 3)) {
    for (nn in c(50, 250, 1000, 6400, 12800, 1e5)) {
      expect_equal(wafc_J_theory(nn, s = s), as.integer(brute(nn, s)))
    }
  }
  ## the value the numerical check of E1.6 reports: with s' = 1/4 the rule
  ## is at J = 7 at the largest sample size of that check, n = 12800, where
  ## the design has d = p q (2^7 - 1) columns
  expect_equal(wafc_J_theory(12800, s = 1/4), 7L)
  ## vectorised in n, non-decreasing, and never below j0 + 1
  Jn <- wafc_J_theory(c(100, 1000, 10000), s = 1)
  expect_equal(length(Jn), 3L)
  expect_false(is.unsorted(Jn))
  expect_gte(min(wafc_J_theory(c(3, 10, 100), s = 5)), 1L)
  expect_error(wafc_J_theory(2, s = 1), "at least 3")
  expect_error(wafc_J_theory(100, s = 0), "positive")
})

test_that("wafc_lambda_theory is the penalty level of Corollary 2 of E1.5", {
  d <- wafc_design(x0, u0, J = 3L)
  pen <- seq_len(d[["nvars"]])[-d[["unpenalized"]]]
  Z <- as.matrix(d[["Z"]][, pen, drop = FALSE])
  smax <- sqrt(max(colSums(Z^2) / n))
  alpha <- 0.05
  expect_equal(wafc_lambda_theory(d, sigma = 0.7, alpha = alpha),
               2 * 0.7 * smax * sqrt(2 * log(2 * length(pen) / alpha) / n))
  ## the reading of ||Z||_max that E1.5 calibrated: the largest column norm
  ## is O(1) in J, while the largest entry grows like 2^{J/2}
  norms <- entries <- numeric(0)
  for (J in 2:5) {
    dJ <- wafc_design(x0, u0, J = J)
    penJ <- seq_len(dJ[["nvars"]])[-dJ[["unpenalized"]]]
    ZJ <- as.matrix(dJ[["Z"]][, penJ, drop = FALSE])
    norms <- c(norms, sqrt(max(colSums(ZJ^2) / n)))
    entries <- c(entries, max(abs(ZJ)))
  }
  expect_lt(max(norms) / min(norms), 2)
  expect_gt(max(entries) / min(entries), 2)
  expect_error(wafc_lambda_theory(d), "Supply 'sigma'")
  expect_error(wafc_lambda_theory(d, sigma = -1), "non-negative")
  expect_error(wafc_lambda_theory(d[["Z"]], sigma = 1), "wafc_design")
})

test_that("wafc_sigma returns a scale near the true one at a rich enough sieve", {
  d <- wafc_design(x0, u0, J = 3L)
  est <- wafc_sigma(d, y0, foldid = folds)
  expect_equal(est[["method"]], "cv")
  expect_gt(est[["sigma"]], 0.5 * dgp[["sigma"]])
  expect_lt(est[["sigma"]], 2 * dgp[["sigma"]])
  expect_equal(est[["lambda"]],
               wafc_lambda_theory(d, sigma = est[["sigma"]]))
  ## the fixed point is the documented alternative, and at this sample size
  ## it stops at the null model, which is exactly what it warns about
  expect_warning(fp <- wafc_sigma(d, y0, method = "fixed.point"),
                 "no non-zero wavelet coefficient")
  expect_equal(fp[["method"]], "fixed.point")
  expect_true(fp[["converged"]])
  expect_gt(fp[["sigma"]], 0)
})

## ---------------------------------------------------------------------------
## The common entry point
## ---------------------------------------------------------------------------

test_that("every rule returns a pair on the path whose optimality conditions hold", {
  rules <- c("cv.min", "cv.1se", "bic", "ebic", "theory")
  for (rule in rules) {
    tn <- wafc_tune(x0, u0, y0, rule = rule, J = 2:4, foldid = folds,
                    s = 3/2, sigma = dgp[["sigma"]])
    expect_s3_class(tn, "wafc_tune")
    expect_equal(tn[["rule"]], rule)
    expect_true(tn[["lambda"]] %in% tn[["fit"]][["lambda"]])
    expect_equal(tn[["fit"]][["design"]][["J"]], rep(tn[["J"]], q))
    expect_true(all(wafc_kkt(tn[["fit"]], s = tn[["lambda"]])[["ok"]]))
    expect_equal(tn[["nzero"]], wafc_blocks(tn[["fit"]], s = tn[["lambda"]]) |>
                   (\(b) sum(b[["nonzero"]]))())
    expect_output(print(tn), "WAFC tuning by rule")
    expect_equal(dim(coef(tn)), c(1L + tn[["fit"]][["nvars"]], 1L))
    expect_equal(nrow(predict(tn, x0, u0)), n)
  }
})

test_that("the theory rule is exactly wafc_J_theory and wafc_lambda_theory", {
  tn <- wafc_tune(x0, u0, y0, rule = "theory", s = 3/2, sigma = dgp[["sigma"]])
  expect_equal(tn[["J"]], wafc_J_theory(n, s = 3/2))
  expect_equal(tn[["lambda"]],
               wafc_lambda_theory(tn[["fit"]][["design"]],
                                  sigma = dgp[["sigma"]]))
  expect_equal(tn[["sigma"]], dgp[["sigma"]])
  ## the path built for it starts where the engine would start, with every
  ## penalized coefficient at zero, and ends exactly at the chosen level
  lam <- tn[["fit"]][["lambda"]]
  expect_equal(min(lam), tn[["lambda"]])
  expect_equal(tn[["fit"]][["nzero"]][1L], 0L)
  ## and with sigma unknown it is estimated instead of being assumed
  tu <- wafc_tune(x0, u0, y0, rule = "theory", s = 3/2, foldid = folds)
  expect_equal(tu[["J"]], tn[["J"]])
  expect_gt(tu[["sigma"]], 0)
})

test_that("the cross-validation rules agree with cv.wafc on the same folds", {
  cv <- cv.wafc(x0, u0, y0, J = 2:4, foldid = folds)
  a <- wafc_tune(x0, u0, y0, rule = "cv.min", J = 2:4, foldid = folds)
  b <- wafc_tune(x0, u0, y0, rule = "cv.1se", J = 2:4, foldid = folds)
  expect_equal(a[["J"]], cv[["J.min"]])
  expect_equal(a[["lambda"]], cv[["lambda.min"]])
  expect_equal(b[["lambda"]], cv[["lambda.1se"]])
  expect_equal(a[["tab"]], cv[["cvtab"]])
})

test_that("the ebic rule is at least as sparse as the bic rule here", {
  a <- wafc_tune(x0, u0, y0, rule = "bic", J = 2:4)
  b <- wafc_tune(x0, u0, y0, rule = "ebic", J = 2:4)
  expect_lte(b[["nzero"]], a[["nzero"]])
  expect_equal(names(a[["tab"]]), c("J", "lambda", "nzero", "df", "bic"))
  expect_equal(names(b[["tab"]]), c("J", "lambda", "nzero", "df", "ebic"))
  expect_equal(a[["tab"]][["J"]], 2:4)
})

test_that("wafc_path_to ends at the level asked for and starts at the entry point", {
  d <- wafc_design(x0, u0, J = 3L)
  lam <- wafc_lambda_theory(d, sigma = dgp[["sigma"]])
  path <- wafc_path_to(lam, d, y0)
  expect_equal(min(path), lam)
  expect_true(!is.unsorted(rev(path)))
  fit <- wafc(design = d, y = y0, lambda = path)
  expect_equal(fit[["nzero"]][1L], 0L)
  ## the entry point is the level at which the first coefficient enters: a
  ## shade below it, the fit is no longer empty
  top <- wafc_lambda_max(d, y0)
  expect_gt(wafc(design = d, y = y0,
                 lambda = c(top, 0.99 * top))[["nzero"]][2L], 0L)
  expect_error(wafc_path_to(-1, d, y0), "positive")
})

test_that("wafc_tune validates its arguments", {
  ## the message of match.arg() is translated, so only the failure is checked
  expect_error(wafc_tune(x0, u0, y0, rule = "aic"))
  expect_error(wafc_tune(x0, u0, y0, rule = "ebic", gamma = -1, J = 2L),
               "in \\[0, 1\\]")
  expect_error(wafc_tune(x0, u0, y0, rule = "theory", s = -1), "positive")
  expect_error(wafc_tune(x0, u0, y0, rule = "bic", J = 0L), "larger than j0")
})
