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
  cv <- cv.wafc(x0, u0, y0, J = 2:4, foldid = folds,
                penalty = "lasso", threshold = "none")
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
  a <- cv.wafc(x0, u0, y0, J = 2:3, foldid = folds,
               penalty = "lasso", threshold = "none")
  b <- cv.wafc(x0, u0, y0, J = 2:3, foldid = folds,
               penalty = "lasso", threshold = "none")
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
                lambda = c(0.1, 0.01, 1e-4, 1e-8),
                penalty = "lasso", threshold = "none")
  expect_equal(cv[["J.min"]], 3L)
  lam <- cv[["wafc.fit"]][["lambda"]]
  expect_equal(cv[["lambda.min"]], min(lam))
  expect_lt(cv[["cvm.min"]], 1e-6)
  cf <- coef(cv, s = "lambda.min")
  expect_equal(max(abs(unname(cf[-1L, 1L]) - theta)), 0, tolerance = 1e-3)
})

test_that("coef and predict of a cv.wafc read the selected pair", {
  cv <- cv.wafc(x0, u0, y0, J = 2:3, foldid = folds,
                penalty = "lasso", threshold = "none")
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
  expect_error(cv.wafc(x0, u0, y0, J = 2.5,
                       penalty = "lasso", threshold = "none"), "integer")
  expect_error(cv.wafc(x0, u0, y0, J = 0L,
                       penalty = "lasso", threshold = "none"), "larger than j0")
  expect_error(cv.wafc(x0, u0, y0, J = 2L, nfolds = 2L,
                       penalty = "lasso", threshold = "none"), "at least 3")
  expect_error(cv.wafc(x0, u0, y0, J = 2L, foldid = folds[-1L],
                       penalty = "lasso", threshold = "none"),
               "one entry per observation")
  expect_error(cv.wafc(x0, u0, y0, J = 2L, foldid = rep(1L, n),
                       penalty = "lasso", threshold = "none"),
               "at least 3 folds")
  expect_error(cv.wafc(x0, u0, y0[-1L], J = 2L,
                       penalty = "lasso", threshold = "none"), "length")
})

## ---------------------------------------------------------------------------
## The information criteria
## ---------------------------------------------------------------------------

test_that("wafc_bic and wafc_ebic are the criteria they claim to be", {
  fit <- wafc(x0, u0, y0, J = 3L, penalty = "lasso")
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
  fit <- wafc(x0, u0, y0, J = 3L, penalty = "lasso")
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
  rules <- c("cv.min", "cv.1se", "bic", "ebic", "theory", "qut", "gcv")
  for (rule in rules) {
    tn <- wafc_tune(x0, u0, y0, rule = rule, J = 2:4, foldid = folds,
                    s = 3/2, sigma = dgp[["sigma"]], nsim = 100L,
                    qut.seed = 5L)
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

test_that("the qut rule is the cross-validated J with the penalty level of the QUT", {
  skip_slow()
  ## Step E2.4b moved wafc_lambda_qut() here from competitors.R and made
  ## it the sixth rule. The resolution is not part of the rule, so it is
  ## the one the cross-validation of the other two rules selects, and only
  ## the penalty level changes: that is what the pilot did by hand.
  cv <- cv.wafc(x0, u0, y0, J = 2:4, foldid = folds,
                penalty = "lasso", threshold = "none")
  tn <- wafc_tune(x0, u0, y0, rule = "qut", J = 2:4, foldid = folds,
                  nsim = 200L, qut.seed = 11L)
  expect_equal(tn[["J"]], cv[["J.min"]])
  expect_equal(tn[["lambda"]],
               wafc_lambda_qut(cv[["wafc.fit"]][["design"]], y0, nsim = 200L,
                               seed = 11L))
  expect_equal(tn[["tab"]], cv[["cvtab"]])
  ## and it is the conservative rule the pilot reports: calibrated on the
  ## null, so never below the one calibrated on prediction
  expect_gt(tn[["lambda"]], cv[["lambda.min"]])
  expect_lte(tn[["nzero"]], wafc_nzero_at(cv[["wafc.fit"]], cv[["lambda.min"]]))
})

test_that("the dots of the tuning functions reach wafc() and not only wafc_design()", {
  ## Step E2.3 routed the whole '...' to wafc_design(), so an argument of
  ## the fit stopped with "unused argument" and the convergence threshold
  ## could not be loosened through the interface at all (docs/ESTADO.md,
  ## 2026-09-21). Both callees now get what is theirs, and a name that
  ## belongs to neither is an error instead of a silent default.
  a <- cv.wafc(x0, u0, y0, J = 2:3, foldid = folds, thresh = 1e-10,
               penalty = "lasso", threshold = "none")
  b <- cv.wafc(x0, u0, y0, J = 2:3, foldid = folds, rescale = FALSE,
               penalty = "lasso", threshold = "none")
  expect_s3_class(a, "cv.wafc")
  expect_equal(a[["J.min"]],
               cv.wafc(x0, u0, y0, J = 2:3, foldid = folds,
                       penalty = "lasso", threshold = "none")[["J.min"]])
  expect_false(b[["wafc.fit"]][["design"]][["rescale"]])
  expect_error(cv.wafc(x0, u0, y0, J = 2L, foldid = folds, thrsh = 1e-7,
                       penalty = "lasso", threshold = "none"),
               "unused argument")
  expect_error(wafc_tune(x0, u0, y0, rule = "bic", J = 2L, nope = 1),
               "unused argument")
  ## the two kinds of argument travel together, and to the right callee
  d <- wafc_tune(x0, u0, y0, rule = "bic", J = 2:3, rescale = FALSE,
                 thresh = 1e-10)
  expect_false(d[["fit"]][["design"]][["rescale"]])
  ## the tighter tolerance is what the optimality check wants, and it is
  ## reachable now
  expect_true(all(wafc_kkt(d[["fit"]], s = d[["lambda"]])[["ok"]]))
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
  cv <- cv.wafc(x0, u0, y0, J = 2:4, foldid = folds,
                penalty = "lasso", threshold = "none")
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
  fit <- wafc(design = d, y = y0, lambda = path, penalty = "lasso")
  expect_equal(fit[["nzero"]][1L], 0L)
  ## the entry point is the level at which the first coefficient enters: a
  ## shade below it, the fit is no longer empty
  top <- wafc_lambda_max(d, y0)
  expect_gt(wafc(design = d, y = y0,
                 lambda = c(top, 0.99 * top),
                 penalty = "lasso")[["nzero"]][2L], 0L)
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

test_that("cv.wafc keeps the fit of the selected level and only that one", {
  ## Only the best fit so far is held during the loop over J; the one
  ## returned has to be the fit at J.min, identical to fitting it directly.
  cv <- cv.wafc(x0, u0, y0, J = 2:4, foldid = folds,
                penalty = "lasso", threshold = "none")
  direct <- wafc(x0, u0, y0, J = cv[["J.min"]], penalty = "lasso")
  expect_equal(cv[["J.min"]], cv[["J"]][which.min(cv[["cvtab"]][["mse"]])])
  expect_equal(cv[["wafc.fit"]][["lambda"]], direct[["lambda"]])
  expect_equal(as.matrix(cv[["wafc.fit"]][["beta"]]),
               as.matrix(direct[["beta"]]))
  expect_equal(cv[["wafc.fit"]][["a0"]], direct[["a0"]])
})

test_that("wafc_lambda_max is the entry point of the engine with and without an intercept", {
  for (ic in c(TRUE, FALSE)) {
    d1 <- simulate_wafc(n, p = p, q = q, scenario = "smooth", seed = 4L,
                        intercept = ic)
    y1 <- d1[["y"]] + 5
    des <- wafc_design(d1[["x"]], d1[["u"]], J = 3L)
    fit <- wafc(design = des, y = y1, penalty = "lasso")
    expect_equal(fit[["intercept"]], ic)
    expect_equal(wafc_lambda_max(des, y1), fit[["lambda"]][1L],
                 tolerance = 1e-8)
  }
})

## ---------------------------------------------------------------------------
## Generalized cross-validation (step E2.5h)
## ---------------------------------------------------------------------------

test_that("wafc_gcv is n RSS / (n - df)^2 recomputed by hand", {
  fit <- wafc(x0, u0, y0, J = 3L, penalty = "lasso")
  gc <- wafc_gcv(fit)
  expect_equal(gc[["lambda"]], fit[["lambda"]])
  ## the degrees of freedom of wafc_bic(): the non-zero wavelet
  ## coefficients plus the p levels (decision D19)
  expect_equal(gc[["df"]], fit[["nzero"]] + p)
  rss <- unname(colSums((y0 - predict(fit, x0, u0))^2))
  expect_equal(unname(gc[["rss"]]), rss)
  expect_equal(gc[["gcv"]], n * rss / (n - (fit[["nzero"]] + p))^2)
  k <- attr(gc, "which.min")
  expect_equal(gc[["gcv"]][k], min(gc[["gcv"]][!gc[["excluded"]]]))
  expect_equal(attr(gc, "lambda.min"), gc[["lambda"]][k])
  ## and at a chosen penalty level
  g1 <- wafc_gcv(fit, s = fit[["lambda"]][10L])
  expect_equal(g1[["gcv"]], gc[["gcv"]][10L])
  expect_error(wafc_gcv(fit, guard = 0), "in \\(0, 1\\]")
})

test_that("the guard of the GCV leaves out the points with df >= n/2", {
  skip_slow()
  ## J = 6 has 378 penalized columns for n = 250, so the end of the path
  ## goes past n/2 and, for some points, past n
  fit <- wafc(x0, u0, y0, J = 6L, lambda.min.ratio = 1e-4, penalty = "lasso")
  gc <- wafc_gcv(fit)
  df <- fit[["nzero"]] + p
  expect_true(any(df >= n / 2))
  expect_identical(gc[["excluded"]], df >= n / 2)
  expect_true(all(is.na(gc[["gcv"]][df >= n])))
  k <- attr(gc, "which.min")
  expect_false(gc[["excluded"]][k])
  ku <- attr(gc, "which.min.unguarded")
  expect_equal(ku, which.min(gc[["gcv"]]))
  ## guard = 1 keeps every point with a criterion, and is the unguarded
  ## choice
  g1 <- wafc_gcv(fit, guard = 1)
  expect_identical(g1[["excluded"]], df >= n)
  expect_equal(attr(g1, "which.min"), ku)
  ## the rule over the grid: the minimum over J of the guarded minima,
  ## the points left out counted, and the guard said to decide exactly when
  ## the unguarded minimum over the grid is another pair
  tn <- wafc_tune(x0, u0, y0, rule = "gcv", J = c(3L, 6L),
                  lambda.min.ratio = 1e-4)
  per <- lapply(c(3L, 6L), function(Ji) {
    wafc_gcv(wafc(x0, u0, y0, J = Ji, lambda.min.ratio = 1e-4,
                  penalty = "lasso"))
  })
  mins <- vapply(per, function(g) g[["gcv"]][attr(g, "which.min")], 0)
  free <- vapply(per, function(g) g[["gcv"]][attr(g, "which.min.unguarded")],
                 0)
  expect_equal(tn[["tab"]][["gcv"]], mins)
  expect_equal(tn[["tab"]][["gcv.unguarded"]], free)
  expect_equal(tn[["J"]], c(3L, 6L)[which.min(mins)])
  b <- per[[which.min(mins)]]
  expect_equal(tn[["lambda"]], attr(b, "lambda.min"))
  expect_equal(tn[["tab"]][["excluded"]],
               vapply(per, function(g) sum(g[["excluded"]]), 0L))
  expect_equal(tn[["guard"]][["excluded"]], sum(tn[["tab"]][["excluded"]]))
  expect_equal(tn[["guard"]][["points"]], sum(tn[["tab"]][["points"]]))
  jf <- c(3L, 6L)[which.min(free)]
  gf <- per[[which.min(free)]]
  expect_equal(tn[["guard"]][["J.unguarded"]], jf)
  expect_equal(tn[["guard"]][["lambda.unguarded"]],
               gf[["lambda"]][attr(gf, "which.min.unguarded")])
  expect_identical(tn[["guard"]][["decides"]],
                   !(jf == tn[["J"]] &&
                       tn[["guard"]][["lambda.unguarded"]] == tn[["lambda"]]))
})

test_that("with theta* in the basis and no noise, the GCV recovers it", {
  ## the residual sum of squares goes to zero at the small end of the path
  ## while df stays below n/2 (the 42 columns of J = 3 plus the p levels:
  ## at 1e-8 the LASSO leaves the null coefficients at rounding size, not at
  ## zero, and the count of D19 counts them); an explicit
  ## path, as in the test of cv.wafc above, because the engine stops its own
  ## path early on a noiseless response (and wafc_tune() has no 'lambda':
  ## the name would match 'lambda.min.ratio')
  fit <- wafc(x0, u0, y_exact, J = 3L, rescale = FALSE,
              lambda = c(0.1, 0.01, 1e-4, 1e-8), penalty = "lasso")
  gc <- wafc_gcv(fit)
  expect_equal(attr(gc, "lambda.min"), 1e-8)
  expect_lt(gc[["df"]][4L], n / 2)
  expect_equal(attr(gc, "which.min.unguarded"), attr(gc, "which.min"))
  cf <- coef(fit, s = attr(gc, "lambda.min"))
  expect_equal(max(abs(unname(cf[-1L, 1L]) - theta)), 0, tolerance = 1e-3)
})

## ---------------------------------------------------------------------------
## The convergence along the cross-validation (step E2.5j)
## ---------------------------------------------------------------------------

test_that("wafc_cv_convergence reads the fits of every J and fold", {
  cv <- cv.wafc(x0, u0, y0, J = 2:4, foldid = folds,
                penalty = "lasso", threshold = "none")
  cc <- wafc_cv_convergence(cv)
  expect_identical(cc[["J"]], 2:4)
  expect_identical(cc[["chosen"]], 2:4 == cv[["J.min"]])
  expect_identical(cc[["nlambda"]], rep(100L, 3L))
  expect_identical(cc[["nreturned"]],
                   vapply(cv[["cv"]], function(z) length(z[["lambda"]]), 0L))
  expect_identical(cc[["lambda.min"]],
                   vapply(cv[["cv"]], `[[`, 0, "lambda.min"))
  expect_true(all(cc[["jerr"]] == 0L))
  expect_false(any(cc[["cut"]] | cc[["cut.at.min"]]))
  expect_identical(cc[["folds.cut"]], rep(0L, 3L))
  expect_true(all(is.na(cc[["fold.lambda.cut"]])))
  ## the record changes nothing: the table of the cross-validation is the
  ## one of a fit whose records are removed
  expect_identical(cv[["cvtab"]][["lambda.min"]], cc[["lambda.min"]])
})

test_that("a cut in the folds is seen, and where it falls", {
  ## a budget of passes small enough for the folds to stop before the end of
  ## the path of the whole sample, which wafc_raw_coef() then reads at the
  ## last point a fold returned; the budget goes to the fold fits only, as
  ## wafc_cv_design() passes it
  des <- wafc_design(x0, u0, J = 3L)
  full <- wafc(design = des, y = y0, penalty = "lasso")
  z <- suppressWarnings(wafc_cv_design(des, y0, full, folds,
                                       function(e) e^2, "lasso", maxit = 200L))
  z[["J"]] <- 3L
  z[["conv"]] <- full[["conv"]]
  obj <- structure(list(cv = list(z), J.min = 3L), class = "cv.wafc")
  cc <- wafc_cv_convergence(obj)
  fo <- z[["conv.folds"]]
  nl <- length(full[["lambda"]])
  expect_identical(nrow(fo), 5L)
  cut <- fo[["nreturned"]] < nl | fo[["jerr"]] != 0L
  expect_true(all(cut))
  expect_identical(cc[["folds.cut"]], 5L)
  expect_true(all(fo[["jerr"]] < 0L))
  expect_identical(fo[["nreturned"]], -fo[["jerr"]] - 1L)
  expect_identical(cc[["fold.jerr"]], fo[["jerr"]][1L])
  expect_equal(cc[["fold.lambda.cut"]], max(fo[["lambda.last"]]))
  expect_identical(cc[["fold.cut.above.min"]],
                   cc[["fold.lambda.cut"]] > z[["lambda.min"]])
  expect_false(cc[["cut"]])
  ## a cut of the fit on the whole sample: jerr = -k, the path stops at
  ## k - 1, and the folds, asked for that shorter path, are not cut
  cv <- suppressWarnings(cv.wafc(x0, u0, y0, J = 3L, foldid = folds,
                                 maxit = 200L,
                                 penalty = "lasso", threshold = "none"))
  c2 <- wafc_cv_convergence(cv)
  expect_true(c2[["cut"]])
  expect_identical(c2[["nreturned"]], -c2[["jerr"]] - 1L)
  expect_identical(c2[["cut.at.min"]],
                   c2[["lambda.min"]] == c2[["lambda.last"]])
})

test_that("wafc_lambda_qut is unchanged by the pivot it now shares", {
  ## a frozen copy of the function as step E2.4b left it
  old_qut <- function(design, y, alpha = 0.05, nsim = 200L, seed = NULL) {
    n <- design[["n"]]
    if (!is.null(seed)) set.seed(seed)
    unp <- design[["unpenalized"]]
    pen <- seq_len(design[["nvars"]])[-unp]
    W <- as.matrix(design[["Z"]][, unp, drop = FALSE])
    Zp <- design[["Z"]][, pen, drop = FALSE]
    qrW <- qr(W)
    resid_of <- function(v) as.numeric(v - W %*% wafc_qr_coef(qrW, v))
    E <- matrix(stats::rnorm(n * nsim), n, nsim)
    R <- E - W %*% wafc_qr_coef(qrW, E)
    G <- as.matrix(Matrix::crossprod(Zp, R))
    ratio <- apply(abs(G), 2L, max) / sqrt(colSums(R^2))
    r0 <- resid_of(y)
    as.numeric(stats::quantile(ratio, 1 - alpha, names = FALSE)) *
      sqrt(sum(r0^2)) / n
  }
  for (J in 2:4) {
    des <- wafc_design(x0, u0, J = J)
    for (sd in c(1L, 11L)) {
      expect_identical(wafc_lambda_qut(des, y0, nsim = 150L, seed = sd),
                       old_qut(des, y0, nsim = 150L, seed = sd))
    }
  }
})
