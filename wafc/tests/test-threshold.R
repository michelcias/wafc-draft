## wafc/tests/test-threshold.R -- tests of the estimation followed by a
## threshold of step E2.5g (wafc/R/threshold.R). Run from the root of the
## repository:
##
##     Rscript -e 'testthat::test_dir("wafc/tests")'
##
## Four tests carry the step, the ones its catalogue entry names. A
## threshold of zero returns the fit it was given, to 1e-12, for the WAFC
## and for the block LASSO; a threshold above the largest norm zeroes every
## block and leaves the levels alone; the sandwich of Lemma 13,
## {N > t + D} in S-hat(t) in {N > t - D}, holds on a case built by hand and
## on a fit to a truth in the basis whose norms are known; and the rule
## "cv" is reproducible with fixed folds. A fifth ties the folds of the rule
## "cv" to the cross-validation that chose the fit: at a threshold of zero
## its error is the one cv.wafc() reported at lambda.min. The others check
## the rule "max" against its definition, the oracle against its own
## candidates and the truth, and the two refits against lm().
##
## Step E2.5j adds the tests its catalogue entry names: the rules "cvrel"
## and "cv1se" against the search redone by hand on the fold fits, for the
## WAFC and for the block LASSO; "cvrel" invariant to the scale of the
## response (y times 10 keeps the same blocks); and the gate of the QUT,
## which zeroes a pure null at a large n and does not fire on a strong
## signal, with its decision checked against the QUT penalty level and its
## seed leaving the random stream of the caller where it was.

library(testthat)

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

n <- 250L
p <- 3L
q <- 2L
dgp <- simulate_wafc(n, p = p, q = q, scenario = "smooth", seed = 20261001L,
                     snr = 4)
x0 <- dgp[["x"]]
u0 <- dgp[["u"]]
y0 <- dgp[["y"]]
active0 <- nzchar(dgp[["structure"]])
folds <- rep_len(1:5, n)
test0 <- simulate_wafc(400L, p = p, q = q, scenario = "smooth",
                       seed = 20261002L, snr = 4)
grid0 <- cbind(seq(min(u0[, 1L]), max(u0[, 1L]), length.out = 64),
               seq(min(u0[, 2L]), max(u0[, 2L]), length.out = 64))

cv0 <- cv.wafc(x0, u0, y0, J = 2:4, foldid = folds)
f0 <- cv0[["wafc.fit"]]
s0 <- cv0[["lambda.min"]]

klopp0 <- if (has("grpreg")) {
  wafc_competitor("klopp", x0, u0, y0, J = 3:4, foldid = folds,
                  penalize.levels = FALSE, balanced = TRUE)
} else NULL

## ---------------------------------------------------------------------------
## The two ends of the threshold
## ---------------------------------------------------------------------------

test_that("t = 0 returns the WAFC fit it was given, to 1e-12", {
  th <- wafc_threshold(cv0, t = 0)
  expect_s3_class(th, "wafc_threshold")
  expect_s3_class(th, "wafc_competitor")
  expect_identical(th[["extra"]][["rule"]], "fixed")
  expect_equal(th[["extra"]][["lambda"]], s0)
  expect_equal(th[["extra"]][["J"]], cv0[["J.min"]])
  ## prediction, beta(u), components and levels against the fit
  d_t <- wafc_design(test0[["x"]], test0[["u"]], spec = f0[["design"]])
  cf <- wafc_raw_coef(f0, s = s0)
  fh <- as.numeric(d_t[["Z"]] %*% cf[-1L, 1L]) + cf[1L, 1L]
  expect_equal(predict(th, test0[["x"]], test0[["u"]]), fh, tolerance = 1e-12)
  expect_equal(unname(th[["beta"]](test0[["u"]])),
               unname(predict(f0, newu = test0[["u"]], s = s0, type = "beta")),
               tolerance = 1e-12)
  ga <- wafc_grid_components(th, grid0)
  gb <- wafc_grid_components(f0, grid0, s = s0)
  for (i in seq_along(ga)) expect_equal(ga[[i]], gb[[i]], tolerance = 1e-12)
  expect_equal(unname(th[["cc"]]), unname(coef(f0, s = s0)[1L + f0[["design"]][["unpenalized"]], 1L]),
               tolerance = 1e-12)
  blk <- wafc_blocks(f0, s = s0)
  expect_equal(unname(th[["blocks"]]), unname(blk[["nonzero"]] > 0L))
  expect_equal(unname(th[["extra"]][["norm"]]), unname(blk[["norm"]]),
               tolerance = 1e-12)
  expect_equal(th[["fitted"]],
               as.numeric(f0[["design"]][["Z"]] %*% cf[-1L, 1L]) + cf[1L, 1L],
               tolerance = 1e-12)
})

test_that("t = 0 returns the block LASSO fit it was given, to 1e-12", {
  skip_if_not(has("grpreg"))
  th <- wafc_threshold(klopp0, t = 0)
  expect_equal(predict(th, test0[["x"]], test0[["u"]]),
               predict(klopp0, test0[["x"]], test0[["u"]]), tolerance = 1e-12)
  expect_equal(th[["beta"]](test0[["u"]]), klopp0[["beta"]](test0[["u"]]),
               tolerance = 1e-12)
  ga <- wafc_grid_components(th, grid0)
  gb <- wafc_grid_components(klopp0, grid0)
  for (i in seq_along(ga)) expect_equal(ga[[i]], gb[[i]], tolerance = 1e-12)
  expect_equal(th[["blocks"]], klopp0[["blocks"]])
  expect_equal(th[["fitted"]], predict(klopp0), tolerance = 1e-12)
  expect_equal(th[["extra"]][["J"]], klopp0[["extra"]][["J"]])
  expect_match(th[["method"]], "^klopp\\+fixed$")
})

test_that("a t above the largest norm zeroes every block and keeps the levels", {
  th0 <- wafc_threshold(cv0, t = 0)
  top <- max(th0[["extra"]][["norm"]])
  for (tt in c(top, 2 * top)) {
    th <- wafc_threshold(cv0, t = tt)
    expect_false(any(th[["blocks"]]))
    expect_identical(th[["extra"]][["nzero"]], 0L)
    g <- th[["g"]](grid0)
    for (i in seq_along(g)) expect_true(all(g[[i]] == 0))
    expect_equal(th[["cc"]], th0[["cc"]], tolerance = 1e-12)
  }
  skip_if_not(has("grpreg"))
  thk <- wafc_threshold(klopp0, t = 1e6)
  expect_false(any(thk[["blocks"]]))
  expect_equal(thk[["cc"]], klopp0[["cc"]], tolerance = 1e-12)
})

## ---------------------------------------------------------------------------
## Lemma 13
## ---------------------------------------------------------------------------

test_that("the sandwich of Lemma 13 holds on a case built by hand", {
  d <- wafc_design(x0, u0, J = 3L, rescale = FALSE)
  ## norms 0.9, 0.5, 0.05 and 0 on four of the six blocks, 0.3 and 0.2 on
  ## the other two, spread over the columns so that no single coefficient
  ## is the block
  target <- matrix(c(0.9, 0.5, 0.05, 0, 0.3, 0.2), p, q)
  b <- numeric(d[["nvars"]])
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      idx <- d[["blocks"]][[wafc_block_name(d, l, m)]]
      v <- seq_along(idx) %% 3 - 1
      if (sum(v^2) > 0 && target[l, m] > 0) {
        b[idx] <- v / sqrt(sum(v^2)) * target[l, m]
      }
    }
  }
  nrm <- wafc_thr_norms(b, d)
  expect_equal(unname(nrm), target, tolerance = 1e-14)
  ## D = 0: the estimate is the truth, and S-hat(t) = {N > t} exactly
  for (tt in c(0, 0.04, 0.05, 0.2, 0.25, 0.5, 0.89, 0.9)) {
    kept <- wafc_thr_apply(0, b, d, nrm > tt)
    expect_equal(unname(wafc_thr_norms(kept[["b"]], d) > 0), target > tt)
  }
  ## a perturbation of norm at most D in every block moves S-hat(t) only
  ## inside the sandwich
  set.seed(1)
  lower <- upper <- logical(0)
  for (r in 1:200) {
    D <- runif(1, 0, 0.3)
    nh <- pmax(target + matrix(runif(p * q, -D, D), p, q), 0)
    Dr <- max(abs(nh - target))
    for (tt in seq(0, 1, by = 0.05)) {
      sh <- nh > tt
      lower <- c(lower, all(sh[target > tt + Dr]))
      upper <- c(upper, all(target[sh] > tt - Dr))
    }
  }
  expect_true(all(lower))
  expect_true(all(upper))
})

test_that("the sandwich of Lemma 13 holds on a fit to a truth in the basis", {
  d <- wafc_design(x0, u0, J = 3L, rescale = FALSE)
  theta <- numeric(d[["nvars"]])
  theta[d[["unpenalized"]]] <- c(1, 2, -1.5)
  theta[d[["blocks"]][["x1:u1"]][c(1L, 4L)]] <- c(1.5, -1)
  theta[d[["blocks"]][["x2:u2"]][2L]] <- 0.8
  theta[d[["blocks"]][["x1:u2"]][3L]] <- 0.4
  N <- wafc_thr_norms(theta, d)
  set.seed(2)
  y <- as.numeric(d[["Z"]] %*% theta) + rnorm(n, sd = 0.5)
  cv <- cv.wafc(x0, u0, y, J = 3L, foldid = folds, rescale = FALSE)
  th0 <- wafc_threshold(cv, t = 0)
  Nh <- th0[["extra"]][["norm"]]
  D <- max(abs(Nh - N))
  lower <- upper <- logical(0)
  for (tt in seq(0, 2, by = 0.01)) {
    sh <- wafc_threshold(cv, t = tt)[["blocks"]]
    lower <- c(lower, all(sh[N > tt + D]))
    upper <- c(upper, all(N[sh] > tt - D))
  }
  expect_true(all(lower))
  expect_true(all(upper))
  ## part (iii): when the separation exceeds 2D, every t in [D, delta - D)
  ## recovers the structure
  delta <- min(N[N > 0])
  if (delta > 2 * D) {
    for (tt in seq(D, delta - D, length.out = 7)[-7]) {
      expect_equal(wafc_threshold(cv, t = tt)[["blocks"]], N > 0)
    }
  }
})

## ---------------------------------------------------------------------------
## The rules
## ---------------------------------------------------------------------------

test_that("rule = \"max\" thresholds at c times the largest norm", {
  th <- wafc_threshold(cv0, rule = "max")
  nrm <- th[["extra"]][["norm"]]
  expect_equal(th[["extra"]][["t"]], 0.15 * max(nrm))
  expect_equal(th[["extra"]][["c"]], 0.15)
  expect_equal(th[["blocks"]], nrm > 0.15 * max(nrm))
  th2 <- wafc_threshold(cv0, rule = "max", c = 0.5)
  expect_equal(th2[["blocks"]], nrm > 0.5 * max(nrm))
  expect_identical(th2[["method"]], "wafc.lasso+max")
  expect_error(wafc_threshold(cv0, rule = "max", c = -1), "'c'")
})

test_that("rule = \"cv\" is reproducible with fixed folds", {
  a <- wafc_threshold(cv0, rule = "cv")
  b <- wafc_threshold(cv0, rule = "cv")
  expect_identical(a[["extra"]][["t"]], b[["extra"]][["t"]])
  expect_identical(a[["extra"]][["candidates"]], b[["extra"]][["candidates"]])
  expect_identical(a[["blocks"]], b[["blocks"]])
  ## the same folds given by hand to the "wafc" object, and the fold fits
  ## computed once and passed, give the same answer
  ff <- wafc_threshold_folds(f0, s = s0, foldid = folds)
  c1 <- wafc_threshold(f0, rule = "cv", s = s0, foldid = folds)
  c2 <- wafc_threshold(cv0, rule = "cv", fold.fits = ff)
  expect_identical(c1[["extra"]][["candidates"]], a[["extra"]][["candidates"]])
  expect_identical(c2[["extra"]][["candidates"]], a[["extra"]][["candidates"]])
  ## the chosen t is inside the best interval, or the largest norm when
  ## that interval zeroes everything
  cand <- a[["extra"]][["candidates"]]
  i <- which.min(cand[["error"]])
  tt <- a[["extra"]][["t"]]
  expect_true(tt >= cand[["lower"]][i])
  expect_true(tt < cand[["upper"]][i] || is.infinite(cand[["upper"]][i]))
  ## the cv rule needs folds
  expect_error(wafc_threshold(f0, rule = "cv", s = s0), "folds")
})

test_that("the folds of rule = \"cv\" are the ones cv.wafc() scored", {
  ## at the interval of t = 0 nothing is zeroed, so the error is the
  ## cross-validated error of the fit at lambda.min, mean of the fold means
  a <- wafc_threshold(cv0, rule = "cv")
  cand <- a[["extra"]][["candidates"]]
  expect_equal(cand[["lower"]][1L], 0)
  expect_equal(cand[["error"]][1L], cv0[["cvm.min"]], tolerance = 1e-10)
  ## and the candidates are exhaustive: the breakpoints are the norms of
  ## the fold fits
  ff <- wafc_threshold_folds(cv0)
  br <- sort(unique(c(0, unlist(lapply(ff[["fits"]], function(f)
    wafc_thr_norms(f[["b"]], f0[["design"]]))))))
  expect_equal(cand[["lower"]], br)
})

test_that("rule = \"cv\" works on the block LASSO and is reproducible", {
  skip_if_not(has("grpreg"))
  a <- wafc_threshold(klopp0, rule = "cv")
  b <- wafc_threshold(klopp0, rule = "cv", foldid = folds, y = y0)
  ## the response grpreg stores is the one given up to rounding
  expect_equal(a[["extra"]][["candidates"]], b[["extra"]][["candidates"]],
               tolerance = 1e-12)
  expect_equal(a[["extra"]][["t"]], b[["extra"]][["t"]], tolerance = 1e-12)
  expect_identical(wafc_threshold(klopp0, rule = "cv")[["extra"]],
                   a[["extra"]])
  expect_identical(a[["method"]], "klopp+cv")
  ## the fold fit is the grpreg fit of the fold at the selected lambda
  ff <- wafc_threshold_folds(klopp0)
  Z <- as.matrix(klopp0[["design"]][["Z"]])
  ex <- klopp0[["extra"]]
  grp <- wafc_kp_groups(klopp0[["design"]], ex[["block.size"]],
                        ex[["penalize.levels"]], ex[["merge.coarse"]],
                        ex[["balanced"]], ex[["free.coarse"]])
  rows <- which(folds != 1L)
  fi <- grpreg::grpreg(Z[rows, ], y0[rows], group = grp, penalty = "grLasso",
                       lambda = klopp0[["fit"]][["lambda"]])
  v <- as.numeric(stats::coef(fi, lambda = ex[["lambda"]]))
  expect_equal(c(ff[["fits"]][[1L]][["a0"]], ff[["fits"]][[1L]][["b"]]), v,
               tolerance = 1e-8)
})

test_that("rule = \"oracle\" picks the candidate closest to the truth", {
  truth <- list(x = test0[["x"]], u = test0[["u"]], f = test0[["f"]])
  th <- wafc_threshold(cv0, rule = "oracle", truth = truth)
  cand <- th[["extra"]][["candidates"]]
  err <- sqrt(mean((predict(th, truth[["x"]], truth[["u"]]) - truth[["f"]])^2))
  expect_equal(err, min(cand[["error"]]), tolerance = 1e-10)
  ## every candidate is scored by the threshold at its lower end
  for (i in seq_len(nrow(cand))) {
    e <- wafc_threshold(cv0, t = cand[["lower"]][i])
    expect_equal(sqrt(mean((predict(e, truth[["x"]], truth[["u"]]) -
                              truth[["f"]])^2)),
                 cand[["error"]][i], tolerance = 1e-10)
  }
  expect_lte(err, sqrt(mean((predict(wafc_threshold(cv0, t = 0), truth[["x"]],
                                     truth[["u"]]) - truth[["f"]])^2)))
  expect_error(wafc_threshold(cv0, rule = "oracle"), "truth")
})

## ---------------------------------------------------------------------------
## The refits
## ---------------------------------------------------------------------------

test_that("the two refits are least squares on what is kept", {
  th <- wafc_threshold(cv0, rule = "max")
  des <- f0[["design"]]
  b <- wafc_raw_coef(f0, s = s0)[-1L, 1L]
  kept <- th[["extra"]][["kept"]]
  for (rf in c("support", "block")) {
    r <- wafc_threshold(cv0, rule = "max", refit = rf)
    expect_equal(r[["blocks"]] | !kept, matrix(TRUE, p, q,
                                               dimnames = dimnames(kept)))
    cols <- setdiff(des[["unpenalized"]], des[["constant"]])
    for (l in seq_len(p)) {
      for (m in seq_len(q)) {
        if (!kept[l, m]) next
        idx <- des[["blocks"]][[wafc_block_name(des, l, m)]]
        cols <- c(cols, if (rf == "support") idx[b[idx] != 0] else idx)
      }
    }
    W <- as.matrix(des[["Z"]][, cols])
    ref <- stats::lm(y0 ~ W)
    expect_equal(r[["fitted"]], unname(stats::fitted(ref)), tolerance = 1e-8)
    expect_identical(r[["method"]], paste0("wafc.lasso+max+", rf))
    ## the residual is orthogonal to the columns refitted
    expect_lt(max(abs(crossprod(cbind(1, W), y0 - r[["fitted"]]))), 1e-8)
  }
  ## the support refit on the fit at t = 0 is least squares on its support
  r0 <- wafc_threshold(cv0, t = 0, refit = "support")
  sup <- c(setdiff(des[["unpenalized"]], des[["constant"]]),
           setdiff(which(b != 0), des[["unpenalized"]]))
  W <- as.matrix(des[["Z"]][, sup])
  expect_equal(r0[["fitted"]], unname(stats::fitted(stats::lm(y0 ~ W))),
               tolerance = 1e-8)
})

test_that("the rules score the refitted estimator when there is a refit", {
  truth <- list(x = test0[["x"]], u = test0[["u"]], f = test0[["f"]])
  th <- wafc_threshold(cv0, rule = "oracle", refit = "support", truth = truth)
  err <- sqrt(mean((predict(th, truth[["x"]], truth[["u"]]) - truth[["f"]])^2))
  expect_equal(err, min(th[["extra"]][["candidates"]][["error"]]),
               tolerance = 1e-10)
  a <- wafc_threshold(cv0, rule = "cv", refit = "block")
  b <- wafc_threshold(cv0, rule = "cv", refit = "block")
  expect_identical(a[["extra"]][["candidates"]], b[["extra"]][["candidates"]])
  expect_false(isTRUE(all.equal(
    a[["extra"]][["candidates"]][["error"]],
    wafc_threshold(cv0, rule = "cv")[["extra"]][["candidates"]][["error"]])))
})

test_that("the block LASSO refit keeps grpreg's intercept convention", {
  skip_if_not(has("grpreg"))
  r <- wafc_threshold(klopp0, rule = "max", refit = "block")
  expect_equal(predict(r), r[["fitted"]], tolerance = 1e-10)
  expect_equal(predict(r, x0, u0), r[["fitted"]], tolerance = 1e-10)
})

## ---------------------------------------------------------------------------
## Step E2.5j: the rules "cvrel" and "cv1se", and the gate
## ---------------------------------------------------------------------------

## The search of the two rules redone by hand, from the fold fits and with
## nothing of threshold.R but the norms: every candidate is a breakpoint,
## each scored by the mean over folds of the squared error of the fold fit
## thresholded there on the observations the fold left out.
hand_search <- function(ff, des, y, relative) {
  fid <- ff[["foldid"]]
  K <- length(ff[["fits"]])
  nk <- lapply(ff[["fits"]], function(f) wafc_thr_norms(f[["b"]], des))
  if (relative) nk <- lapply(nk, function(v) if (max(v) > 0) v / max(v) else v)
  cand <- sort(unique(c(0, unlist(nk))))
  E <- matrix(NA_real_, K, length(cand))
  for (k in seq_len(K)) {
    out <- fid == k
    f <- ff[["fits"]][[k]]
    for (i in seq_along(cand)) {
      b <- f[["b"]]
      for (l in seq_len(des[["p"]])) {
        for (m in seq_len(des[["q"]])) {
          if (nk[[k]][l, m] > cand[i]) next
          b[des[["blocks"]][[wafc_block_name(des, l, m)]]] <- 0
        }
      }
      eta <- as.numeric(des[["Z"]][out, , drop = FALSE] %*% b) + f[["a0"]]
      E[k, i] <- mean((y[out] - eta)^2)
    }
  }
  list(cand = cand, err = colMeans(E), se = apply(E, 2L, sd) / sqrt(K))
}
hand_mid <- function(cand, i, top) {
  if (i == length(cand)) top else (cand[i] + cand[i + 1L]) / 2
}

test_that("rule = \"cvrel\" is the search over c redone by hand", {
  ff <- wafc_threshold_folds(cv0)
  des <- f0[["design"]]
  h <- hand_search(ff, des, y0, relative = TRUE)
  i <- which.min(h[["err"]])
  th <- wafc_threshold(cv0, rule = "cvrel", fold.fits = ff)
  nrm <- th[["extra"]][["norm"]]
  chat <- hand_mid(h[["cand"]], i, 1)
  expect_equal(th[["extra"]][["c"]], chat, tolerance = 1e-12)
  expect_equal(th[["extra"]][["t"]], chat * max(nrm), tolerance = 1e-12)
  expect_equal(unname(th[["blocks"]]), unname(nrm > chat * max(nrm)))
  cand <- th[["extra"]][["candidates"]]
  expect_equal(cand[["lower"]], h[["cand"]])
  expect_equal(cand[["error"]], h[["err"]], tolerance = 1e-10)
  ## every candidate is a fraction of the largest norm of some fold, and the
  ## last one, 1, zeroes every block
  expect_equal(max(cand[["lower"]]), 1)
  expect_true(all(cand[["lower"]] >= 0 & cand[["lower"]] <= 1))
  expect_identical(th[["method"]], "wafc.lasso+cvrel")
  ## the folds fitted inside give the same answer as the ones passed
  expect_identical(wafc_threshold(cv0, rule = "cvrel")[["extra"]][["t"]],
                   th[["extra"]][["t"]])
})

test_that("rule = \"cv1se\" is the one-standard-error search redone by hand", {
  ff <- wafc_threshold_folds(cv0)
  des <- f0[["design"]]
  h <- hand_search(ff, des, y0, relative = FALSE)
  imin <- which.min(h[["err"]])
  i <- max(which(h[["err"]] <= h[["err"]][imin] + h[["se"]][imin]))
  th <- wafc_threshold(cv0, rule = "cv1se", fold.fits = ff)
  nrm <- th[["extra"]][["norm"]]
  tt <- hand_mid(h[["cand"]], i, max(c(nrm, h[["cand"]])))
  expect_equal(th[["extra"]][["t"]], tt, tolerance = 1e-12)
  expect_equal(unname(th[["blocks"]]), unname(nrm > tt))
  cand <- th[["extra"]][["candidates"]]
  expect_equal(cand[["se"]], h[["se"]], tolerance = 1e-10)
  ## the scores are the ones of "cv"; only the pick moves, towards the
  ## larger t, so it never keeps more blocks than "cv"
  cv <- wafc_threshold(cv0, rule = "cv", fold.fits = ff)
  expect_identical(cand, cv[["extra"]][["candidates"]])
  expect_gte(th[["extra"]][["t"]], cv[["extra"]][["t"]])
  expect_true(all(cv[["blocks"]] | !th[["blocks"]]))
  expect_identical(th[["method"]], "wafc.lasso+cv1se")
})

test_that("the two rules redone by hand on the block LASSO", {
  skip_if_not(has("grpreg"))
  ff <- wafc_threshold_folds(klopp0)
  des <- klopp0[["design"]]
  for (rel in c(TRUE, FALSE)) {
    h <- hand_search(ff, des, y0, relative = rel)
    th <- wafc_threshold(klopp0, rule = if (rel) "cvrel" else "cv1se",
                         fold.fits = ff)
    nrm <- th[["extra"]][["norm"]]
    imin <- which.min(h[["err"]])
    if (rel) {
      tt <- hand_mid(h[["cand"]], imin, 1) * max(nrm)
    } else {
      i <- max(which(h[["err"]] <= h[["err"]][imin] + h[["se"]][imin]))
      tt <- hand_mid(h[["cand"]], i, max(c(nrm, h[["cand"]])))
    }
    expect_equal(th[["extra"]][["t"]], tt, tolerance = 1e-10)
    expect_equal(th[["blocks"]], nrm > tt)
  }
})

test_that("rule = \"cvrel\" is invariant to the scale of the response", {
  ## the same J and the same folds; the path of the LASSO scales with y, so
  ## the fits at the same point of it scale too, and the fractions c do not
  cv1 <- cv.wafc(x0, u0, y0, J = 3L, foldid = folds)
  cv10 <- cv.wafc(x0, u0, 10 * y0, J = 3L, foldid = folds)
  expect_equal(cv10[["lambda.min"]], 10 * cv1[["lambda.min"]],
               tolerance = 1e-6)
  a <- wafc_threshold(cv1, rule = "cvrel")
  b <- wafc_threshold(cv10, rule = "cvrel")
  expect_identical(b[["blocks"]], a[["blocks"]])
  expect_equal(b[["extra"]][["c"]], a[["extra"]][["c"]], tolerance = 1e-6)
  expect_equal(b[["extra"]][["t"]], 10 * a[["extra"]][["t"]],
               tolerance = 1e-6)
  skip_if_not(has("grpreg"))
  k1 <- wafc_competitor("klopp", x0, u0, y0, J = 3L, foldid = folds,
                        penalize.levels = FALSE, balanced = TRUE)
  k10 <- wafc_competitor("klopp", x0, u0, 10 * y0, J = 3L, foldid = folds,
                         penalize.levels = FALSE, balanced = TRUE)
  expect_identical(wafc_threshold(k10, rule = "cvrel")[["blocks"]],
                   wafc_threshold(k1, rule = "cvrel")[["blocks"]])
})

test_that("the gate zeroes a pure null at a large n", {
  nb <- 1000L
  d <- simulate_wafc(nb, p = p, q = q, scenario = "null", seed = 20261003L,
                     sigma = 0.62)
  cv <- cv.wafc(d[["x"]], d[["u"]], d[["y"]], J = 3:4,
                foldid = rep_len(1:5, nb))
  g <- wafc_threshold_gate(cv, seed = 1L)
  expect_false(g[["reject"]])
  expect_lt(g[["stat"]], g[["quantile"]])
  for (rule in c("max", "cv", "cvrel", "cv1se")) {
    th <- wafc_threshold(cv, rule = rule, gate = "qut", gate.seed = 1L)
    expect_false(any(th[["blocks"]]))
    expect_identical(th[["extra"]][["nzero"]], 0L)
    expect_identical(th[["method"]], paste0("wafc.lasso+", rule, "+qut"))
    ## the levels are the ones of the fit
    expect_equal(th[["cc"]], wafc_threshold(cv, t = 0)[["cc"]],
                 tolerance = 1e-12)
  }
  ## over 40 null responses on the same design the test rejects about
  ## alpha of the time (P(more than 8 of 40) < 1e-3 at alpha = 0.05)
  des <- cv[["wafc.fit"]][["design"]]
  rej <- vapply(1:40, function(i) {
    set.seed(500L + i)
    yy <- 1 + stats::rnorm(nb, sd = 0.62)
    wafc_thr_gate(list(design = des, y = yy, J = des[["J"]][1L]),
                  alpha = 0.05, nsim = 200L, seed = i)[["reject"]]
  }, TRUE)
  expect_lte(sum(rej), 8L)
})

test_that("the gate does not fire on a strong signal, and leaves the rule", {
  g <- wafc_threshold_gate(cv0, seed = 1L)
  expect_true(g[["reject"]])
  expect_identical(g[["J"]], cv0[["J.min"]])
  for (rule in c("max", "cv", "cvrel", "cv1se")) {
    a <- wafc_threshold(cv0, rule = rule)
    b <- wafc_threshold(cv0, rule = rule, gate = "qut", gate.test = g)
    expect_identical(b[["blocks"]], a[["blocks"]])
    expect_identical(b[["extra"]][["t"]], a[["extra"]][["t"]])
    expect_identical(b[["extra"]][["gate"]], g)
  }
  skip_if_not(has("grpreg"))
  gk <- wafc_threshold_gate(klopp0, seed = 1L)
  expect_true(gk[["reject"]])
  ## the statistic is coordinatewise and does not depend on the fit, so on
  ## the same design and response the two fits are tested alike
  if (identical(klopp0[["extra"]][["J"]], cv0[["J.min"]])) {
    expect_equal(gk[["stat"]], g[["stat"]], tolerance = 1e-12)
    expect_equal(gk[["quantile"]], g[["quantile"]], tolerance = 1e-12)
  }
})

test_that("the gate rejects exactly when the QUT leaves a coefficient", {
  des <- f0[["design"]]
  for (i in 1:6) {
    set.seed(700L + i)
    yy <- if (i <= 3) y0 else 1 + stats::rnorm(n, sd = 0.62)
    g <- wafc_thr_gate(list(design = des, y = yy, J = des[["J"]][1L]),
                       alpha = 0.05, nsim = 100L, seed = i)
    lam <- wafc_lambda_qut(des, yy, alpha = 0.05, nsim = 100L, seed = i)
    expect_identical(g[["reject"]], wafc_lambda_max(des, yy) > lam)
  }
})

test_that("the seed of the gate leaves the random stream where it was", {
  set.seed(42L)
  before <- .Random.seed
  wafc_threshold_gate(cv0, seed = 7L)
  expect_identical(.Random.seed, before)
  set.seed(42L)
  a <- stats::runif(3)
  set.seed(42L)
  wafc_threshold(cv0, rule = "max", gate = "qut", gate.seed = 7L)
  expect_identical(stats::runif(3), a)
  ## and the same seed gives the same quantile
  expect_identical(wafc_threshold_gate(cv0, seed = 7L)[["quantile"]],
                   wafc_threshold_gate(cv0, seed = 7L)[["quantile"]])
})
