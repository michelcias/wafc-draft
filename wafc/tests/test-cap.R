## wafc/tests/test-cap.R -- tests of step E3.4: the level of each modulating
## covariate capped by the number of distinct points at which its basis is
## evaluated, J_m = min(J, floor(log2(d_m))) (decision D57, lesson 3 of step
## E6.1b). Run from the root of the repository:
##
##     Rscript -e 'testthat::test_dir("wafc/tests")'
##
## Three things carry the step. The rule itself, with d_m counted on the
## circle of the periodized basis, where the smallest and the largest value
## of a covariate are one point when eps = 0. The proof that nothing moved
## where it should not: on a continuous design with 2^J < n the cap does not
## act, and every fit with and without it agrees to 1e-12. And the reason
## for the step: with the hour of the day (24 values) a block of level 6 has
## 63 columns on 23 points, most of its norm lies along directions no
## observation sees, and with the cap the block is identified and its norm
## is a function of the fitted values.

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

## A continuous design: n - 1 = 259 points on the circle, so the cap does
## not act below J = 9.
dgp <- simulate_wafc(260L, p = 3L, q = 2L, scenario = "smooth",
                     seed = 20261005L, snr = 4)
fold_c <- rep_len(1:5, 260L)

## The hour of the day, ten days of it, as the only modulator, and next to a
## continuous one.
n_h <- 240L
hour <- rep(0:23, length.out = n_h)
set.seed(20261005L)
x_h <- cbind(1, rnorm(n_h))
u_c <- runif(n_h)
u_ch <- cbind(u1 = u_c, u2 = hour)
y_h <- 1 + 0.8 * sin(2 * pi * hour / 24) +
  x_h[, 2] * (1 + cos(2 * pi * hour / 24)) + rnorm(n_h, sd = 0.3)
y_h2 <- y_h + 0.5 * cos(2 * pi * u_c)
fold_h <- rep_len(1:5, n_h)

## Rank of a matrix by its singular values, which, unlike qr() without
## pivoting on the norms, does not depend on the order of the columns.
svd_rank <- function(M) {
  d <- svd(as.matrix(M), nu = 0L, nv = 0L)[["d"]]
  sum(d > 1e-9 * d[1L])
}

## The coefficients of a block that its fitted contribution determines: the
## minimum-norm least squares of Z_b theta_b on Z_b, which is theta_b itself
## exactly when Z_b has no direction without data.
seen_part <- function(Zb, b) {
  sv <- svd(as.matrix(Zb))
  r <- sv[["d"]] > 1e-9 * sv[["d"]][1L]
  as.numeric(sv[["v"]][, r, drop = FALSE] %*%
               (crossprod(sv[["u"]][, r, drop = FALSE], Zb %*% b) /
                  sv[["d"]][r]))
}

test_that("the cap is the floor of log2 of the distinct points on the circle", {
  d <- wafc_design(dgp[["x"]], dgp[["u"]], J = 4)
  expect_identical(d[["ndistinct"]], c(259L, 259L))
  expect_identical(d[["J"]], c(4L, 4L))
  expect_identical(d[["J.requested"]], c(4L, 4L))
  expect_true(d[["cap.J"]])
  ## a continuous covariate is capped only when 2^J >= n: 130 values are
  ## 129 points, and 2^7 = 128 <= 129 < 2^8
  expect_identical(wafc_design(matrix(1, 130L), runif(130L), J = 8)[["J"]],
                   7L)

  ## the hour: with eps = 0 the hours 0 and 23 are the point 0 = 1 of the
  ## periodized basis, so 23 points; with a margin, 24; J = 4 in both
  xx <- x_h
  h0 <- wafc_design(xx, hour, J = 6)
  expect_identical(h0[["ndistinct"]], 23L)
  expect_identical(h0[["J"]], 4L)
  expect_identical(h0[["J.requested"]], 6L)
  h5 <- wafc_design(xx, hour, J = 6, eps = 0.05)
  expect_identical(h5[["ndistinct"]], 24L)
  expect_identical(h5[["J"]], 4L)

  ## the case the literal count of distinct values gets wrong: 16 values
  ## are 15 points at eps = 0, which identify J = 3 and not 4
  v16 <- rep(0:15, length.out = n_h)
  expect_identical(wafc_design(xx, v16, J = 6)[["J"]], 3L)
  expect_identical(wafc_design(xx, v16, J = 6, eps = 0.05)[["J"]], 4L)
  v17 <- rep(0:16, length.out = n_h)
  expect_identical(wafc_design(xx, v17, J = 6)[["J"]], 4L)

  ## the cap only lowers, one covariate at a time
  hc <- wafc_design(xx, u_ch, J = c(6, 2))
  expect_identical(hc[["J"]], c(6L, 2L))
  hc <- wafc_design(xx, u_ch, J = 6)
  expect_identical(hc[["J"]], c(6L, 4L))
  expect_identical(hc[["NJ"]], c(63L, 15L))
  expect_identical(lengths(hc[["blocks"]]),
                   c(`x1:u1` = 63L, `x1:u2` = 15L, `x2:u1` = 63L,
                     `x2:u2` = 15L))
  expect_identical(colnames(hc[["Z"]])[hc[["blocks"]][["x1:u2"]]],
                   wafc_colnames("x1", "u2", 0L, 4L))

  ## switched off, every block at the J asked for, the counts still recorded
  off <- wafc_design(xx, u_ch, J = 6, cap.J = FALSE)
  expect_identical(off[["J"]], c(6L, 6L))
  expect_identical(off[["ndistinct"]], hc[["ndistinct"]])
  expect_false(off[["cap.J"]])
})

test_that("with the cap a block of the hour is identified, without it it is not", {
  for (cap in c(TRUE, FALSE)) {
    d <- wafc_design(x_h, hour, J = 6, cap.J = cap)
    for (l in 1:2) {
      idx <- d[["blocks"]][[l]]
      W <- cbind(d[["Z"]][, l], d[["Z"]][, idx])
      ## [X_l, X_l psi(hour)] spans at most the 23 functions of the points
      expect_identical(svd_rank(W), if (cap) 16L else 23L)
      expect_identical(ncol(W), if (cap) 16L else 64L)
    }
  }
  ## 16 values: one direction without data at J = 4, none at the capped 3
  v16 <- rep(0:15, length.out = n_h)
  for (cap in c(TRUE, FALSE)) {
    d <- wafc_design(x_h, v16, J = 4, cap.J = cap)
    W <- cbind(1, d[["Z"]][, d[["blocks"]][[1L]]])
    expect_identical(ncol(W) - svd_rank(W), if (cap) 0L else 1L)
  }
})

test_that("the norm of a capped block does not depend on directions without data", {
  skip_no_grpreg()
  for (cap in c(TRUE, FALSE)) {
    f <- wafc(x_h, hour, y_h, J = 6, cap.J = cap)
    s <- f[["lambda"]][60L]
    b <- wafc_raw_coef(f, s = s)[-1L, 1L]
    nrm <- wafc_blocks(f, s = s)[["norm"]]
    d <- f[["design"]]
    for (l in 1:2) {
      idx <- d[["blocks"]][[l]]
      Zb <- as.matrix(d[["Z"]][, idx])
      seen <- seen_part(Zb, b[idx])
      off <- sqrt(sum((b[idx] - seen)^2))
      expect_equal(nrm[l, 1L], sqrt(sum(b[idx]^2)), tolerance = 1e-12)
      if (cap) {
        ## the coefficients are the ones the fitted values determine
        expect_lt(off, 1e-10 * nrm[l, 1L])
      } else {
        ## most of the norm lies where no observation is
        expect_gt(off, 0.5 * nrm[l, 1L])
      }
    }
    ## the optimality conditions hold on both designs
    expect_true(all(wafc_kkt(f, s = f[["lambda"]][c(20L, 60L)])[["ok"]]))
  }
})

test_that("on a continuous design the fits with and without the cap agree to 1e-12", {
  skip_no_grpreg()
  x0 <- dgp[["x"]]
  u0 <- dgp[["u"]]
  y0 <- dgp[["y"]]
  a <- wafc_design(x0, u0, J = 4)
  b <- wafc_design(x0, u0, J = 4, cap.J = FALSE)
  expect_identical(a[["Z"]], b[["Z"]])
  expect_identical(a[["blocks"]], b[["blocks"]])
  expect_identical(a[["J"]], b[["J"]])

  for (pen in c("block", "lasso")) {
    fa <- wafc(x0, u0, y0, J = 4, penalty = pen)
    fb <- wafc(x0, u0, y0, J = 4, penalty = pen, cap.J = FALSE)
    expect_equal(fa[["lambda"]], fb[["lambda"]], tolerance = 1e-12)
    expect_equal(as.matrix(fa[["beta"]]), as.matrix(fb[["beta"]]),
                 tolerance = 1e-12)
    expect_equal(fa[["cc"]], fb[["cc"]], tolerance = 1e-12)
  }

  ca <- cv.wafc(x0, u0, y0, J = 3, foldid = fold_c)
  cb <- cv.wafc(x0, u0, y0, J = 3, foldid = fold_c, cap.J = FALSE)
  expect_identical(ca[["J.min"]], cb[["J.min"]])
  expect_equal(ca[["cvtab"]], cb[["cvtab"]], tolerance = 1e-12)
  expect_equal(ca[["lambda.min"]], cb[["lambda.min"]], tolerance = 1e-12)
  expect_identical(ca[["J.eff"]], cb[["J.eff"]])
  expect_identical(ca[["threshold"]][["kept"]], cb[["threshold"]][["kept"]])
  expect_equal(ca[["threshold"]][["t"]], cb[["threshold"]][["t"]],
               tolerance = 1e-12)
  expect_equal(coef(ca), coef(cb), tolerance = 1e-12)
  expect_equal(predict(ca, x0, u0), predict(cb, x0, u0), tolerance = 1e-12)
  expect_null(wafc_cv_cap_text(ca))
})

test_that("cross-validation over J takes a capped candidate from the one it repeats", {
  skip_no_grpreg()
  cv <- cv.wafc(x_h, hour, y_h, J = 2:6, foldid = fold_h)
  expect_identical(cv[["J"]], 2:6)
  expect_identical(unname(cv[["J.eff"]][, 1L]), c(2L, 3L, 4L, 4L, 4L))
  for (i in 4:5) {
    expect_identical(cv[["cv"]][[i]][["cvm"]], cv[["cv"]][[3L]][["cvm"]])
    expect_identical(cv[["cv"]][[i]][["lambda"]], cv[["cv"]][[3L]][["lambda"]])
    expect_identical(cv[["cv"]][[i]][["J"]], cv[["J"]][i])
  }
  expect_true(cv[["J.min"]] <= 4L)
  expect_identical(wafc_cv_convergence(cv)[["J"]], 2:6)

  ## the candidate taken from J = 4 is what fitting J = 5 by itself gives
  c5 <- cv.wafc(x_h, hour, y_h, J = 5, foldid = fold_h, threshold = "none")
  expect_identical(c5[["J.eff"]][1L, 1L], 4L)
  expect_identical(c5[["cv"]][[1L]][["cvm"]], cv[["cv"]][[4L]][["cvm"]])
  expect_identical(c5[["cv"]][[1L]][["lambda"]], cv[["cv"]][[4L]][["lambda"]])

  ## the cap is said, and the J asked for is the one reported
  out <- capture.output(print(cv))
  expect_true(any(grepl(paste0("Levels capped by the distinct values: ",
                               "u1 at 4 \\(23 points\\) for J = 5, 6"),
                        out)))
  expect_identical(cv[["wafc.fit"]][["design"]][["J.requested"]],
                   cv[["J.min"]])

  ## switched off through '...': every candidate is its own design
  c6 <- cv.wafc(x_h, hour, y_h, J = 6, foldid = fold_h, threshold = "none",
                cap.J = FALSE)
  expect_identical(c6[["J.eff"]][1L, 1L], 6L)
  expect_identical(c6[["wafc.fit"]][["npen"]], 2L * 63L)
})

test_that("what reads the structure of the design reads the capped levels", {
  skip_no_grpreg()
  f <- wafc(x_h, u_ch, y_h2, J = 6)
  d <- f[["design"]]
  s <- f[["lambda"]][50L]
  expect_identical(d[["J"]], c(6L, 4L))

  ## the spec route: the design of new data, with fewer hours, keeps the
  ## levels of the training sample
  new_u <- cbind(runif(30), rep(0:4, 6))
  dn <- wafc_design(x_h[1:30, ], new_u, spec = d)
  expect_identical(dn[["J"]], d[["J"]])
  expect_identical(dn[["J.requested"]], d[["J.requested"]])
  expect_identical(colnames(dn[["Z"]]), colnames(d[["Z"]]))
  expect_equal(predict(f, x_h, u_ch, s = s), predict(f, s = s),
               tolerance = 1e-12)

  ## the component of the hour at the 24 hours is the fitted contribution
  ## of its block where X_1 = 1
  fn <- wafc_functions(f, s = s, grid = cbind(0.5, 0:23))
  b <- wafc_raw_coef(f, s = s)[-1L, 1L]
  idx <- d[["blocks"]][["x1:u2"]]
  rows <- match(0:23, hour)
  expect_equal(fn[["g"]][[1L, 2L]],
               as.numeric(d[["Z"]][rows, idx] %*% b[idx]), tolerance = 1e-10)
  expect_equal(wafc_blocks(f, s = s)[["norm"]][1L, 2L], sqrt(sum(b[idx]^2)),
               tolerance = 1e-12)

  ## the chunks of the block LASSO follow the 15 columns of the capped block
  grp <- f[["group"]][["group"]]
  expect_identical(length(unique(grp[idx])), 2L)
  expect_true(all(grp[idx] > 0L))

  ## print and the plots say it
  expect_true(any(grepl("levels capped by the distinct values of the modulators: u2 at 4 \\(23 points\\)",
                        capture.output(print(f)))))
  expect_identical(wafc_J_text(d), "6, 6 (capped: u2 at 4 (23 points))")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  pc <- plot(f, s = s)
  expect_equal(pc[["components"]][["g"]][[1L, 2L]],
               wafc_functions(f, s = s)[["g"]][[1L, 2L]], tolerance = 1e-12)

  ## the threshold reads the norms of the capped blocks
  th <- wafc_threshold(f, s = s, t = 0)
  expect_equal(th[["extra"]][["norm"]], wafc_blocks(f, s = s)[["norm"]],
               tolerance = 1e-12)
  expect_identical(th[["extra"]][["J"]], 6L)
})

test_that("the cap validates its argument and warns on a covariate of one point", {
  expect_error(wafc_design(x_h, hour, J = 4, cap.J = NA), "cap.J")
  expect_error(wafc_design(x_h, hour, J = 4, cap.J = c(TRUE, FALSE)), "cap.J")
  bin <- rep(0:1, length.out = n_h)
  expect_warning(d <- wafc_design(x_h, bin, J = 4), "single point")
  expect_identical(d[["J"]], 1L)
  expect_identical(d[["ndistinct"]], 1L)
  expect_silent(d <- wafc_design(x_h, bin, J = 4, eps = 0.05))
  expect_identical(d[["J"]], 1L)
  expect_identical(d[["ndistinct"]], 2L)

  ## a spec built before the cap existed has the levels it was built at
  old <- wafc_design(x_h, hour, J = 3)
  old[c("J.requested", "ndistinct", "cap.J")] <- NULL
  dn <- wafc_design(x_h, hour, spec = old)
  expect_identical(dn[["J"]], 3L)
  expect_identical(dn[["J.requested"]], 3L)
  expect_identical(colnames(dn[["Z"]]), colnames(old[["Z"]]))
})
