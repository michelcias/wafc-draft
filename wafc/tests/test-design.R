## wafc/tests/test-design.R -- tests of the design matrix and of the data
## generating processes of step E2.1. Run from the root of the repository:
##
##     Rscript -e 'testthat::test_dir("wafc/tests")'
##
## The first test is the exact recovery required by the plan: with the
## components in the basis (theta* known) and no noise, the LASSO with
## lambda -> 0 returns theta*.

library(testthat)

## Loads wafc/R/ whether the tests are run from the root of the repository
## (source) or from wafc/tests/ (testthat::test_dir, which sets the working
## directory to the test folder).
local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

## Reference construction of derivations/check/03-desenho-produtos.R: the
## design as the row-wise Kronecker product of X and (1, Psi(U_1), ...,
## Psi(U_q)), permuted to the order of D12. Independent of design.R and used
## here to check the column order.
ref_kron_rows <- function(X, P) {
  X[, rep(seq_len(ncol(X)), each = ncol(P)), drop = FALSE] *
    P[, rep(seq_len(ncol(P)), times = ncol(X)), drop = FALSE]
}
ref_perm_D12 <- function(p, q, NJ) {
  K <- 1L + q * NJ
  kron_idx <- function(l, col) as.integer((l - 1L) * K + col)
  c(vapply(seq_len(p), function(l) kron_idx(l, 1L), 1L),
    unlist(lapply(seq_len(p), function(l)
      lapply(seq_len(q), function(m)
        kron_idx(l, 1L + (m - 1L) * NJ + seq_len(NJ))))))
}
ref_design <- function(X, U, J, filter.size = 8L) {
  NJ <- 2L^J - 1L
  P <- cbind(1, do.call(cbind, lapply(seq_len(ncol(U)), function(m)
    WaveBased::wbasis(U[, m], j0 = 0, J = J,
                      filter.size = filter.size)[, -1L, drop = FALSE])))
  ref_kron_rows(X, P)[, ref_perm_D12(ncol(X), ncol(U), NJ), drop = FALSE]
}

set.seed(20260919)
n <- 200L
p <- 2L
q <- 2L
J <- 3L
NJ <- 2L^J - 1L
x0 <- cbind(1, stats::rnorm(n))
u0 <- matrix(stats::runif(n * q), n, q)

## ---------------------------------------------------------------------------

test_that("the loader declares the dependencies it uses", {
  expect_true(all(c("WaveBased", "glmnet", "Matrix") %in% wafc_depends))
  expect_true(all(vapply(wafc_depends, requireNamespace, TRUE, quietly = TRUE)))
  expect_silent(wafc_attach())
})

test_that("wafc_rescale maps onto [eps, 1 - eps] and reuses the transformation", {
  r <- wafc_rescale(u0, eps = 0.05)
  expect_equal(apply(r[["u"]], 2L, min), rep(0.05, q))
  expect_equal(apply(r[["u"]], 2L, max), rep(0.95, q))
  ## the transformation is affine and invertible
  back <- sweep(sweep(r[["u"]], 2L, r[["scale"]], "*"), 2L, r[["location"]], "+")
  expect_equal(back, u0)
  ## reusing location and scale on values outside the original range
  unew <- rbind(c(-1, 2), c(0.5, 0.5))
  r2 <- wafc_rescale(unew, eps = r[["eps"]], location = r[["location"]],
                     scale = r[["scale"]], clip = TRUE)
  expect_true(all(r2[["u"]] >= 0.05 - 1e-12 & r2[["u"]] <= 0.95 + 1e-12))
  expect_error(wafc_rescale(u0, eps = 0.6), "\\[0, 0.5\\)")
  expect_error(wafc_rescale(cbind(rep(1, 5)), eps = 0), "constant")
  ## rescale = FALSE does not rescale, and accepts a constant column
  d <- wafc_design(x0, cbind(u0[, 1L], rep(0.5, n)), J = J, rescale = FALSE)
  expect_equal(d[["location"]], rep(0, q))
  expect_equal(d[["scale"]], rep(1, q))
  expect_warning(wafc_design(x0, u0 + 2, J = J, rescale = FALSE), "outside")
})

test_that("the design has the columns of D12, named and in order", {
  d <- wafc_design(x0, u0, J = J, rescale = FALSE)
  expect_s3_class(d, "wafc_design")
  expect_equal(dim(d[["Z"]]), c(n, p + p * q * NJ))
  expect_equal(d[["nvars"]], 2L + 4L * 7L)
  ## first the unpenalized columns, then the blocks (l, m) lexicographic
  expect_equal(colnames(d[["Z"]])[seq_len(p)], c("x1", "x2"))
  expect_equal(names(d[["blocks"]]), c("x1:u1", "x1:u2", "x2:u1", "x2:u2"))
  expect_equal(unlist(d[["blocks"]], use.names = FALSE), (p + 1L):d[["nvars"]])
  ## inside a block, increasing j and, within a level, increasing k
  expect_equal(colnames(d[["Z"]])[d[["blocks"]][["x2:u1"]]],
               c("x2:u1:psi0.0", "x2:u1:psi1.0", "x2:u1:psi1.1",
                 "x2:u1:psi2.0", "x2:u1:psi2.1", "x2:u1:psi2.2",
                 "x2:u1:psi2.3"))
  ## penalty factors: zero on the p levels c_l, one on every wavelet
  expect_equal(d[["penalty.factor"]], c(0, 0, rep(1, p * q * NJ)))
  expect_equal(d[["unpenalized"]], seq_len(p))
  ## column names are unique, which coef() and the block norms rely on
  expect_false(anyDuplicated(colnames(d[["Z"]])) > 0L)
})

test_that("the design equals the reference Kronecker construction", {
  d <- wafc_design(x0, u0, J = J, rescale = FALSE)
  expect_equal(unname(as.matrix(d[["Z"]])), unname(ref_design(x0, u0, J)),
               tolerance = 1e-12)
})

test_that("the constant scaling function phi_00 is discarded", {
  d <- wafc_design(x0, u0, J = J, rescale = FALSE)
  B <- WaveBased::wbasis(u0[, 1L], j0 = 0, J = J, filter.size = 8L)
  expect_equal(max(abs(B[, 1L] - 1)), 0, tolerance = 1e-12)
  ## the block has 2^J - 1 columns, and none of them is the constant
  expect_equal(length(d[["blocks"]][["x1:u1"]]), NJ)
  expect_equal(unname(as.matrix(d[["Z"]])[, d[["blocks"]][["x1:u1"]]]),
               unname(B[, -1L, drop = FALSE]), tolerance = 1e-12)
  ## keeping it would make the Gram singular with nullity p q (Lemma 1(iii))
  Zfull <- cbind(as.matrix(d[["Z"]]), x0[, 1L], x0[, 1L], x0[, 2L], x0[, 2L])
  expect_equal(qr(Zfull)$rank, d[["nvars"]])
})

test_that("the Gram matrix of the design has full rank", {
  d <- wafc_design(x0, u0, J = J, rescale = FALSE)
  Z <- as.matrix(d[["Z"]])
  expect_equal(qr(Z)$rank, d[["nvars"]])
  ev <- eigen(crossprod(Z) / n, symmetric = TRUE, only.values = TRUE)[["values"]]
  expect_gt(min(ev), 1e-4)
})

test_that("with theta* in the basis and no noise, the LASSO with lambda -> 0 recovers theta*", {
  d <- wafc_design(x0, u0, J = J, rescale = FALSE)
  Z <- as.matrix(d[["Z"]])
  ## sparse truth: the levels c_l and a few wavelet coefficients
  b_true <- numeric(d[["nvars"]])
  b_true[d[["unpenalized"]]] <- c(1, 2)
  active <- c(d[["blocks"]][["x1:u1"]][c(1L, 4L)],
              d[["blocks"]][["x2:u2"]][c(2L, 7L)])
  b_true[active] <- c(1.5, -1, 0.8, 1.2)
  y <- as.numeric(Z %*% b_true)

  ## least squares recovers it exactly (identifiability of the sieve)
  b_ls <- qr.solve(Z, y)
  expect_equal(max(abs(b_ls - b_true)), 0, tolerance = 1e-10)

  ## and so does the LASSO along a path ending near zero
  ## and so does the LASSO along a path ending near zero. glmnet drops the
  ## constant column x1 from the fit and its intercept carries c_1, which is
  ## why the fit uses intercept = TRUE and the design reports the column in
  ## d[["constant"]].
  expect_equal(d[["constant"]], 1L)
  fit <- glmnet::glmnet(Z, y, family = "gaussian", intercept = TRUE,
                        standardize = FALSE,
                        penalty.factor = d[["penalty.factor"]],
                        lambda = c(0.1, 0.01, 1e-4, 1e-8))
  cf <- as.numeric(stats::coef(fit, s = 1e-8))
  expect_equal(cf[1L], b_true[1L], tolerance = 1e-5)
  expect_equal(cf[2L], 0)
  expect_equal(max(abs(cf[-(1:2)] - b_true[-1L])), 0, tolerance = 1e-4)
  ## with no constant column there is nothing for the intercept to absorb
  xv <- cbind(stats::rnorm(n), x0[, 2L])
  dv <- wafc_design(xv, u0, J = J, rescale = FALSE)
  expect_equal(dv[["constant"]], integer(0))
  yv <- as.numeric(as.matrix(dv[["Z"]]) %*% b_true)
  fv <- glmnet::glmnet(as.matrix(dv[["Z"]]), yv, intercept = FALSE,
                       standardize = FALSE,
                       penalty.factor = dv[["penalty.factor"]],
                       lambda = c(0.1, 0.01, 1e-4, 1e-8))
  expect_equal(max(abs(as.numeric(stats::coef(fv, s = 1e-8))[-1L] - b_true)), 0,
               tolerance = 1e-4)
  ## at a lambda of order one the penalized part is shrunk, the levels are not
  b_mid <- as.numeric(stats::coef(fv, s = 0.1))[-1L]
  expect_true(sum(b_mid[-dv[["unpenalized"]]] != 0) < length(active) + 10L)
  expect_true(all(b_mid[dv[["unpenalized"]]] != 0))
})

test_that("the sparse design is the dense one", {
  dd <- wafc_design(x0, u0, J = 4L, rescale = FALSE, sparse = "never")
  ds <- wafc_design(x0, u0, J = 4L, rescale = FALSE, sparse = "always")
  expect_false(dd[["sparse"]])
  expect_true(ds[["sparse"]])
  expect_true(inherits(ds[["Z"]], "Matrix"))
  expect_equal(as.matrix(ds[["Z"]]), dd[["Z"]], tolerance = 1e-12)
  expect_equal(colnames(ds[["Z"]]), colnames(dd[["Z"]]))
  expect_equal(ds[["penalty.factor"]], dd[["penalty.factor"]])
})

test_that("the lookup table reproduces the exact basis evaluation", {
  de <- wafc_design(x0, u0, J = J, rescale = FALSE, use.table = "never")
  dt <- wafc_design(x0, u0, J = J, rescale = FALSE, use.table = "always")
  expect_null(de[["wavelet.table"]])
  expect_false(is.null(dt[["wavelet.table"]]))
  expect_equal(as.matrix(dt[["Z"]]), as.matrix(de[["Z"]]), tolerance = 1e-4)
  ## a table supplied by the caller is used as it is
  tb <- WaveBased::wtable(filter.size = 8L, check = FALSE)
  d2 <- wafc_design(x0, u0, J = J, rescale = FALSE, wavelet.table = tb,
                    use.table = "never")
  expect_equal(as.matrix(d2[["Z"]]), as.matrix(dt[["Z"]]), tolerance = 1e-12)
})

test_that("the design works with one linear and one modulating covariate", {
  d <- wafc_design(x0[, 2L], u0[, 1L], J = J, rescale = FALSE)
  expect_equal(d[["nvars"]], 1L + NJ)
  expect_equal(names(d[["blocks"]]), "x1:u1")
  expect_equal(d[["penalty.factor"]], c(0, rep(1, NJ)))
  expect_equal(qr(as.matrix(d[["Z"]]))$rank, d[["nvars"]])
})

test_that("a design built with spec repeats the basis and the rescaling", {
  d <- wafc_design(x0, u0, J = J, rescale = TRUE)
  expect_equal(unname(d[["eps"]]), rep(wafc_eps_periodic, q))
  ## the same data give the same design
  d2 <- wafc_design(x0, u0, spec = d)
  expect_equal(as.matrix(d2[["Z"]]), as.matrix(d[["Z"]]), tolerance = 1e-12)
  expect_equal(d2[["location"]], d[["location"]])
  ## new data outside the training range are clipped, not extrapolated
  unew <- rbind(c(-5, 5), c(0.4, 0.6))
  d3 <- wafc_design(x0[1:2, ], unew, spec = d)
  expect_equal(nrow(d3[["Z"]]), 2L)
  expect_equal(colnames(d3[["Z"]]), colnames(d[["Z"]]))
  expect_true(all(is.finite(as.matrix(d3[["Z"]]))))
  expect_error(wafc_design(x0, u0[, 1L, drop = FALSE], spec = d), "column")
})

## The margin of the rescaling is a feature of the target, not of the sieve
## (decision D23, step E2.1b): it used to be 1.9^(-J), inherited from
## WaveBased::wall, which moved the rescaled support with J and was not even
## admissible at J = 1.
test_that("the default margin does not depend on the resolution level", {
  eps_of <- function(Jm) unname(wafc_design(x0, u0, J = Jm)[["eps"]])
  expect_equal(eps_of(2L), rep(wafc_eps_periodic, q))
  expect_equal(eps_of(5L), rep(wafc_eps_periodic, q))
  expect_equal(eps_of(c(2L, 5L)), rep(wafc_eps_periodic, q))
  expect_true(wafc_eps_periodic > 0 && wafc_eps_periodic < 0.5)
  ## the rescaled support, and so the estimated target, is the same for
  ## every candidate of the grid of cv.wafc()
  d2 <- wafc_design(x0, u0, J = 2L)
  d5 <- wafc_design(x0, u0, J = 5L)
  expect_equal(d2[["location"]], d5[["location"]])
  expect_equal(d2[["scale"]], d5[["scale"]])
})

test_that("a design at J = 1 is built with the defaults", {
  d <- wafc_design(x0, u0, J = 1L)
  expect_equal(d[["nvars"]], p + p * q * 1L)
  expect_equal(unname(d[["NJ"]]), rep(1L, q))
  expect_true(all(is.finite(as.matrix(d[["Z"]]))))
  expect_equal(qr(as.matrix(d[["Z"]]))$rank, d[["nvars"]])
})

test_that("the margin asked for by the caller is the one used", {
  d <- wafc_design(x0, u0, J = J, eps = 0.2)
  expect_equal(unname(d[["eps"]]), rep(0.2, q))
  expect_equal(range(wafc_rescale(u0, eps = 0.2)[["u"]]), c(0.2, 0.8))
  d1 <- wafc_design(x0, u0, J = J, eps = c(0.1, 0.3))
  expect_equal(unname(d1[["eps"]]), c(0.1, 0.3))
  expect_equal(unname(wafc_design(x0, u0, J = J, rescale = FALSE)[["eps"]]),
               rep(0, q))
})

## boundary = "interval" still raises an error in wafc_design(), so the
## branch is read where it lives: with no periodization there is nothing to
## keep the data away from and the default margin is zero.
test_that("the default margin of the interval basis is zero", {
  expect_equal(wafc_eps(NULL, q = 2L, rescale = TRUE, boundary = "interval"),
               c(0, 0))
  expect_equal(wafc_eps(NULL, q = 2L, rescale = TRUE, boundary = "periodic"),
               rep(wafc_eps_periodic, 2L))
  expect_equal(wafc_eps(0.1, q = 3L, rescale = TRUE, boundary = "interval"),
               rep(0.1, 3L))
  expect_error(wafc_eps(0.5, q = 1L, rescale = TRUE), "\\[0, 0.5\\)")
})

test_that("the design validates its arguments", {
  expect_error(wafc_design(x0, u0), "'J' must be provided")
  expect_error(wafc_design(x0, u0[-1L, ], J = J), "same number of rows")
  expect_error(wafc_design(x0, u0, J = 0L), "larger than 'j0'")
  expect_error(wafc_design(x0, u0, J = J, j0 = 1L), "Only j0 = 0")
  expect_error(wafc_design(x0, u0, J = J, boundary = "interval"),
               "not implemented")
  expect_error(wafc_design(x0, u0, J = c(3L, 3L, 3L)), "length 1 or one entry")
  expect_error(wafc_design(x0, u0, J = J, eps = 0.7), "\\[0, 0.5\\)")
  ## one J per modulating covariate is allowed
  d <- wafc_design(x0, u0, J = c(2L, 4L), rescale = FALSE)
  expect_equal(d[["nvars"]], p + p * ((2^2 - 1) + (2^4 - 1)))
  expect_equal(lengths(d[["blocks"]], use.names = FALSE), c(3L, 15L, 3L, 15L))
})

## ---------------------------------------------------------------------------
## Data generating processes
## ---------------------------------------------------------------------------

test_that("the components are centred and normalised in [0,1]", {
  grid <- (seq_len(2^14) - 0.5) / 2^14
  for (nm in c("sine", "cosine", "cubic", "bumps", "blocks", "heavisine")) {
    g <- wafc_component(nm)
    v <- g(grid)
    expect_lt(abs(mean(v)), 1e-3)
    expect_equal(sqrt(mean(v^2)), 1, tolerance = 1e-3)
  }
  expect_equal(wafc_component("zero")(grid), rep(0, length(grid)))
})

test_that("simulate_wafc returns the model it claims", {
  d <- simulate_wafc(300L, p = 3L, q = 2L, scenario = "smooth", seed = 1L)
  expect_s3_class(d, "wafc_dgp")
  expect_equal(dim(d[["x"]]), c(300L, 3L))
  expect_equal(dim(d[["u"]]), c(300L, 2L))
  expect_true(all(d[["x"]][, 1L] == 1))
  expect_true(all(d[["u"]] >= 0 & d[["u"]] <= 1))
  ## beta_1 is additive in the two modulating covariates, beta_3 is constant
  expect_equal(d[["structure"]][1L, ], c("sine", "cubic"))
  expect_equal(d[["structure"]][3L, ], c("", ""))
  expect_equal(d[["beta"]][, 3L], rep(d[["cc"]][3L], 300L))
  ## f is the regression function and y - f has the right scale
  expect_equal(d[["f"]], as.numeric(rowSums(d[["x"]] * d[["beta"]])))
  expect_equal(stats::sd(d[["y"]] - d[["f"]]), d[["sigma"]], tolerance = 0.2)
  expect_equal(d[["sigma"]], stats::sd(d[["f"]]) / 4, tolerance = 1e-10)
  ## the seed makes it reproducible
  expect_equal(simulate_wafc(300L, seed = 1L)[["y"]],
               simulate_wafc(300L, seed = 1L)[["y"]])
  ## wafc_beta reproduces the coefficients of the sample
  expect_equal(wafc_beta(d, d[["u"]]), d[["beta"]])
})

test_that("the null scenario has no functional part", {
  d <- simulate_wafc(200L, p = 3L, q = 2L, scenario = "null", seed = 2L,
                     sigma = 0.5)
  expect_true(all(d[["structure"]] == ""))
  expect_equal(d[["beta"]], matrix(rep(d[["cc"]], each = 200L), 200L, 3L,
                                   dimnames = list(NULL, colnames(d[["x"]]))))
  expect_equal(d[["sigma"]], 0.5)
  ## with p = 1 and no components the regression function is constant and
  ## the signal to noise ratio cannot set sigma
  expect_error(simulate_wafc(50L, p = 1L, scenario = "null"), "supply 'sigma'")
})

test_that("the inhomogeneous scenario and the optional laws of x and u", {
  d <- simulate_wafc(400L, p = 2L, q = 2L, scenario = "inhomogeneous",
                     seed = 3L, x_dist = "uniform", u_dist = "beta",
                     u_rho = 0.6)
  expect_equal(d[["structure"]][1L, ], c("bumps", "blocks"))
  expect_equal(d[["structure"]][2L, ], c("heavisine", ""))
  expect_true(all(abs(d[["x"]][, 2L]) <= 1))
  ## u_rho couples the modulating covariates without changing the support
  expect_gt(stats::cor(d[["u"]][, 1L], d[["u"]][, 2L]), 0.3)
  expect_true(all(d[["u"]] > 0 & d[["u"]] < 1))
})
