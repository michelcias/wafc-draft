## wafc/tests/test-dgp.R -- tests of the data generating processes of step
## E2.1 and of the two things step E2.4b added to them. Run from the root of
## the repository:
##
##     Rscript -e 'testthat::test_dir("wafc/tests")'
##
## The two tests that carry the step are the one on the seam, which is what
## makes the scenario "uneven" measure uneven curvature and not a corner in
## the periodic extension, and the one on wafc_sprime(), which is the
## number decision D27 fixed and which step E2.4 had to keep in a script
## because dgp.R was outside its catalogue.

library(testthat)

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

grid <- (seq_len(2^16) - 0.5) / 2^16
comps <- c("sine", "cosine", "cubic", "gaussians", "chirp", "bumps",
           "blocks", "heavisine")

## ---------------------------------------------------------------------------
## The components
## ---------------------------------------------------------------------------

test_that("every component is centred and normalised on the unit interval", {
  for (nm in comps) {
    g <- wafc_component(nm)
    v <- g(grid)
    expect_equal(mean(v), 0, tolerance = 1e-8, label = nm)
    expect_equal(mean(v^2), 1, tolerance = 1e-6, label = nm)
  }
  z <- wafc_component("zero")
  expect_equal(z(grid), rep(0, length(grid)))
  expect_error(wafc_component("doppler"))
})

test_that("the window vanishes to every order at the endpoints and is one inside", {
  expect_equal(wafc_window(c(0, 1)), c(0, 0))
  expect_equal(wafc_window(c(0.2, 0.5, 0.8)), c(1, 1, 1))
  ## the transition is C^infinity: the first difference quotient at the
  ## endpoint is below anything a corner would produce
  h <- 10^-(3:6)
  expect_true(all(wafc_window(h) / h < 1e-8))
  expect_true(all(diff(wafc_window(seq(0, 0.08, length.out = 50))) >= 0))
  expect_error(wafc_window(0.5, d = 0.6), "in \\(0, 0.5\\)")
})

test_that("the two components of decision D30 have no seam, and the cubic has one", {
  ## A component whose periodic extension is C^infinity matches value and
  ## derivative where 1 joins 0. That is what separates the uneven scenario
  ## from the smooth one, where the corner of the cubic is the reason the
  ## effective regularity reads 3/2 and not 4 (Lemma 12 of E1.8).
  h <- 1e-6
  seam <- function(nm) {
    g <- wafc_component(nm)
    c(value = abs(g(1) - g(0)),
      deriv = abs((g(1) - g(1 - h)) / h - (g(h) - g(0)) / h))
  }
  for (nm in c("gaussians", "chirp")) {
    s <- seam(nm)
    expect_lt(unname(s[["value"]]), 1e-8)
    expect_lt(unname(s[["deriv"]]), 1e-6)
  }
  sc <- seam("cubic")
  expect_lt(unname(sc[["value"]]), 1e-8)
  expect_gt(unname(sc[["deriv"]]), 0.5)
})

test_that("the uneven components vary their scale and the smooth ones do not", {
  ## The point of decision D30 is that a single smoothing parameter has to
  ## compromise, so what has to be shown is a comparison and not a
  ## threshold: the local curvature of the two new components is far less
  ## uniform along the domain than that of any component of the smooth
  ## scenario. Measured on sixteen windows as the largest local curvature
  ## over the typical one; an absolute cut would only say where the
  ## constants of these particular functions happen to fall.
  rough <- function(nm) {
    g <- wafc_component(nm)
    v <- g(seq(0, 1, length.out = 4097))
    d2 <- abs(diff(v, differences = 2L))
    m <- vapply(split(d2, cut(seq_along(d2), 16L)), max, 0)
    max(m) / stats::median(m)
  }
  uneven <- vapply(c("gaussians", "chirp"), rough, 0)
  smooth <- vapply(c("sine", "cosine", "cubic"), rough, 0)
  expect_gt(min(uneven), 2 * max(smooth))
})

## ---------------------------------------------------------------------------
## The scenarios and the regularity they declare
## ---------------------------------------------------------------------------

test_that("the scenarios activate three blocks and name their components", {
  for (sc in c("smooth", "uneven", "inhomogeneous")) {
    st <- wafc_scenario(sc, 3L, 2L)
    expect_equal(dim(st), c(3L, 2L))
    expect_equal(sum(nzchar(st)), 3L)
    expect_true(all(st[nzchar(st)] %in% comps))
    expect_equal(which(nzchar(st)), c(1L, 2L, 4L))
  }
  expect_equal(sum(nzchar(wafc_scenario("null", 3L, 2L))), 0L)
  expect_equal(wafc_scenario("uneven", 3L, 2L)[1L, ], c("gaussians", "chirp"))
  ## q = 1 drops the second component, p = 1 the third
  expect_equal(sum(nzchar(wafc_scenario("uneven", 3L, 1L))), 2L)
  expect_equal(sum(nzchar(wafc_scenario("uneven", 1L, 2L))), 2L)
})

test_that("wafc_sprime is the number of decision D27, and carries its regime", {
  expect_equal(as.numeric(wafc_sprime("smooth")), 3/2)
  expect_equal(as.numeric(wafc_sprime("inhomogeneous")), 1/2)
  expect_equal(as.numeric(wafc_sprime("uneven")), 4)
  expect_equal(as.numeric(wafc_sprime("null")), 3/2)
  ## the corner of the cubic is bought away by the margin and absent from
  ## the boundary corrected basis, so only the smooth scenario moves
  for (rg in c("margin", "interval")) {
    expect_equal(as.numeric(wafc_sprime("smooth", rg)), 4)
    expect_equal(as.numeric(wafc_sprime("inhomogeneous", rg)), 1/2)
    expect_equal(as.numeric(wafc_sprime("uneven", rg)), 4)
    expect_equal(attr(wafc_sprime("smooth", rg), "regime"), rg)
  }
  expect_equal(attr(wafc_sprime("smooth"), "regime"), "periodic")
  expect_error(wafc_sprime("bumps"))
  expect_error(wafc_sprime("smooth", "cdv"))
})

test_that("the scenario and the sample carry the declared regularity", {
  st <- wafc_scenario("inhomogeneous", 3L, 2L, regime = "margin")
  expect_equal(as.numeric(attr(st, "sprime")), 1/2)
  expect_equal(attr(st, "regime"), "margin")
  d <- simulate_wafc(200L, p = 3L, q = 2L, scenario = "uneven", seed = 1L)
  expect_equal(d[["sprime"]], 4)
  expect_equal(d[["regime"]], "periodic")
  dm <- simulate_wafc(200L, p = 3L, q = 2L, scenario = "smooth", seed = 1L,
                      regime = "margin")
  expect_equal(dm[["sprime"]], 4)
  ## the regime labels a number and changes nothing about the sample
  d0 <- simulate_wafc(200L, p = 3L, q = 2L, scenario = "smooth", seed = 1L)
  expect_equal(dm[["y"]], d0[["y"]])
  expect_equal(d0[["sprime"]], 3/2)
})

## ---------------------------------------------------------------------------
## The sample
## ---------------------------------------------------------------------------

test_that("the uneven scenario simulates the model it declares", {
  d <- simulate_wafc(300L, p = 3L, q = 2L, scenario = "uneven", seed = 7L,
                     snr = 4)
  expect_s3_class(d, "wafc_dgp")
  expect_equal(d[["scenario"]], "uneven")
  ## beta_l(U) is the level plus the components, and the regression
  ## function is rowSums(x * beta)
  b <- wafc_beta(d, d[["u"]])
  expect_equal(b, d[["beta"]])
  expect_equal(d[["f"]], as.numeric(rowSums(d[["x"]] * d[["beta"]])))
  expect_equal(d[["beta"]][, 3L], rep(d[["cc"]][3L], 300L))
  ## the two components of beta_1 are far from proportional, unlike the
  ## sine and the cubic of the smooth scenario (step E1.7a measured 0.969),
  ## so the uneven scenario is not a gratuitous worst case for any measure
  ## of structure
  gr <- (seq_len(4096L) - 0.5) / 4096L
  ip <- function(a, b) mean(a(gr) * b(gr))
  expect_lt(abs(ip(wafc_component("gaussians"), wafc_component("chirp"))), 0.3)
  expect_gt(abs(ip(wafc_component("sine"), wafc_component("cubic"))), 0.9)
})

test_that("the reproducibility of a draw does not depend on the scenario", {
  a <- simulate_wafc(150L, scenario = "uneven", seed = 3L)
  b <- simulate_wafc(150L, scenario = "uneven", seed = 3L)
  expect_equal(a[["y"]], b[["y"]])
  expect_equal(a[["u"]], b[["u"]])
})

## ---------------------------------------------------------------------------
## Linear covariates dependent on the modulating ones (step E4.1b)
## ---------------------------------------------------------------------------

test_that("the default draws the samples drawn before 'x_u_rho' existed", {
  ## Fingerprints printed by the code of commit 35b1ac3, before the
  ## argument existed: a weighted sum of y, of x and of u, so that a change
  ## in the order of the draws, and not only in their values, moves them.
  cfgs <- list(
    list(n = 200L, p = 3L, q = 2L, scenario = "smooth", seed = 1L),
    list(n = 300L, p = 3L, q = 2L, scenario = "uneven", seed = 7L),
    list(n = 250L, p = 4L, q = 4L, scenario = "inhomogeneous", seed = 2026L,
         snr = 3),
    list(n = 200L, p = 3L, q = 2L, scenario = "null", seed = 11L,
         sigma = 0.62),
    list(n = 200L, p = 3L, q = 3L, scenario = "smooth", seed = 3L,
         x_dist = "uniform", u_dist = "beta", u_rho = 0.4),
    list(n = 150L, p = 3L, q = 2L, scenario = "inhomogeneous", seed = 5L,
         intercept = FALSE, amplitude = 2)
  )
  want <- rbind(c(24006.790944074142, 8562.0987401374896, 38808.819022007519),
                c(37534.638510872632, 29338.894679189467, 90081.430520175258),
                c(22989.429750043091, 15373.521886492303, 252506.25148645299),
                c(20932.220559466303, 32083.591502384064, 40518.953004339244),
                c(21967.580523552428, 25896.257231129333, 84249.95270883858),
                c(5425.1663801643826, 770.37371897757487, 21732.94920101529))
  fp <- function(d) {
    c(sum(d[["y"]] * seq_along(d[["y"]])), sum(d[["x"]] * seq_along(d[["x"]])),
      sum(d[["u"]] * seq_along(d[["u"]])))
  }
  ## every element but the call, which records the argument, and the
  ## components, which are closures with an environment of their own
  keep <- setdiff(names(simulate_wafc(10L, seed = 1L)), c("call", "g"))
  for (i in seq_along(cfgs)) {
    d <- do.call(simulate_wafc, cfgs[[i]])
    expect_equal(fp(d), want[i, ], tolerance = 1e-10, label = paste("cfg", i))
    ## the explicit zero is the default, down to the state of the generator
    r_default <- .Random.seed
    d0 <- do.call(simulate_wafc, c(cfgs[[i]], list(x_u_rho = 0)))
    expect_identical(d0[keep], d[keep])
    expect_identical(.Random.seed, r_default)
  }
})

test_that("a positive 'x_u_rho' adds no draw and keeps u, Z and the errors", {
  for (xd in c("gaussian", "uniform")) {
    a <- list(n = 300L, p = 4L, q = 2L, scenario = "inhomogeneous",
              seed = 13L, x_dist = xd)
    d0 <- do.call(simulate_wafc, a)
    r0 <- .Random.seed
    d1 <- do.call(simulate_wafc, c(a, list(x_u_rho = 0.6)))
    expect_identical(.Random.seed, r0)
    expect_identical(d1[["u"]], d0[["u"]])
    expect_equal((d1[["y"]] - d1[["f"]]) / d1[["sigma"]],
                 (d0[["y"]] - d0[["f"]]) / d0[["sigma"]], tolerance = 1e-12)
    ## the constant column is untouched, and Z is the old draw at unit
    ## variance once the score of the paired modulating covariate is removed
    expect_identical(d1[["x"]][, 1L], rep(1, 300L))
    expect_equal(unname(d1[["x_u"]]), c(NA, 1L, 2L, 1L))
    h <- sqrt(12) * (d1[["u"]][, c(1L, 2L, 1L)] - 0.5)
    z <- (d1[["x"]][, 2:4] - 0.6 * h) / sqrt(1 - 0.6^2)
    s <- if (xd == "uniform") sqrt(3) else 1
    expect_equal(z, s * d0[["x"]][, 2:4], tolerance = 1e-12,
                 ignore_attr = TRUE)
  }
  ## without the intercept the first covariate is paired too, and with no
  ## non-constant covariate nothing moves
  d <- simulate_wafc(50L, p = 3L, q = 2L, intercept = FALSE, seed = 1L,
                     x_u_rho = 0.5)
  expect_equal(unname(d[["x_u"]]), c(1L, 2L, 1L))
  expect_equal(names(d[["x_u"]]), c("x1", "x2", "x3"))
  d <- simulate_wafc(50L, p = 1L, q = 2L, scenario = "smooth", seed = 1L,
                     x_u_rho = 0.5)
  expect_identical(d[["x"]][, 1L], rep(1, 50L))
  expect_true(all(is.na(d[["x_u"]])))
})

test_that("a paired covariate has unit variance and correlation rho with its modulating covariate", {
  ## n = 1e5, so the standard error of a correlation is at most 0.0032 and
  ## that of a variance about 0.005; the cuts are five of them or more.
  ## u_dist = "beta" makes U_2 a Beta(2,3), which the score has to centre.
  for (xd in c("gaussian", "uniform")) for (ud in c("uniform", "beta")) {
    for (rho in c(0.3, 0.7)) {
      d <- simulate_wafc(1e5L, p = 3L, q = 2L, scenario = "smooth",
                         seed = 99L, x_dist = xd, u_dist = ud,
                         x_u_rho = rho)
      x <- d[["x"]]
      u <- d[["u"]]
      lab <- paste(xd, ud, rho)
      for (l in 2:3) {
        expect_lt(abs(mean(x[, l])), 0.02, label = lab)
        expect_lt(abs(stats::var(x[, l]) - 1), 0.03, label = lab)
        expect_lt(abs(stats::cor(x[, l], u[, l - 1L]) - rho), 0.015,
                  label = lab)
      }
      ## the pairing is one modulating covariate each
      expect_lt(abs(stats::cor(x[, 2L], u[, 2L])), 0.015, label = lab)
      expect_lt(abs(stats::cor(x[, 3L], u[, 1L])), 0.015, label = lab)
    }
  }
})

test_that("E(XX'|U) keeps its smallest eigenvalue at 1 - rho^2 once the intercept is removed", {
  ## Estimated on a 4 by 4 grid of cells of (U_1, U_2), about 12 500 points
  ## each. Inside a cell the covariance of (X_2, X_3) is (1 - rho^2) I plus
  ## rho^2 times the covariance of the scores over the cell, rho^2/16 on the
  ## diagonal, so its smallest eigenvalue sits just above 1 - rho^2: it is
  ## bounded below by it, and close to it, which says that the dependence
  ## is there. With the intercept, the second moment matrix stays above
  ## (1 - rho^2)/(2 - rho^2 + 6 rho^2), the bound of the roxygen at the
  ## largest scores, |h| <= sqrt(3) for each of the two.
  for (rho in c(0, 0.5, 0.8)) {
    d <- simulate_wafc(2e5L, p = 3L, q = 2L, scenario = "smooth",
                       seed = 2026L, x_u_rho = rho)
    x <- d[["x"]]
    cell <- interaction(cut(d[["u"]][, 1L], 0:4 / 4, include.lowest = TRUE),
                        cut(d[["u"]][, 2L], 0:4 / 4, include.lowest = TRUE))
    lo <- 1 - rho^2
    bound <- lo / (2 - rho^2 + 6 * rho^2)
    ev <- vapply(levels(cell), function(k) {
      xs <- x[cell == k, , drop = FALSE]
      lmin <- function(a) {
        min(eigen(a, symmetric = TRUE, only.values = TRUE)[["values"]])
      }
      c(schur = lmin(stats::cov(xs[, 2:3])),
        full = lmin(crossprod(xs) / nrow(xs)))
    }, numeric(2L))
    expect_gt(min(ev["schur", ]), lo - 0.05, label = paste("rho", rho))
    expect_lt(max(ev["schur", ]), lo + rho^2 / 16 + 0.05,
              label = paste("rho", rho))
    expect_gt(min(ev["full", ]), bound - 0.02, label = paste("rho", rho))
  }
})

test_that("'x_u_rho' outside [0, 1) is an informative error", {
  for (bad in list(1, -0.1, NA_real_, Inf, c(0.1, 0.2), "0.5", TRUE)) {
    expect_error(simulate_wafc(50L, seed = 1L, x_u_rho = bad),
                 "'x_u_rho' must be a single value in \\[0, 1\\)")
  }
})

test_that("the signal to noise ratio is read on the regression function of the sample", {
  d <- simulate_wafc(400L, p = 3L, q = 2L, scenario = "uneven", seed = 4L,
                     snr = 2, x_u_rho = 0.5)
  expect_identical(d[["sigma"]], stats::sd(d[["f"]]) / 2)
  expect_equal(d[["f"]], as.numeric(rowSums(d[["x"]] * d[["beta"]])))
  expect_equal(wafc_beta(d, d[["u"]]), d[["beta"]])
  expect_equal(d[["x_u_rho"]], 0.5)
})

## ---------------------------------------------------------------------------
## A structure of its own (step E4.1c)
## ---------------------------------------------------------------------------

test_that("the default draws the samples drawn before 'structure' existed", {
  ## Fingerprints printed by the code of commit 901fac7, before the
  ## argument existed, as in the test of 'x_u_rho' above; the fifth draw
  ## has 'x_u_rho' and the sixth no intercept and the copula.
  cfgs <- list(
    list(n = 200L, p = 3L, q = 2L, scenario = "smooth", seed = 1L),
    list(n = 300L, p = 3L, q = 2L, scenario = "uneven", seed = 7L),
    list(n = 250L, p = 4L, q = 4L, scenario = "inhomogeneous", seed = 2026L,
         snr = 3),
    list(n = 200L, p = 3L, q = 2L, scenario = "null", seed = 11L,
         sigma = 0.62),
    list(n = 250L, p = 6L, q = 4L, scenario = "inhomogeneous", seed = 31L,
         snr = 3, x_u_rho = 0.5),
    list(n = 150L, p = 3L, q = 2L, scenario = "inhomogeneous", seed = 5L,
         intercept = FALSE, amplitude = 2, u_rho = 0.4)
  )
  want <- rbind(c(24006.790944074142, 8562.0987401374896, 38808.819022007519),
                c(37534.638510872632, 29338.894679189467, 90081.430520175258),
                c(22989.429750043091, 15373.521886492303, 252506.25148645299),
                c(20932.220559466303, 32083.591502384064, 40518.953004339244),
                c(34277.880721259746, 18597.019907742815, 254045.81913873018),
                c(-9184.3354015473051, 6750.654140973762, 22455.129269914039))
  fp <- function(d) {
    c(sum(d[["y"]] * seq_along(d[["y"]])), sum(d[["x"]] * seq_along(d[["x"]])),
      sum(d[["u"]] * seq_along(d[["u"]])))
  }
  keep <- setdiff(names(simulate_wafc(10L, seed = 1L)), c("call", "g"))
  for (i in seq_along(cfgs)) {
    d <- do.call(simulate_wafc, cfgs[[i]])
    expect_equal(fp(d), want[i, ], tolerance = 1e-10, label = paste("cfg", i))
    r_default <- .Random.seed
    d0 <- do.call(simulate_wafc, c(cfgs[[i]], list(structure = NULL)))
    expect_identical(d0[keep], d[keep])
    expect_identical(.Random.seed, r_default)
  }
})

test_that("the structure of a scenario, given, draws the sample of the scenario", {
  for (sc in c("smooth", "uneven", "inhomogeneous", "null")) {
    for (pq in list(c(3L, 2L), c(4L, 4L))) {
      a <- list(n = 200L, p = pq[1L], q = pq[2L], scenario = sc, seed = 17L)
      if (sc == "null") a[["sigma"]] <- 0.5
      d <- do.call(simulate_wafc, a)
      r <- .Random.seed
      s <- wafc_scenario(sc, pq[1L], pq[2L])
      ## a zero block may also be named "zero"
      s0 <- s
      s0[!nzchar(s0)] <- "zero"
      for (st in list(s, s0)) {
        e <- do.call(simulate_wafc, c(a, list(structure = st)))
        expect_identical(.Random.seed, r)
        for (nm in c("y", "x", "u", "f", "beta", "sigma", "cc", "x_u")) {
          expect_identical(e[[nm]], d[[nm]], label = paste(sc, nm))
        }
        expect_identical(as.vector(e[["structure"]]),
                         as.vector(d[["structure"]]))
        expect_identical(e[["scenario"]], sc)
        expect_true(is.na(e[["sprime"]]))
        expect_identical(attr(e[["structure"]], "regime"), "periodic")
      }
    }
  }
})

test_that("a structure of its own is drawn on the covariates and errors of the seed", {
  ## the larger model of decision D63: the inhomogeneous structure twice
  s <- matrix("", 6L, 4L)
  s[1L, 1:2] <- c("bumps", "blocks")
  s[2L, 1L] <- "heavisine"
  s[3L, 3:4] <- c("bumps", "blocks")
  s[4L, 3L] <- "heavisine"
  for (xu in c(0, 0.5)) {
    a <- list(n = 300L, p = 6L, q = 4L, seed = 23L, x_u_rho = xu)
    d <- do.call(simulate_wafc, c(a, list(scenario = "inhomogeneous", snr = 3,
                                          structure = s)))
    ## the null scenario with unit errors draws the same u, x and errors
    d0 <- do.call(simulate_wafc, c(a, list(scenario = "null", sigma = 1)))
    expect_identical(d[["u"]], d0[["u"]])
    expect_identical(d[["x"]], d0[["x"]])
    expect_equal((d[["y"]] - d[["f"]]) / d[["sigma"]], d0[["y"]] - d0[["f"]],
                 tolerance = 1e-12)
    ## the coefficients are the components named, block by block
    beta <- matrix(rep(d[["cc"]], each = 300L), 300L, 6L)
    for (l in 1:6) for (m in 1:4) {
      if (nzchar(s[l, m])) {
        beta[, l] <- beta[, l] + wafc_component(s[l, m])(d[["u"]][, m])
      }
    }
    expect_equal(unname(d[["beta"]]), beta, tolerance = 1e-14)
    expect_identical(d[["sigma"]], stats::sd(d[["f"]]) / 3)
    expect_identical(!vapply(d[["g"]], is.null, NA), as.vector(nzchar(s)))
    expect_identical(as.vector(d[["structure"]]), as.vector(s))
    expect_true(is.na(d[["sprime"]]))
  }
})

test_that("an ill-formed 'structure' is an informative error", {
  s <- wafc_scenario("smooth", 3L, 2L)
  sim <- function(st) simulate_wafc(50L, p = 3L, q = 2L, seed = 1L,
                                    structure = st)
  expect_error(sim(as.vector(s)), "character matrix of component names")
  expect_error(sim(matrix(0, 3L, 2L)), "character matrix of component names")
  expect_error(sim(matrix("", 2L, 3L)), "must be 3 by 2 \\(p by q\\); it is 2 by 3")
  s_na <- s
  s_na[3L, 2L] <- NA
  expect_error(sim(s_na), "must not contain NA")
  s_bad <- s
  s_bad[3L, 2L] <- "doppler"
  expect_error(sim(s_bad), "unknown component\\(s\\) in 'structure': \"doppler\"")
})
