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
