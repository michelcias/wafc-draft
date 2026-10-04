## wafc/tests/test-plot.R -- tests of the plots of step E3.2
## (wafc/R/plot.R). Run from the root of the repository:
##
##     Rscript -e 'testthat::test_dir("wafc/tests")'
##
## Every plot is drawn on a null device, pdf(NULL), and what each method
## returns is held against the quantities it says it drew, read from the
## object by the accessors and not by the plotting code: the curves of
## cv.wafc() by J, the candidates and the threshold of the rule, the
## components of wafc_functions(), the zero blocks and the norms, and the
## path of the block norms of wafc_blocks(). The truth is checked to be
## the true component shifted by a constant. The layout of the device is
## restored after the grid of components, and the errors are the ones the
## help pages state.

library(testthat)

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

n <- 200L
dgp <- simulate_wafc(n, p = 3L, q = 2L, scenario = "smooth",
                     seed = 20261003L, snr = 4)
folds <- rep_len(1:5, n)
cvb <- if (has("grpreg")) {
  cv.wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], J = 2:3, foldid = folds)
} else NULL
cvl <- cv.wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], J = 2:3, foldid = folds,
               penalty = "lasso", threshold = "none")
fitl <- cvl[["wafc.fit"]]

## Draws on a null device that is closed afterwards, and returns the value
## of the call together with its visibility.
on_null <- function(expr) {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev))
  withVisible(expr)
}

## ---------------------------------------------------------------------------
## plot.cv.wafc()
## ---------------------------------------------------------------------------

test_that("plot.cv.wafc draws the cross-validation and the threshold by default", {
  skip_if_not(has("grpreg"))
  r <- on_null(plot(cvb))
  expect_false(r[["visible"]])
  out <- r[["value"]]
  expect_identical(names(out), c("cv", "threshold"))
  ## the curves are the ones cv.wafc() stored, J by J
  tab <- out[["cv"]][["table"]]
  for (i in seq_along(cvb[["J"]])) {
    z <- cvb[["cv"]][[i]]
    rows <- tab[["J"]] == cvb[["J"]][i]
    expect_identical(tab[["lambda"]][rows], z[["lambda"]])
    expect_identical(tab[["cvm"]][rows], z[["cvm"]])
    expect_identical(tab[["cvsd"]][rows], z[["cvsd"]])
    expect_true(all(tab[["selected"]][rows] == (cvb[["J"]][i] == cvb[["J.min"]])))
  }
  expect_identical(out[["cv"]][["J.min"]], cvb[["J.min"]])
  expect_identical(out[["cv"]][["lambda.min"]], cvb[["lambda.min"]])
  expect_identical(out[["cv"]][["lambda.1se"]], cvb[["lambda.1se"]])
  ## the threshold curve is the one the rule scored, and the band is the
  ## smallest error and one standard error above it
  th <- cvb[["threshold"]]
  expect_identical(out[["threshold"]][["candidates"]], th[["candidates"]])
  expect_identical(out[["threshold"]][["t"]], th[["t"]])
  expect_identical(out[["threshold"]][["rule"]], "cv1se")
  expect_identical(out[["threshold"]][["kept"]], th[["kept"]])
  cand <- th[["candidates"]]
  i <- which.min(cand[["error"]])
  expect_equal(out[["threshold"]][["band"]],
               c(cand[["error"]][i], cand[["error"]][i] + cand[["se"]][i]))
  ## the t chosen is the largest one inside the band, the rule cv1se
  inside <- which(cand[["error"]] <= out[["threshold"]][["band"]][2L])
  expect_gte(th[["t"]], cand[["lower"]][max(inside)])
})

test_that("the components of a cv.wafc are the thresholded ones, with the truth", {
  skip_if_not(has("grpreg"))
  out <- on_null(plot(cvb, which = "components", truth = dgp,
                      n_grid = 64))[["value"]]
  expect_identical(names(out), "components")
  cp <- out[["components"]]
  fn <- wafc_functions(cvb, n_grid = 64)
  expect_identical(cp[["g"]], fn[["g"]])
  expect_identical(cp[["grid"]], fn[["grid"]])
  expect_identical(cp[["s"]], cvb[["lambda.min"]])
  th <- cvb[["threshold"]]
  expect_identical(cp[["zero"]], !th[["kept"]])
  ## the norm in the corner is the one before the threshold
  expect_identical(cp[["norm"]], th[["norm"]])
  expect_identical(cp[["t"]], th[["t"]])
  ## a zeroed block is drawn at zero
  for (i in which(cp[["zero"]])) expect_true(all(cp[["g"]][[i]] == 0))
  ## the truth is the true component shifted to the mean of the estimate
  for (l in 1:3) {
    for (m in 1:2) {
      gm <- dgp[["g"]][[l, m]]
      v <- fn[["grid"]][, m]
      gt <- if (is.null(gm)) numeric(length(v)) else gm(v)
      d <- cp[["truth"]][[l, m]] - gt
      expect_lt(diff(range(d)), 1e-12)
      expect_equal(mean(cp[["truth"]][[l, m]]), mean(fn[["g"]][[l, m]]),
                   tolerance = 1e-12)
    }
  }
  ## a list of functions is the same truth as the simulation object
  out2 <- on_null(plot(cvb, which = "components", truth = dgp[["g"]],
                       n_grid = 64))[["value"]]
  expect_identical(out2[["components"]][["truth"]], cp[["truth"]])
  ## without the threshold the components are the path at lambda.min
  out3 <- on_null(plot(cvb, which = "components", thresholded = FALSE,
                       n_grid = 64))[["value"]][["components"]]
  fn3 <- wafc_functions(cvb, n_grid = 64, thresholded = FALSE)
  expect_identical(out3[["g"]], fn3[["g"]])
  expect_identical(out3[["zero"]], fn3[["nonzero"]] == 0L)
  expect_null(out3[["t"]])
  expect_null(out3[["truth"]])
})

test_that("the path of a cv.wafc is the one of wafc.fit, with the threshold", {
  skip_if_not(has("grpreg"))
  out <- on_null(plot(cvb, which = "path"))[["value"]][["path"]]
  f <- cvb[["wafc.fit"]]
  expect_identical(out[["lambda"]], f[["lambda"]])
  blk <- wafc_blocks(f, s = f[["lambda"]])
  des <- f[["design"]]
  expect_identical(colnames(out[["norm"]]), names(des[["blocks"]]))
  for (l in 1:3) {
    for (m in 1:2) {
      expect_equal(out[["norm"]][, wafc_block_name(des, l, m)],
                   unname(blk[["norm"]][l, m, ]), tolerance = 1e-12)
    }
  }
  expect_identical(out[["marked"]], c(lambda.min = cvb[["lambda.min"]],
                                      lambda.1se = cvb[["lambda.1se"]]))
  expect_identical(out[["t"]], cvb[["threshold"]][["t"]])
})

test_that("without a scored threshold only the cross-validation is drawn", {
  out <- on_null(plot(cvl))[["value"]]
  expect_identical(names(out), "cv")
  expect_error(on_null(plot(cvl, which = "threshold")), "no threshold")
  skip_if_not(has("grpreg"))
  cvm <- cv.wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], J = 2:3,
                 foldid = folds, threshold = "max")
  expect_identical(names(on_null(plot(cvm))[["value"]]), "cv")
  expect_error(on_null(plot(cvm, which = "threshold")),
               "scores no candidate")
  out <- on_null(plot(cvm, which = c("cv", "path")))[["value"]]
  expect_identical(names(out), c("cv", "path"))
  expect_identical(out[["path"]][["t"]], cvm[["threshold"]][["t"]])
})

## ---------------------------------------------------------------------------
## plot.wafc()
## ---------------------------------------------------------------------------

test_that("plot.wafc draws the components of wafc_functions() at s", {
  s <- cvl[["lambda.min"]]
  r <- on_null(plot(fitl, s = s, truth = dgp, n_grid = 64))
  expect_false(r[["visible"]])
  cp <- r[["value"]][["components"]]
  fn <- wafc_functions(fitl, s = s, n_grid = 64)
  expect_identical(cp[["g"]], fn[["g"]])
  expect_identical(cp[["s"]], s)
  expect_identical(cp[["zero"]], fn[["nonzero"]] == 0L)
  expect_identical(cp[["norm"]], fn[["norm"]])
  expect_null(cp[["t"]])
  ## s = NULL is the smallest level of the path, as in wafc_functions()
  cp0 <- on_null(plot(fitl, n_grid = 32))[["value"]][["components"]]
  expect_identical(cp0[["s"]], min(fitl[["lambda"]]))
})

test_that("plot.wafc draws the path of the block norms", {
  out <- on_null(plot(fitl, which = "path"))[["value"]]
  expect_identical(names(out), "path")
  blk <- wafc_blocks(fitl, s = fitl[["lambda"]])
  des <- fitl[["design"]]
  for (l in 1:3) {
    for (m in 1:2) {
      expect_equal(out[["path"]][["norm"]][, wafc_block_name(des, l, m)],
                   unname(blk[["norm"]][l, m, ]), tolerance = 1e-12)
    }
  }
  expect_null(out[["path"]][["marked"]])
  s <- fitl[["lambda"]][20]
  both <- on_null(plot(fitl, which = c("path", "components"), s = s,
                       n_grid = 32))[["value"]]
  expect_identical(names(both), c("path", "components"))
  expect_identical(both[["path"]][["marked"]], c(s = s))
})

test_that("the grid of components restores the layout of the device", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev))
  graphics::par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)
  before <- graphics::par(c("mfrow", "mar", "oma"))
  plot(fitl, s = cvl[["lambda.min"]], n_grid = 32)
  expect_identical(graphics::par(c("mfrow", "mar", "oma")), before)
})

test_that("graphical parameters reach the panels, and bad arguments stop", {
  expect_silent(on_null(plot(fitl, which = "path", main = "my title",
                             col.axis = "grey30")))
  expect_error(on_null(plot(fitl, which = "cv")), "unknown panel")
  expect_error(on_null(plot(cvl, which = "everything")), "unknown panel")
  expect_error(on_null(plot(fitl, s = cvl[["lambda.min"]],
                            truth = dgp[["g"]][1:2, ], n_grid = 16)),
               "list of functions")
  bad <- dgp[["g"]]
  bad[[3L, 2L]] <- 1
  expect_error(on_null(plot(fitl, s = cvl[["lambda.min"]], truth = bad,
                            n_grid = 16)), "must be a function or NULL")
  expect_error(on_null(plot(fitl, which = "path", ask = NA)), "'ask'")
})
