## wafc/tests/test-competitors.R -- tests of the competitors of step E2.4
## and of the penalty rule the pilot adds. Run from the root of the
## repository:
##
##     Rscript -e 'testthat::test_dir("wafc/tests")'
##
## Four tests carry the step. The first is the one docs/instrucoes.md,
## section 6, makes the first of every estimator: with the truth inside the
## model and no noise, the method recovers it. The second checks that every
## competitor answers the three questions the pilot asks in the same
## coordinates as wafc(), which is what makes the table a comparison and
## not a list. The third checks the block LASSO of Klopp and Pensky against
## a group vector built by hand, since the whole point of that column is to
## be their penalty and not another one; since step E2.5c it also checks
## the weight of each chunk and the merged coarse levels against vectors
## built by hand, and the default against the call it replaced, to 1e-12;
## since step E2.5e, the balanced chunks against vectors built by hand at
## the chunk sizes 6 and 7 of the pilot, and the three earlier forms
## against a frozen copy of the grouping of before that step, the groups
## bit for bit and the fits to 1e-12; since step E2.5f, the free coarse
## levels against vectors built by hand at 6 and 7, the least squares fit
## of a resolution with nothing to penalize against lm() and against the
## scale of cv.grpreg, and the four earlier forms against a frozen copy of
## the grouping of step E2.5e, the same way.
## Since step E4.1c the oracle of the structure has the block LASSO too,
## tested as the coordinatewise one is (exact recovery), against
## cv.wafc() with every block active, and on the candidates of J that the
## cap of the levels makes equal.
## The fourth checks the quantile
## universal threshold against its own definition: it is the quantile of
## the smallest penalty level that kills the penalized part under the null,
## so simulating from the null and solving for that level has to reproduce
## it.

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
                     snr = 6)
x0 <- dgp[["x"]]
u0 <- dgp[["u"]]
y0 <- dgp[["y"]]
active0 <- nzchar(dgp[["structure"]])
folds <- rep_len(1:5, n)
grid0 <- cbind(seq(min(u0[, 1L]), max(u0[, 1L]), length.out = 64),
               seq(min(u0[, 2L]), max(u0[, 2L]), length.out = 64))

has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

## ---------------------------------------------------------------------------
## Exact recovery
## ---------------------------------------------------------------------------

test_that("every competitor recovers a model of its own form without noise", {
  ## A model that is inside the span of all the sieves at once: the
  ## coefficients vary with u through a single smooth function, so the
  ## spline, the wavelet and the tree method can all reproduce it, and
  ## without noise the recovery is a question about the plumbing and not
  ## about the basis.
  gg <- function(v) sin(2 * pi * v)
  beta <- cbind(1 + gg(u0[, 1L]), 2 + numeric(n), -1.5 + numeric(n))
  f <- as.numeric(rowSums(x0 * beta))
  for (mth in c("gam", "bsgl", "klopp", "aspline", "linear")) {
    if (mth %in% c("bsgl", "klopp") && !has("grpreg")) next
    if (mth == "gam" && !has("mgcv")) next
    ## df fixed so that the test does not run the whole grid of D36
    fit <- wafc_competitor(mth, x0, u0, f, foldid = folds, J = 4L, df = 16L)
    fh <- predict(fit, x0, u0)
    ## the linear fit cannot reproduce the varying part and is the control
    ## that the tolerance is not vacuous
    tol <- if (mth == "linear") 10 else 0.05 * stats::sd(f)
    if (mth == "linear") {
      expect_gt(sqrt(mean((fh - f)^2)), 0.1 * stats::sd(f))
    } else {
      expect_lt(sqrt(mean((fh - f)^2)), tol)
    }
  }
})

test_that("the oracle of the structure recovers a sparse truth in the basis", {
  d0 <- wafc_design(x0, u0, J = 3L, rescale = FALSE)
  theta <- numeric(d0[["nvars"]])
  theta[d0[["unpenalized"]]] <- c(1, 2, -1.5)
  theta[c(d0[["blocks"]][["x1:u1"]][c(1L, 4L)],
          d0[["blocks"]][["x2:u1"]][2L])] <- c(1.5, -1, 0.8)
  y_exact <- as.numeric(d0[["Z"]] %*% theta)
  act <- matrix(FALSE, p, q)
  act[1L, 1L] <- TRUE
  act[2L, 1L] <- TRUE
  fit <- wafc_competitor("oracle", x0, u0, y_exact, active = act,
                         foldid = folds, J = 3L, rescale = FALSE,
                         lambda = c(1e-1, 1e-3, 1e-6, 1e-9))
  expect_lt(sqrt(mean((predict(fit, x0, u0) - y_exact)^2)), 1e-4)
  ## it cannot have used a block it was told is zero
  expect_false(any(fit[["blocks"]][!act]))
  expect_identical(fit[["extra"]][["penalty"]], "lasso")
})

test_that("the block oracle recovers a sparse truth in the basis", {
  skip_if_not(has("grpreg"))
  ## the truth of the test above; the tolerance of grpreg is lowered from
  ## its 1e-4, which stops the path at an error of 4.5e-6
  d0 <- wafc_design(x0, u0, J = 3L, rescale = FALSE)
  theta <- numeric(d0[["nvars"]])
  theta[d0[["unpenalized"]]] <- c(1, 2, -1.5)
  theta[c(d0[["blocks"]][["x1:u1"]][c(1L, 4L)],
          d0[["blocks"]][["x2:u1"]][2L])] <- c(1.5, -1, 0.8)
  y_exact <- as.numeric(d0[["Z"]] %*% theta)
  act <- matrix(FALSE, p, q)
  act[1L, 1L] <- TRUE
  act[2L, 1L] <- TRUE
  fit <- wafc_competitor("oracle", x0, u0, y_exact, active = act,
                         foldid = folds, J = 3L, rescale = FALSE,
                         penalty = "block", thresh = 1e-8,
                         lambda = c(1e-1, 1e-3, 1e-6, 1e-9))
  expect_lt(sqrt(mean((predict(fit, x0, u0) - y_exact)^2)), 1e-6)
  cf <- coef.wafc(fit[["fit"]], s = fit[["extra"]][["lambda"]])[-1L, 1L]
  expect_lt(max(abs(cf - theta[fit[["fit"]][["design"]][["keep"]]])), 1e-6)
  expect_equal(unname(fit[["cc"]]), c(1, 2, -1.5), tolerance = 1e-6)
  expect_identical(unname(fit[["blocks"]]), act)
  ## the constant covariate carries the intercept of grpreg
  expect_identical(fit[["intercept"]], 0)
  expect_identical(fit[["extra"]][["penalty"]], "block")
  expect_identical(fit[["fit"]][["penalty"]], "block")
  expect_identical(fit[["extra"]][["block.size"]],
                   as.integer(ceiling(log(n))))
})

test_that("with every block active the block oracle is cv.wafc without threshold", {
  skip_if_not(has("grpreg"))
  fit <- wafc_fit_oracle(x0, u0, y0, active = matrix(TRUE, p, q), J = 2:4,
                         foldid = folds, penalty = "block")
  cv <- cv.wafc(x0, u0, y0, J = 2:4, foldid = folds, penalty = "block",
                threshold = "none")
  expect_identical(fit[["extra"]][["J"]], cv[["J.min"]])
  expect_identical(fit[["extra"]][["lambda"]], cv[["lambda.min"]])
  expect_identical(unname(fit[["cc"]]), unname(coef(cv)[1L + seq_len(p)]))
  expect_equal(fit[["fitted"]], as.numeric(predict(cv, x0, u0)),
               tolerance = 1e-12)
  expect_identical(unname(fit[["blocks"]]),
                   unname(wafc_blocks(cv)[["nonzero"]] > 0L))
  expect_identical(fit[["extra"]][["nzero"]],
                   as.integer(sum(wafc_blocks(cv)[["nonzero"]])))
  ## one row of the path of grpreg per candidate, the selected one marked
  conv <- fit[["extra"]][["conv"]]
  expect_s3_class(conv, "data.frame")
  expect_identical(conv[["J"]], 2:4)
  expect_identical(conv[["J"]][conv[["chosen"]]], cv[["J.min"]])
})

test_that("a candidate J capped to one already fitted is not refitted", {
  skip_if_not(has("grpreg"))
  ## modulators with 16 values have 15 distinct points on the circle, so
  ## J = 4 is built at 3 (step E3.4) and is the candidate J = 3 again, as
  ## J = 8 is the candidate J = 7 at n = 250 with continuous modulators
  u1 <- round(u0 * 15) / 15
  expect_identical(wafc_design(x0, u1, J = 4L)[["J"]],
                   wafc_design(x0, u1, J = 3L)[["J"]])
  for (pen in c("lasso", "block")) {
    a <- wafc_fit_oracle(x0, u1, y0, active = active0, J = 3L,
                         foldid = folds, penalty = pen)
    b <- wafc_fit_oracle(x0, u1, y0, active = active0, J = 3:4,
                         foldid = folds, penalty = pen)
    expect_identical(b[["extra"]][["J"]], 3L)
    expect_identical(b[["extra"]][["lambda"]], a[["extra"]][["lambda"]])
    expect_identical(b[["fitted"]], a[["fitted"]])
    if (pen == "block") {
      conv <- b[["extra"]][["conv"]]
      expect_identical(conv[["J"]], 3:4)
      expect_identical(conv[["chosen"]], c(TRUE, FALSE))
      same <- setdiff(names(conv), c("J", "chosen"))
      expect_identical(unlist(conv[2L, same]), unlist(conv[1L, same]))
    }
  }
})

test_that("the block oracle keeps the intercept of grpreg when no covariate is constant", {
  skip_if_not(has("grpreg"))
  d1 <- simulate_wafc(n, p = 3L, q = 2L, scenario = "smooth", seed = 3L,
                      intercept = FALSE)
  fit <- wafc_competitor("oracle", d1[["x"]], d1[["u"]], d1[["y"]],
                         active = nzchar(d1[["structure"]]), J = 2:4,
                         foldid = folds, penalty = "block")
  a0 <- fit[["intercept"]]
  expect_true(is.numeric(a0) && length(a0) == 1L && a0 != 0)
  expect_equal(predict(fit, d1[["x"]], d1[["u"]]),
               a0 + rowSums(d1[["x"]] * fit[["beta"]](d1[["u"]])),
               tolerance = 1e-12)
  expect_equal(predict(fit, d1[["x"]], d1[["u"]]), fit[["fitted"]],
               tolerance = 1e-12)
})

test_that("an absent, empty or ill-formed 'active' is an error, and no active block is least squares", {
  for (pen in c("lasso", "block")) {
    ora <- function(a) wafc_fit_oracle(x0, u0, y0, active = a, J = 3L,
                                       foldid = folds, penalty = pen)
    expect_error(ora(NULL), "needs 'active'")
    expect_error(ora(logical(0)), "needs 'active'")
    expect_error(ora(matrix(logical(0), 0L, 0L)), "needs 'active'")
    expect_error(ora(matrix(TRUE, q, p)), "must be 3 by 2 \\(p by q\\); it is 2 by 3")
    expect_error(ora(TRUE), "must have p q = 6 entries; it has 1")
    expect_error(ora(c(TRUE, NA, FALSE, FALSE, FALSE, FALSE)), "must not contain NA")
    expect_error(ora(matrix("a", p, q)), "logical p by q matrix")
    ## a matrix of FALSE is the null cell, where the oracle is the linear fit
    fit <- ora(matrix(FALSE, p, q))
    lin <- wafc_fit_linear(x0, u0, y0)
    expect_identical(fit[["fitted"]], lin[["fitted"]])
    expect_identical(fit[["extra"]][["note"]],
                     "no active block: least squares on x")
  }
})

## ---------------------------------------------------------------------------
## The common interface
## ---------------------------------------------------------------------------

test_that("every competitor answers in the coordinates of the model", {
  mths <- c("gam", "bsgl", "klopp", "aspline", "linear", "oracle",
            "oracle.block")
  if (has("VCBART")) mths <- c(mths, "vcbart")
  for (mth in mths) {
    if (mth %in% c("bsgl", "klopp", "oracle.block") && !has("grpreg")) next
    if (mth == "gam" && !has("mgcv")) next
    ## the oracle in blocks of step E4.1c is the oracle with penalty = "block"
    args <- list(sub("\\.block$", "", mth), x0, u0, y0, active = active0,
                 foldid = folds, J = 3L, df = 8L, burn = 100L, nd = 100L)
    if (mth == "oracle.block") args[["penalty"]] <- "block"
    fit <- do.call(wafc_competitor, args)
    expect_s3_class(fit, "wafc_competitor")
    ## beta(u) has one column per linear covariate, in the order of the design
    b <- fit[["beta"]](u0)
    expect_equal(dim(b), c(n, p))
    expect_equal(colnames(b), colnames(x0))
    ## the prediction is the model's own: rowSums(x * beta(u))
    expect_equal(predict(fit, x0, u0), as.numeric(rowSums(x0 * b)),
                 tolerance = 1e-8)
    ## predict on the training data agrees with the stored fitted values
    expect_equal(predict(fit, x0, u0), fit[["fitted"]], tolerance = 1e-6)
    ## the components, when the method has them, are centred on the grid and
    ## add up to beta minus its level
    g <- wafc_grid_components(fit, grid0)
    if (is.null(g)) {
      expect_true(mth == "vcbart")
      next
    }
    expect_equal(dim(g), c(p, q))
    expect_true(all(vapply(g, function(v) abs(mean(v)) < 1e-8, TRUE)))
    gb <- fit[["beta"]](grid0)
    for (l in seq_len(p)) {
      s <- Reduce(`+`, g[l, ])
      expect_equal(gb[, l] - mean(gb[, l]), s - mean(s), tolerance = 1e-6)
    }
    if (!is.null(fit[["blocks"]])) {
      expect_equal(dim(fit[["blocks"]]), c(p, q))
      expect_type(fit[["blocks"]], "logical")
    }
  }
})

test_that("a competitor ignores an argument meant for another one", {
  ## the pilot calls every method with one list of arguments
  expect_s3_class(wafc_competitor("linear", x0, u0, y0, foldid = folds,
                                  J = 3L, nfolds = 5L),
                  "wafc_competitor")
  expect_error(wafc_competitor("oracle", x0, u0, y0), "needs 'active'")
})

## ---------------------------------------------------------------------------
## The block LASSO of Klopp and Pensky
## ---------------------------------------------------------------------------

test_that("the chunks of the block LASSO are the ones of Klopp and Pensky", {
  skip_if_not(has("grpreg"))
  des <- wafc_design(x0, u0, J = 4L)
  bs <- 3L
  grp <- wafc_kp_groups(des, bs, penalize.levels = TRUE)
  ## every column belongs to exactly one group, and the p level terms are a
  ## group of their own, which is their norm (3.1) and not decision D3
  expect_length(grp, des[["nvars"]])
  expect_true(all(grp > 0L))
  expect_equal(length(unique(grp[des[["unpenalized"]]])), 1L)
  ## a chunk never straddles two resolution levels: inside a block the
  ## columns are ordered by increasing j (D12), so the sizes are
  ## min(2^j, bs) repeated as the level is cut
  idx <- des[["blocks"]][["x1:u1"]]
  sizes <- as.integer(table(grp[idx])[as.character(unique(grp[idx]))])
  want <- unlist(lapply(0:3, function(j) {
    nj <- 2^j
    diff(c(seq(0L, nj - 1L, by = bs), nj))
  }))
  expect_equal(sizes, as.integer(want))
  ## and the unpenalized version leaves the level terms out
  g0 <- wafc_kp_groups(des, bs, penalize.levels = FALSE)
  expect_true(all(g0[des[["unpenalized"]]] == 0L))
  expect_equal(max(g0), max(grp) - 1L)
})

test_that("the block LASSO zeroes a whole chunk at a time", {
  skip_if_not(has("grpreg"))
  fit <- wafc_competitor("klopp", x0, u0, y0, foldid = folds, J = 3L)
  des <- fit[["design"]]
  b <- as.numeric(stats::coef(fit[["fit"]]))[-1L]
  grp <- wafc_kp_groups(des, fit[["extra"]][["block.size"]],
                        fit[["extra"]][["penalize.levels"]])
  ## Only the wavelet chunks: the group of the level terms contains the
  ## constant covariate, which grpreg drops as a column of zero variance and
  ## leaves at zero whatever the rest of the group does. Columns that are
  ## identically zero in the design (they appear from J = 8 with the default
  ## margin) are out for the same reason.
  alive <- as.numeric(Matrix::colSums(abs(des[["Z"]]))) > 0
  pen <- setdiff(seq_len(des[["nvars"]]), des[["unpenalized"]])
  for (g in unique(grp[pen])) {
    v <- b[grp == g & alive]
    if (length(v) == 0L) next
    expect_true(all(v == 0) || all(v != 0))
  }
})

test_that("the weights and the merged coarse levels are the ones asked for", {
  skip_slow()
  ## Step E2.5c, open question 34 of docs/ESTADO.md. At J = 4 a block has
  ## the levels 0 to 3, with 1, 2, 4 and 8 columns in that order (D12).
  skip_if_not(has("grpreg"))
  des <- wafc_design(x0, u0, J = 4L)
  nb <- length(des[["blocks"]])
  by_hand <- function(one, levels.group) {
    grp <- integer(des[["nvars"]])
    off <- as.integer(levels.group)
    if (levels.group) grp[des[["unpenalized"]]] <- 1L
    for (b in seq_len(nb)) {
      grp[des[["blocks"]][[b]]] <- one + off + (b - 1L) * max(one)
    }
    grp
  }
  ## chunks of 7: level 3 is cut into 7 + 1; merged, the levels 0 to 2
  ## (1 + 2 + 4 = 7 columns) are one chunk and level 3 is cut as before
  sep7 <- c(1L, 2L, 2L, 3L, 3L, 3L, 3L, rep(4L, 7L), 5L)
  mrg7 <- c(rep(1L, 7L), rep(2L, 7L), 3L)
  expect_identical(wafc_kp_groups(des, 7L, FALSE), by_hand(sep7, FALSE))
  expect_identical(wafc_kp_groups(des, 7L, FALSE, merge.coarse = TRUE),
                   by_hand(mrg7, FALSE))
  expect_identical(wafc_kp_groups(des, 7L, TRUE, merge.coarse = TRUE),
                   by_hand(mrg7, TRUE))
  ## the merge takes 2^j strictly below the chunk size: with chunks of 4,
  ## level 2 is a full chunk and stays one, and level 3 is two of them
  mrg4 <- c(1L, 1L, 1L, rep(2L, 4L), rep(3L, 4L), rep(4L, 4L))
  expect_identical(wafc_kp_groups(des, 4L, FALSE, merge.coarse = TRUE),
                   by_hand(mrg4, FALSE))
  ## a chunk size of 1 has no level below it, and the merge changes nothing
  expect_identical(wafc_kp_groups(des, 1L, FALSE, merge.coarse = TRUE),
                   wafc_kp_groups(des, 1L, FALSE))

  ## The weights, read from the multiplier grpreg records it used: sqrt(|G|)
  ## by default, one with chunk.weights = "unit", and sqrt(|G|) of the
  ## merged chunks with merge.coarse = TRUE.
  kp <- function(...) {
    wafc_fit_klopp(x0, u0, y0, J = 4L, block.size = 7L,
                   penalize.levels = FALSE, foldid = folds, ...)
  }
  mult <- function(fit) unname(as.numeric(fit[["fit"]][["fit"]][["group.multiplier"]]))
  f_sqrt <- kp()
  f_unit <- kp(chunk.weights = "unit")
  f_mrg <- kp(merge.coarse = TRUE)
  expect_equal(mult(f_sqrt), sqrt(rep(c(1, 2, 4, 7, 1), nb)))
  expect_equal(mult(f_unit), rep(1, 5L * nb))
  expect_equal(mult(f_mrg), sqrt(rep(c(7, 7, 1), nb)))
  expect_identical(f_mrg[["extra"]][["ngroups"]], 3L * nb)
  expect_identical(f_unit[["extra"]][["chunk.weights"]], "unit")
  ## match.arg() translates its message, so only the error is checked
  expect_error(kp(chunk.weights = "log"))
  ## different weights are a different estimator
  expect_false(isTRUE(all.equal(as.numeric(stats::coef(f_unit[["fit"]])),
                                as.numeric(stats::coef(f_sqrt[["fit"]])))))

  ## The default is the fit of before the two arguments existed, which
  ## called grpreg with no multiplier on the groups of wafc_kp_groups(); this
  ## is what keeps every number of steps E2.5a and E2.5b.
  for (pl in c(TRUE, FALSE)) {
    fit <- wafc_fit_klopp(x0, u0, y0, J = 3:4, penalize.levels = pl,
                          foldid = folds)
    des <- fit[["design"]]
    ref <- grpreg::cv.grpreg(as.matrix(des[["Z"]]), y0,
                             group = wafc_kp_groups(des, fit[["extra"]][["block.size"]], pl),
                             penalty = "grLasso", fold = folds)
    expect_equal(as.numeric(stats::coef(fit[["fit"]])),
                 as.numeric(stats::coef(ref)), tolerance = 1e-12)
    expect_equal(fit[["fit"]][["cve"]], ref[["cve"]], tolerance = 1e-12)
    expect_equal(fit[["fitted"]],
                 as.numeric(predict(ref, as.matrix(des[["Z"]]))),
                 tolerance = 1e-12)
  }
})

## The grouping of wafc_kp_groups() as it was at the end of step E2.5c,
## frozen here so that the forms measured in steps E2.5a to E2.5c can be
## checked against it once the function has the balanced chunks.
kp_groups_e25c <- function(design, block.size, penalize.levels,
                           merge.coarse = FALSE) {
  nvars <- design[["nvars"]]
  grp <- integer(nvars)
  g <- 0L
  unp <- design[["unpenalized"]]
  if (penalize.levels) {
    g <- g + 1L
    grp[unp] <- g
  }
  j0 <- design[["j0"]]
  for (nm in names(design[["blocks"]])) {
    idx <- design[["blocks"]][[nm]]
    pos <- 0L
    Jm <- j0 + as.integer(round(log2(length(idx) + 2^j0)))
    levs <- j0:(Jm - 1L)
    if (merge.coarse) {
      coarse <- levs[2^levs < block.size]
      if (length(coarse) > 0L) {
        pos <- as.integer(sum(2^coarse))
        g <- g + 1L
        grp[idx[seq_len(pos)]] <- g
        levs <- setdiff(levs, coarse)
      }
    }
    for (j in levs) {
      nj <- 2^j
      lev <- idx[pos + seq_len(nj)]
      pos <- pos + nj
      for (start in seq(1L, nj, by = block.size)) {
        g <- g + 1L
        grp[lev[start:min(start + block.size - 1L, nj)]] <- g
      }
    }
  }
  grp
}

test_that("the balanced chunks are the ones asked for", {
  ## Step E2.5e. At J = 6 a block has the levels 0 to 5, with 1, 2, 4, 8,
  ## 16 and 32 columns in that order (D12).
  des <- wafc_design(x0, u0, J = 6L)
  nb <- length(des[["blocks"]])
  by_hand <- function(sizes, levels.group, d = des) {
    one <- rep(seq_along(sizes), sizes)
    grp <- integer(d[["nvars"]])
    off <- as.integer(levels.group)
    if (levels.group) grp[d[["unpenalized"]]] <- 1L
    for (b in seq_along(d[["blocks"]])) {
      grp[d[["blocks"]][[b]]] <- one + off + (b - 1L) * length(sizes)
    }
    grp
  }
  ## b_n = 6, the chunk size at n = 250: levels 0 to 2 are one chunk of
  ## 1 + 2 + 4 = 7; level 3 is 6 + 2, one chunk of 8 once the 2 is
  ## absorbed; level 4 is 6 + 6 + 4, so 6 and 10; level 5 is 5 x 6 + 2, so
  ## four of 6 and one of 8
  bal6 <- c(7L, 8L, 6L, 10L, 6L, 6L, 6L, 6L, 8L)
  ## b_n = 7, at n = 500 and 1000: level 3 is 7 + 1, so 8; level 4 is
  ## 7 + 7 + 2, so 7 and 9; level 5 is 4 x 7 + 4, so three of 7 and one
  ## of 11
  bal7 <- c(7L, 8L, 7L, 9L, 7L, 7L, 7L, 11L)
  expect_identical(wafc_kp_groups(des, 6L, FALSE, balanced = TRUE),
                   by_hand(bal6, FALSE))
  expect_identical(wafc_kp_groups(des, 7L, FALSE, balanced = TRUE),
                   by_hand(bal7, FALSE))
  expect_identical(wafc_kp_groups(des, 6L, TRUE, balanced = TRUE),
                   by_hand(bal6, TRUE))
  expect_identical(wafc_kp_groups(des, 7L, TRUE, balanced = TRUE),
                   by_hand(bal7, TRUE))
  ## balanced merges the coarse levels whatever merge.coarse says
  expect_identical(wafc_kp_groups(des, 7L, FALSE, merge.coarse = TRUE,
                                  balanced = TRUE),
                   wafc_kp_groups(des, 7L, FALSE, balanced = TRUE))
  ## the two cases in which the coarse chunk falls below the chunk size: a
  ## power of two (4: levels 0 and 1, 3 columns, and every finer level cut
  ## exactly), and 2^J at most the chunk size (J = 2 with 6: one chunk of 3)
  expect_identical(wafc_kp_groups(des, 4L, FALSE, balanced = TRUE),
                   by_hand(c(3L, rep(4L, 15L)), FALSE))
  d2 <- wafc_design(x0, u0, J = 2L)
  expect_identical(wafc_kp_groups(d2, 6L, FALSE, balanced = TRUE),
                   by_hand(3L, FALSE, d2))
  ## when the chunk size divides every finer level there is nothing to
  ## absorb, and with chunks of one there is nothing to merge either
  expect_identical(wafc_kp_groups(des, 8L, FALSE, balanced = TRUE),
                   wafc_kp_groups(des, 8L, FALSE, merge.coarse = TRUE))
  expect_identical(wafc_kp_groups(des, 1L, FALSE, balanced = TRUE),
                   wafc_kp_groups(des, 1L, FALSE))
  ## the sizes the roxygen of wafc_fit_klopp() states, for every chunk size
  ## from 2 to 70 at J = 6: the coarse chunk has 2^(j* + 1) - 1 columns,
  ## j* the finest level with 2^j below the chunk size, and every chunk of
  ## a finer level has between the chunk size and twice it minus one
  idx <- des[["blocks"]][[1L]]
  as_stated <- vapply(2:70, function(bs) {
    grp <- wafc_kp_groups(des, bs, FALSE, balanced = TRUE)[idx]
    sizes <- as.integer(table(grp)[as.character(unique(grp))])
    js <- max(which(2^(0:5) < bs)) - 1L
    sizes[1L] == 2^(js + 1L) - 1L &&
      all(sizes[-1L] >= bs & sizes[-1L] <= 2L * bs - 1L)
  }, logical(1))
  expect_true(all(as_stated))

  ## The weights grpreg uses are sqrt(|G|) of the balanced chunks: at J = 5
  ## and chunks of 7, the coarse 7, level 3 as 8 and level 4 as 7 + 9.
  skip_if_not(has("grpreg"))
  f_bal <- wafc_fit_klopp(x0, u0, y0, J = 5L, block.size = 7L,
                          penalize.levels = FALSE, foldid = folds,
                          balanced = TRUE)
  mult <- unname(as.numeric(f_bal[["fit"]][["fit"]][["group.multiplier"]]))
  expect_equal(mult, sqrt(rep(c(7, 8, 7, 9), nb)))
  expect_identical(f_bal[["extra"]][["ngroups"]], 4L * nb)
  expect_true(f_bal[["extra"]][["balanced"]])
  expect_true(f_bal[["extra"]][["merge.coarse"]])
})

test_that("the forms of before the balanced chunks do not move", {
  skip_slow()
  ## The grouping of the three earlier forms is the frozen one, bit for bit,
  ## over resolutions, chunk sizes and both treatments of the levels.
  same <- logical(0)
  for (J in 2:7) {
    des <- wafc_design(x0, u0, J = J)
    for (bs in c(1:9, 16L, 64L)) {
      for (pl in c(TRUE, FALSE)) {
        for (mc in c(FALSE, TRUE)) {
          same <- c(same,
                    identical(wafc_kp_groups(des, bs, pl, merge.coarse = mc),
                              kp_groups_e25c(des, bs, pl, merge.coarse = mc)))
        }
      }
    }
  }
  expect_length(same, 264L)
  expect_true(all(same))
  ## and the fits of the three forms are grpreg on the frozen groups, to
  ## 1e-12, with the levels penalized (klopp) and free (the forms of steps
  ## E2.5b and E2.5c)
  skip_if_not(has("grpreg"))
  forms <- list(sqrt = list(), unit = list(chunk.weights = "unit"),
                merged = list(merge.coarse = TRUE))
  for (nm in names(forms)) {
    for (pl in c(TRUE, FALSE)) {
      fit <- do.call(wafc_fit_klopp,
                     c(list(x0, u0, y0, J = 3:4, penalize.levels = pl,
                            foldid = folds), forms[[nm]]))
      des <- fit[["design"]]
      grp <- kp_groups_e25c(des, fit[["extra"]][["block.size"]], pl,
                            merge.coarse = isTRUE(forms[[nm]][["merge.coarse"]]))
      args <- list(as.matrix(des[["Z"]]), y0, group = grp,
                   penalty = "grLasso", fold = folds)
      if (nm == "unit") args[["group.multiplier"]] <- rep(1, max(grp))
      ref <- do.call(grpreg::cv.grpreg, args)
      expect_equal(as.numeric(stats::coef(fit[["fit"]])),
                   as.numeric(stats::coef(ref)), tolerance = 1e-12)
      expect_equal(fit[["fit"]][["cve"]], ref[["cve"]], tolerance = 1e-12)
      expect_false(fit[["extra"]][["balanced"]])
    }
  }
})

## The grouping of wafc_kp_groups() as it was at the end of step E2.5e,
## frozen here so that the four forms measured in steps E2.5a to E2.5e can
## be checked against it once the function has the free coarse levels.
kp_groups_e25e <- function(design, block.size, penalize.levels,
                            merge.coarse = FALSE, balanced = FALSE) {
  nvars <- design[["nvars"]]
  grp <- integer(nvars)
  g <- 0L
  unp <- design[["unpenalized"]]
  if (penalize.levels) {
    g <- g + 1L
    grp[unp] <- g
  }
  j0 <- design[["j0"]]
  for (nm in names(design[["blocks"]])) {
    idx <- design[["blocks"]][[nm]]
    pos <- 0L
    Jm <- j0 + as.integer(round(log2(length(idx) + 2^j0)))
    levs <- j0:(Jm - 1L)
    if (merge.coarse || balanced) {
      coarse <- levs[2^levs < block.size]
      if (length(coarse) > 0L) {
        pos <- as.integer(sum(2^coarse))
        g <- g + 1L
        grp[idx[seq_len(pos)]] <- g
        levs <- setdiff(levs, coarse)
      }
    }
    for (j in levs) {
      nj <- 2^j
      lev <- idx[pos + seq_len(nj)]
      pos <- pos + nj
      starts <- seq(1L, nj, by = block.size)
      if (balanced && length(starts) > 1L && nj %% block.size != 0) {
        starts <- starts[-length(starts)]
      }
      ends <- c(starts[-1L] - 1L, nj)
      for (i in seq_along(starts)) {
        g <- g + 1L
        grp[lev[starts[i]:ends[i]]] <- g
      }
    }
  }
  grp
}

test_that("the free coarse levels are the ones asked for", {
  skip_slow()
  ## Step E2.5f, open question 37. At J = 6 a block has the levels 0 to 5,
  ## with 1, 2, 4, 8, 16 and 32 columns in that order (D12); with chunks of
  ## 6 or 7 the levels 0 to 2 (7 columns) are coarse, and with free.coarse
  ## they go to group 0 while the finer levels keep the balanced chunks of
  ## step E2.5e.
  des <- wafc_design(x0, u0, J = 6L)
  nb <- length(des[["blocks"]])
  by_hand <- function(ncoarse, sizes, levels.group, d = des) {
    one <- c(integer(ncoarse), rep(seq_along(sizes), sizes))
    grp <- integer(d[["nvars"]])
    off <- as.integer(levels.group)
    if (levels.group) grp[d[["unpenalized"]]] <- 1L
    for (b in seq_along(d[["blocks"]])) {
      grp[d[["blocks"]][[b]]] <- ifelse(one > 0L,
                                        one + off + (b - 1L) * length(sizes),
                                        0L)
    }
    grp
  }
  ## b_n = 6: the balanced chunks of before, 7 | 8 | 6, 10 | 6, 6, 6, 6, 8,
  ## without the coarse 7; b_n = 7: 7 | 8 | 7, 9 | 7, 7, 7, 11, likewise
  fc6 <- c(8L, 6L, 10L, 6L, 6L, 6L, 6L, 8L)
  fc7 <- c(8L, 7L, 9L, 7L, 7L, 7L, 11L)
  expect_identical(wafc_kp_groups(des, 6L, FALSE, balanced = TRUE,
                                  free.coarse = TRUE),
                   by_hand(7L, fc6, FALSE))
  expect_identical(wafc_kp_groups(des, 7L, FALSE, balanced = TRUE,
                                  free.coarse = TRUE),
                   by_hand(7L, fc7, FALSE))
  expect_identical(wafc_kp_groups(des, 6L, TRUE, balanced = TRUE,
                                  free.coarse = TRUE),
                   by_hand(7L, fc6, TRUE))
  expect_identical(wafc_kp_groups(des, 7L, TRUE, balanced = TRUE,
                                  free.coarse = TRUE),
                   by_hand(7L, fc7, TRUE))
  ## the free coarse levels override the merge, and without 'balanced' the
  ## finer levels keep the short piece at the end: 7 + 1, 7 + 7 + 2, ...
  expect_identical(wafc_kp_groups(des, 7L, FALSE, merge.coarse = TRUE,
                                  balanced = TRUE, free.coarse = TRUE),
                   by_hand(7L, fc7, FALSE))
  expect_identical(wafc_kp_groups(des, 7L, FALSE, free.coarse = TRUE),
                   by_hand(7L, c(7L, 1L, 7L, 7L, 2L, 7L, 7L, 7L, 7L, 4L),
                           FALSE))
  ## the columns in group 0 are the p level terms (when free) and the
  ## p q (2^(j* + 1) - 1) coarse ones, j* the finest level with 2^j < b_n:
  ## that is p_0 of section 12 of derivations/08a-sondagem-blocos.md
  for (bs in c(2L, 3L, 4L, 6L, 7L, 9L, 17L)) {
    js <- max(which(2^(0:5) < bs)) - 1L
    g <- wafc_kp_groups(des, bs, FALSE, balanced = TRUE, free.coarse = TRUE)
    expect_identical(sum(g == 0L), as.integer(p + p * q * (2^(js + 1L) - 1L)))
  }
  ## a block with no level of 2^j at least the chunk size is wholly free: at
  ## J = 3 with 6 or 7 nothing is penalized, and at J = 4 with 7 level 3 is
  ## one chunk of 8
  d3 <- wafc_design(x0, u0, J = 3L)
  expect_true(all(wafc_kp_groups(d3, 6L, FALSE, balanced = TRUE,
                                 free.coarse = TRUE) == 0L))
  d4 <- wafc_design(x0, u0, J = 4L)
  expect_identical(wafc_kp_groups(d4, 7L, FALSE, balanced = TRUE,
                                  free.coarse = TRUE),
                   by_hand(7L, 8L, FALSE, d4))

  skip_if_not(has("grpreg"))
  ## The fit is grpreg on those groups, with the sqrt(|G|) of grpreg on the
  ## fine chunks alone: at J = 5 and chunks of 7, 8 | 7, 9 per block.
  fc <- function(J, bs = 7L) {
    wafc_fit_klopp(x0, u0, y0, J = J, block.size = bs,
                   penalize.levels = FALSE, foldid = folds, balanced = TRUE,
                   free.coarse = TRUE)
  }
  f5 <- fc(5L)
  mult <- unname(as.numeric(f5[["fit"]][["fit"]][["group.multiplier"]]))
  expect_equal(mult, sqrt(rep(c(8, 7, 9), nb)))
  expect_identical(f5[["extra"]][["ngroups"]], 3L * nb)
  expect_true(f5[["extra"]][["free.coarse"]])
  d5 <- f5[["design"]]
  grp5 <- wafc_kp_groups(d5, 7L, FALSE, balanced = TRUE, free.coarse = TRUE)
  ref <- grpreg::cv.grpreg(as.matrix(d5[["Z"]]), y0, group = grp5,
                           penalty = "grLasso", fold = folds)
  expect_equal(as.numeric(stats::coef(f5[["fit"]])),
               as.numeric(stats::coef(ref)), tolerance = 1e-12)
  ## the coarse coefficients are never zero, so every block is in the
  ## estimate; the selection is the penalized part, and nzero counts it
  b5 <- as.numeric(stats::coef(f5[["fit"]]))[-1L]
  expect_true(all(f5[["blocks"]]))
  for (nm in names(d5[["blocks"]])) {
    idx <- d5[["blocks"]][[nm]]
    expect_true(all(b5[idx[1:7]] != 0))
  }
  fine <- unlist(lapply(d5[["blocks"]], function(idx) idx[-(1:7)]))
  expect_identical(f5[["extra"]][["nzero"]], sum(b5[fine] != 0))
  bf <- vapply(d5[["blocks"]], function(idx) any(b5[idx[-(1:7)]] != 0), TRUE)
  expect_identical(as.logical(t(f5[["extra"]][["blocks.fine"]])),
                   unname(bf))

  ## With nothing to penalize the fit is least squares, scored by the
  ## squared error cv.grpreg reports: exactly the least squares of each
  ## fold, and what grpreg gives at a vanishing penalty (the last of the
  ## two levels) on one extra column of noise.
  f3 <- fc(3L)
  d3 <- f3[["design"]]
  Z3 <- as.matrix(d3[["Z"]])
  expect_true(is.na(f3[["extra"]][["lambda"]]))
  expect_identical(f3[["extra"]][["ngroups"]], 0L)
  expect_identical(f3[["extra"]][["nzero"]], 0L)
  ls <- stats::lm(y0 ~ Z3)
  expect_equal(f3[["fitted"]], unname(stats::fitted(ls)), tolerance = 1e-10)
  err <- numeric(n)
  for (k in unique(folds)) {
    o <- folds == k
    fk <- stats::lm(y0[!o] ~ Z3[!o, ])
    bk <- stats::coef(fk)
    bk[is.na(bk)] <- 0
    err[o] <- (y0[o] - cbind(1, Z3[o, ]) %*% bk)^2
  }
  expect_equal(f3[["extra"]][["cve"]], mean(err), tolerance = 1e-10)
  set.seed(5)
  zn <- stats::rnorm(n)
  g0 <- grpreg::cv.grpreg(cbind(Z3, zn), y0, group = c(integer(ncol(Z3)), 1L),
                          penalty = "grLasso", fold = folds,
                          lambda = c(1e-6, 1e-9), eps = 1e-10)
  expect_equal(g0[["cve"]][2L],
               wafc_kp_cv_ols(cbind(Z3, zn), y0, folds)[["cve"]],
               tolerance = 1e-8)
  ## and the search over J compares the two kinds of candidate on that one
  ## scale: J = 3:5 keeps the one with the smaller error
  f35 <- fc(3:5)
  expect_identical(f35[["extra"]][["J"]],
                   c(3L, 5L)[which.min(c(f3[["extra"]][["cve"]],
                                         f5[["extra"]][["cve"]]))])
})

test_that("the forms of before the free coarse levels do not move", {
  skip_slow()
  ## The four earlier forms group as the frozen copy of step E2.5e, bit for
  ## bit, over resolutions, chunk sizes, both treatments of the levels, the
  ## merge and the balance; and their fits are grpreg on the frozen groups
  ## to 1e-12.
  same <- logical(0)
  for (J in 2:7) {
    des <- wafc_design(x0, u0, J = J)
    for (bs in c(1:9, 16L, 64L)) {
      for (pl in c(TRUE, FALSE)) {
        for (mc in c(FALSE, TRUE)) {
          for (bal in c(FALSE, TRUE)) {
            same <- c(same,
                      identical(wafc_kp_groups(des, bs, pl, merge.coarse = mc,
                                               balanced = bal),
                                kp_groups_e25e(des, bs, pl, merge.coarse = mc,
                                               balanced = bal)),
                      identical(wafc_kp_groups(des, bs, pl, mc, bal, FALSE),
                                kp_groups_e25e(des, bs, pl, mc, bal)))
          }
        }
      }
    }
  }
  expect_length(same, 1056L)
  expect_true(all(same))
  skip_if_not(has("grpreg"))
  forms <- list(sqrt = list(), unit = list(chunk.weights = "unit"),
                merged = list(merge.coarse = TRUE),
                balanced = list(balanced = TRUE))
  for (nm in names(forms)) {
    for (pl in c(TRUE, FALSE)) {
      fit <- do.call(wafc_fit_klopp,
                     c(list(x0, u0, y0, J = 4:5, penalize.levels = pl,
                            foldid = folds), forms[[nm]]))
      des <- fit[["design"]]
      grp <- kp_groups_e25e(des, fit[["extra"]][["block.size"]], pl,
                            merge.coarse = isTRUE(forms[[nm]][["merge.coarse"]]),
                            balanced = isTRUE(forms[[nm]][["balanced"]]))
      args <- list(as.matrix(des[["Z"]]), y0, group = grp,
                   penalty = "grLasso", fold = folds)
      if (nm == "unit") args[["group.multiplier"]] <- rep(1, max(grp))
      ref <- do.call(grpreg::cv.grpreg, args)
      expect_equal(as.numeric(stats::coef(fit[["fit"]])),
                   as.numeric(stats::coef(ref)), tolerance = 1e-12)
      expect_equal(fit[["fit"]][["cve"]], ref[["cve"]], tolerance = 1e-12)
      expect_equal(fit[["fitted"]],
                   as.numeric(predict(ref, as.matrix(des[["Z"]]))),
                   tolerance = 1e-12)
      expect_false(fit[["extra"]][["free.coarse"]])
      ## with no column of a block in group 0, the two readings of the
      ## selection coincide
      expect_identical(fit[["extra"]][["blocks.fine"]], fit[["blocks"]])
    }
  }
})

## ---------------------------------------------------------------------------
## The quantile universal threshold
## ---------------------------------------------------------------------------

test_that("the QUT is the quantile of the penalty level that kills the fit", {
  skip_slow()
  des <- wafc_design(x0, u0, J = 3L)
  lam <- wafc_lambda_qut(des, y0, alpha = 0.05, nsim = 400L, seed = 11L)
  expect_true(is.finite(lam) && lam > 0)
  expect_equal(lam, wafc_lambda_qut(des, y0, alpha = 0.05, nsim = 400L,
                                    seed = 11L))
  ## the rule is free of sigma: multiplying the response by a constant (on
  ## the residual of the unpenalized part, which is what carries the scale)
  ## multiplies the penalty level by the same constant
  unp <- des[["unpenalized"]]
  W <- as.matrix(des[["Z"]][, unp, drop = FALSE])
  r0 <- as.numeric(y0 - W %*% qr.coef(qr(W), y0))
  y2 <- y0 + 4 * r0
  expect_equal(wafc_lambda_qut(des, y2, nsim = 400L, seed = 11L) / lam, 5,
               tolerance = 1e-8)
  ## and it is the quantile of a level that really kills the penalized part:
  ## fitting a null sample at the level the rule returns for that sample
  ## leaves no non-zero wavelet coefficient in at least 1 - alpha of the
  ## draws, by construction of the statistic
  kill <- logical(40L)
  for (i in seq_along(kill)) {
    yy <- as.numeric(W %*% c(1, 2, -1.5)) + stats::rnorm(n, sd = 0.5)
    li <- wafc_lambda_qut(des, yy, alpha = 0.05, nsim = 200L, seed = 100L + i)
    fi <- wafc(design = des, y = yy,
               lambda = wafc_path_to(li, des, yy), penalty = "lasso")
    kill[i] <- fi[["nzero"]][length(fi[["lambda"]])] == 0L
  }
  expect_gt(mean(kill), 0.85)
})

test_that("the QUT is more conservative than the cross-validated penalty", {
  ## the property the pilot reports: the rule is calibrated on the null, so
  ## it penalizes more than the rule calibrated on prediction
  cv <- cv.wafc(x0, u0, y0, J = 3L, foldid = folds,
                penalty = "lasso", threshold = "none")
  lam <- wafc_lambda_qut(cv[["wafc.fit"]][["design"]], y0, nsim = 400L,
                         seed = 11L)
  expect_gt(lam, cv[["lambda.min"]])
})

## ---------------------------------------------------------------------------
## The metric the pilot reads
## ---------------------------------------------------------------------------

test_that("wafc_grid_components gives the same components as wafc_functions", {
  fit <- wafc(x0, u0, y0, J = 3L, penalty = "lasso")
  s <- fit[["lambda"]][40L]
  fn <- wafc_functions(fit, s = s, grid = grid0)
  g <- wafc_grid_components(fit, grid0, s = s)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      a <- fn[["g"]][[l, m]]
      expect_equal(g[[l, m]], a - mean(a), tolerance = 1e-10)
    }
  }
})

test_that("the components of a fit reproduce a truth that is in the basis", {
  ## the integrated squared error the pilot reports is zero when the fit is
  ## the truth, which is what says the grid, the centring and the range are
  ## the same on both sides
  d0 <- wafc_design(x0, u0, J = 3L, rescale = FALSE)
  theta <- numeric(d0[["nvars"]])
  theta[d0[["unpenalized"]]] <- c(1, 2, -1.5)
  theta[d0[["blocks"]][["x1:u1"]][2L]] <- 1.5
  y_exact <- as.numeric(d0[["Z"]] %*% theta)
  fit <- wafc(design = d0, y = y_exact, lambda = c(1e-3, 1e-8),
              penalty = "lasso")
  g <- wafc_grid_components(fit, grid0, s = 1e-8)
  d_grid <- wafc_design(matrix(1, nrow(grid0), p), grid0, spec = d0)
  truth <- as.numeric(d_grid[["Z"]][, d0[["blocks"]][["x1:u1"]][2L]]) * 1.5
  expect_equal(g[[1L, 1L]], truth - mean(truth), tolerance = 1e-4)
  expect_lt(max(abs(unlist(g[-1L]))), 1e-4)
})

## ---------------------------------------------------------------------------
## Corrections of the review of 2026-09-28
## ---------------------------------------------------------------------------

test_that("an argument of the basis reaches only the fitters that build one", {
  skip_slow()
  ## The pilot passes 'wavelet.table' to every method (decision D31). It
  ## used to reach the engine of vcbart, which stopped with "unused
  ## argument" and removed the column from the table.
  tb <- WaveBased::wtable(family = "Daublets", filter.size = 8L,
                          prec.wavelet = 30L, check = FALSE)
  if (has("VCBART")) {
    expect_s3_class(wafc_competitor("vcbart", x0, u0, y0, burn = 20L,
                                    nd = 20L, wavelet.table = tb),
                    "wafc_competitor")
  }
  if (has("grpreg")) {
    kp <- wafc_competitor("klopp", x0, u0, y0, foldid = folds, J = 3L,
                          wavelet.table = tb)
    expect_false(is.null(kp[["design"]][["wavelet.table"]]))
    expect_s3_class(wafc_competitor("bsgl", x0, u0, y0, foldid = folds,
                                    df = 8L, wavelet.table = tb),
                    "wafc_competitor")
  }
  if (has("mgcv")) {
    expect_s3_class(wafc_competitor("gam", x0, u0, y0, wavelet.table = tb),
                    "wafc_competitor")
  }
})

test_that("bsgl and klopp keep the engine intercept out of beta", {
  ## Without a constant covariate grpreg still fits an intercept, which no
  ## level can carry. It used to be added to beta_1, so every prediction
  ## multiplied it by X_1.
  skip_if_not(has("grpreg"))
  d1 <- simulate_wafc(n, p = p, q = q, scenario = "smooth", seed = 11L,
                      snr = 6, intercept = FALSE)
  y1 <- d1[["y"]] + 5
  for (mth in c("bsgl", "klopp")) {
    fit <- wafc_competitor(mth, d1[["x"]], d1[["u"]], y1, foldid = folds,
                           df = 8L, J = 3L)
    expect_equal(fit[["intercept"]], 5, tolerance = 0.2)
    expect_equal(predict(fit, d1[["x"]], d1[["u"]]), fit[["fitted"]],
                 tolerance = 1e-8)
    expect_equal(predict(fit, d1[["x"]], d1[["u"]]),
                 as.numeric(fit[["intercept"]] +
                              rowSums(d1[["x"]] * fit[["beta"]](d1[["u"]]))),
                 tolerance = 1e-8)
    expect_lt(sqrt(mean((fit[["fitted"]] - d1[["f"]] - 5)^2)),
              0.5 * stats::sd(d1[["f"]]))
  }
  ## with a constant covariate the intercept is folded into its level
  fit <- wafc_competitor("bsgl", x0, u0, y0, foldid = folds, df = 8L)
  expect_equal(fit[["intercept"]], 0)
})

test_that("the knot search scores the candidates as a refit would", {
  ## Reference: the search as it was written first, one B-spline basis and
  ## one least squares fit per candidate knot.
  refit_search <- function(r, w, v, rng, cand, max.knots, degree, n) {
    basis <- function(kn) {
      as.matrix(splines::bs(pmin(pmax(v, rng[1L]), rng[2L]), knots = kn,
                            degree = degree, Boundary.knots = rng,
                            intercept = TRUE))
    }
    score <- function(kn) {
      fk <- stats::.lm.fit(basis(kn) * w, r)
      rss <- sum(fk[["residuals"]]^2)
      list(bic = n * log(max(rss, .Machine[["double.eps"]]) / n) +
             fk[["rank"]] * log(n), rss = rss, knots = kn)
    }
    kn <- numeric(0)
    path <- list(score(kn))
    while (length(kn) < max.knots) {
      pool <- setdiff(cand, kn)
      if (length(pool) == 0L) break
      best <- NULL
      for (k in pool) {
        sc <- score(sort(c(kn, k)))
        if (is.null(best) || sc[["rss"]] < best[["rss"]]) best <- sc
      }
      kn <- best[["knots"]]
      path[[length(path) + 1L]] <- best
    }
    path[[which.min(vapply(path, `[[`, 0, "bic"))]][["knots"]]
  }
  for (sc in c("inhomogeneous", "uneven")) {
    d1 <- simulate_wafc(n, p = p, q = q, scenario = sc, seed = 5L, snr = 3)
    r <- as.numeric(d1[["y"]] - d1[["x"]] %*%
                      qr.coef(qr(d1[["x"]]), d1[["y"]]))
    for (l in c(1L, 2L)) {
      v <- d1[["u"]][, 1L]
      rng <- range(v)
      cand <- unique(stats::quantile(v, seq_len(20L) / 21, names = FALSE))
      a <- wafc_knot_search(r, d1[["x"]][, l], v, rng, cand, 8L, 3L, n)
      expect_equal(a[["knots"]],
                   refit_search(r, d1[["x"]][, l], v, rng, cand, 8L, 3L, n))
    }
  }
})

test_that("bsgl searches the grid of the WAFC and respects discrete modulators", {
  ## decision D36: the candidate dimensions are 2^J for J = 2, ..., 8, the
  ## grid of decision D34
  expect_equal(eval(formals(wafc_fit_bsgl)[["df"]]), 2L^(2:8))
  skip_if_not(has("grpreg"))
  ## a modulating covariate with ten distinct values admits at most nine
  ## basis functions per block, so the larger candidates are dropped
  ud <- u0
  ud[, 1L] <- round(ud[, 1L] * 9) / 9
  fit <- wafc_competitor("bsgl", x0, ud, y0, foldid = folds)
  expect_lte(fit[["extra"]][["df"]], 9L)
})

test_that("bsgl records the dimension it chose as a level of the WAFC", {
  ## step E2.5b: the pilot read the J of bsgl from 'extra' and found none
  skip_if_not(has("grpreg"))
  fit <- wafc_competitor("bsgl", x0, u0, y0, foldid = folds, df = c(4L, 16L))
  expect_true(fit[["extra"]][["df"]] %in% c(4L, 16L))
  expect_identical(fit[["extra"]][["J"]],
                   as.integer(log2(fit[["extra"]][["df"]])))
  ## df - 1 columns per block, the 2^J - 1 of a WAFC block at that J
  des <- wafc_design(x0, u0, J = fit[["extra"]][["J"]])
  expect_equal(length(des[["blocks"]][[1L]]), fit[["extra"]][["df"]] - 1L)
  ## a candidate that is not a power of two has no level
  fit <- wafc_competitor("bsgl", x0, u0, y0, foldid = folds, df = 6L)
  expect_identical(fit[["extra"]][["J"]], NA_integer_)
})

## ---------------------------------------------------------------------------
## The basis dimension of the spline chosen by a criterion (step E2.5h)
## ---------------------------------------------------------------------------

## A small replicate, so that the search by hand below costs seconds.
dk <- simulate_wafc(150L, p = 2L, q = 2L, scenario = "smooth", seed = 7L,
                    snr = 4)
xk <- dk[["x"]]
uk <- dk[["u"]]
yk <- dk[["y"]]
gridk <- c(5L, 10L, 20L)

## The model of wafc_fit_gam() at a common k, written out by hand.
gam_by_hand <- function(k, method, engine = "gam", u = uk) {
  dat <- data.frame(x1 = xk[, 1L], x2 = xk[, 2L], u1 = u[, 1L], u2 = u[, 2L],
                    y = yk)
  fo <- stats::as.formula(sprintf(paste(
    "y ~ 0 + x1 + x2 + s(u1, k = %d, by = x1) + s(u2, k = %d, by = x1) +",
    "s(u1, k = %d, by = x2) + s(u2, k = %d, by = x2)"), k[1L], k[2L], k[1L],
    k[2L]))
  if (engine == "gam") {
    mgcv::gam(fo, data = dat, method = method, select = TRUE)
  } else {
    mgcv::bam(fo, data = dat, method = "fREML", discrete = TRUE,
              select = TRUE)
  }
}

## The Gaussian restricted negative log-likelihood of y with the parametric
## terms as fixed effects, written with the n by n covariance
## V = I + Z S^{-1} Z' of the mixed model (the form of Harville, and of
## Kauermann and Opsomer, 2011, with REML in place of ML), at the smoothing
## parameters of a fit and at a given scale.
reml_dense <- function(fit, y, s2) {
  X <- stats::predict(fit, type = "lpmatrix")
  np <- fit[["nsdf"]]
  P <- ncol(X)
  S <- matrix(0, P, P)
  j <- 0L
  for (sm in fit[["smooth"]]) {
    ii <- sm[["first.para"]]:sm[["last.para"]]
    for (Sj in sm[["S"]]) {
      j <- j + 1L
      S[ii, ii] <- S[ii, ii] + fit[["sp"]][j] * Sj
    }
  }
  Xp <- X[, seq_len(np), drop = FALSE]
  Z <- X[, -seq_len(np), drop = FALSE]
  V <- diag(nrow(X)) + Z %*% solve(S[-seq_len(np), -seq_len(np)], t(Z))
  Vi <- solve(V)
  A <- t(Xp) %*% Vi %*% Xp
  r <- y - Xp %*% solve(A, t(Xp) %*% Vi %*% y)
  ld <- function(M) as.numeric(determinant(M, logarithm = TRUE)[["modulus"]])
  n <- nrow(X)
  0.5 * ((n - np) * log(2 * pi * s2) + ld(V) + ld(A) +
           sum(r * (Vi %*% r)) / s2)
}

test_that("the REML score is one likelihood of the same data at every k", {
  skip_slow()
  skip_if_not(has("mgcv"))
  for (k in gridk) {
    ## with the exact engine the score is the one mgcv optimizes
    fg <- gam_by_hand(c(k, k), "REML")
    sg <- wafc_gam_reml(fg, yk)
    expect_equal(as.numeric(sg), as.numeric(fg[["gcv.ubre"]]), tolerance = 1e-6)
    ## and with both engines it is the restricted likelihood of y with the
    ## same two level terms as fixed effects, so the scores of two values
    ## of k are values of one function of the data: comparable
    for (eng in c("gam", "bam")) {
      fk <- if (eng == "gam") fg else gam_by_hand(c(k, k), "REML", "bam")
      sk <- wafc_gam_reml(fk, yk)
      expect_equal(as.numeric(sk),
                   reml_dense(fk, yk, attr(sk, "sigma2")),
                   tolerance = 1e-8)
    }
  }
})

test_that("the choice of k in each criterion is the search redone by hand", {
  skip_slow()
  skip_if_not(has("mgcv"))
  ## REML with the exact engine: the score of each candidate is the one
  ## mgcv reports, so the search by hand reads it from mgcv alone
  fr <- wafc_competitor("gam", xk, uk, yk, k = gridk, k.select = "reml")
  hand <- vapply(gridk, function(k) {
    as.numeric(gam_by_hand(c(k, k), "REML")[["gcv.ubre"]])
  }, 0)
  tab <- fr[["extra"]][["k.table"]]
  expect_equal(tab[["k"]], gridk)
  expect_equal(tab[["score"]], hand, tolerance = 1e-6)
  expect_identical(fr[["extra"]][["k"]], rep(gridk[which.min(hand)], 2L))
  expect_identical(fr[["extra"]][["k.top"]], which.min(hand) == 3L)
  expect_identical(fr[["extra"]][["smooth.method"]], "REML")
  ## the fit kept is the one at the chosen k
  fh <- gam_by_hand(fr[["extra"]][["k"]], "REML")
  expect_equal(fr[["fitted"]], as.numeric(stats::fitted(fh)),
               tolerance = 1e-8)
  ## GCV: the score of a candidate is n RSS / (n - tau)^2 of its fit, with
  ## tau the total effective degrees of freedom
  fv <- wafc_competitor("gam", xk, uk, yk, k = gridk, k.select = "gcv")
  hv <- vapply(gridk, function(k) {
    f <- gam_by_hand(c(k, k), "GCV.Cp")
    nn <- length(yk)
    nn * sum(stats::residuals(f)^2) / (nn - sum(f[["edf"]]))^2
  }, 0)
  tv <- fv[["extra"]][["k.table"]]
  expect_equal(tv[["score"]], hv, tolerance = 1e-6)
  expect_identical(fv[["extra"]][["k"]], rep(gridk[which.min(hv)], 2L))
  expect_identical(fv[["extra"]][["smooth.method"]], "GCV.Cp")
  ## with the bam engine, REML searches by the score recomputed here and
  ## keeps the one bam reports beside it
  fb <- wafc_competitor("gam", xk, uk, yk, k = gridk, k.select = "reml",
                        engine = "bam")
  hb <- vapply(gridk, function(k) {
    f <- gam_by_hand(c(k, k), "REML", "bam")
    c(as.numeric(wafc_gam_reml(f, yk)), as.numeric(f[["gcv.ubre"]]))
  }, c(0, 0))
  tb <- fb[["extra"]][["k.table"]]
  expect_equal(tb[["score"]], hb[1L, ], tolerance = 1e-8)
  expect_equal(tb[["score.engine"]], hb[2L, ], tolerance = 1e-8)
  expect_identical(fb[["extra"]][["k"]], rep(gridk[which.min(hb[1L, ])], 2L))
})

test_that("the grid of k is the one of D41, truncated at the distinct values", {
  skip_slow()
  ## the default grid of the search: Ruppert (2002, section 3) without his
  ## 120, one grid for the three criteria (decision D41)
  expect_identical(wafc_k_grid, c(5L, 10L, 20L, 40L, 80L))
  ## a modulator with 25 distinct values admits at most 24 basis functions:
  ## 40 and 80 become 24 on it, and the candidates the truncation makes
  ## equal are fitted once
  ud <- uk
  ud[, 2L] <- round(ud[, 2L] * 24) / 24
  cd <- wafc_k_candidates(ud, wafc_k_grid)
  expect_identical(attr(cd, "grid"), wafc_k_grid)
  expect_identical(cd[[4L]], c(40L, 24L))
  expect_identical(cd[[5L]], c(80L, 24L))
  ud[, 1L] <- round(ud[, 1L] * 24) / 24
  cd <- wafc_k_candidates(ud, wafc_k_grid)
  expect_identical(attr(cd, "grid"), c(5L, 10L, 20L, 40L))
  expect_identical(cd, structure(list(c(5L, 5L), c(10L, 10L), c(20L, 20L),
                                      c(24L, 24L)),
                                 grid = c(5L, 10L, 20L, 40L)))
  ## and the search fits exactly those, with the top read on them
  skip_if_not(has("mgcv"))
  ud <- uk
  ud[, 1L] <- round(ud[, 1L] * 11) / 11
  ud[, 2L] <- round(ud[, 2L] * 11) / 11
  fit <- wafc_competitor("gam", xk, ud, yk, k = c(5L, 20L, 40L),
                         k.select = "reml")
  tab <- fit[["extra"]][["k.table"]]
  expect_identical(tab[["k"]], c(5L, 20L))
  expect_identical(tab[["k.used"]], c("5,5", "11,11"))
  hand <- vapply(list(c(5L, 5L), c(11L, 11L)), function(k) {
    as.numeric(gam_by_hand(k, "REML", u = ud)[["gcv.ubre"]])
  }, 0)
  expect_equal(tab[["score"]], hand, tolerance = 1e-6)
  expect_identical(fit[["extra"]][["k.top"]], which.min(hand) == 2L)
  ## without a search nothing changes: k = 10, no table
  f0 <- wafc_competitor("gam", xk, uk, yk)
  expect_identical(f0[["extra"]][["k"]], c(10L, 10L))
  expect_null(f0[["extra"]][["k.table"]])
  expect_identical(f0[["extra"]][["k.select"]], "none")
  expect_error(wafc_competitor("gam", xk, uk, yk, k = c(2L, 10L),
                               k.select = "gcv"), "at least 3")
})

test_that("the choice of k by cross-validation is the search redone by hand", {
  skip_slow()
  skip_if_not(has("mgcv"))
  ## the folds of the replicate, REML and bam inside each fold, and the
  ## loss of cv.wafc(): the mean over folds of the mean squared error on the
  ## observations the fold leaves out
  fk <- rep_len(1:5, length(yk))
  fit <- wafc_competitor("gam", xk, uk, yk, k = gridk, k.select = "cv",
                         engine = "bam", foldid = fk)
  dat <- data.frame(x1 = xk[, 1L], x2 = xk[, 2L], u1 = uk[, 1L],
                    u2 = uk[, 2L], y = yk)
  hand <- t(vapply(gridk, function(k) {
    fm <- vapply(1:5, function(f) {
      out <- fk == f
      fo <- stats::as.formula(sprintf(paste(
        "y ~ 0 + x1 + x2 + s(u1, k = %d, by = x1) + s(u2, k = %d, by = x1) +",
        "s(u1, k = %d, by = x2) + s(u2, k = %d, by = x2)"), k, k, k, k))
      ff <- mgcv::bam(fo, data = dat[!out, ], method = "fREML",
                      discrete = TRUE, select = TRUE)
      mean((yk[out] - as.numeric(stats::predict(ff, newdata = dat[out, ])))^2)
    }, 0)
    c(mean(fm), stats::sd(fm) / sqrt(5))
  }, c(0, 0)))
  tab <- fit[["extra"]][["k.table"]]
  expect_equal(tab[["k"]], gridk)
  expect_equal(tab[["score"]], hand[, 1L], tolerance = 1e-10)
  expect_equal(tab[["cvsd"]], hand[, 2L], tolerance = 1e-10)
  expect_true(all(is.na(tab[["score.engine"]])))
  kh <- gridk[which.min(hand[, 1L])]
  expect_identical(fit[["extra"]][["k"]], c(kh, kh))
  expect_identical(fit[["extra"]][["k.top"]], which.min(hand[, 1L]) == 3L)
  expect_identical(fit[["extra"]][["foldid"]], fk)
  expect_identical(fit[["extra"]][["smooth.method"]], "REML")
  ## the fit returned is the one on the whole sample at the k chosen
  fh <- gam_by_hand(c(kh, kh), "REML", "bam")
  expect_equal(fit[["fitted"]], as.numeric(stats::fitted(fh)),
               tolerance = 1e-10)
  ## with the folds fixed the search reproduces exactly
  again <- wafc_competitor("gam", xk, uk, yk, k = gridk, k.select = "cv",
                           engine = "bam", foldid = fk)
  expect_identical(again[["extra"]][["k.table"]][["score"]], tab[["score"]])
  expect_identical(again[["fitted"]], fit[["fitted"]])
  ## the folds are validated as the ones of cv.wafc()
  expect_error(wafc_competitor("gam", xk, uk, yk, k = gridk, k.select = "cv",
                               foldid = rep_len(1:2, length(yk))),
               "at least 3 folds")
})

test_that("a candidate of k that fails leaves the search (step E4.3b)", {
  skip_if_not(has("mgcv"))
  ## The failure the probe of D78 met, in the cell "mixed" at n = 1000 with
  ## k up to 240, does not happen at a size a test can afford, so it is
  ## planted: the score of the REML search stops, with the message of
  ## wafc_gam_reml(), on the fit at the top of the grid, or the fit itself
  ## stops there. Either way the search is the one on the other candidates,
  ## and the table keeps the candidate left out with its reason.
  rank_msg <- "the penalty of a smooth has rank below the one mgcv declares."
  top_of <- function(fit) {
    any(vapply(fit[["smooth"]], function(sm) sm[["bs.dim"]], 0) == 20)
  }
  real_reml <- wafc_gam_reml
  plant <- function(fails) {
    assign("wafc_gam_reml", function(fit, y) {
      if (fails(fit)) stop(rank_msg, call. = FALSE)
      real_reml(fit, y)
    }, envir = globalenv())
  }
  withr::defer(assign("wafc_gam_reml", real_reml, envir = globalenv()))
  ## the search on the two candidates that do not fail
  ref <- wafc_competitor("gam", xk, uk, yk, k = c(5L, 10L),
                         k.select = "reml", engine = "bam")
  rt <- ref[["extra"]][["k.table"]]
  expect_true(all(is.na(rt[["error"]])))
  plant(top_of)
  fs <- wafc_competitor("gam", xk, uk, yk, k = gridk, k.select = "reml",
                        engine = "bam")
  ts <- fs[["extra"]][["k.table"]]
  expect_identical(ts[["k"]], gridk)
  expect_identical(ts[["error"]], c(NA, NA, rank_msg))
  expect_identical(ts[["score"]][1:2], rt[["score"]])
  expect_true(is.na(ts[["score"]][3L]))
  ## the fit was there, so its size is recorded
  expect_equal(ts[["ncoef"]][3L], 4 * 20)
  expect_identical(fs[["extra"]][["k"]], ref[["extra"]][["k"]])
  expect_identical(fs[["fitted"]], ref[["fitted"]])
  ## the top is read on the candidates fitted and scored
  expect_identical(fs[["extra"]][["k.top"]], ref[["extra"]][["k.top"]])
  ## a fit that fails is left out the same way, with nothing of it in the
  ## table but the reason (inside the folds of the cross-validation too, in
  ## the slow test below)
  assign("wafc_gam_reml", real_reml, envir = globalenv())
  real_bam <- mgcv::bam
  local_mocked_bindings(bam = function(formula, ...) {
    if (grepl("k = 20", paste(deparse(formula), collapse = ""))) {
      stop("planted failure of the fit", call. = FALSE)
    }
    real_bam(formula, ...)
  }, .package = "mgcv")
  ff <- wafc_competitor("gam", xk, uk, yk, k = gridk, k.select = "reml",
                        engine = "bam")
  tf <- ff[["extra"]][["k.table"]]
  expect_identical(tf[["error"]], c(NA, NA, "planted failure of the fit"))
  expect_true(all(is.na(unlist(tf[3L, c("score", "edf", "ncoef")]))))
  expect_identical(tf[["score"]][1:2], rt[["score"]])
  expect_identical(ff[["fitted"]], ref[["fitted"]])
  ## when every candidate fails, the search stops with the error of the
  ## first, as it did before the step
  plant(function(fit) TRUE)
  expect_error(wafc_competitor("gam", xk, uk, yk, k = gridk,
                               k.select = "reml", engine = "bam"),
               rank_msg, fixed = TRUE)
  expect_error(wafc_competitor("gam", xk, uk, yk, k = 20L, engine = "bam"),
               "planted failure of the fit", fixed = TRUE)
})

test_that("a fold that fails or does not converge leaves its k out (E4.3b)", {
  skip_slow()
  skip_if_not(has("mgcv"))
  fk <- rep_len(1:5, length(yk))
  ## a fit that fails inside a fold of the cross-validation
  real_bam <- mgcv::bam
  local_mocked_bindings(bam = function(formula, ...) {
    if (grepl("k = 20", paste(deparse(formula), collapse = ""))) {
      stop("planted failure of the fit", call. = FALSE)
    }
    real_bam(formula, ...)
  }, .package = "mgcv")
  rc <- wafc_competitor("gam", xk, uk, yk, k = c(5L, 10L), k.select = "cv",
                        engine = "bam", foldid = fk)
  fc <- wafc_competitor("gam", xk, uk, yk, k = gridk, k.select = "cv",
                        engine = "bam", foldid = fk)
  tc <- fc[["extra"]][["k.table"]]
  expect_identical(tc[["error"]], c(NA, NA, "planted failure of the fit"))
  expect_true(all(is.na(unlist(tc[3L, c("score", "cvsd", "edf", "ncoef")]))))
  expect_identical(tc[["score"]][1:2], rc[["extra"]][["k.table"]][["score"]])
  expect_identical(fc[["extra"]][["k"]], rc[["extra"]][["k"]])
  expect_identical(fc[["fitted"]], rc[["fitted"]])
  ## a fold that does not converge within the limit (decision D80): on
  ## this replicate the folds take at most 35, 30 and 30 outer iterations at
  ## k = 5, 10 and 20, so 31 stops k = 5 only (folds at 33 and 35)
  local_mocked_bindings(bam = real_bam, .package = "mgcv")
  msg <- "not converged in 31 iterations of the smoothing parameters"
  rc <- wafc_competitor("gam", xk, uk, yk, k = c(10L, 20L), k.select = "cv",
                        engine = "bam", foldid = fk)
  fc <- wafc_competitor("gam", xk, uk, yk, k = gridk, k.select = "cv",
                        engine = "bam", foldid = fk, k.maxit = 31L)
  tc <- fc[["extra"]][["k.table"]]
  expect_identical(tc[["error"]], c(msg, NA, NA))
  expect_identical(tc[["iter"]], c(NA, 30, 30))
  expect_identical(tc[["score"]][2:3], rc[["extra"]][["k.table"]][["score"]])
  expect_identical(fc[["fitted"]], rc[["fitted"]])
})

test_that("the size rule of D79 leaves a k out before it is fitted", {
  skip_if_not(has("mgcv"))
  ## the count: one smooth of k_m columns per block, p sum(k_m), and the
  ## bound 2n is kept, only what passes it is left out
  expect_identical(wafc_k_size(list(c(75L, 75L), c(76L, 75L)), 2L),
                   c(300, 302))
  expect_identical(wafc_k_size_ratio, 2)
  ## and only above the common grid of D41 (decision D80): 80 stays even
  ## past the bound, 90 does not
  expect_identical(wafc_k_left_out(c(80L, 90L, 90L, NA), c(320, 360, 300, 999),
                                   150L), c(FALSE, TRUE, FALSE, FALSE))
  ## n = 150 with p = q = 2: 4 smooths, so k = 90 has 360 > 300
  ## coefficients and is never fitted; the search is the one without it
  real_bam <- mgcv::bam
  fitted_k <- integer(0)
  local_mocked_bindings(bam = function(formula, ...) {
    fo <- paste(deparse(formula), collapse = "")
    fitted_k <<- c(fitted_k, as.integer(sub(".*k = ([0-9]+).*", "\\1", fo)))
    real_bam(formula, ...)
  }, .package = "mgcv")
  ref <- wafc_competitor("gam", xk, uk, yk, k = c(5L, 10L),
                         k.select = "reml", engine = "bam")
  fitted_k <- integer(0)
  fs <- wafc_competitor("gam", xk, uk, yk, k = c(5L, 10L, 90L),
                        k.select = "reml", engine = "bam")
  expect_false(90L %in% fitted_k)
  ts <- fs[["extra"]][["k.table"]]
  expect_identical(ts[["k"]], c(5L, 10L, 90L))
  expect_identical(ts[["error"]][3L], paste("not fitted: 360 coefficients in",
                                            "the smooths, more than 2 n = 300"))
  expect_true(all(is.na(unlist(ts[3L, c("score", "edf", "ncoef")]))))
  expect_identical(ts[["score"]][1:2], ref[["extra"]][["k.table"]][["score"]])
  expect_identical(fs[["extra"]][["k"]], ref[["extra"]][["k"]])
  expect_identical(fs[["fitted"]], ref[["fitted"]])
  expect_identical(fs[["extra"]][["k.top"]], ref[["extra"]][["k.top"]])
  ## the same in the cross-validation, where nothing of 80 is fitted in
  ## any fold either
  fitted_k <- integer(0)
  fc <- wafc_competitor("gam", xk, uk, yk, k = c(5L, 90L), k.select = "cv",
                        engine = "bam", foldid = rep_len(1:5, length(yk)))
  expect_false(90L %in% fitted_k)
  expect_identical(fc[["extra"]][["k"]], c(5L, 5L))
  ## a search with every candidate past the bound stops on the rule
  expect_error(wafc_competitor("gam", xk, uk, yk, k = c(90L, 160L),
                               k.select = "reml", engine = "bam"),
               "every candidate of k has more than 2 n = 300", fixed = TRUE)
})

test_that("a candidate of k that does not converge leaves the search (D80)", {
  skip_if_not(has("mgcv"))
  ## On this replicate bam converges in 35, 24 and 16 outer iterations at
  ## k = 5, 10 and 20, and in at most 35, 30 and 30 in the folds. A limit
  ## of 20 stops the first two; the search is then the one on k = 20 alone,
  ## and that fit, within the limit, is the fit without it.
  f20 <- wafc_competitor("gam", xk, uk, yk, k = 20L, k.select = "reml",
                         engine = "bam", k.maxit = 1000L)
  fs <- wafc_competitor("gam", xk, uk, yk, k = gridk, k.select = "reml",
                        engine = "bam", k.maxit = 20L)
  ts <- fs[["extra"]][["k.table"]]
  msg <- "not converged in 20 iterations of the smoothing parameters"
  expect_identical(ts[["error"]], c(msg, msg, NA))
  expect_identical(ts[["iter"]], c(NA, NA, 16))
  expect_identical(ts[["score"]][3L], f20[["extra"]][["k.table"]][["score"]])
  expect_identical(fs[["extra"]][["k"]], c(20L, 20L))
  expect_identical(fs[["fitted"]], f20[["fitted"]])
  ## with a limit no fit reaches, the search is the one of mgcv's default
  ## of 200 iterations, and the default limit (80) is one of them here
  fa <- wafc_competitor("gam", xk, uk, yk, k = gridk, k.select = "reml",
                        engine = "bam")
  fb <- wafc_competitor("gam", xk, uk, yk, k = gridk, k.select = "reml",
                        engine = "bam", k.maxit = 200L)
  expect_identical(fa[["extra"]][["k.table"]][["iter"]], c(35, 24, 16))
  expect_identical(fa[["fitted"]], fb[["fitted"]])
  expect_identical(fa[["extra"]][["k.table"]][["score"]],
                   fb[["extra"]][["k.table"]][["score"]])
  ## (in the folds of the cross-validation too, in the slow test above)
  ## no candidate converges: the error of the first; the limit is checked
  expect_error(wafc_competitor("gam", xk, uk, yk, k = gridk,
                               k.select = "reml", engine = "bam", k.maxit = 2L),
               "not converged in 2 iterations", fixed = TRUE)
  expect_error(wafc_competitor("gam", xk, uk, yk, k = gridk,
                               k.select = "reml", engine = "bam", k.maxit = 0),
               "'k.maxit' must be one positive number", fixed = TRUE)
})

## ---------------------------------------------------------------------------
## What grpreg returned of the path of the block LASSO (step E2.5j)
## ---------------------------------------------------------------------------

test_that("the block LASSO records the path of grpreg at every J", {
  skip_slow()
  skip_if_not(has("grpreg"))
  k <- wafc_competitor("klopp", x0, u0, y0, J = 2:4, foldid = folds,
                       penalize.levels = FALSE, balanced = TRUE)
  cc <- k[["extra"]][["conv"]]
  expect_identical(cc[["J"]], 2:4)
  expect_identical(cc[["chosen"]], 2:4 == k[["extra"]][["J"]])
  expect_identical(cc[["nlambda"]], rep(100L, 3L))
  expect_identical(cc[["max.iter"]], rep(10000L, 3L))
  ## the row of the J kept is the cv.grpreg object the fit carries
  ch <- cc[cc[["chosen"]], ]
  cvo <- k[["fit"]]
  expect_identical(ch[["nreturned"]], length(cvo[["fit"]][["lambda"]]))
  expect_identical(ch[["ncv"]], length(cvo[["lambda"]]))
  expect_identical(ch[["iter.total"]], as.integer(sum(cvo[["fit"]][["iter"]])))
  expect_identical(ch[["lambda.min"]], k[["extra"]][["lambda"]])
  expect_identical(ch[["lambda.cv.last"]], min(cvo[["lambda"]]))
  expect_false(any(cc[["budget"]] | cc[["cv.cut"]] | cc[["cut.at.min"]]))
  ## the free coarse levels at a coarse J are least squares, with no path
  kf <- wafc_competitor("klopp", x0, u0, y0, J = 2:3, foldid = folds,
                        penalize.levels = FALSE, balanced = TRUE,
                        free.coarse = TRUE)
  cf <- kf[["extra"]][["conv"]]
  expect_identical(nrow(cf), 2L)
  expect_true(is.na(cf[["nreturned"]][1L]))
})

test_that("the flags of the table read a cut as grpreg makes it", {
  ## a budget spent and a cross-validated path shorter than the fit's, built
  ## by hand on the rows wafc_kp_conv() returns
  r <- data.frame(J = 3:4, nlambda = 100L, nreturned = c(100L, 60L),
                  ncv = c(90L, 60L), iter.total = c(500L, 10000L),
                  max.iter = 10000L, n.maxiter = 0L, lambda.maxiter = NA_real_,
                  lambda.last = c(1e-4, 0.01), lambda.cv.last = c(2e-4, 0.01),
                  lambda.min = c(2e-4, 0.05))
  d <- wafc_kp_conv_table(list(r[1L, ], r[2L, ]), J.min = 4L)
  expect_identical(d[["chosen"]], c(FALSE, TRUE))
  expect_identical(d[["budget"]], c(FALSE, TRUE))
  expect_identical(d[["cv.cut"]], c(TRUE, FALSE))
  expect_identical(d[["cut.at.min"]], c(TRUE, FALSE))
})

## ---------------------------------------------------------------------------
## The two public functions that had no test of their own (step E3.1)
## ---------------------------------------------------------------------------

test_that("wafc_k_matched is 2^J truncated at the distinct values", {
  expect_identical(wafc_k_matched(u0, J = 4L), c(16L, 16L))
  expect_identical(wafc_k_matched(u0, J = c(3L, 5L)), c(8L, 32L))
  ## 250 distinct values admit at most 249 basis functions
  expect_identical(wafc_k_matched(u0, J = 8L), c(249L, 249L))
  ud <- cbind(rep_len(1:5, n), u0[, 2L])
  expect_identical(wafc_k_matched(ud, J = 4L), c(4L, 16L))
  expect_identical(wafc_k_matched(ud, J = 1L, kmin = 3L), c(3L, 3L))
  expect_error(wafc_k_matched(u0, J = 0L), "at least 1")
})

test_that("a competitor prints its method, its levels and its blocks", {
  fit <- wafc_competitor("linear", x0, u0, y0)
  out <- utils::capture.output(print(fit))
  expect_match(out[1L], "WAFC competitor \"linear\": n = 250, p = 3, q = 2")
  expect_true(any(grepl("levels c:", out)))
  if (has("grpreg")) {
    kp <- wafc_competitor("klopp", x0, u0, y0, J = 3L, foldid = folds,
                          penalize.levels = FALSE, balanced = TRUE)
    expect_output(print(kp), "blocks kept: [0-6] of 6")
  }
})
