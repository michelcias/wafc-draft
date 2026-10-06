## R/data.R -- the data of one replicate.
##
## A replicate is a training sample of size n, a test sample of size
## `test_size` from the same law, the grid on which the components are
## compared, and the partition in folds every method of the replicate uses.
## All four depend on the cell, n and the replicate only, never on the
## method, so the comparison between methods is paired.

#' Draw a sample of a cell
#'
#' The entries of the cell other than the reserved ones are passed to
#' simulate_wafc(). A cell with an explicit `structure` (a p by q matrix of
#' component names, "" for a zero block) is drawn as simulate_wafc() draws
#' a scenario: the covariates and the standard normal errors come from
#' simulate_wafc() in the null scenario with unit error standard deviation,
#' and the regression function is built from the components of
#' wafc_component(), each multiplied by its amplitude, before the error
#' scale is set from `snr` or `sigma`.
#'
#' @return An object of class "wafc_dgp".
draw_cell <- function(cell, n, seed) {
  a <- cell[setdiff(names(cell), c(cell_reserved, "name"))]
  a[["n"]] <- n
  a[["seed"]] <- seed
  if (is.null(cell[["structure"]])) return(do.call(simulate_wafc, a))
  snr <- a[["snr"]]
  sigma <- a[["sigma"]]
  a[["snr"]] <- NULL
  a[["scenario"]] <- "null"
  a[["sigma"]] <- 1
  d <- do.call(simulate_wafc, a)
  z <- d[["y"]] - d[["f"]]
  s <- cell[["structure"]]
  p <- d[["p"]]
  q <- d[["q"]]
  amp <- matrix(rep_len(if (is.null(a[["amplitude"]])) 1 else
    unlist(a[["amplitude"]]), p * q), p, q)
  g <- vector("list", p * q)
  dim(g) <- c(p, q)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      if (!nzchar(s[l, m])) next
      g[[l, m]] <- local({
        gg <- wafc_component(s[l, m])
        aa <- amp[l, m]
        function(v) aa * gg(v)
      })
    }
  }
  d[["g"]] <- g
  d[["structure"]] <- s
  d[["scenario"]] <- cell[["scenario"]]
  d[["sprime"]] <- NA_real_
  d[["beta"]] <- wafc_beta(d, d[["u"]])
  d[["f"]] <- as.numeric(rowSums(d[["x"]] * d[["beta"]]))
  d[["sigma"]] <- if (is.null(sigma)) stats::sd(d[["f"]]) / snr else sigma
  d[["y"]] <- d[["f"]] + d[["sigma"]] * z
  d
}

#' The test sample of a replicate
#'
#' Drawn from the law of the cell with its own seed; the errors of the test
#' responses are drawn after it at the error scale of the training sample,
#' and not at the scale the signal to noise ratio would set on the test
#' sample itself: the two differ slightly, and the difference would end up
#' inside the out-of-sample error.
draw_test <- function(cell, train, size, seed) {
  test <- draw_cell(cell, size, seed)
  test[["y"]] <- test[["f"]] + stats::rnorm(size, sd = train[["sigma"]])
  test[["sigma"]] <- train[["sigma"]]
  test
}

#' The grid on which the components are compared
#'
#' The observed range of each modulating covariate in the training sample,
#' which is where every method is identified.
grid_of <- function(train, size) {
  q <- train[["q"]]
  g <- matrix(0, size, q)
  for (m in seq_len(q)) {
    rg <- range(train[["u"]][, m])
    g[, m] <- seq(rg[1L], rg[2L], length.out = size)
  }
  g
}

#' Everything a method needs from one replicate
#'
#' Draws the training sample, the test sample and the folds from the
#' seeds of the replicate, in that order, and leaves the random stream
#' where the folds leave it: every method starts from there.
#'
#' @param seeds The list of unit_seeds().
#' @param tuning The tuning of the sample size (method_tuning()), for the
#'   number of folds.
#' @param table The wavelet table of the sample size.
replicate_context <- function(cfg, cell, n, rep, seeds, tuning, table) {
  train <- draw_cell(cell, n, seeds[["data"]])
  test <- draw_test(cell, train, as.integer(cfg[["test_size"]]),
                    seeds[["test"]])
  set.seed(seeds[["folds"]])
  foldid <- sample(rep_len(seq_len(tuning[["nfolds"]]), n))
  list(cell = cell, n = n, rep = rep, train = train, test = test,
       grid = grid_of(train, as.integer(cfg[["grid_size"]])),
       active = nzchar(train[["structure"]]), foldid = foldid,
       table = table, rng = get(".Random.seed", envir = globalenv()))
}

## Puts the random stream back where the folds left it.
restore_stream <- function(ctx) {
  assign(".Random.seed", ctx[["rng"]], envir = globalenv())
}
