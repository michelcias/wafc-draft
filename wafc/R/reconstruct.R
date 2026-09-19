## wafc/R/reconstruct.R -- reading a WAFC fit as functions (step E2.2):
## the levels c_l, the additive components g_{lm} on a grid, the functional
## coefficients beta_l(u) and the structure selected block by block.
##
## Everything here goes through the 'spec' route of wafc_design(): a design
## built with spec = the fitted design has the same basis, the same
## rescaling and the same columns in the same order, so a block of
## coefficients can be multiplied by the block of a design evaluated
## anywhere. Evaluating it with X_l = 1 is what turns the block (l, m) into
## the function g_{lm}.
##
## One caveat that E2.1 recorded and that the theory does not see (open
## question 10 of docs/ESTADO.md): with rescale = TRUE the basis is
## evaluated on [eps, 1 - eps] of the rescaled covariate, so the estimated
## component integrates to zero over the observed range of U_m and not over
## [0,1]. It is a shift of level, not of shape: the difference belongs to
## c_l, and beta_l(u) = c_l + sum_m g_{lm}(u_m) does not see it.

#' Components of a WAFC fit as functions
#'
#' Reconstructs, at one penalty level, the levels \eqn{\hat c_\ell} and the
#' additive components \eqn{\hat g_{\ell m}} of a fit on a grid of values of
#' each modulating covariate.
#'
#' @param object An object of class \code{"wafc"}.
#' @param s The penalty level, a single value on the scale of the objective
#'   of \code{\link{wafc}}. \code{NULL} takes the smallest value of the
#'   path, which is the least penalized fit.
#' @param grid The values at which the components are evaluated, on the
#'   scale of the modulating covariates as they were given to
#'   \code{\link{wafc}}: a matrix with \eqn{q} columns, a vector used for
#'   every covariate, or \code{NULL} (the default) for \code{n_grid} equally
#'   spaced points over the range each covariate had in the training
#'   sample.
#' @param n_grid Number of grid points when \code{grid} is \code{NULL}.
#'
#' @return An object of class \code{"wafc_functions"}: a list with the
#'   penalty level \code{s}, the \code{grid} (\code{n_grid} by \eqn{q}), the
#'   levels \code{cc}, the \eqn{p} by \eqn{q} list \code{g} of evaluated
#'   components, the \eqn{p} by \eqn{q} matrices \code{nonzero} and
#'   \code{norm} of \code{\link{wafc_blocks}}, and the names of the
#'   covariates.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' fit <- wafc(d$x, d$u, d$y, J = 4)
#' fn <- wafc_functions(fit, s = fit$lambda[40], n_grid = 128)
#' fn$nonzero
#' plot(fn$grid[, 1], fn$g[[1, 1]], type = "l")
#'
#' @export
wafc_functions <- function(object, s = NULL, grid = NULL, n_grid = 512L) {
  if (!inherits(object, "wafc")) {
    stop("'object' must be an object returned by wafc().", call. = FALSE)
  }
  design <- object[["design"]]
  p <- design[["p"]]
  q <- design[["q"]]
  s <- wafc_single_s(object, s)
  grid <- wafc_grid(design, grid, n_grid)
  cf <- coef.wafc(object, s = s)[, 1L]
  b <- cf[-1L]

  ## One design on the grid gives every block at once: with X_l = 1 the
  ## block (l, m) is the basis of U_m evaluated at grid[, m].
  d_grid <- wafc_design(matrix(1, nrow(grid), p), grid, spec = design)
  g <- vector("list", p * q)
  dim(g) <- c(p, q)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      idx <- design[["blocks"]][[wafc_block_name(design, l, m)]]
      g[[l, m]] <- as.numeric(d_grid[["Z"]][, idx, drop = FALSE] %*% b[idx])
    }
  }
  dimnames(g) <- list(design[["xnames"]], design[["unames"]])
  blk <- wafc_blocks(object, s = s)
  out <- list(s = s, grid = grid, cc = cf[1L + design[["unpenalized"]]],
              g = g, nonzero = blk[["nonzero"]], norm = blk[["norm"]],
              xnames = design[["xnames"]], unames = design[["unames"]],
              intercept = cf[[1L]])
  names(out[["cc"]]) <- design[["xnames"]]
  class(out) <- "wafc_functions"
  out
}

#' Structure selected by a WAFC fit, block by block
#'
#' The number of non-zero wavelet coefficients and the Euclidean norm of
#' \eqn{\hat\theta_{\ell m}} in each block \eqn{(\ell, m)}. A block that is
#' entirely zero says that the fit left \eqn{\hat\beta_\ell} free of the
#' modulating covariate \eqn{U_m}, which is the selection of structure the
#' sparse group LASSO variant is meant to do.
#'
#' @param object An object of class \code{"wafc"}.
#' @param s Penalty levels, on the scale of the objective of
#'   \code{\link{wafc}}. \code{NULL} takes the smallest value of the path.
#'
#' @return With a single \code{s}, a list of two \eqn{p} by \eqn{q}
#'   matrices, \code{nonzero} and \code{norm}. With several, the two entries
#'   are arrays of dimension \eqn{p \times q \times} \code{length(s)}.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' fit <- wafc(d$x, d$u, d$y, J = 4, penalty = "sglasso")
#' wafc_blocks(fit, s = fit$lambda[30])$nonzero
#'
#' @export
wafc_blocks <- function(object, s = NULL) {
  if (!inherits(object, "wafc")) {
    stop("'object' must be an object returned by wafc().", call. = FALSE)
  }
  design <- object[["design"]]
  p <- design[["p"]]
  q <- design[["q"]]
  if (is.null(s)) s <- min(object[["lambda"]])
  cf <- wafc_raw_coef(object, s = s)[-1L, , drop = FALSE]
  dn <- list(design[["xnames"]], design[["unames"]])
  nz <- array(0L, c(p, q, ncol(cf)))
  nm <- array(0, c(p, q, ncol(cf)))
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      idx <- design[["blocks"]][[wafc_block_name(design, l, m)]]
      nz[l, m, ] <- as.integer(colSums(cf[idx, , drop = FALSE] != 0))
      nm[l, m, ] <- sqrt(colSums(cf[idx, , drop = FALSE]^2))
    }
  }
  if (ncol(cf) == 1L) {
    return(list(nonzero = matrix(nz[, , 1L], p, q, dimnames = dn),
                norm = matrix(nm[, , 1L], p, q, dimnames = dn)))
  }
  dimnames(nz) <- dimnames(nm) <- c(dn, list(colnames(cf)))
  list(nonzero = nz, norm = nm)
}

#' @export
print.wafc_functions <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 ...) {
  cat(sprintf("WAFC components at lambda = %s\n",
              format(x[["s"]], digits = digits)))
  cat("  levels c:", paste(sprintf("%s = %s", names(x[["cc"]]),
                                   format(x[["cc"]], digits = digits)),
                           collapse = ", "), "\n")
  cat("  non-zero wavelet coefficients per block (l, m):\n")
  print(x[["nonzero"]])
  invisible(x)
}

## ---------------------------------------------------------------------------
## Internals
## ---------------------------------------------------------------------------

## The functional coefficients beta_l(u) = c_l + sum_m g_{lm}(u_m) at one
## penalty level, evaluated at a matrix of modulating covariates (the
## training one when 'u' is NULL). Used by predict(type = "beta").
wafc_beta_hat <- function(object, u = NULL, s = NULL) {
  design <- object[["design"]]
  p <- design[["p"]]
  q <- design[["q"]]
  s <- wafc_single_s(object, s)
  if (is.null(u)) {
    stop("type = \"beta\" needs 'newu', the modulating covariates at which ",
         "the functional coefficients are wanted.", call. = FALSE)
  }
  u <- wafc_as_matrix(u, "newu")
  if (ncol(u) != q) {
    stop("'newu' must have ", q, " column(s).", call. = FALSE)
  }
  cf <- coef.wafc(object, s = s)[, 1L]
  b <- cf[-1L]
  d_u <- wafc_design(matrix(1, nrow(u), p), u, spec = design)
  beta <- matrix(rep(cf[1L + design[["unpenalized"]]], each = nrow(u)),
                 nrow(u), p)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      idx <- design[["blocks"]][[wafc_block_name(design, l, m)]]
      beta[, l] <- beta[, l] +
        as.numeric(d_u[["Z"]][, idx, drop = FALSE] %*% b[idx])
    }
  }
  colnames(beta) <- design[["xnames"]]
  beta
}

wafc_block_name <- function(design, l, m) {
  paste0(design[["xnames"]][l], ":", design[["unames"]][m])
}

wafc_single_s <- function(object, s) {
  if (is.null(s)) return(min(object[["lambda"]]))
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s < 0) {
    stop("'s' must be a single non-negative value.", call. = FALSE)
  }
  as.numeric(s)
}

## Grid of values of the modulating covariates on their original scale. The
## default is the range of the training sample, which the rescaling of the
## design stores implicitly: a rescaled value lies in [eps, 1 - eps], so the
## observed range was [location + eps * scale, location + (1 - eps) * scale]
## (and [0, 1] when rescale = FALSE).
wafc_grid <- function(design, grid, n_grid) {
  q <- design[["q"]]
  if (is.null(grid)) {
    if (length(n_grid) != 1L || !is.finite(n_grid) || n_grid < 2) {
      stop("'n_grid' must be a single integer of at least 2.", call. = FALSE)
    }
    n_grid <- as.integer(n_grid)
    grid <- matrix(0, n_grid, q)
    for (m in seq_len(q)) {
      lo <- design[["location"]][m] + design[["eps"]][m] * design[["scale"]][m]
      hi <- design[["location"]][m] +
        (1 - design[["eps"]][m]) * design[["scale"]][m]
      grid[, m] <- seq(lo, hi, length.out = n_grid)
    }
  } else {
    grid <- wafc_as_matrix(grid, "grid")
    if (ncol(grid) == 1L && q > 1L) {
      grid <- matrix(grid[, 1L], nrow(grid), q)
    }
    if (ncol(grid) != q) {
      stop("'grid' must have 1 or ", q, " column(s).", call. = FALSE)
    }
  }
  colnames(grid) <- design[["unames"]]
  grid
}
