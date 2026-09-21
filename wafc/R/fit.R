## wafc/R/fit.R -- the WAFC estimator (step E2.2): the LASSO of decision D3
## and the sparse group LASSO variant, both on the design of wafc_design().
##
## Two facts about the engines are built into this file, because both cost
## half a day to find (docs/ESTADO.md, E2.1 and E1.5):
##
##  (i)  glmnet drops columns of zero variance from the fit even with
##       standardize = FALSE, so the constant covariate X_1 = 1 gets
##       coefficient zero and its level c_1 is carried by the intercept of
##       the engine. The fit therefore uses intercept = TRUE whenever the
##       design has a constant column, and the level reported for that
##       covariate is the coefficient plus the intercept (divided by the
##       value of the column, which need not be one).
##  (ii) glmnet rescales penalty.factor so that it sums to nvars. With the
##       penalty factors of wafc_design() (zero on the p levels, one on the
##       wavelets) the penalty actually applied to a wavelet coefficient is
##       therefore lambda * nvars / npen and not lambda. The lambda of a
##       "wafc" object is always the effective one, the constant that
##       multiplies the penalty as it is written in the theory (E1.5), and
##       wafc_kkt() verifies it on the stationarity conditions instead of
##       trusting the argument. sparsegl does not rescale its penalty
##       factors, so for the group variant the two scales coincide.

#' Fit the WAFC model
#'
#' Penalized least squares on the design of \code{\link{wafc_design}}: the
#' wavelet coefficients of every additive component are penalized and the
#' \eqn{p} level terms \eqn{c_\ell} are not (decision D3). Two penalties are
#' available, the LASSO and the sparse group LASSO with one group per block
#' \eqn{(\ell, m)}, the unit that decision D12 keeps contiguous in the
#' columns of the design.
#'
#' The objective is
#' \deqn{\frac{1}{2n}\|y - Z\theta\|_2^2 + \lambda P(\theta),}
#' with \eqn{P(\theta) = \sum_{\ell m jk} |\theta_{\ell m, jk}|} for
#' \code{penalty = "lasso"} and
#' \deqn{P(\theta) = \alpha \sum_{\ell m jk} |\theta_{\ell m, jk}| +
#'   (1 - \alpha) \sum_{\ell m} \sqrt{N_J}\,\|\theta_{\ell m}\|_2}
#' for \code{penalty = "sglasso"}, where \eqn{\alpha} is \code{asparse}. In
#' both cases the level terms are outside the penalty. The \code{lambda} of
#' the returned object is the \eqn{\lambda} of this objective, which for the
#' LASSO is not the \eqn{\lambda} reported by \pkg{glmnet}; see
#' \code{\link{wafc_kkt}}.
#'
#' @param x Matrix (or data frame, or vector) of linear covariates, with
#'   \eqn{n} rows and \eqn{p} columns; a constant column, the usual
#'   \eqn{X_1 \equiv 1}, gives the additive intercept of the model. Ignored,
#'   and allowed to be missing, when \code{design} is supplied.
#' @param u Matrix (or data frame, or vector) of modulating covariates, with
#'   \eqn{n} rows and \eqn{q} columns. Ignored, and allowed to be missing,
#'   when \code{design} is supplied.
#' @param y Numeric response of length \eqn{n}.
#' @param J Resolution level of the sieve, passed to
#'   \code{\link{wafc_design}}.
#' @param penalty \code{"lasso"} (the default, decision D3) or
#'   \code{"sglasso"}, the sparse group LASSO by block \eqn{(\ell, m)},
#'   which needs the package \pkg{sparsegl}.
#' @param lambda Optional decreasing sequence of penalty levels, on the
#'   scale of the objective above. When \code{NULL} (the default) the
#'   engine builds its own path.
#' @param nlambda Length of the path built by the engine.
#' @param lambda.min.ratio Ratio between the smallest and the largest
#'   penalty level of that path. \code{NULL} uses the engine default,
#'   \eqn{0.01} when \eqn{n < } \code{nvars} and \eqn{10^{-4}} otherwise.
#' @param asparse Weight of the \eqn{\ell_1} part of the sparse group
#'   LASSO, in \eqn{[0, 1]}; ignored when \code{penalty = "lasso"}.
#'   \code{asparse = 0} is the pure group LASSO. \code{asparse = 1} is the
#'   LASSO in principle, but \pkg{sparsegl} does not converge reliably
#'   there and warns: the LASSO is \code{penalty = "lasso"}.
#' @param intercept Logical, or \code{NULL} (the default) for
#'   \code{TRUE} when the design has a constant column and \code{FALSE}
#'   otherwise. It cannot be \code{FALSE} when there is a constant column:
#'   the engines drop that column and its level would be lost.
#' @param standardize Passed to the engine. The default \code{FALSE} is the
#'   right one here, as in \code{WaveBased::wall}: the wavelet basis is
#'   orthonormal and standardizing the columns would change the penalty
#'   from level to level.
#' @param thresh Convergence threshold of the coordinate descent,
#'   \eqn{10^{-9}} by default. \code{NULL} uses the engine default,
#'   \eqn{10^{-7}} for \pkg{glmnet} and \eqn{10^{-8}} for \pkg{sparsegl}.
#'   The default is not the engine one because the design is not
#'   standardized (decision D17), so the relative criterion of the engine
#'   is read on a scale the penalty does not share: step E2.4b swept
#'   \eqn{10^{-7}} to \eqn{10^{-10}} over two scenarios, two sample sizes
#'   and three resolution levels, and \code{\link{wafc_kkt}} flagged
#'   between 19 and 73 points of the path in every cell at \eqn{10^{-7}}
#'   and in ten of twelve cells at \eqn{10^{-8}}, against none at
#'   \eqn{10^{-9}} but one cell with three. The selected penalty level is
#'   the same at every tolerance except one cell at \eqn{10^{-7}}, so what
#'   a looser value buys is speed and what it costs is the verification,
#'   not the tuning. Loosening it deliberately, for a simulation that
#'   checks the optimality conditions elsewhere, is a legitimate choice and
#'   is what the argument is for.
#' @param maxit Maximum number of passes; \code{NULL} uses the engine
#'   default.
#' @param design An object returned by \code{\link{wafc_design}}, to be
#'   used instead of building one from \code{x} and \code{u}. This is how
#'   several fits share one design.
#' @param ... Further arguments passed to \code{\link{wafc_design}}
#'   (\code{j0}, \code{family}, \code{filter.size}, \code{boundary},
#'   \code{rescale}, \code{eps}, \code{use.table}, \code{sparse}, ...).
#'
#' @return An object of class \code{"wafc"}: a list with the fitted
#'   \code{design}, the response \code{y}, the \code{penalty}, the path
#'   \code{lambda} on the scale of the objective, the raw coefficients
#'   \code{beta} (\code{nvars} by \code{nlambda}, as the engine returns
#'   them, so the coefficient of a constant covariate is zero) and
#'   \code{a0}, the levels \code{cc} (\eqn{p} by \code{nlambda}, with the
#'   intercept already folded into the constant covariate), the number
#'   \code{nzero} of non-zero wavelet coefficients, the fitted engine
#'   object \code{fit}, the factor \code{lambda.factor} between the two
#'   scales of \eqn{\lambda}, and the \code{group} vector of the sparse
#'   group LASSO.
#'
#' @seealso \code{\link{wafc_kkt}} for the verification of the
#'   optimality conditions, \code{\link{wafc_functions}} for the
#'   reconstruction of the components, \code{\link{predict.wafc}}.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' fit <- wafc(d$x, d$u, d$y, J = 3)
#' print(fit)
#' round(coef(fit, s = fit$lambda[30])[1:5, ], 4)
#'
#' @export
wafc <- function(x, u, y, J = 4L, penalty = c("lasso", "sglasso"),
                 lambda = NULL, nlambda = 100L, lambda.min.ratio = NULL,
                 asparse = 0.05, intercept = NULL, standardize = FALSE,
                 thresh = 1e-9, maxit = NULL, design = NULL, ...) {

  this_call <- match.call()
  penalty <- match.arg(penalty)

  if (is.null(design)) {
    if (missing(x) || missing(u)) {
      stop("Supply 'x' and 'u', or a 'design' built by wafc_design().",
           call. = FALSE)
    }
    design <- wafc_design(x, u, J = J, ...)
  } else if (!inherits(design, "wafc_design")) {
    stop("'design' must be an object returned by wafc_design().", call. = FALSE)
  }

  n <- design[["n"]]
  y <- wafc_check_y(y, n)
  nvars <- design[["nvars"]]
  npen <- nvars - length(design[["unpenalized"]])
  if (npen == 0L) {
    stop("The design has no penalized column; there is nothing to fit.",
         call. = FALSE)
  }

  ## The carrier of the intercept: at most one constant linear covariate,
  ## whose level the engine reads into its own intercept (note (i) above).
  carrier <- wafc_carrier(design)
  if (is.null(intercept)) intercept <- !is.na(carrier[["index"]])
  if (!isTRUE(intercept) && !is.na(carrier[["index"]])) {
    stop("intercept = FALSE with the constant linear covariate '",
         design[["xnames"]][carrier[["index"]]],
         "' in the design: the engine drops that column from the fit and ",
         "its level c would be lost. Use intercept = TRUE, which is the ",
         "default, or remove the column from 'x'.", call. = FALSE)
  }

  ## Scale of lambda. glmnet rescales penalty.factor to sum to nvars, so
  ## what multiplies the L1 norm of the wavelet coefficients is lambda
  ## times nvars/npen; sparsegl does not rescale (note (ii) above).
  lambda.factor <- if (penalty == "lasso") nvars / npen else 1
  lambda <- wafc_check_lambda(lambda)
  group <- wafc_groups(design)

  args <- list(x = design[["Z"]], y = y, intercept = intercept,
               standardize = standardize)
  if (!is.null(lambda)) {
    args[["lambda"]] <- lambda / lambda.factor
  } else {
    args[["nlambda"]] <- as.integer(nlambda)
  }

  if (penalty == "lasso") {
    args[["family"]] <- "gaussian"
    args[["penalty.factor"]] <- design[["penalty.factor"]]
    ctrl <- list()
    if (!is.null(thresh)) ctrl[["thresh"]] <- thresh
    if (!is.null(maxit)) ctrl[["maxit"]] <- maxit
    if (length(ctrl) > 0L) args[["control"]] <- ctrl
    if (!is.null(lambda.min.ratio) && is.null(lambda)) {
      args[["lambda.min.ratio"]] <- lambda.min.ratio
    }
    fit <- do.call(glmnet::glmnet, args)
    a0 <- as.numeric(fit[["a0"]])
    beta <- fit[["beta"]]
    lam <- fit[["lambda"]] * lambda.factor
  } else {
    if (!requireNamespace("sparsegl", quietly = TRUE)) {
      stop("penalty = \"sglasso\" needs the package 'sparsegl' ",
           "(see docs/CONTINUAR.md, section 2).", call. = FALSE)
    }
    if (length(asparse) != 1L || !is.finite(asparse) || asparse < 0 ||
        asparse > 1) {
      stop("'asparse' must be a single value in [0, 1].", call. = FALSE)
    }
    if (asparse == 1) {
      warning("asparse = 1 removes the group term, and sparsegl does not ",
              "converge reliably on that corner: the optimality conditions ",
              "are violated by a few percent of lambda (see wafc_kkt()). ",
              "Use penalty = \"lasso\" for the LASSO.", call. = FALSE)
    }
    args[["group"]] <- group[["group"]]
    args[["pf_group"]] <- group[["pf_group"]]
    args[["pf_sparse"]] <- group[["pf_sparse"]]
    args[["asparse"]] <- asparse
    args[["family"]] <- "gaussian"
    if (!is.null(thresh)) args[["eps"]] <- thresh
    if (!is.null(maxit)) args[["maxit"]] <- maxit
    if (!is.null(lambda.min.ratio) && is.null(lambda)) {
      args[["lambda.factor"]] <- lambda.min.ratio
    }
    fit <- do.call(sparsegl::sparsegl, args)
    a0 <- as.numeric(fit[["b0"]])
    beta <- fit[["beta"]]
    lam <- fit[["lambda"]]
  }

  dimnames(beta) <- list(colnames(design[["Z"]]), NULL)
  pen_idx <- seq_len(nvars)[-design[["unpenalized"]]]
  ## At a penalty level that kills every penalized coefficient the solution
  ## is the least squares fit on the level terms alone. sparsegl leaves the
  ## whole vector at zero at the first point of its own path, a warm start
  ## artifact that the KKT check catches; glmnet already returns the right
  ## values and the correction leaves them alone.
  fixed <- wafc_fix_null_point(beta, a0, design, carrier, y, intercept, pen_idx)
  beta <- fixed[["beta"]]
  a0 <- fixed[["a0"]]
  out <- list(design = design, y = y, penalty = penalty, lambda = lam,
              beta = beta, a0 = a0,
              cc = wafc_levels(beta, a0, design, carrier),
              nzero = as.integer(Matrix::colSums(beta[pen_idx, , drop = FALSE] != 0)),
              intercept = intercept, asparse = if (penalty == "sglasso") asparse else NA_real_,
              lambda.factor = lambda.factor, group = group,
              carrier = carrier, fit = fit, n = n, p = design[["p"]],
              q = design[["q"]], nvars = nvars, npen = npen,
              call = this_call)
  class(out) <- "wafc"
  out
}

#' Coefficients of a WAFC fit
#'
#' The coefficients in the parametrisation of the model, that is with the
#' intercept of the engine already folded into the level of the constant
#' linear covariate: the row \code{"(Intercept)"} is zero whenever the
#' design has such a covariate, and the linear predictor is always the
#' intercept plus the design times the remaining rows.
#'
#' @param object An object of class \code{"wafc"}.
#' @param s Penalty levels at which the coefficients are wanted, on the
#'   scale of the objective of \code{\link{wafc}}. \code{NULL} returns the
#'   whole path. Values between two points of the path are interpolated
#'   linearly, as in \code{\link[glmnet]{coef.glmnet}}.
#' @param ... Ignored.
#'
#' @return A matrix with \code{1 + nvars} rows, the first one the
#'   intercept, and one column per value of \code{s}.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' fit <- wafc(d$x, d$u, d$y, J = 3)
#' cf <- coef(fit, s = fit$lambda[40])
#' cf[1:4, ]
#'
#' @export
coef.wafc <- function(object, s = NULL, ...) {
  cf <- wafc_raw_coef(object, s)
  a0 <- cf[1L, , drop = TRUE]
  carrier <- object[["carrier"]]
  if (!is.na(carrier[["index"]])) {
    ## the engine left this coefficient at zero and read the level into its
    ## intercept; the column is constant equal to carrier$value, so the
    ## level is the intercept divided by that value
    cf[1L + carrier[["index"]], ] <- cf[1L + carrier[["index"]], ] +
      a0 / carrier[["value"]]
    cf[1L, ] <- 0
  }
  cf
}

#' Predictions and derived quantities of a WAFC fit
#'
#' @param object An object of class \code{"wafc"}.
#' @param newx,newu New linear and modulating covariates, with the same
#'   numbers of columns as the data the model was fitted on. The design is
#'   rebuilt with the \code{spec} of the fitted one, so the basis and the
#'   rescaling are the ones of the training sample and modulating
#'   covariates outside its range are truncated to it. Needed by
#'   \code{type = "response"}; \code{type = "beta"} needs only \code{newu}.
#'   When both are missing the training sample is used.
#' @param s Penalty levels, on the scale of the objective of
#'   \code{\link{wafc}}. \code{NULL} means the whole path, except for
#'   \code{type = "beta"}, which needs a single value.
#' @param type \code{"response"} (or \code{"link"}, the same thing in the
#'   Gaussian case) for the fitted values, \code{"beta"} for the matrix of
#'   functional coefficients \eqn{\hat\beta_\ell(U_i)}, \code{"coefficients"}
#'   for \code{\link{coef.wafc}}, \code{"nonzero"} for the indices of the
#'   non-zero wavelet coefficients.
#' @param ... Ignored.
#'
#' @return A matrix of fitted values (\eqn{n} by the number of penalty
#'   levels), a matrix of functional coefficients (\eqn{n} by \eqn{p}), a
#'   matrix of coefficients, or a list of index vectors.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' fit <- wafc(d$x, d$u, d$y, J = 3)
#' yhat <- predict(fit, d$x, d$u, s = fit$lambda[40])
#' b <- predict(fit, newu = d$u, s = fit$lambda[40], type = "beta")
#' head(round(b, 3))
#'
#' @export
predict.wafc <- function(object, newx, newu, s = NULL,
                         type = c("response", "link", "coefficients",
                                  "beta", "nonzero"), ...) {
  type <- match.arg(type)
  if (type == "coefficients") return(coef.wafc(object, s = s))
  if (type == "nonzero") {
    cf <- wafc_raw_coef(object, s)[-1L, , drop = FALSE]
    pen <- seq_len(object[["nvars"]])[-object[["design"]][["unpenalized"]]]
    return(lapply(seq_len(ncol(cf)), function(k) pen[cf[pen, k] != 0]))
  }
  if (type == "beta") {
    if (missing(newu)) newu <- NULL
    return(wafc_beta_hat(object, newu, s))
  }
  design <- object[["design"]]
  if (missing(newx) && missing(newu)) {
    Z <- design[["Z"]]
  } else {
    if (missing(newx) || missing(newu)) {
      stop("Supply both 'newx' and 'newu', or neither.", call. = FALSE)
    }
    Z <- wafc_design(newx, newu, spec = design)[["Z"]]
  }
  ## raw coefficients: the engine intercept plus the design, which is
  ## correct whatever the value of a constant covariate in the new data
  cf <- wafc_raw_coef(object, s)
  eta <- as.matrix(Z %*% cf[-1L, , drop = FALSE])
  eta <- sweep(eta, 2L, cf[1L, ], "+")
  dimnames(eta) <- list(rownames(Z), colnames(cf))
  eta
}

#' Print a WAFC fit
#'
#' @param x An object of class \code{"wafc"}.
#' @param digits Number of significant digits printed.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' print(wafc(d$x, d$u, d$y, J = 3))
#'
#' @export
print.wafc <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  d <- x[["design"]]
  cat("WAFC fit:", switch(x[["penalty"]], lasso = "LASSO",
                          sglasso = sprintf("sparse group LASSO (asparse = %s)",
                                            format(x[["asparse"]], digits = digits))),
      "\n")
  cat(sprintf("  n = %d, p = %d, q = %d, J = %s, %d columns (%d penalized, %d level term(s))\n",
              x[["n"]], x[["p"]], x[["q"]], paste(d[["J"]], collapse = ", "),
              x[["nvars"]], x[["npen"]], length(d[["unpenalized"]])))
  cat(sprintf("  %d penalty level(s), lambda from %s down to %s\n",
              length(x[["lambda"]]), format(max(x[["lambda"]]), digits = digits),
              format(min(x[["lambda"]]), digits = digits)))
  cat(sprintf("  non-zero wavelet coefficients: %d to %d of %d\n",
              min(x[["nzero"]]), max(x[["nzero"]]), x[["npen"]]))
  if (!is.na(x[["carrier"]][["index"]])) {
    cat(sprintf("  the level of '%s' is carried by the intercept of the engine\n",
                d[["xnames"]][x[["carrier"]][["index"]]]))
  }
  invisible(x)
}

#' Verify the optimality conditions of a WAFC fit
#'
#' Checks the Karush-Kuhn-Tucker conditions of the objective of
#' \code{\link{wafc}} at one or more points of the path. This is how the
#' scale of \eqn{\lambda} is established rather than assumed: \pkg{glmnet}
#' rescales \code{penalty.factor} so that it sums to \code{nvars}, so the
#' penalty it applies to a wavelet coefficient is its own \eqn{\lambda}
#' times \code{nvars/npen}, and it is the stationarity conditions that say
#' so.
#'
#' For the LASSO, writing \eqn{r} for the residual and
#' \eqn{g = Z' r / n}, the conditions are \eqn{|g_a| = \lambda} on the
#' non-zero wavelet coefficients, \eqn{|g_a| \le \lambda} on the zero ones
#' and \eqn{g_a = 0} on the level terms. For the sparse group LASSO they
#' are the usual ones of the block-separable subdifferential: on a block
#' that is entirely zero,
#' \eqn{\|S(g_{\ell m}, \lambda\alpha)\|_2 \le \lambda(1-\alpha)\sqrt{N_J}}
#' with \eqn{S} the soft threshold, and on an active block the
#' stationarity equation holds coordinate by coordinate.
#'
#' @param object An object of class \code{"wafc"}.
#' @param s Penalty levels to check. \code{NULL} checks the whole path.
#' @param tol Tolerance, relative to \eqn{\lambda}, for the stationarity
#'   residual and for the excess of a subgradient over its bound.
#' @param tol.abs Absolute tolerance for the same two residuals and for the
#'   gradient on the unpenalized terms. A penalty level passes when each
#'   residual is below \code{tol.abs} or below \code{tol} times
#'   \eqn{\lambda}.
#'
#' @return A data frame with one row per penalty level: \code{lambda}, the
#'   engine's \code{lambda.engine}, the number \code{nzero} of non-zero
#'   wavelet coefficients, the largest stationarity residual
#'   \code{stationarity} and the same residual \code{stationarity.rel}
#'   relative to \eqn{\lambda}, the largest subgradient ratio
#'   \code{subgradient} (at most one at a solution) and the largest
#'   absolute excess \code{subgradient.abs} over the corresponding bound,
#'   the largest gradient \code{unpenalized} on the level terms, the
#'   gradient \code{intercept} of the intercept, and \code{ok}.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' fit <- wafc(d$x, d$u, d$y, J = 3)
#' wafc_kkt(fit, s = fit$lambda[c(20, 40, 60)])
#'
#' @export
wafc_kkt <- function(object, s = NULL, tol = 1e-2, tol.abs = 1e-4) {
  if (!inherits(object, "wafc")) {
    stop("'object' must be an object returned by wafc().", call. = FALSE)
  }
  design <- object[["design"]]
  Z <- design[["Z"]]
  y <- object[["y"]]
  n <- object[["n"]]
  unp <- design[["unpenalized"]]
  carrier <- object[["carrier"]]
  ## the constant covariate is out: its gradient is the one of the
  ## intercept, which is the equation sum(r) = 0
  unp_free <- setdiff(unp, carrier[["index"]])
  cf <- wafc_raw_coef(object, s)
  lam <- if (is.null(s)) object[["lambda"]] else as.numeric(s)
  soft <- function(z, t) sign(z) * pmax(abs(z) - t, 0)
  out <- data.frame(lambda = lam, lambda.engine = lam / object[["lambda.factor"]],
                    nzero = NA_integer_, stationarity = NA_real_,
                    stationarity.rel = NA_real_, subgradient = NA_real_,
                    subgradient.abs = NA_real_, unpenalized = NA_real_,
                    intercept = NA_real_)
  for (k in seq_along(lam)) {
    b <- cf[-1L, k]
    r <- as.numeric(y - cf[1L, k] - Z %*% b)
    g <- as.numeric(Matrix::crossprod(Z, r)) / n
    le <- lam[k]
    stat <- 0
    sub <- 0
    sub_abs <- 0
    if (object[["penalty"]] == "lasso") {
      pen <- seq_along(b)[-unp]
      act <- pen[b[pen] != 0]
      inact <- pen[b[pen] == 0]
      if (length(act)) stat <- max(abs(abs(g[act]) - le))
      if (length(inact)) {
        sub <- max(abs(g[inact])) / le
        sub_abs <- max(0, max(abs(g[inact])) - le)
      }
    } else {
      grp <- object[["group"]]
      a <- object[["asparse"]]
      for (gg in grp[["penalized"]]) {
        idx <- which(grp[["group"]] == gg)
        pg <- grp[["pf_group"]][gg]
        ps <- grp[["pf_sparse"]][idx]
        bg <- b[idx]
        if (all(bg == 0)) {
          nrm <- sqrt(sum(soft(g[idx], le * a * ps)^2))
          bound <- le * (1 - a) * pg
          if (bound > 0) {
            sub <- max(sub, nrm / bound)
            sub_abs <- max(sub_abs, nrm - bound)
          } else {
            ## asparse = 1: no group term, so the condition is the one of
            ## the LASSO coordinate by coordinate
            sub <- max(sub, max(abs(g[idx]) / (le * a * ps)))
            sub_abs <- max(sub_abs, max(abs(g[idx]) - le * a * ps))
          }
        } else {
          nz <- bg != 0
          v <- g[idx][nz] - le * a * ps[nz] * sign(bg[nz]) -
            le * (1 - a) * pg * bg[nz] / sqrt(sum(bg^2))
          stat <- max(stat, max(abs(v)))
          if (any(!nz)) {
            bound <- le * a * ps[!nz]
            sub <- max(sub, max(abs(g[idx][!nz]) / bound))
            sub_abs <- max(sub_abs, max(abs(g[idx][!nz]) - bound))
          }
        }
      }
    }
    out[k, "nzero"] <- sum(b[-unp] != 0)
    out[k, "stationarity"] <- stat
    out[k, "stationarity.rel"] <- if (le > 0) stat / le else NA_real_
    out[k, "subgradient"] <- sub
    out[k, "subgradient.abs"] <- max(0, sub_abs)
    out[k, "unpenalized"] <- if (length(unp_free)) max(abs(g[unp_free])) else 0
    out[k, "intercept"] <- if (object[["intercept"]]) abs(sum(r)) / n else 0
  }
  ## A residual passes when it is below the absolute tolerance or below
  ## 'tol' times lambda: the coordinate descent stops on the objective, so
  ## its gradient error has a fixed size and a purely relative rule would
  ## fail at the small end of the path by construction.
  slack <- pmax(tol.abs, tol * out[["lambda"]])
  out[["ok"]] <- out[["stationarity"]] <= slack &
    out[["subgradient.abs"]] <= slack &
    out[["unpenalized"]] <= tol.abs &
    out[["intercept"]] <= tol.abs
  out
}

## ---------------------------------------------------------------------------
## Internals
## ---------------------------------------------------------------------------

wafc_check_y <- function(y, n) {
  if (is.matrix(y) && ncol(y) == 1L) y <- as.numeric(y)
  if (!is.numeric(y) || length(y) != n) {
    stop("'y' must be numeric of length ", n, ".", call. = FALSE)
  }
  if (anyNA(y) || any(!is.finite(y))) {
    stop("'y' must not contain NA or non-finite values.", call. = FALSE)
  }
  as.numeric(y)
}

wafc_check_lambda <- function(lambda) {
  if (is.null(lambda)) return(NULL)
  if (!is.numeric(lambda) || length(lambda) == 0L || any(!is.finite(lambda)) ||
      any(lambda < 0)) {
    stop("'lambda' must be a vector of non-negative finite values.",
         call. = FALSE)
  }
  sort(as.numeric(lambda), decreasing = TRUE)
}

## The linear covariate whose level the engine reads into its intercept:
## the single constant column of the design, if there is one. More than one
## makes the level terms collinear and is refused here rather than silently
## fitted.
wafc_carrier <- function(design) {
  idx <- design[["constant"]]
  if (length(idx) == 0L) return(list(index = NA_integer_, value = NA_real_))
  if (length(idx) > 1L) {
    stop("The design has ", length(idx),
         " constant linear covariates (", paste(design[["xnames"]][idx],
                                                collapse = ", "),
         "): their level terms are not identified. Keep one of them.",
         call. = FALSE)
  }
  value <- as.numeric(design[["Z"]][1L, idx])
  if (value == 0) {
    stop("The linear covariate '", design[["xnames"]][idx],
         "' is identically zero; remove it from 'x'.", call. = FALSE)
  }
  list(index = as.integer(idx), value = value)
}

## Groups of the sparse group LASSO: one group per block (l, m), which D12
## keeps contiguous, plus the group of the level terms, which is left
## unpenalized by zero penalty factors (sparsegl does not rescale them).
wafc_groups <- function(design) {
  nvars <- design[["nvars"]]
  grp <- integer(nvars)
  grp[design[["unpenalized"]]] <- 1L
  blocks <- design[["blocks"]]
  for (b in seq_along(blocks)) grp[blocks[[b]]] <- 1L + b
  sizes <- as.integer(tabulate(grp, nbins = 1L + length(blocks)))
  pf_group <- sqrt(sizes)
  pf_group[1L] <- 0
  pf_sparse <- rep_len(1, nvars)
  pf_sparse[design[["unpenalized"]]] <- 0
  list(group = grp, pf_group = pf_group, pf_sparse = pf_sparse,
       sizes = sizes, penalized = 1L + seq_along(blocks),
       names = c("(levels)", names(blocks)))
}

## Level terms c_l of the model along the path, with the intercept of the
## engine folded into the constant covariate.
wafc_levels <- function(beta, a0, design, carrier) {
  unp <- design[["unpenalized"]]
  cc <- as.matrix(beta[unp, , drop = FALSE])
  if (!is.na(carrier[["index"]])) {
    row <- match(carrier[["index"]], unp)
    cc[row, ] <- cc[row, ] + a0 / carrier[["value"]]
  }
  rownames(cc) <- design[["xnames"]]
  cc
}

## Least squares on the level terms alone, written into the path points at
## which every penalized coefficient is zero. The constant covariate stays
## at zero and the intercept carries its level, which is the convention of
## both engines (note (i) at the top of the file).
wafc_fix_null_point <- function(beta, a0, design, carrier, y, intercept,
                                pen_idx) {
  null_pt <- which(Matrix::colSums(beta[pen_idx, , drop = FALSE] != 0) == 0)
  if (length(null_pt) == 0L) return(list(beta = beta, a0 = a0))
  unp <- design[["unpenalized"]]
  free <- setdiff(unp, carrier[["index"]])
  W <- as.matrix(design[["Z"]][, free, drop = FALSE])
  if (intercept) W <- cbind(`(Intercept)` = 1, W)
  cf <- if (ncol(W) == 0L) numeric(0) else {
    qrW <- qr(W)
    v <- qr.coef(qrW, y)
    v[is.na(v)] <- 0
    v
  }
  a0_new <- if (intercept) cf[[1L]] else 0
  b_new <- if (intercept) cf[-1L] else cf
  for (k in null_pt) {
    beta[unp, k] <- 0
    if (length(free)) beta[free, k] <- b_new
    a0[k] <- a0_new
  }
  list(beta = beta, a0 = a0)
}

## Coefficients in the parametrisation of the engines: first row the
## intercept, then the nvars columns of the design. The path is the
## corrected one stored in the object, so the interpolation cannot be
## delegated to the engine; the rule is the linear interpolation in lambda
## of glmnet's internal lambda.interp, read and rewritten here (the same
## rule this project applies to wall(), docs/instrucoes.md section 6).
wafc_raw_coef <- function(object, s = NULL) {
  if (!inherits(object, "wafc")) {
    stop("'object' must be an object returned by wafc().", call. = FALSE)
  }
  cf <- rbind(object[["a0"]], as.matrix(object[["beta"]]))
  rownames(cf) <- c("(Intercept)", colnames(object[["design"]][["Z"]]))
  if (is.null(s)) {
    colnames(cf) <- paste0("s", seq_along(object[["lambda"]]) - 1L)
    return(cf)
  }
  if (!is.numeric(s) || length(s) == 0L || any(!is.finite(s)) || any(s < 0)) {
    stop("'s' must be a vector of non-negative finite values.", call. = FALSE)
  }
  ip <- wafc_lambda_interp(object[["lambda"]], as.numeric(s))
  out <- cf[, ip[["left"]], drop = FALSE] * rep(ip[["frac"]], each = nrow(cf)) +
    cf[, ip[["right"]], drop = FALSE] * rep(1 - ip[["frac"]], each = nrow(cf))
  dimnames(out) <- list(rownames(cf), paste0("s", seq_along(s) - 1L))
  out
}

## Neighbours and weights of the linear interpolation of the path at the
## penalty levels 's', which are truncated to the range of the path.
wafc_lambda_interp <- function(lambda, s) {
  k <- length(lambda)
  if (k == 1L) {
    return(list(left = rep(1L, length(s)), right = rep(1L, length(s)),
                frac = rep(1, length(s))))
  }
  s <- pmin(pmax(s, min(lambda)), max(lambda))
  ## both mapped to [0, 1] with the path decreasing, as in lambda.interp
  span <- lambda[1L] - lambda[k]
  sfrac <- (lambda[1L] - s) / span
  lfrac <- (lambda[1L] - lambda) / span
  coord <- stats::approx(lfrac, seq_along(lfrac), sfrac)[["y"]]
  left <- as.integer(floor(coord))
  right <- as.integer(ceiling(coord))
  frac <- (sfrac - lfrac[right]) / (lfrac[left] - lfrac[right])
  frac[left == right] <- 1
  frac[abs(lfrac[left] - lfrac[right]) < .Machine[["double.eps"]]] <- 1
  list(left = left, right = right, frac = frac)
}
