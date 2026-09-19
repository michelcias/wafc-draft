## wafc/R/tune.R -- choosing the pair (J, lambda) of the WAFC fit (step
## E2.3): cross-validation over the grid, the information criteria, and the
## resolution rule that the theory of E1.6 prescribes.
##
## Four selection rules live here, and they answer different questions.
##
##  (i)   cv.wafc() cross-validates over (J, lambda) with folds that are
##        fixed once for the whole grid, as cv.wall() does in WaveBased: one
##        design per candidate J, one lambda path per candidate J, and the
##        pair (J, lambda) that minimises the cross-validated loss.
##        lambda.min and lambda.1se are read inside the selected J.
##  (ii)  wafc_bic() and wafc_ebic() score a fitted path with degrees of
##        freedom equal to the number of non-zero coefficients, which is the
##        unbiased count for the LASSO in an orthonormal-like design
##        (Zou, Hastie and Tibshirani, 2007); wafc_tune() minimises them
##        over the same grid of J, so the three rules are compared on the
##        same candidates.
##  (iii) wafc_J_theory() and wafc_lambda_theory() implement the pair that
##        Theorem 2 of E1.6 balances: the resolution
##        J_n = min{J : 2^J >= (n/log n)^{1/(2s'+1)}} of equation (Jn) and
##        the penalty level lambda_n of Corollary 2 of E1.5. This rule uses
##        no data-driven search at all; it needs the effective regularity
##        s' of the components and the error scale sigma, and what it costs
##        to use it is the fourth column of the comparison of E2.3.
##
## Two things this file inherits from E2.2 and does not renegotiate
## (decision D17): the lambda of every object is the lambda of the
## objective of E1.5, not the one of the engine, and wafc_kkt() is what
## verifies it. Every lambda that goes in or comes out of the functions
## below is therefore on the scale of the objective.

#' Cross-validation of a WAFC fit over the resolution level and the penalty
#'
#' \eqn{k}-fold cross-validation over the pair \eqn{(J, \lambda)}: for each
#' candidate resolution level the design is built once on the whole sample,
#' the penalty path is the one of the fit on the whole sample, and the folds
#' are the same for every candidate, so that the levels are compared on the
#' same partition. This is the scheme of \code{\link[WaveBased]{cv.wall}},
#' read and rewritten for the design of \code{\link{wafc_design}}
#' (\file{docs/instrucoes.md}, section 6).
#'
#' The fold fits reuse the columns of the design of the whole sample rather
#' than rebuilding a basis on each training set. As in \code{cv.wall}, this
#' means that the rescaling of the modulating covariates
#' (\code{\link{wafc_rescale}}) is the one of the whole sample: it is a
#' range, not a fitted quantity, and keeping it fixed is what makes the
#' columns of the folds comparable.
#'
#' @param x,u,y The data, as in \code{\link{wafc}}.
#' @param J The grid of candidate resolution levels. \code{NULL} (the
#'   default) uses \code{2:max(2, ceiling(log2(n)/2))}, the rule of
#'   \code{cv.wall} truncated below at \eqn{J = 2}: the default rescaling
#'   constant \eqn{\epsilon = 1.9^{-J}} of \code{\link{wafc_design}} is
#'   not admissible at \eqn{J = 1}.
#' @param penalty \code{"lasso"} or \code{"sglasso"}, as in
#'   \code{\link{wafc}}.
#' @param nfolds Number of folds, at least 3.
#' @param foldid Optional vector of length \eqn{n} with the fold of each
#'   observation, which overrides \code{nfolds}. Supplying it is how two
#'   calls are made to use the same partition.
#' @param lambda Optional penalty path, on the scale of the objective,
#'   used for every candidate \eqn{J}. \code{NULL} (the default) lets each
#'   candidate have the path its own fit on the whole sample produces,
#'   which is what \code{cv.glmnet} does at a fixed design.
#' @param nlambda,lambda.min.ratio Length and range of those paths.
#' @param type.measure \code{"mse"} (the default) or \code{"mae"}: the loss
#'   averaged over the held-out observations.
#' @param trace If \code{TRUE}, prints one line per candidate \eqn{J}.
#' @param ... Further arguments passed to \code{\link{wafc}} and through it
#'   to \code{\link{wafc_design}}.
#'
#' @return An object of class \code{"cv.wafc"}: a list with the grid
#'   \code{J}, the list \code{cv} of per-candidate results (\code{lambda},
#'   \code{cvm}, \code{cvsd}, \code{cvup}, \code{cvlo}, \code{nzero},
#'   \code{lambda.min}, \code{lambda.1se} and the values attained at
#'   \code{lambda.min}), the summary \code{cvtab}, the selected
#'   \code{J.min}, \code{lambda.min} and \code{lambda.1se}, the
#'   \code{foldid} used, and \code{wafc.fit}, the fit on the whole sample at
#'   \code{J.min}, which is the fit the selected pair refers to.
#'
#' @seealso \code{\link{wafc_bic}}, \code{\link{wafc_tune}},
#'   \code{\link{wafc_J_theory}}.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' cvfit <- cv.wafc(d$x, d$u, d$y, J = 2:4, nfolds = 5)
#' cvfit
#' cf <- coef(cvfit, s = "lambda.1se")
#'
#' @export
cv.wafc <- function(x, u, y, J = NULL, penalty = c("lasso", "sglasso"),
                    nfolds = 10L, foldid = NULL, lambda = NULL,
                    nlambda = 100L, lambda.min.ratio = NULL,
                    type.measure = c("mse", "mae"), trace = FALSE, ...) {

  this_call <- match.call()
  penalty <- match.arg(penalty)
  type.measure <- match.arg(type.measure)
  x <- wafc_as_matrix(x, "x")
  u <- wafc_as_matrix(u, "u")
  n <- nrow(x)
  y <- wafc_check_y(y, n)
  J <- wafc_J_grid(J, n)
  foldid <- wafc_foldid(foldid, n, nfolds)
  nfolds <- max(foldid)
  lambda <- wafc_check_lambda(lambda)

  loss <- switch(type.measure,
                 mse = function(e) e^2,
                 mae = function(e) abs(e))

  run_one <- function(Ji) {
    design <- wafc_design(x, u, J = Ji, ...)
    full <- wafc(design = design, y = y, penalty = penalty, lambda = lambda,
                 nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, ...)
    z <- wafc_cv_design(design, y, full, foldid, loss, penalty, ...)
    z[["J"]] <- Ji
    list(cv = z, fit = full)
  }

  runs <- vector("list", length(J))
  for (i in seq_along(J)) {
    runs[[i]] <- run_one(J[i])
    if (trace) {
      z <- runs[[i]][["cv"]]
      cat(sprintf("J = %d: %s = %.5f at lambda = %.5g (%d nonzero of %d)\n",
                  J[i], type.measure, z[["cvm.min"]], z[["lambda.min"]],
                  z[["nzero.min"]], runs[[i]][["fit"]][["npen"]]))
    }
  }

  cvlist <- lapply(runs, `[[`, "cv")
  cvm_all <- vapply(cvlist, `[[`, 0, "cvm.min")
  ## ties go to the smallest J: the grid is increasing and which.min takes
  ## the first minimum, which is the parsimonious reading
  best <- which.min(cvm_all)
  z <- cvlist[[best]]

  cvtab <- data.frame(J = J,
                      lambda.min = vapply(cvlist, `[[`, 0, "lambda.min"),
                      measure = cvm_all,
                      sd = vapply(cvlist, `[[`, 0, "cvsd.min"),
                      nzero = vapply(cvlist,
                                     function(w) as.integer(w[["nzero.min"]]),
                                     0L))
  names(cvtab)[3L] <- type.measure

  out <- list(call = this_call, J = J, cv = cvlist, cvtab = cvtab,
              J.min = J[best], lambda.min = z[["lambda.min"]],
              lambda.1se = z[["lambda.1se"]], cvm.min = z[["cvm.min"]],
              cvsd.min = z[["cvsd.min"]], nzero.min = z[["nzero.min"]],
              type.measure = type.measure, nfolds = nfolds, foldid = foldid,
              penalty = penalty, wafc.fit = runs[[best]][["fit"]])
  class(out) <- "cv.wafc"
  out
}

#' @rdname cv.wafc
#' @param digits Number of significant digits printed.
#' @export
print.cv.wafc <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCross-validated WAFC fit\n\n")
  cat("Call:", deparse(x[["call"]]), "\n\n")
  cat("Measure:", x[["type.measure"]], "(", x[["nfolds"]], "folds )\n\n")
  tab <- x[["cvtab"]]
  tab[["lambda.min"]] <- signif(tab[["lambda.min"]], digits)
  tab[[x[["type.measure"]]]] <- signif(tab[[x[["type.measure"]]]], digits)
  tab[["sd"]] <- signif(tab[["sd"]], digits)
  print(tab, row.names = FALSE)
  cat("\nSelected: J =", x[["J.min"]],
      "with lambda.min =", signif(x[["lambda.min"]], digits),
      "( lambda.1se =", signif(x[["lambda.1se"]], digits), ")\n")
  invisible(x)
}

#' Coefficients and predictions at the cross-validated pair
#'
#' @param object An object of class \code{"cv.wafc"}.
#' @param s \code{"lambda.min"}, \code{"lambda.1se"}, or a numeric penalty
#'   level on the scale of the objective of \code{\link{wafc}}.
#' @param newx,newu New covariates, as in \code{\link{predict.wafc}}.
#' @param ... Passed to the method of the underlying \code{"wafc"} object.
#'
#' @return As \code{\link{coef.wafc}} and \code{\link{predict.wafc}}, at the
#'   resolution level \code{J.min} selected by the cross-validation.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' cvfit <- cv.wafc(d$x, d$u, d$y, J = 2:3, nfolds = 5)
#' dim(coef(cvfit))
#' head(predict(cvfit, d$x, d$u, s = "lambda.min"))
#'
#' @export
coef.cv.wafc <- function(object, s = c("lambda.min", "lambda.1se"), ...) {
  coef(object[["wafc.fit"]], s = wafc_cv_s(object, s), ...)
}

#' @rdname coef.cv.wafc
#' @export
predict.cv.wafc <- function(object, newx, newu,
                            s = c("lambda.min", "lambda.1se"), ...) {
  s <- wafc_cv_s(object, s)
  if (missing(newx) && missing(newu)) {
    return(predict(object[["wafc.fit"]], s = s, ...))
  }
  if (missing(newx)) return(predict(object[["wafc.fit"]], newu = newu, s = s, ...))
  predict(object[["wafc.fit"]], newx, newu, s = s, ...)
}

#' Information criteria of a WAFC fit
#'
#' The Bayesian information criterion and its extended version (Chen and
#' Chen, 2008) along the penalty path of a fit, with degrees of freedom
#' equal to the number of non-zero coefficients: the non-zero wavelet
#' coefficients plus the \eqn{p} level terms \eqn{c_\ell}, which are never
#' penalized and always count. Writing \eqn{k} for the number of non-zero
#' wavelet coefficients, \eqn{df = k + p}, \eqn{d} for the number of
#' penalized columns and \eqn{RSS} for the residual sum of squares,
#' \deqn{BIC = n\log(RSS/n) + df\log n, \qquad
#'       EBIC = BIC + 2\gamma\log\binom{d}{k},}
#' so that \code{gamma = 0} returns the BIC. The scale \eqn{\sigma} is
#' profiled out, which is the usual form when it is unknown.
#'
#' @param object An object of class \code{"wafc"}.
#' @param s Penalty levels at which the criterion is wanted, on the scale of
#'   the objective of \code{\link{wafc}}. \code{NULL} (the default) is the
#'   whole path of the fit.
#' @param gamma The parameter of the extended criterion, in \eqn{[0, 1]}.
#'   The default \eqn{1} is the most conservative member of the family and
#'   the usual choice when the number of columns is of the order of the
#'   sample size.
#'
#' @return A data frame with one row per penalty level: \code{lambda},
#'   \code{nzero} (non-zero wavelet coefficients), \code{df}, \code{rss} and
#'   the column \code{bic} or \code{ebic}. The index of the minimising row
#'   is the attribute \code{"which.min"} and the penalty level attaining it
#'   the attribute \code{"lambda.min"}.
#'
#' @references Chen, J. and Chen, Z. (2008). Extended Bayesian information
#'   criteria for model selection with large model spaces. \emph{Biometrika}
#'   95(3), 759-771.
#'
#'   Zou, H., Hastie, T. and Tibshirani, R. (2007). On the degrees of
#'   freedom of the lasso. \emph{The Annals of Statistics} 35(5), 2173-2192.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' fit <- wafc(d$x, d$u, d$y, J = 3)
#' ic <- wafc_bic(fit)
#' attr(ic, "lambda.min")
#'
#' @export
wafc_bic <- function(object, s = NULL) {
  wafc_ic(object, s = s, gamma = 0, name = "bic")
}

#' @rdname wafc_bic
#' @export
wafc_ebic <- function(object, s = NULL, gamma = 1) {
  if (length(gamma) != 1L || !is.finite(gamma) || gamma < 0 || gamma > 1) {
    stop("'gamma' must be a single value in [0, 1].", call. = FALSE)
  }
  wafc_ic(object, s = s, gamma = gamma, name = "ebic")
}

#' Resolution level and penalty level prescribed by the theory
#'
#' The pair that Theorem 2 of \file{derivations/05-taxas.tex} balances. The
#' resolution is
#' \deqn{J_n = \min\{J \in \mathbb{N} : 2^J \ge (n/\log n)^{1/(2s'+1)}\},}
#' equation (Jn) of that file, which equates the estimation term and the
#' bias term of the oracle inequality of E1.5 and yields the rate
#' \eqn{(\log n/n)^{2s'/(2s'+1)}}. The penalty is the one of Corollary 2 of
#' E1.5,
#' \deqn{\lambda_n = 2\,\sigma\,\hat\sigma_{\max}
#'        \sqrt{2\log(2d/\alpha)/n},}
#' twice the level \eqn{\lambda_0} that dominates
#' \eqn{\|\tilde B'\varepsilon/n\|_\infty} with probability
#' \eqn{1-\alpha}, where \eqn{d} is the number of penalized columns and
#' \eqn{\hat\sigma_{\max} = \max_a (\hat\Sigma_{aa})^{1/2}} is the largest
#' empirical column norm of the penalized part. That reading of
#' \eqn{\|Z\|_{\max}} is the one E1.5 calibrated: the literal
#' \eqn{\max_{i,a}|Z_{ia}|} is of order \eqn{2^{J/2}} and would destroy the
#' rate of E1.6.
#'
#' Neither quantity is estimated from the data: \eqn{s'} is the effective
#' regularity \eqn{s - (1/\pi - 1/2)_+} of the components, an assumption,
#' and \eqn{\sigma} is the error scale. What it costs to follow the rule
#' rather than to search is measured in step E2.3.
#'
#' @param n Sample size.
#' @param s The effective regularity \eqn{s'} of the components,
#'   \eqn{s - (1/\pi - 1/2)_+} in the notation of E1.3. For reference, a
#'   periodic function with a corner has \eqn{s' = 3/2} and a piecewise
#'   smooth function with a jump inside the interval has \eqn{s' = 1/2}.
#' @param design An object of class \code{"wafc_design"}.
#' @param sigma The standard deviation of the error. \code{NULL} uses
#'   \code{\link{wafc_sigma}}.
#' @param alpha The confidence level of the event \eqn{\mathcal{T}} of E1.5.
#' @param y The response, needed only when \code{sigma} is \code{NULL}.
#'
#' @return \code{wafc_J_theory} returns an integer, \code{wafc_lambda_theory}
#'   a single penalty level on the scale of the objective of
#'   \code{\link{wafc}}.
#'
#' @examples
#' wafc_J_theory(c(250, 1000, 12800), s = 1/4)
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' des <- wafc_design(d$x, d$u, J = wafc_J_theory(300, s = 3/2))
#' wafc_lambda_theory(des, sigma = d$sigma)
#'
#' @export
wafc_J_theory <- function(n, s = 1) {
  if (!is.numeric(n) || length(n) == 0L || any(!is.finite(n)) || any(n < 3)) {
    stop("'n' must be numeric and at least 3.", call. = FALSE)
  }
  if (length(s) != 1L || !is.finite(s) || s <= 0) {
    stop("'s' must be a single positive value.", call. = FALSE)
  }
  ## min{J : 2^J >= t} = ceiling(log2(t)), and J > j0 = 0 always
  pmax(1L, as.integer(ceiling(log2((n / log(n))^(1 / (2 * s + 1))))))
}

#' @rdname wafc_J_theory
#' @export
wafc_lambda_theory <- function(design, sigma = NULL, alpha = 0.05, y = NULL) {
  if (!inherits(design, "wafc_design")) {
    stop("'design' must be an object returned by wafc_design().", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single value in (0, 1).", call. = FALSE)
  }
  if (is.null(sigma)) {
    if (is.null(y)) {
      stop("Supply 'sigma', or 'y' so that it can be estimated by wafc_sigma().",
           call. = FALSE)
    }
    sigma <- wafc_sigma(design, y, alpha = alpha)[["sigma"]]
  }
  if (length(sigma) != 1L || !is.finite(sigma) || sigma < 0) {
    stop("'sigma' must be a single non-negative value.", call. = FALSE)
  }
  n <- design[["n"]]
  pen <- seq_len(design[["nvars"]])[-design[["unpenalized"]]]
  d <- length(pen)
  Z <- design[["Z"]][, pen, drop = FALSE]
  smax <- sqrt(max(Matrix::colSums(Z^2) / n))
  2 * sigma * smax * sqrt(2 * log(2 * d / alpha) / n)
}

#' Plug-in estimate of the error scale
#'
#' The penalty level of \code{\link{wafc_lambda_theory}} needs
#' \eqn{\sigma}, which the theory of E1.5 and E1.6 treats as known. This
#' function supplies it when it is not, by either of two devices, and
#' nothing is claimed here about the fit that follows: the rates of E1.6 are
#' stated for \eqn{\sigma} known, and a random penalty level is listed in
#' "what this does not cover" of both E1.5 and E1.6.
#'
#' \code{method = "cv"}, the default, cross-validates the penalty level at
#' the given design and returns
#' \eqn{\hat\sigma^2 = RSS/(n - df)} at \code{lambda.min}, with the
#' degrees of freedom of \code{\link{wafc_bic}}. It uses the data twice,
#' but only to fix a scale.
#'
#' \code{method = "fixed.point"} iterates
#' \eqn{\sigma \mapsto \hat\sigma(\lambda_n(\sigma))}: fit at the
#' penalty level the current \eqn{\sigma} prescribes, re-estimate
#' \eqn{\sigma} from the residuals, repeat. The map is increasing, so it
#' can have several fixed points, and in the regime measured in step E2.3 it
#' has a bad one: \eqn{\lambda_n} is conservative enough at
#' \eqn{n \le 1000} to leave the wavelet part nearly empty, the residual
#' scale is then the standard deviation of the varying part of the
#' regression function, and that value reproduces itself. The iteration
#' therefore starts from below, at the residual scale of the least penalized
#' fit whose degrees of freedom do not exceed \eqn{n/2}, and a warning is
#' issued when it lands on a fit with no non-zero wavelet coefficient.
#'
#' @param design An object of class \code{"wafc_design"}.
#' @param y The response.
#' @param method \code{"cv"} or \code{"fixed.point"}; see above.
#' @param penalty \code{"lasso"} or \code{"sglasso"}, as in
#'   \code{\link{wafc}}.
#' @param alpha The confidence level of the penalty level of E1.5.
#' @param nfolds,foldid Folds used by \code{method = "cv"}.
#' @param sigma0 Starting value of the iteration; \code{NULL} is the
#'   default described above.
#' @param maxit,tol Maximum number of iterations and relative tolerance of
#'   the fixed point.
#' @param ... Passed to \code{\link{wafc}}.
#'
#' @return A list with the estimate \code{sigma}, the penalty level
#'   \code{lambda} it prescribes, the \code{method}, the number
#'   \code{iter} of iterations, the sequence \code{trace} of estimates and
#'   the flag \code{converged}.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' des <- wafc_design(d$x, d$u, J = 3)
#' c(wafc_sigma(des, d$y)$sigma, d$sigma)
#'
#' @export
wafc_sigma <- function(design, y, method = c("cv", "fixed.point"),
                       penalty = c("lasso", "sglasso"), alpha = 0.05,
                       nfolds = 10L, foldid = NULL, sigma0 = NULL,
                       maxit = 20L, tol = 1e-4, ...) {
  if (!inherits(design, "wafc_design")) {
    stop("'design' must be an object returned by wafc_design().", call. = FALSE)
  }
  method <- match.arg(method)
  penalty <- match.arg(penalty)
  n <- design[["n"]]
  y <- wafc_check_y(y, n)

  if (method == "cv") {
    foldid <- wafc_foldid(foldid, n, nfolds)
    full <- wafc(design = design, y = y, penalty = penalty, ...)
    z <- wafc_cv_design(design, y, full, foldid, function(e) e^2, penalty, ...)
    ic <- wafc_ic(full, s = z[["lambda.min"]], gamma = 0)
    sigma <- sqrt(ic[["rss"]][1L] / max(n - ic[["df"]][1L], 1))
    return(list(sigma = sigma,
                lambda = wafc_lambda_theory(design, sigma = sigma,
                                            alpha = alpha),
                method = method, iter = 1L, trace = sigma, converged = TRUE))
  }

  if (is.null(sigma0)) {
    fit0 <- wafc(design = design, y = y, penalty = penalty, ...)
    ic0 <- wafc_ic(fit0, gamma = 0)
    ok <- which(ic0[["df"]] <= max(1, floor(n / 2)))
    k <- if (length(ok) > 0L) max(ok) else 1L
    sigma <- sqrt(ic0[["rss"]][k] / max(n - ic0[["df"]][k], 1))
  } else {
    sigma <- sigma0
  }
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    stop("'sigma0' must be a single positive value.", call. = FALSE)
  }
  trace <- sigma
  converged <- FALSE
  nz <- NA_integer_
  it <- 0L
  for (it in seq_len(maxit)) {
    lam <- wafc_lambda_theory(design, sigma = sigma, alpha = alpha)
    fit <- wafc(design = design, y = y, penalty = penalty,
                lambda = wafc_path_to(lam, design, y), ...)
    ic <- wafc_ic(fit, s = lam, gamma = 0)
    nz <- ic[["nzero"]][1L]
    new <- sqrt(ic[["rss"]][1L] / max(n - ic[["df"]][1L], 1))
    trace <- c(trace, new)
    if (abs(new - sigma) <= tol * sigma) {
      sigma <- new
      converged <- TRUE
      break
    }
    sigma <- new
  }
  if (isTRUE(nz == 0L)) {
    warning("The fixed point of wafc_sigma() stopped at a fit with no ",
            "non-zero wavelet coefficient, so the estimate is the scale of ",
            "the whole varying part and not of the error. Supply 'sigma' or ",
            "use method = \"cv\".", call. = FALSE)
  }
  list(sigma = sigma,
       lambda = wafc_lambda_theory(design, sigma = sigma, alpha = alpha),
       method = method, iter = it, trace = trace, converged = converged)
}

#' Select the pair (J, lambda) of a WAFC fit by one rule
#'
#' One entry point for the four selection rules of step E2.3, so that they
#' are applied to the same data through the same interface and their cost
#' can be compared. \code{"cv.min"} and \code{"cv.1se"} call
#' \code{\link{cv.wafc}}; \code{"bic"} and \code{"ebic"} minimise the
#' criterion of \code{\link{wafc_bic}} over the same grid of \eqn{J} and
#' over the path of each candidate; \code{"theory"} takes the pair of
#' \code{\link{wafc_J_theory}} and \code{\link{wafc_lambda_theory}} without
#' looking at any loss.
#'
#' @param x,u,y The data, as in \code{\link{wafc}}.
#' @param rule The selection rule.
#' @param J The grid of candidate resolution levels, as in
#'   \code{\link{cv.wafc}}; ignored by \code{rule = "theory"}, which
#'   computes its own.
#' @param penalty \code{"lasso"} or \code{"sglasso"}.
#' @param nfolds,foldid Folds of the cross-validation rules.
#' @param nlambda,lambda.min.ratio The paths of the candidates.
#' @param gamma The parameter of \code{\link{wafc_ebic}}.
#' @param s The effective regularity used by \code{rule = "theory"}.
#' @param sigma The error scale used by \code{rule = "theory"};
#'   \code{NULL} estimates it with \code{\link{wafc_sigma}}.
#' @param alpha The confidence level of the penalty level of E1.5.
#' @param ... Passed to \code{\link{wafc}} and to
#'   \code{\link{wafc_design}}.
#'
#' @return An object of class \code{"wafc_tune"}: a list with the
#'   \code{rule}, the selected \code{J} and \code{lambda}, the fit
#'   \code{fit} on the whole sample at that \eqn{J} (whose path contains the
#'   selected \eqn{\lambda}), the number \code{nzero} of non-zero wavelet
#'   coefficients there, the table \code{tab} of the rule over the grid, the
#'   \code{sigma} used by \code{rule = "theory"} and, for the
#'   cross-validation rules, the whole \code{"cv.wafc"} object in \code{cv}.
#'
#' @seealso \code{\link{cv.wafc}}, \code{\link{wafc_bic}},
#'   \code{\link{wafc_J_theory}}.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' t1 <- wafc_tune(d$x, d$u, d$y, rule = "ebic", J = 2:4)
#' t1
#' t2 <- wafc_tune(d$x, d$u, d$y, rule = "theory", s = 3/2, sigma = d$sigma)
#' c(t2$J, t2$lambda)
#'
#' @export
wafc_tune <- function(x, u, y,
                      rule = c("cv.min", "cv.1se", "bic", "ebic", "theory"),
                      J = NULL, penalty = c("lasso", "sglasso"),
                      nfolds = 10L, foldid = NULL, nlambda = 100L,
                      lambda.min.ratio = NULL, gamma = 1, s = 1,
                      sigma = NULL, alpha = 0.05, ...) {

  this_call <- match.call()
  rule <- match.arg(rule)
  penalty <- match.arg(penalty)
  x <- wafc_as_matrix(x, "x")
  u <- wafc_as_matrix(u, "u")
  n <- nrow(x)
  y <- wafc_check_y(y, n)

  if (rule == "theory") {
    Jn <- wafc_J_theory(n, s = s)
    if (Jn < 2L && is.null(list(...)[["eps"]])) {
      stop("The rule of the theory gives J_n = ", Jn, " at n = ", n,
           " with s' = ", s, ", and wafc_design() cannot build a design at ",
           "J = 1 with its default rescaling constant eps = 1.9^{-J} = ",
           signif(1.9^(-Jn), 3), ", which is outside [0, 0.5). Supply 'eps' ",
           "or use a larger sample.", call. = FALSE)
    }
    design <- wafc_design(x, u, J = Jn, ...)
    if (is.null(sigma)) {
      est <- wafc_sigma(design, y, alpha = alpha, nfolds = nfolds,
                        foldid = foldid, penalty = penalty, ...)
      sigma <- est[["sigma"]]
    }
    lam <- wafc_lambda_theory(design, sigma = sigma, alpha = alpha)
    fit <- wafc(design = design, y = y, penalty = penalty,
                lambda = wafc_path_to(lam, design, y), ...)
    tab <- data.frame(J = Jn, lambda = lam, sigma = sigma,
                      nzero = fit[["nzero"]][length(fit[["lambda"]])])
    out <- list(rule = rule, J = Jn, lambda = lam, fit = fit,
                nzero = tab[["nzero"]], tab = tab, sigma = sigma, cv = NULL,
                call = this_call)
    class(out) <- "wafc_tune"
    return(out)
  }

  if (rule %in% c("cv.min", "cv.1se")) {
    cv <- cv.wafc(x, u, y, J = J, penalty = penalty, nfolds = nfolds,
                  foldid = foldid, nlambda = nlambda,
                  lambda.min.ratio = lambda.min.ratio, ...)
    lam <- if (rule == "cv.min") cv[["lambda.min"]] else cv[["lambda.1se"]]
    fit <- cv[["wafc.fit"]]
    nz <- wafc_nzero_at(fit, lam)
    out <- list(rule = rule, J = cv[["J.min"]], lambda = lam, fit = fit,
                nzero = nz, tab = cv[["cvtab"]], sigma = NA_real_, cv = cv,
                call = this_call)
    class(out) <- "wafc_tune"
    return(out)
  }

  ## bic and ebic: the same grid of J, the criterion minimised over the
  ## path of each candidate and then over the grid.
  J <- wafc_J_grid(J, n)
  name <- rule
  if (length(gamma) != 1L || !is.finite(gamma) || gamma < 0 || gamma > 1) {
    stop("'gamma' must be a single value in [0, 1].", call. = FALSE)
  }
  gam <- if (rule == "bic") 0 else gamma
  best <- NULL
  rows <- vector("list", length(J))
  for (i in seq_along(J)) {
    fit <- wafc(x, u, y, J = J[i], penalty = penalty, nlambda = nlambda,
                lambda.min.ratio = lambda.min.ratio, ...)
    ic <- wafc_ic(fit, s = NULL, gamma = gam, name = name)
    k <- attr(ic, "which.min")
    rows[[i]] <- data.frame(J = J[i], lambda = ic[["lambda"]][k],
                            nzero = ic[["nzero"]][k], df = ic[["df"]][k],
                            value = ic[[name]][k])
    if (is.null(best) || rows[[i]][["value"]] < best[["value"]]) {
      best <- list(value = rows[[i]][["value"]], J = J[i],
                   lambda = ic[["lambda"]][k], fit = fit,
                   nzero = ic[["nzero"]][k])
    }
  }
  tab <- do.call(rbind, rows)
  names(tab)[names(tab) == "value"] <- name
  out <- list(rule = rule, J = best[["J"]], lambda = best[["lambda"]],
              fit = best[["fit"]], nzero = best[["nzero"]], tab = tab,
              sigma = NA_real_, cv = NULL, call = this_call)
  class(out) <- "wafc_tune"
  out
}

#' @rdname wafc_tune
#' @param object,digits An object of class \code{"wafc_tune"} and the number
#'   of significant digits printed.
#' @export
print.wafc_tune <- function(x, digits = max(3L, getOption("digits") - 3L),
                            ...) {
  cat(sprintf("WAFC tuning by rule \"%s\": J = %d, lambda = %s, %d non-zero wavelet coefficient(s) of %d\n",
              x[["rule"]], x[["J"]], format(x[["lambda"]], digits = digits),
              x[["nzero"]], x[["fit"]][["npen"]]))
  if (!is.na(x[["sigma"]])) {
    cat(sprintf("  sigma = %s\n", format(x[["sigma"]], digits = digits)))
  }
  print(x[["tab"]], row.names = FALSE, digits = digits)
  invisible(x)
}

#' @rdname wafc_tune
#' @export
coef.wafc_tune <- function(object, ...) {
  coef(object[["fit"]], s = object[["lambda"]], ...)
}

#' @rdname wafc_tune
#' @param newx,newu New covariates, as in \code{\link{predict.wafc}}.
#' @export
predict.wafc_tune <- function(object, newx, newu, ...) {
  if (missing(newx) && missing(newu)) {
    return(predict(object[["fit"]], s = object[["lambda"]], ...))
  }
  if (missing(newx)) {
    return(predict(object[["fit"]], newu = newu, s = object[["lambda"]], ...))
  }
  predict(object[["fit"]], newx, newu, s = object[["lambda"]], ...)
}

## ---------------------------------------------------------------------------
## Internals
## ---------------------------------------------------------------------------

## Grid of candidate resolution levels: the rule of cv.wall, with j0 = 0,
## starting at J = 2. It does not start at J = 1 for two reasons: the
## default rescaling constant of wafc_design() is eps = 1.9^{-J}, which at
## J = 1 is 0.526 and outside the admissible [0, 0.5), so the design cannot
## be built with the defaults; and a block of a single wavelet column is not
## a sieve. A user who wants J = 1 has to supply 'eps' as well.
wafc_J_grid <- function(J, n) {
  if (is.null(J)) return(seq(2L, max(2L, ceiling(log2(n) / 2))))
  if (!is.numeric(J) || length(J) == 0L || any(!is.finite(J)) ||
      any(J != round(J))) {
    stop("'J' must contain only integer values.", call. = FALSE)
  }
  if (any(J < 1)) {
    stop("'J' must be larger than j0 (= 0).", call. = FALSE)
  }
  sort(unique(as.integer(J)))
}

wafc_foldid <- function(foldid, n, nfolds) {
  if (is.null(foldid)) {
    if (length(nfolds) != 1L || !is.finite(nfolds) || nfolds < 3 ||
        nfolds != round(nfolds)) {
      stop("'nfolds' must be a single integer of at least 3.", call. = FALSE)
    }
    if (nfolds > n) {
      stop("'nfolds' cannot exceed the sample size (", n, ").", call. = FALSE)
    }
    return(sample(rep_len(seq_len(as.integer(nfolds)), n)))
  }
  if (length(foldid) != n) {
    stop("'foldid' must have one entry per observation (", n, ").",
         call. = FALSE)
  }
  foldid <- as.integer(foldid)
  if (anyNA(foldid) || min(foldid) < 1L) {
    stop("'foldid' must contain positive integers.", call. = FALSE)
  }
  if (max(foldid) < 3L) stop("'foldid' must define at least 3 folds.",
                             call. = FALSE)
  if (!all(seq_len(max(foldid)) %in% foldid)) {
    stop("'foldid' must use every fold from 1 to ", max(foldid), ".",
         call. = FALSE)
  }
  foldid
}

## The design restricted to a set of rows, used by the folds. Everything
## that describes the basis and the rescaling is kept, which is the point:
## the columns of a fold are the columns of the whole sample, so the fits
## are comparable and no basis is evaluated twice. 'constant' is kept as
## well, so that the carrier of the intercept cannot change from fold to
## fold (it is the constant covariate of the whole sample).
wafc_subset_design <- function(design, rows) {
  design[["Z"]] <- design[["Z"]][rows, , drop = FALSE]
  design[["n"]] <- nrow(design[["Z"]])
  design
}

## The cross-validation of one design: the fold fits at the penalty path of
## the fit on the whole sample, and the held-out loss averaged over folds.
## Shared by cv.wafc(), which calls it once per candidate J, and by
## wafc_sigma(), which calls it once.
wafc_cv_design <- function(design, y, full, foldid, loss, penalty, ...) {
  n <- design[["n"]]
  nfolds <- max(foldid)
  lam <- full[["lambda"]]
  nl <- length(lam)
  err <- matrix(NA_real_, n, nl)
  for (k in seq_len(nfolds)) {
    out <- which(foldid == k)
    fit <- wafc(design = wafc_subset_design(design, -out), y = y[-out],
                penalty = penalty, lambda = lam, ...)
    cf <- wafc_raw_coef(fit, s = lam)
    eta <- as.matrix(design[["Z"]][out, , drop = FALSE] %*%
                       cf[-1L, , drop = FALSE])
    eta <- sweep(eta, 2L, cf[1L, ], "+")
    err[out, ] <- loss(y[out] - eta)
  }
  ## Mean of the fold means, and the standard error of that mean: the folds
  ## are the replicates, which is the reading that makes cvsd the standard
  ## error of a mean of nfolds numbers.
  fold_mean <- matrix(NA_real_, nfolds, nl)
  for (k in seq_len(nfolds)) {
    fold_mean[k, ] <- colMeans(err[foldid == k, , drop = FALSE])
  }
  cvm <- colMeans(fold_mean)
  cvsd <- apply(fold_mean, 2L, stats::sd) / sqrt(nfolds)
  imin <- which.min(cvm)
  i1se <- min(which(cvm <= cvm[imin] + cvsd[imin]))
  list(lambda = lam, cvm = cvm, cvsd = cvsd, cvup = cvm + cvsd,
       cvlo = cvm - cvsd, nzero = full[["nzero"]],
       lambda.min = lam[imin], lambda.1se = lam[i1se],
       cvm.min = cvm[imin], cvsd.min = cvsd[imin],
       nzero.min = full[["nzero"]][imin])
}

wafc_cv_s <- function(object, s) {
  if (is.character(s)) {
    s <- match.arg(s[1L], c("lambda.min", "lambda.1se"))
    return(object[[s]])
  }
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s < 0) {
    stop("'s' must be \"lambda.min\", \"lambda.1se\" or a single ",
         "non-negative value.", call. = FALSE)
  }
  as.numeric(s)
}

## Residual sum of squares, degrees of freedom and criterion along a path.
wafc_ic <- function(object, s = NULL, gamma = 0, name = "bic") {
  if (!inherits(object, "wafc")) {
    stop("'object' must be an object returned by wafc().", call. = FALSE)
  }
  design <- object[["design"]]
  n <- object[["n"]]
  y <- object[["y"]]
  cf <- wafc_raw_coef(object, s = s)
  lam <- if (is.null(s)) object[["lambda"]] else as.numeric(s)
  b <- cf[-1L, , drop = FALSE]
  eta <- as.matrix(design[["Z"]] %*% b)
  eta <- sweep(eta, 2L, cf[1L, ], "+")
  rss <- colSums((y - eta)^2)
  pen <- seq_len(object[["nvars"]])[-design[["unpenalized"]]]
  npen <- object[["npen"]]
  nz <- as.integer(colSums(b[pen, , drop = FALSE] != 0))
  ## the p level terms are never penalized and always count
  df <- nz + length(design[["unpenalized"]])
  value <- n * log(pmax(rss, .Machine[["double.eps"]]) / n) + df * log(n)
  if (gamma > 0) value <- value + 2 * gamma * lchoose(npen, nz)
  out <- data.frame(lambda = lam, nzero = nz, df = df, rss = rss,
                    value = value)
  names(out)[5L] <- name
  k <- which.min(value)
  attr(out, "which.min") <- k
  attr(out, "lambda.min") <- lam[k]
  out
}

## The entry point of the path: the smallest penalty level at which every
## penalized coefficient is zero, which is max_a |B_a' r_0| / n with r_0 the
## residual of the least squares fit on the unpenalized columns alone. On
## the scale of the objective, this is the lambda at which glmnet starts.
wafc_lambda_max <- function(design, y) {
  n <- design[["n"]]
  unp <- design[["unpenalized"]]
  pen <- seq_len(design[["nvars"]])[-unp]
  W <- cbind(`(Intercept)` = 1, as.matrix(design[["Z"]][, unp, drop = FALSE]))
  r <- qr.resid(qr(W), y)
  max(abs(as.numeric(Matrix::crossprod(design[["Z"]][, pen, drop = FALSE], r)))) / n
}

## A decreasing path ending exactly at 'lambda', so that the engine reaches
## it by warm starts instead of being asked for a single penalty level. It
## starts at the entry point of the path when that is above 'lambda', which
## is the sequence the engine would have built by itself.
wafc_path_to <- function(lambda, design, y, length.out = 50L) {
  if (length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
    stop("'lambda' must be a single positive value.", call. = FALSE)
  }
  top <- max(wafc_lambda_max(design, y), 1.1 * lambda)
  exp(seq(log(top), log(lambda), length.out = length.out))
}

## Number of non-zero wavelet coefficients of a fit at one penalty level,
## read from the path rather than interpolated: interpolation of the
## coefficients would count a coefficient that is zero at both neighbours
## as zero, but one that changes sign as non-zero.
wafc_nzero_at <- function(object, lambda) {
  k <- which.min(abs(object[["lambda"]] - lambda))
  object[["nzero"]][k]
}
