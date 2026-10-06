## wafc/R/tune.R -- choosing the pair (J, lambda) of the WAFC fit (step
## E2.3): cross-validation over the grid, the information criteria, and the
## resolution rule that the theory of E1.6 prescribes.
##
## Five selection rules live here, and they answer different questions.
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
##  (iv)  wafc_lambda_qut() is the quantile universal threshold of
##        Giacobino et al. (2017), added by step E2.4 and moved here by
##        step E2.4b, which is where a penalty rule belongs. It is the only
##        rule that needs neither sigma, nor s', nor folds, and the pilot
##        measured it as the only one that never turns on a false block.
##        The resolution is not part of it, so wafc_tune(rule = "qut")
##        takes J from the cross-validation and changes only lambda, which
##        is what the pilot did.
##  (v)   wafc_gcv() scores a fitted path by generalized cross-validation,
##        GCV = n RSS / (n - df)^2, with the degrees of freedom of
##        wafc_bic(), which Zou, Hastie and Tibshirani (2007) justify for
##        the LASSO; wafc_tune(rule = "gcv") minimises it over the same grid
##        of J, so that the WAFC can be compared with a spline whose basis
##        dimension is chosen by the same criterion (step E2.5h, after
##        Ruppert, 2002). Points with df >= n/2 are left out, a guard
##        declared before measuring: the denominator vanishes at df = n,
##        and near it the criterion rewards an interpolating fit.
##
## Two things this file inherits from E2.2 and does not renegotiate
## (decision D17): the lambda of every object is the lambda of the
## objective of E1.5, not the one of the engine, and wafc_kkt() is what
## verifies it. Every lambda that goes in or comes out of the functions
## below is therefore on the scale of the objective.

#' Cross-validation of a WAFC fit over the resolution level and the penalty,
#' followed by a threshold on the blocks
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
#' With \code{penalty = "block"}, the default since step E3.1 (decision
#' D44), the cross-validation of each candidate \eqn{J} is the one of the
#' engine, \code{\link[grpreg]{cv.grpreg}} on the same folds, which is how
#' \code{\link{wafc_fit_klopp}} chose \eqn{(J, \lambda)} in every
#' measurement of step E2.5: the fold fits use the standardization of the
#' whole sample, the path is the one of the fit on the whole sample cut at
#' the last penalty level every fold reached, \code{cvm} is the loss
#' averaged over the \eqn{n} held-out observations (the mean of the fold
#' means of the other two penalties when the folds have equal sizes) and
#' \code{cvsd} its standard error over the observations,
#' \code{sd/sqrt(n)}, where the other two penalties take the standard
#' error over the folds. The pair selected is therefore the one of
#' \code{wafc_fit_klopp(penalize.levels = FALSE, balanced = TRUE)} on the
#' same folds, coefficient for coefficient.
#'
#' The fit at the selected pair is then thresholded on the blocks
#' \eqn{(\ell, m)} (decision D45): \code{\link{wafc_threshold}} with the
#' rule \code{threshold}, on the same folds, without refit. The blocks
#' whose coefficient norm is not above the chosen \eqn{t} are set to zero,
#' and the result is the WAFC estimator: \code{\link[=coef.cv.wafc]{coef}},
#' \code{\link[=predict.cv.wafc]{predict}}, \code{\link{wafc_functions}}
#' and \code{\link{wafc_blocks}} read it at \code{s = "lambda.min"}. The
#' path itself is untouched in \code{wafc.fit}, and every other penalty
#' level is read there, without threshold.
#'
#' @param x,u,y The data, as in \code{\link{wafc}}. For \code{print},
#'   \code{x} is an object of class \code{"cv.wafc"}.
#' @param J The grid of candidate resolution levels. \code{NULL} (the
#'   default) uses \code{2:8} (decision D34). The rule of \code{cv.wall},
#'   \code{2:ceiling(log2(n)/2)}, used before, truncated the expansion
#'   exactly where the choice was being made: in the pilot of step E2.4 its
#'   top was the oracle of the grid in every replicate with a component to
#'   resolve, and going to 8 cut the prediction error by 17.6\% and the
#'   error of the components by 32\% in the inhomogeneous scenario at
#'   \eqn{n = 1000}. The grid starts at 2 because step E2.4 measured that
#'   \eqn{J = 1} is never selected; passing \code{J = 1:8} is allowed.
#'   Each candidate is the level asked of \code{\link{wafc_design}}, which
#'   caps the level of a modulating covariate at what its distinct values
#'   identify (\code{cap.J}, step E3.4): two candidates whose designs have
#'   the same levels are the same design, and the second one takes the
#'   cross-validation of the first instead of refitting it, so the tie goes
#'   to the smaller \eqn{J}, as every tie does.
#' @param penalty \code{"block"} (the default), \code{"lasso"} or
#'   \code{"sglasso"}, as in \code{\link{wafc}}.
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
#' @param threshold The rule that chooses the threshold \eqn{t} of the
#'   blocks at the selected pair: \code{"cv1se"} (the default, decision
#'   D45: the largest \eqn{t} whose cross-validated error is within one
#'   standard error between folds of the smallest), \code{"cv"} (the
#'   \eqn{t} of the smallest cross-validated error) or \code{"max"}
#'   (\eqn{t = 0.15 \max_{\ell m} \hat\nu_{\ell m}}); see
#'   \code{\link{wafc_threshold}}. \code{"none"} returns the fit without
#'   threshold, which is what \code{cv.wafc()} returned before step E3.1.
#'   It comes after \code{...}, so it is matched only by its full name, and
#'   the \code{thresh} of \code{\link{wafc}} cannot be taken for it.
#'
#' @return An object of class \code{"cv.wafc"}: a list with the grid
#'   \code{J}, the list \code{cv} of per-candidate results (\code{lambda},
#'   \code{cvm}, \code{cvsd}, \code{cvup}, \code{cvlo}, \code{nzero},
#'   \code{lambda.min}, \code{lambda.1se} and the values attained at
#'   \code{lambda.min}), the summary \code{cvtab}, the selected
#'   \code{J.min}, \code{lambda.min} and \code{lambda.1se}, the
#'   matrix \code{J.eff} of the levels each candidate was built at (one row
#'   per candidate, one column per modulating covariate; equal to \code{J}
#'   wherever the cap of \code{\link{wafc_design}} did not act), the
#'   \code{foldid} used, the \code{penalty}, \code{wafc.fit}, the fit on
#'   the whole sample at \code{J.min}, which is the fit the selected pair
#'   refers to, and \code{threshold}, \code{NULL} with
#'   \code{threshold = "none"} and otherwise a list with the \code{rule},
#'   the threshold \code{t}, the fraction \code{c} of the largest norm,
#'   the \eqn{p \times q} matrices \code{norm} (before the threshold) and
#'   \code{kept}, the thresholded coefficients \code{a0} and \code{b} in
#'   the coordinates of the design, their number \code{nzero} of non-zero
#'   wavelet coefficients, the \code{candidates} the rule scored and the
#'   seconds \code{time} it took (the fold fits included). \code{print}
#'   returns \code{x} invisibly.
#'
#' @seealso \code{\link{wafc_bic}}, \code{\link{wafc_tune}},
#'   \code{\link{wafc_J_theory}}, \code{\link{wafc_threshold}},
#'   \code{\link{plot.cv.wafc}}.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' cvfit <- cv.wafc(d$x, d$u, d$y, J = 2:4, nfolds = 5)
#' cvfit
#' cvfit$threshold$kept
#' cf <- coef(cvfit)
#' cl <- cv.wafc(d$x, d$u, d$y, J = 2:4, nfolds = 5, penalty = "lasso",
#'               threshold = "none")
#' cf1 <- coef(cl, s = "lambda.1se")
#'
#' @export
cv.wafc <- function(x, u, y, J = NULL,
                    penalty = c("block", "lasso", "sglasso"),
                    nfolds = 10L, foldid = NULL, lambda = NULL,
                    nlambda = 100L, lambda.min.ratio = NULL,
                    type.measure = c("mse", "mae"), trace = FALSE, ...,
                    threshold = c("cv1se", "cv", "max", "none")) {

  this_call <- wafc_compact_call(match.call(), "cv.wafc")
  penalty <- match.arg(penalty)
  type.measure <- match.arg(type.measure)
  threshold <- match.arg(threshold)
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

  dots <- wafc_split_dots(list(...))

  run_one <- function(Ji) {
    design <- do.call(wafc_design,
                      c(list(x = x, u = u, J = Ji), dots[["design"]]))
    ## a candidate whose levels, after the cap of wafc_design(), are the
    ## ones of a candidate already run is the same design on the same folds:
    ## its cross-validation is that one, and refitting it would only repeat
    ## it (step E3.4)
    prev <- match(paste(design[["J"]], collapse = ","), J_seen)
    if (!is.na(prev)) {
      z <- cvlist[[prev]]
      z[["J"]] <- Ji
      ## the convergence table of the block LASSO names its J
      if (is.data.frame(z[["conv"]])) z[["conv"]][["J"]] <- Ji
      return(list(cv = z, fit = NULL, J.eff = design[["J"]], same = prev))
    }
    if (penalty == "block") {
      return(c(wafc_cv_block(design, y, foldid, lambda, nlambda,
                             lambda.min.ratio, type.measure, dots[["fit"]],
                             Ji),
               list(J.eff = design[["J"]], same = NA_integer_)))
    }
    full <- do.call(wafc, c(list(design = design, y = y, penalty = penalty,
                                 lambda = lambda, nlambda = nlambda,
                                 lambda.min.ratio = lambda.min.ratio),
                            dots[["fit"]]))
    z <- do.call(wafc_cv_design,
                 c(list(design, y, full, foldid, loss, penalty),
                   dots[["fit"]]))
    z[["J"]] <- Ji
    z[["conv"]] <- full[["conv"]]
    list(cv = z, fit = full, J.eff = design[["J"]], same = NA_integer_)
  }

  ## Only the fit of the best candidate so far is kept: every fit carries its
  ## design, and holding all of them until the end made the memory peak the
  ## sum over the grid instead of its largest term. The candidate kept is the
  ## one which.min() over the whole grid selects, since which.min() over the
  ## first i values points at i exactly when i is a new first minimum.
  cvlist <- vector("list", length(J))
  J_seen <- character(length(J))
  J.eff <- matrix(NA_integer_, length(J), ncol(u),
                  dimnames = list(J, wafc_names(u, ncol(u), "u")))
  best_fit <- NULL
  npen <- NA_integer_
  for (i in seq_along(J)) {
    run <- run_one(J[i])
    cvlist[[i]] <- run[["cv"]]
    J.eff[i, ] <- run[["J.eff"]]
    J_seen[i] <- paste(run[["J.eff"]], collapse = ",")
    ## a repeated design ties with its first occurrence, which which.min()
    ## keeps, so its fit (NULL here) is never the one kept
    cvm_sofar <- vapply(cvlist[seq_len(i)], `[[`, 0, "cvm.min")
    if (identical(which.min(cvm_sofar), i)) best_fit <- run[["fit"]]
    if (!is.null(run[["fit"]])) npen <- run[["fit"]][["npen"]]
    if (trace) {
      z <- run[["cv"]]
      cat(sprintf("J = %d: %s = %.5f at lambda = %.5g (%d nonzero of %d)%s\n",
                  J[i], type.measure, z[["cvm.min"]], z[["lambda.min"]],
                  z[["nzero.min"]], npen,
                  if (is.na(run[["same"]])) "" else
                    sprintf(", the design of J = %d", J[run[["same"]]])))
    }
    rm(run)
  }

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
              J.min = J[best], J.eff = J.eff, lambda.min = z[["lambda.min"]],
              lambda.1se = z[["lambda.1se"]], cvm.min = z[["cvm.min"]],
              cvsd.min = z[["cvsd.min"]], nzero.min = z[["nzero.min"]],
              type.measure = type.measure, nfolds = nfolds, foldid = foldid,
              penalty = penalty, wafc.fit = best_fit, threshold = NULL)
  class(out) <- "cv.wafc"

  ## The threshold of decision D45, on the folds of the cross-validation and
  ## at the pair it selected. The fold fits of wafc_threshold() refit the
  ## engine, so the controls of the engine the caller passed go with them
  ## (the ones the object records, the penalty, the chunks and the
  ## intercept, are set there).
  if (threshold != "none") {
    eng <- dots[["fit"]][names(dots[["fit"]]) %in%
                           c("thresh", "maxit", "standardize")]
    th <- do.call(wafc_threshold, c(list(out, rule = threshold), eng))
    out[["threshold"]] <- wafc_cv_thr_slot(out, th)
  }
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
  cap <- wafc_cv_cap_text(x)
  if (!is.null(cap)) cat(cap, "\n")
  th <- x[["threshold"]]
  if (is.null(th)) {
    cat("No threshold on the blocks.\n")
  } else {
    kept <- th[["kept"]]
    cat(sprintf("Threshold \"%s\" at lambda.min: t = %s, %d of %d block(s) kept",
                th[["rule"]], format(th[["t"]], digits = digits), sum(kept),
                length(kept)))
    if (any(kept)) {
      cat(":", paste(outer(rownames(kept), colnames(kept), paste,
                           sep = ":")[kept], collapse = ", "))
    }
    cat("\n")
  }
  invisible(x)
}

#' Coefficients and predictions at the cross-validated pair
#'
#' At \code{s = "lambda.min"}, the default, these read the thresholded fit
#' when the object has one (decision D45; \code{\link{cv.wafc}}, argument
#' \code{threshold}), which is the WAFC estimator; at any other penalty
#' level they read the path of \code{wafc.fit}, without threshold, since
#' the threshold was chosen at \code{lambda.min}.
#'
#' @param object An object of class \code{"cv.wafc"}.
#' @param s \code{"lambda.min"}, \code{"lambda.1se"}, or a numeric penalty
#'   level on the scale of the objective of \code{\link{wafc}}.
#' @param thresholded \code{NULL} (the default) reads the thresholded fit
#'   exactly when the object has one and \code{s} is \code{lambda.min};
#'   \code{FALSE} reads the path at \code{s} in every case; \code{TRUE}
#'   asks for the thresholded fit and is an error where there is none.
#' @param newx,newu New covariates, as in \code{\link{predict.wafc}}.
#' @param ... Passed to the method of the underlying \code{"wafc"} object
#'   (\code{type} of \code{\link{predict.wafc}}).
#'
#' @return As \code{\link{coef.wafc}} and \code{\link{predict.wafc}}, at the
#'   resolution level \code{J.min} selected by the cross-validation.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' cvfit <- cv.wafc(d$x, d$u, d$y, J = 2:3, nfolds = 5)
#' dim(coef(cvfit))
#' head(predict(cvfit, d$x, d$u))
#' head(predict(cvfit, d$x, d$u, thresholded = FALSE))
#' head(predict(cvfit, newu = d$u, type = "beta"))
#'
#' @export
coef.cv.wafc <- function(object, s = c("lambda.min", "lambda.1se"),
                         thresholded = NULL, ...) {
  f <- wafc_cv_fit(object, s, thresholded)
  coef(f[["fit"]], s = f[["s"]], ...)
}

#' @rdname coef.cv.wafc
#' @export
predict.cv.wafc <- function(object, newx, newu,
                            s = c("lambda.min", "lambda.1se"),
                            thresholded = NULL, ...) {
  f <- wafc_cv_fit(object, s, thresholded)
  fit <- f[["fit"]]
  s <- f[["s"]]
  if (missing(newx) && missing(newu)) {
    return(predict(fit, s = s, ...))
  }
  if (missing(newx)) return(predict(fit, newu = newu, s = s, ...))
  predict(fit, newx, newu, s = s, ...)
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

#' Generalized cross-validation of a WAFC fit
#'
#' The criterion of generalized cross-validation along the penalty path of
#' a fit,
#' \deqn{GCV = \frac{n\,RSS}{(n - df)^2},}
#' with the degrees of freedom of \code{\link{wafc_bic}}: the non-zero
#' wavelet coefficients plus the \eqn{p} level terms. For the LASSO the
#' number of non-zero coefficients is an unbiased estimate of the degrees of
#' freedom (Zou, Hastie and Tibshirani, 2007), which is what makes it the
#' trace that the criterion of generalized cross-validation needs; it is the count of
#' decision D19, and the one the rule for the sparse group LASSO borrows
#' without that justification. This is the criterion Ruppert (2002) uses to
#' choose the number of knots of a penalized spline, so minimising it over
#' \eqn{(J, \lambda)} (\code{wafc_tune(rule = "gcv")}) puts the WAFC and
#' a spline tuned by GCV on the same footing (step E2.5h).
#'
#' The guard is declared here, before any measurement: points of the path
#' with \eqn{df \ge} \code{guard}\eqn{\cdot n} are left out of the
#' minimisation. The denominator vanishes at \eqn{df = n}, and near it a
#' fit that interpolates the sample has a small criterion for the wrong
#' reason; the same \eqn{n/2} bounds the starting fit of
#' \code{\link{wafc_sigma}}. Points with \eqn{df \ge n} have no criterion
#' (\code{NA}) and are left out with or without the guard.
#'
#' @param object An object of class \code{"wafc"}.
#' @param s Penalty levels at which the criterion is wanted, on the scale of
#'   the objective of \code{\link{wafc}}. \code{NULL} (the default) is the
#'   whole path of the fit.
#' @param guard Fraction of \eqn{n} at or above which a point is left out.
#'   \code{1} keeps every point with a criterion.
#'
#' @return A data frame with one row per penalty level: \code{lambda},
#'   \code{nzero}, \code{df}, \code{rss}, \code{gcv} and the flag
#'   \code{excluded}. The index of the minimising row among the ones kept
#'   is the attribute \code{"which.min"}, and the penalty level attaining
#'   it \code{"lambda.min"}; the index of the minimising row among every
#'   row with a criterion is \code{"which.min.unguarded"}.
#'
#' @references Ruppert, D. (2002). Selecting the number of knots for penalized
#'   splines. \emph{Journal of Computational and Graphical Statistics}
#'   11(4), 735-757.
#'
#'   Zou, H., Hastie, T. and Tibshirani, R. (2007). On the degrees of
#'   freedom of the lasso. \emph{The Annals of Statistics} 35(5), 2173-2192.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' fit <- wafc(d$x, d$u, d$y, J = 3)
#' gc <- wafc_gcv(fit)
#' attr(gc, "lambda.min")
#'
#' @export
wafc_gcv <- function(object, s = NULL, guard = 0.5) {
  if (length(guard) != 1L || !is.finite(guard) || guard <= 0 || guard > 1) {
    stop("'guard' must be a single value in (0, 1].", call. = FALSE)
  }
  ic <- wafc_ic(object, s = s, gamma = 0)
  n <- object[["n"]]
  df <- ic[["df"]]
  gcv <- ifelse(df < n, n * ic[["rss"]] / (n - df)^2, NA_real_)
  excluded <- df >= guard * n
  out <- data.frame(lambda = ic[["lambda"]], nzero = ic[["nzero"]], df = df,
                    rss = ic[["rss"]], gcv = gcv, excluded = excluded)
  ## the first point of a path has df = p and is never left out, so both
  ## minima exist
  keep <- which(!excluded & !is.na(gcv))
  k <- keep[which.min(gcv[keep])]
  ku <- which.min(gcv)
  attr(out, "which.min") <- k
  attr(out, "lambda.min") <- out[["lambda"]][k]
  attr(out, "which.min.unguarded") <- ku
  out
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
#' One entry point for the seven selection rules, so that they are applied
#' to the same data through the same interface and their cost can be
#' compared. \code{"cv.min"} and \code{"cv.1se"} call
#' \code{\link{cv.wafc}}; \code{"bic"}, \code{"ebic"} and \code{"gcv"}
#' minimise the criterion of \code{\link{wafc_bic}} or of
#' \code{\link{wafc_gcv}} over the same grid of \eqn{J} and over the path
#' of each candidate, \code{"gcv"} with its guard (points with
#' \eqn{df \ge n/2} left out); \code{"theory"} takes the pair of
#' \code{\link{wafc_J_theory}} and \code{\link{wafc_lambda_theory}}
#' without looking at any loss; and \code{"qut"} takes \eqn{J} from the
#' same cross-validation the first two use and \eqn{\lambda} from
#' \code{\link{wafc_lambda_qut}}, since the quantile universal threshold
#' is a rule for the penalty level alone.
#'
#' @param x,u,y The data, as in \code{\link{wafc}}. For \code{print},
#'   \code{x} is an object of class \code{"wafc_tune"}.
#' @param rule The selection rule.
#' @param J The grid of candidate resolution levels, as in
#'   \code{\link{cv.wafc}}; ignored by \code{rule = "theory"}, which
#'   computes its own.
#' @param penalty \code{"lasso"} or \code{"sglasso"}.
#' @param nfolds,foldid Folds of the cross-validation rules.
#' @param nlambda,lambda.min.ratio The paths of the candidates.
#' @param gamma The parameter of \code{\link{wafc_ebic}}.
#' @param s The effective regularity used by \code{rule = "theory"}; the
#'   value a scenario declares is \code{\link{wafc_sprime}}.
#' @param sigma The error scale used by \code{rule = "theory"};
#'   \code{NULL} estimates it with \code{\link{wafc_sigma}}.
#' @param alpha The confidence level of the penalty level of E1.5, and of
#'   the quantile of \code{rule = "qut"}.
#' @param nsim Number of null samples simulated by \code{rule = "qut"}.
#' @param qut.seed Optional seed of that simulation, so that two calls on
#'   the same design give the same penalty level.
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
#'   For \code{rule = "gcv"}, \code{tab} also has, for each \eqn{J}, the
#'   number of points of the path the guard left out (\code{excluded}) and
#'   the minimum without the guard (\code{gcv.unguarded}), and the list
#'   \code{guard} says how many points of the whole grid were left out and
#'   whether the guard decided, that is, whether the pair chosen without it
#'   (among the points with \eqn{df < n}) would have been another one.
#'   \code{coef} and \code{predict} read \code{fit} at the selected
#'   \eqn{\lambda}, as \code{\link{coef.wafc}} and
#'   \code{\link{predict.wafc}} do; \code{print} returns \code{x}
#'   invisibly.
#'
#' @seealso \code{\link{cv.wafc}}, \code{\link{wafc_bic}},
#'   \code{\link{wafc_J_theory}}, \code{\link{wafc_lambda_qut}}.
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
                      rule = c("cv.min", "cv.1se", "bic", "ebic", "theory",
                               "qut", "gcv"),
                      J = NULL, penalty = c("lasso", "sglasso"),
                      nfolds = 10L, foldid = NULL, nlambda = 100L,
                      lambda.min.ratio = NULL, gamma = 1, s = 1,
                      sigma = NULL, alpha = 0.05, nsim = 200L,
                      qut.seed = NULL, ...) {

  this_call <- wafc_compact_call(match.call(), "wafc_tune")
  rule <- match.arg(rule)
  penalty <- match.arg(penalty)
  x <- wafc_as_matrix(x, "x")
  u <- wafc_as_matrix(u, "u")
  n <- nrow(x)
  y <- wafc_check_y(y, n)
  dots <- wafc_split_dots(list(...))

  if (rule == "theory") {
    Jn <- wafc_J_theory(n, s = s)
    ## J_n = 1 used to be unreachable, because the inherited margin
    ## eps = 1.9^{-J} was 0.526 there, outside the [0, 0.5) that
    ## wafc_rescale() requires. Since step E2.1b the margin is a fixed
    ## constant and the design exists at every J > j0, so the rule is
    ## followed wherever it leads and there is nothing to guard against.
    design <- do.call(wafc_design,
                      c(list(x = x, u = u, J = Jn), dots[["design"]]))
    if (is.null(sigma)) {
      ## 'maxit' is held back: it is the number of iterations of the fixed
      ## point of wafc_sigma() and the number of passes of the engine, and
      ## the two are not the same number.
      est <- do.call(wafc_sigma,
                     c(list(design, y, alpha = alpha, nfolds = nfolds,
                            foldid = foldid, penalty = penalty),
                       dots[["fit"]][names(dots[["fit"]]) != "maxit"]))
      sigma <- est[["sigma"]]
    }
    lam <- wafc_lambda_theory(design, sigma = sigma, alpha = alpha)
    fit <- do.call(wafc, c(list(design = design, y = y, penalty = penalty,
                                lambda = wafc_path_to(lam, design, y)),
                           dots[["fit"]]))
    tab <- data.frame(J = Jn, lambda = lam, sigma = sigma,
                      nzero = fit[["nzero"]][length(fit[["lambda"]])])
    out <- list(rule = rule, J = Jn, lambda = lam, fit = fit,
                nzero = tab[["nzero"]], tab = tab, sigma = sigma, cv = NULL,
                call = this_call)
    class(out) <- "wafc_tune"
    return(out)
  }

  if (rule %in% c("cv.min", "cv.1se", "qut")) {
    ## no threshold: the rules here choose (J, lambda), and the fit they
    ## return is read on the path (step E3.1 made the threshold the default
    ## of cv.wafc(), which this function did not ask for before)
    cv <- cv.wafc(x, u, y, J = J, penalty = penalty, nfolds = nfolds,
                  foldid = foldid, nlambda = nlambda,
                  lambda.min.ratio = lambda.min.ratio, ...,
                  threshold = "none")
    ## The quantile universal threshold says nothing about the
    ## resolution, so J is the cross-validated one and only lambda
    ## changes. The fit is refitted down to that lambda, because it need
    ## not belong to the path the cross-validation built.
    if (rule == "qut") {
      fit0 <- cv[["wafc.fit"]]
      lam <- wafc_lambda_qut(fit0[["design"]], y, alpha = alpha,
                             nsim = nsim, seed = qut.seed)
      fit <- do.call(wafc,
                     c(list(design = fit0[["design"]], y = y,
                            penalty = penalty,
                            lambda = wafc_path_to(lam, fit0[["design"]], y)),
                       dots[["fit"]]))
    } else {
      lam <- if (rule == "cv.min") cv[["lambda.min"]] else cv[["lambda.1se"]]
      fit <- cv[["wafc.fit"]]
    }
    nz <- wafc_nzero_at(fit, lam)
    out <- list(rule = rule, J = cv[["J.min"]], lambda = lam, fit = fit,
                nzero = nz, tab = cv[["cvtab"]], sigma = NA_real_, cv = cv,
                call = this_call)
    class(out) <- "wafc_tune"
    return(out)
  }

  ## gcv: the same grid of J, the criterion minimised over the points of
  ## each path the guard keeps, and then over the grid. The minimum without
  ## the guard is carried beside it, so that how often the guard decides is
  ## read from the object and not rerun.
  if (rule == "gcv") {
    J <- wafc_J_grid(J, n)
    best <- NULL
    free <- NULL
    rows <- vector("list", length(J))
    nexcl <- 0L
    npts <- 0L
    for (i in seq_along(J)) {
      fit <- do.call(wafc, c(list(x = x, u = u, y = y, J = J[i],
                                  penalty = penalty, nlambda = nlambda,
                                  lambda.min.ratio = lambda.min.ratio),
                             dots[["design"]], dots[["fit"]]))
      gc <- wafc_gcv(fit)
      k <- attr(gc, "which.min")
      ku <- attr(gc, "which.min.unguarded")
      nexcl <- nexcl + sum(gc[["excluded"]])
      npts <- npts + nrow(gc)
      rows[[i]] <- data.frame(J = J[i], lambda = gc[["lambda"]][k],
                              nzero = gc[["nzero"]][k], df = gc[["df"]][k],
                              gcv = gc[["gcv"]][k],
                              excluded = sum(gc[["excluded"]]),
                              points = nrow(gc),
                              gcv.unguarded = gc[["gcv"]][ku],
                              df.unguarded = gc[["df"]][ku])
      if (is.null(best) || rows[[i]][["gcv"]] < best[["value"]]) {
        best <- list(value = rows[[i]][["gcv"]], J = J[i],
                     lambda = gc[["lambda"]][k], fit = fit,
                     nzero = gc[["nzero"]][k])
      }
      if (is.null(free) || gc[["gcv"]][ku] < free[["value"]]) {
        free <- list(value = gc[["gcv"]][ku], J = J[i],
                     lambda = gc[["lambda"]][ku], df = gc[["df"]][ku])
      }
    }
    tab <- do.call(rbind, rows)
    guard <- list(excluded = nexcl, points = npts,
                  decides = !(free[["J"]] == best[["J"]] &&
                                free[["lambda"]] == best[["lambda"]]),
                  J.unguarded = free[["J"]],
                  lambda.unguarded = free[["lambda"]],
                  df.unguarded = free[["df"]])
    out <- list(rule = rule, J = best[["J"]], lambda = best[["lambda"]],
                fit = best[["fit"]], nzero = best[["nzero"]], tab = tab,
                sigma = NA_real_, cv = NULL, guard = guard, call = this_call)
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
    fit <- do.call(wafc, c(list(x = x, u = u, y = y, J = J[i],
                                penalty = penalty, nlambda = nlambda,
                                lambda.min.ratio = lambda.min.ratio),
                           dots[["design"]], dots[["fit"]]))
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
#' @param object An object of class \code{"wafc_tune"}.
#' @param digits Number of significant digits printed.
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

#' Quantile universal threshold for the WAFC
#'
#' The penalty rule of Giacobino, Sardy, Diaz-Rodriguez and Hengartner
#' (2017), the one Sardy and Ma (2024) use, transcribed to the objective of
#' \code{\link{wafc}}: the smallest \eqn{\lambda} that leaves every
#' penalized coefficient at zero is
#' \deqn{\lambda_0 = \|Z_{pen}' r_0\|_\infty / n,}
#' with \eqn{r_0} the residual of the least squares fit on the unpenalized
#' columns alone, and the rule takes the \eqn{1-\alpha} quantile of
#' \eqn{\lambda_0} under the null model \eqn{y = Z_{unp}c + \varepsilon}.
#'
#' It is computed here in the pivotal form, which is what makes the rule
#' free of \eqn{\sigma}: under the null, \eqn{r_0 = (I - P)\varepsilon}, so
#' the ratio \eqn{\|Z_{pen}'(I-P)\varepsilon\|_\infty /
#' \|(I-P)\varepsilon\|_2} does not depend on the error scale, and
#' multiplying its simulated quantile by the observed \eqn{\|r_0\|_2 / n}
#' gives a penalty level that needs neither \eqn{\sigma} nor
#' cross-validation. This is the only rule of the pilot with that property:
#' the rule of the theory needs \eqn{\sigma} and \eqn{s'}, the information
#' criteria need a count of degrees of freedom, and the cross-validation
#' needs the folds.
#'
#' @param design An object of class \code{"wafc_design"}.
#' @param y The response.
#' @param alpha One minus the level of the quantile. The default
#'   \eqn{0.05} is the one of Giacobino et al.
#' @param nsim Number of null samples simulated.
#' @param seed Optional seed, so that two calls on the same design give the
#'   same penalty level.
#'
#' @return A single penalty level, on the scale of the objective of
#'   \code{\link{wafc}}.
#'
#' @references Giacobino, C., Sardy, S., Diaz-Rodriguez, J. and Hengartner,
#'   N. (2017). Quantile universal threshold. \emph{Electronic Journal of
#'   Statistics} 11(2), 4701-4722.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' des <- wafc_design(d$x, d$u, J = 3)
#' wafc_lambda_qut(des, d$y, nsim = 50, seed = 1)
#'
#' @export
wafc_lambda_qut <- function(design, y, alpha = 0.05, nsim = 200L,
                            seed = NULL) {
  z <- wafc_qut_pivot(design, y, alpha = alpha, nsim = nsim, seed = seed)
  z[["quantile"]] * z[["r0norm"]] / design[["n"]]
}

#' Convergence of the engine along the cross-validation of a WAFC fit
#'
#' One row per candidate resolution level of a \code{\link{cv.wafc}} fit:
#' the length of the penalty path asked of the engine and returned by it,
#' the error code of the engine (\code{jerr} of \pkg{glmnet} and
#' \pkg{sparsegl}: \eqn{-k} when the \eqn{k}-th penalty level did not
#' converge, in which case the path stops at the one before), the smallest
#' penalty level of the path, and the same for the fold fits, which are
#' asked for the path of the whole sample. It reads what the fits recorded
#' and changes nothing (step E2.5j, open question 41(d) of
#' \file{docs/ESTADO.md}).
#'
#' A path of \pkg{glmnet} shorter than asked with \code{jerr = 0} is its
#' own early stop (the deviance stopped changing), not a failure. A fold
#' path cut above the \code{lambda.min} of its \eqn{J} is read there at its
#' last point, because \code{\link{wafc_raw_coef}} truncates \eqn{s} to the
#' range of the path.
#'
#' For \code{penalty = "block"} the folds are fitted inside
#' \code{\link[grpreg]{cv.grpreg}}, which does not report them one by one;
#' the table is then the one \code{\link{wafc_fit_klopp}} keeps in
#' \code{extra$conv}, read from the same objects: the path asked, returned
#' and cross-validated (\code{cv.grpreg} drops the penalty levels some
#' fold did not reach), the iterations against the budget of \pkg{grpreg},
#' which it counts over the whole path, and whether the selected level is
#' the last one of a path that was cut (\code{cut.at.min}).
#'
#' @param object An object of class \code{"cv.wafc"}.
#'
#' @return For \code{penalty = "lasso"} and \code{"sglasso"}, a data frame
#'   with \code{J}, \code{chosen}, \code{nlambda}
#'   (asked), \code{nreturned}, \code{jerr}, \code{lambda.last},
#'   \code{lambda.min} (of that \eqn{J}), \code{cut} (\code{jerr != 0}),
#'   \code{cut.at.min} (the path was cut and \code{lambda.min} is its last
#'   point, so the minimum may lie past the cut), \code{folds.cut} (folds
#'   whose path is shorter than asked or whose \code{jerr} is not zero),
#'   \code{fold.jerr} (the first nonzero code among them),
#'   \code{fold.lambda.cut} (the largest last point among the folds cut)
#'   and \code{fold.cut.above.min} (that point is above \code{lambda.min}).
#'   For \code{penalty = "block"}, the table of \code{\link{wafc_fit_klopp}}:
#'   \code{J}, \code{nlambda}, \code{nreturned}, the length \code{ncv} of
#'   the cross-validated path, \code{iter.total}, \code{max.iter},
#'   \code{n.maxiter}, \code{lambda.maxiter}, \code{lambda.last},
#'   \code{lambda.cv.last}, \code{lambda.min}, \code{chosen}, and the flags
#'   \code{budget}, \code{cv.cut} and \code{cut.at.min}.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' cv <- cv.wafc(d$x, d$u, d$y, J = 2:4, nfolds = 5, threshold = "none")
#' wafc_cv_convergence(cv)
#' cl <- cv.wafc(d$x, d$u, d$y, J = 2:4, nfolds = 5, penalty = "lasso",
#'               threshold = "none")
#' wafc_cv_convergence(cl)
#'
#' @export
wafc_cv_convergence <- function(object) {
  if (!inherits(object, "cv.wafc")) {
    stop("'object' must be a \"cv.wafc\" object.", call. = FALSE)
  }
  if (identical(object[["penalty"]], "block")) {
    return(wafc_kp_conv_table(lapply(object[["cv"]], `[[`, "conv"),
                              object[["J.min"]]))
  }
  rows <- lapply(object[["cv"]], function(z) {
    cv <- z[["conv"]]
    fo <- z[["conv.folds"]]
    if (is.null(cv)) {
      stop("This \"cv.wafc\" object does not record the convergence of ",
           "its fits; refit it.", call. = FALSE)
    }
    nl <- length(z[["lambda"]])
    fcut <- fo[["nreturned"]] < nl | fo[["jerr"]] != 0L
    fl <- if (any(fcut)) max(fo[["lambda.last"]][fcut]) else NA_real_
    fj <- fo[["jerr"]][fo[["jerr"]] != 0L]
    data.frame(J = z[["J"]], chosen = z[["J"]] == object[["J.min"]],
               nlambda = cv[["nlambda"]], nreturned = cv[["nreturned"]],
               jerr = cv[["jerr"]], lambda.last = min(z[["lambda"]]),
               lambda.min = z[["lambda.min"]], cut = cv[["jerr"]] != 0L,
               cut.at.min = cv[["jerr"]] != 0L &&
                 z[["lambda.min"]] == min(z[["lambda"]]),
               folds.cut = sum(fcut),
               fold.jerr = if (length(fj)) fj[1L] else 0L,
               fold.lambda.cut = fl,
               fold.cut.above.min = isTRUE(fl > z[["lambda.min"]]))
  })
  do.call(rbind, rows)
}

## The pivotal statistic of the quantile universal threshold and its
## simulated quantile under the null (wafc_lambda_qut()); it is shared with
## the gate of wafc_threshold() (step E2.5j), which compares the two:
## stat > quantile exactly when the QUT penalty level leaves some
## coefficient nonzero.
wafc_qut_pivot <- function(design, y, alpha = 0.05, nsim = 200L,
                           seed = NULL) {
  if (!inherits(design, "wafc_design")) {
    stop("'design' must be an object returned by wafc_design().", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single value in (0, 1).", call. = FALSE)
  }
  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1) {
    stop("'nsim' must be a single positive integer.", call. = FALSE)
  }
  nsim <- as.integer(nsim)
  n <- design[["n"]]
  y <- wafc_check_y(y, n)
  if (!is.null(seed)) set.seed(seed)
  unp <- design[["unpenalized"]]
  pen <- seq_len(design[["nvars"]])[-unp]
  W <- as.matrix(design[["Z"]][, unp, drop = FALSE])
  Zp <- design[["Z"]][, pen, drop = FALSE]
  ## wafc_qr_coef() is the least squares solve of wafc/R/competitors.R; it
  ## is shared and not duplicated.
  qrW <- qr(W)
  resid_of <- function(v) as.numeric(v - W %*% wafc_qr_coef(qrW, v))
  ## Pivotal statistic under the null: the ratio does not depend on sigma,
  ## and the observed residual norm carries the scale.
  E <- matrix(stats::rnorm(n * nsim), n, nsim)
  R <- E - W %*% wafc_qr_coef(qrW, E)
  G <- as.matrix(Matrix::crossprod(Zp, R))
  ratio <- apply(abs(G), 2L, max) / sqrt(colSums(R^2))
  r0 <- resid_of(y)
  r0norm <- sqrt(sum(r0^2))
  list(quantile = as.numeric(stats::quantile(ratio, 1 - alpha,
                                             names = FALSE)),
       r0norm = r0norm,
       stat = max(abs(as.numeric(Matrix::crossprod(Zp, r0)))) / r0norm)
}

## ---------------------------------------------------------------------------
## Internals
## ---------------------------------------------------------------------------

## Routing of the '...' of the tuning functions. It carries arguments for
## two callees: wafc_design(), which takes the basis, the margin and the
## rescaling, and wafc(), which takes the engine (thresh, maxit, asparse,
## standardize, intercept). Step E2.3 sent the whole lot to wafc_design(),
## so cv.wafc(x, u, y, thresh = 1e-7) stopped with "unused argument" and
## the tolerance of the fit could not be reached through the interface at
## all; that is why the default of 1e-10 survived three steps unmeasured
## (docs/ESTADO.md, 2026-09-21). Split by name, and refuse a name that is
## neither, so that a typo is an error and not a silent default.
wafc_split_dots <- function(dots) {
  if (length(dots) == 0L) return(list(design = list(), fit = list()))
  if (is.null(names(dots)) || any(!nzchar(names(dots)))) {
    stop("every argument passed through '...' must be named.", call. = FALSE)
  }
  dnames <- setdiff(names(formals(wafc_design)), c("x", "u", "J", "..."))
  fnames <- setdiff(names(formals(wafc)),
                    c("x", "u", "y", "J", "design", "...", dnames))
  unknown <- setdiff(names(dots), c(dnames, fnames))
  if (length(unknown) > 0L) {
    stop("unused argument(s) in '...': ", paste(unknown, collapse = ", "),
         ". They belong to neither wafc_design() nor wafc().", call. = FALSE)
  }
  list(design = dots[names(dots) %in% dnames],
       fit = dots[names(dots) %in% fnames])
}

## The line print.cv.wafc() adds when the cap of wafc_design() lowered the
## level of some modulating covariate in some candidate (step E3.4): which
## covariate, to what, its number of distinct points, and the candidates it
## acted on. NULL when it acted on none, or the object predates E3.4.
wafc_cv_cap_text <- function(x) {
  Je <- x[["J.eff"]]
  if (is.null(Je)) return(NULL)
  J <- x[["J"]]
  low <- which(colSums(Je < J) > 0L)
  if (length(low) == 0L) return(NULL)
  nd <- x[["wafc.fit"]][["design"]][["ndistinct"]]
  parts <- vapply(low, function(m) {
    sprintf("%s at %d (%d points) for J = %s", colnames(Je)[m],
            max(Je[, m]), nd[m], paste(J[Je[, m] < J], collapse = ", "))
  }, "")
  paste0("Levels capped by the distinct values: ",
         paste(parts, collapse = "; "))
}

## Grid of candidate resolution levels, 2:8 by default (decision D34,
## proposal P1 of step E2.4, ratified on 2026-09-28). It used to be the rule
## of cv.wall, 2:ceiling(log2(n)/2), which is 2:5 at n = 1000; the oracle of
## that grid sat on its top in every replicate of every cell of the pilot
## with a component to resolve, and widening it to 2:8 cut the prediction
## error by 17.6% and the error of the components by 32% in the
## inhomogeneous scenario at n = 1000, in 50 replicates of 50, while
## changing nothing in the smooth one (E1.3 had measured that bumps only
## resolves from J = 9). The top does not grow with n: 8 is the depth that
## was measured, at n from 250 to 1000, and a larger sample that wants more
## asks for it through 'J'. The grid starts at 2 because a block of a single
## wavelet column is not an approximation space, and step E2.4 measured that
## carrying J = 1 buys nothing; the argument 'n' is kept for the callers.
wafc_J_top <- 8L
wafc_J_grid <- function(J, n) {
  if (is.null(J)) return(seq(2L, wafc_J_top))
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
  ## what each fold fit returned of the path it was asked for (step E2.5j)
  conv <- data.frame(fold = seq_len(nfolds), nreturned = NA_integer_,
                     jerr = NA_integer_, lambda.last = NA_real_)
  for (k in seq_len(nfolds)) {
    out <- which(foldid == k)
    fit <- wafc(design = wafc_subset_design(design, -out), y = y[-out],
                penalty = penalty, lambda = lam, ...)
    conv[k, -1L] <- list(fit[["conv"]][["nreturned"]],
                         fit[["conv"]][["jerr"]], min(fit[["lambda"]]))
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
       nzero.min = full[["nzero"]][imin], conv.folds = conv)
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

## The fit a method of "cv.wafc" reads at 's': the thresholded fit when
## there is one and s is lambda.min (unless thresholded = FALSE), the path
## of wafc.fit otherwise. The thresholded fit is the "wafc" object of
## wafc_cv_thr_view(), a path of one point, so every method of "wafc" reads
## it unchanged.
wafc_cv_fit <- function(object, s, thresholded = NULL) {
  s <- wafc_cv_s(object, s)
  if (!is.null(thresholded) &&
      (!is.logical(thresholded) || length(thresholded) != 1L ||
       is.na(thresholded))) {
    stop("'thresholded' must be NULL, TRUE or FALSE.", call. = FALSE)
  }
  has <- !is.null(object[["threshold"]])
  at_min <- isTRUE(s == object[["lambda.min"]])
  use <- if (is.null(thresholded)) has && at_min else thresholded
  if (isTRUE(use)) {
    if (!has) {
      stop("This \"cv.wafc\" object has no threshold (threshold = ",
           "\"none\"); use thresholded = FALSE or NULL.", call. = FALSE)
    }
    if (!at_min) {
      stop("The threshold was chosen at lambda.min; at another penalty ",
           "level use thresholded = FALSE or NULL.", call. = FALSE)
    }
    return(list(fit = wafc_cv_thr_view(object), s = s))
  }
  list(fit = object[["wafc.fit"]], s = s)
}

## What cv.wafc() keeps of the threshold: the rule, the threshold and the
## blocks kept, and the coefficients of the thresholded fit in the
## coordinates of the design (the coefficients of wafc.fit at lambda.min
## with the blocks not kept set to zero, the estimator of Corollary 8 with
## no refit). The object of wafc_threshold() itself is not kept: its
## functions carry the design in their environment, and a saved "cv.wafc"
## would carry it twice.
wafc_cv_thr_slot <- function(object, th) {
  ex <- th[["extra"]]
  fit <- object[["wafc.fit"]]
  cf <- wafc_raw_coef(fit, s = object[["lambda.min"]])[, 1L]
  est <- wafc_thr_apply(cf[[1L]], unname(cf[-1L]), fit[["design"]],
                        ex[["kept"]])
  pen <- seq_len(fit[["nvars"]])[-fit[["design"]][["unpenalized"]]]
  list(rule = ex[["rule"]], t = ex[["t"]], c = ex[["c"]],
       norm = ex[["norm"]], kept = ex[["kept"]], a0 = est[["a0"]],
       b = est[["b"]], nzero = sum(est[["b"]][pen] != 0),
       candidates = ex[["candidates"]], time = th[["time"]])
}

## The thresholded fit of a "cv.wafc" as a "wafc" object: wafc.fit with its
## path replaced by the single point lambda.min and the coefficients there
## by the thresholded ones. coef(), predict(), wafc_functions() and
## wafc_blocks() read it as they read any fit.
wafc_cv_thr_view <- function(object) {
  th <- object[["threshold"]]
  fit <- object[["wafc.fit"]]
  des <- fit[["design"]]
  beta <- Matrix::Matrix(matrix(th[["b"]], ncol = 1L), sparse = TRUE)
  beta <- methods::as(beta, "CsparseMatrix")
  dimnames(beta) <- list(colnames(des[["Z"]]), NULL)
  fit[["lambda"]] <- object[["lambda.min"]]
  fit[["beta"]] <- beta
  fit[["a0"]] <- th[["a0"]]
  fit[["cc"]] <- wafc_levels(beta, th[["a0"]], des, fit[["carrier"]])
  fit[["nzero"]] <- as.integer(th[["nzero"]])
  fit[["threshold"]] <- th[c("rule", "t", "kept")]
  fit
}

## The cross-validation of one design for the block LASSO (step E3.1): the
## one of the engine, grpreg::cv.grpreg() on the folds of cv.wafc(), called
## as wafc_fit_klopp() calls it, so that the pair selected and the fit at it
## are the ones of wafc_fit_klopp(penalize.levels = FALSE, balanced = TRUE).
## The fit on the whole sample is the one cv.grpreg() returns, turned into
## the "wafc" object wafc() returns for the same design. The loss is the
## one of cv.grpreg() for "mse", read as it reports it, and is recomputed
## from its held-out predictions for "mae"; in both the mean is over the n
## observations and the standard error is sd / sqrt(n), the convention of
## cv.grpreg().
wafc_cv_block <- function(design, y, foldid, lambda, nlambda,
                          lambda.min.ratio, type.measure, fit_args, Ji) {
  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop("penalty = \"block\" needs the package 'grpreg' ",
         "(see docs/CONTINUAR.md, section 2).", call. = FALSE)
  }
  wafc_block_check(fit_args[["intercept"]],
                   if ("standardize" %in% names(fit_args))
                     fit_args[["standardize"]] else FALSE)
  ## the default tolerance is the one wafc() has for this penalty, read from
  ## its formals so that the two cannot drift apart
  thresh <- if ("thresh" %in% names(fit_args)) fit_args[["thresh"]] else
    eval(formals(wafc)[["thresh"]], list(penalty = "block"))
  maxit <- fit_args[["maxit"]]
  a <- wafc_block_args(design, y, lambda, nlambda, lambda.min.ratio,
                       fit_args[["block.size"]], thresh, maxit)
  cvg <- do.call(grpreg::cv.grpreg,
                 c(a[["args"]], list(fold = foldid,
                                     returnY = type.measure == "mae")))
  carrier <- wafc_carrier(design)
  full <- wafc_block_new(design, y, cvg[["fit"]], a[["blk"]], carrier,
                         call = as.call(list(as.symbol("wafc"),
                                             design = quote(design),
                                             y = quote(y), penalty = "block")),
                         nlambda = if (is.null(lambda)) as.integer(nlambda)
                           else length(lambda),
                         max.iter = if (is.null(maxit))
                           wafc_grpreg_default("max.iter") else
                             as.integer(maxit))
  lam <- as.numeric(cvg[["lambda"]])
  n <- length(y)
  if (type.measure == "mse") {
    cvm <- as.numeric(cvg[["cve"]])
    cvsd <- as.numeric(cvg[["cvse"]])
  } else {
    E <- abs(y - cvg[["Y"]])
    cvm <- colMeans(E)
    cvsd <- apply(E, 2L, stats::sd) / sqrt(n)
  }
  imin <- which.min(cvm)
  i1se <- min(which(cvm <= cvm[imin] + cvsd[imin]))
  nz <- full[["nzero"]][match(lam, full[["lambda"]])]
  z <- list(lambda = lam, cvm = cvm, cvsd = cvsd, cvup = cvm + cvsd,
            cvlo = cvm - cvsd, nzero = nz, lambda.min = lam[imin],
            lambda.1se = lam[i1se], cvm.min = cvm[imin],
            cvsd.min = cvsd[imin], nzero.min = nz[imin], conv.folds = NULL,
            J = Ji,
            conv = wafc_kp_conv(Ji, cvg,
                                nlambda = full[["conv"]][["nlambda"]],
                                max.iter = full[["conv"]][["max.iter"]]))
  list(cv = z, fit = full)
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
## The intercept column enters only when the fit has one, which by default
## (wafc(), argument 'intercept') is when the design has a constant
## covariate; always adding it, as this function did, gave the entry point
## of a fit without an intercept the wrong residual.
wafc_lambda_max <- function(design, y,
                            intercept = length(design[["constant"]]) > 0L) {
  n <- design[["n"]]
  unp <- design[["unpenalized"]]
  pen <- seq_len(design[["nvars"]])[-unp]
  W <- as.matrix(design[["Z"]][, unp, drop = FALSE])
  if (isTRUE(intercept)) W <- cbind(`(Intercept)` = 1, W)
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
