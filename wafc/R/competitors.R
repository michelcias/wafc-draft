## wafc/R/competitors.R -- the competitors of the pilot (step E2.4).
##
## The list of competitors is the one proposal L2d fixed (docs/ESTADO.md,
## table of decisions), and every one of them is asked the same three
## questions the WAFC is asked: what is beta_l(u), what is g_{lm}, and which
## blocks (l, m) did you keep. Seven methods answer, five of them
## competitors and two of them references:
##
##   "gam"      mgcv::gam with s(u_m, by = x_l), the additive varying
##              coefficient model with smoothing parameters by REML. It is
##              the method a referee will call the obvious baseline.
##   "bsgl"     B-splines plus group LASSO by block (l, m) (grpreg), the
##              same sieve-plus-selection architecture as the WAFC with the
##              wavelet basis replaced by a spline one. It is the comparison
##              that isolates the basis.
##   "klopp"    the block LASSO of Klopp and Pensky (2015) on the WAFC
##              design: the wavelet coefficients of each block (l, m) are
##              grouped in consecutive chunks of about log n, and the level
##              terms are penalized as blocks of their own, which is what
##              their norm (3.1) does and what decision D3 does not do. It
##              is the estimator their theory covers, run on the design the
##              WAFC uses, which is the question decision D18 invites.
##              The weight of a chunk (the sqrt(|G|) of grpreg or their
##              one) and the merging of the coarse levels into one chunk
##              are options, measured in step E2.5c; the balanced chunks,
##              in which every chunk of a finer level has between the
##              chunk size and twice it, are a third, measured in step
##              E2.5e; the coarse levels left unpenalized beside those
##              balanced chunks are a fourth, measured in step E2.5f.
##   "aspline"  a spline with knots chosen adaptively per block, in the
##              spirit of Wang, Jiang and Liu (2024). Their knot search is
##              an exact dynamic program; the one here is greedy forward
##              insertion scored by BIC, inside a block coordinate descent
##              over the pq terms. See the note on the function.
##   "vcbart"   VCBART (Deshpande et al., 2026), Bayesian trees for varying
##              coefficients: several modulating covariates and no additive
##              restriction, so it is the competitor that is not handicapped
##              by a structure the WAFC assumes and the scenarios satisfy.
##   "linear"   ordinary least squares of y on x alone, the model with
##              constant coefficients. It is the oracle of the null
##              scenario and the floor everywhere else.
##   "oracle"   the WAFC restricted to the blocks that are really active,
##              cross-validated over lambda at fixed J. It is not a
##              competitor but a reference: the price of not knowing the
##              structure is the distance from the WAFC to this column.
##              The LASSO by default; since step E4.1c, penalty = "block"
##              makes it the block LASSO of decision D44 (decision D63).
##
## Every fitter returns an object of class "wafc_competitor" with the same
## four accessors, so the pilot loops over methods and not over special
## cases: beta(u) gives the n by p matrix of functional coefficients,
## g(grid) gives the p by q list of additive components (NULL when the
## method does not decompose, which is the case of vcbart), blocks gives
## the p by q matrix of selected blocks (NULL when the method selects
## nothing), and predict() gives rowSums(newx * beta(newu)), plus the
## 'intercept' slot when the method has one that no constant covariate can
## carry (bsgl and klopp, whose engine always fits an intercept; zero when
## the design has a constant covariate, which then carries it).
##
## The penalty rule step E2.4 added, wafc_lambda_qut(), used to live here
## for want of a catalogue entry; step E2.4b moved it to wafc/R/tune.R,
## beside the five rules it belongs with, and exposed it as
## wafc_tune(rule = "qut").

#' Fit a competitor of the WAFC
#'
#' Common entry point to the methods the pilot of step E2.4 compares with
#' \code{\link{wafc}}. Every method is fitted on the same data and answers
#' the same questions, so the comparison is one table and not six.
#'
#' @param method One of \code{"gam"}, \code{"bsgl"}, \code{"klopp"},
#'   \code{"aspline"}, \code{"vcbart"}, \code{"linear"} or \code{"oracle"};
#'   see the file header for what each one is.
#' @param x Matrix of linear covariates, \eqn{n} by \eqn{p}. A constant
#'   column, the usual \eqn{X_1 \equiv 1}, gives the additive intercept.
#'   For \code{print}, an object of class \code{"wafc_competitor"}.
#' @param u Matrix of modulating covariates, \eqn{n} by \eqn{q}.
#' @param y Numeric response of length \eqn{n}.
#' @param active Logical \eqn{p} by \eqn{q} matrix of the blocks that are
#'   really active. Required by \code{method = "oracle"} and ignored by
#'   every other method.
#' @param ... Passed to the fitter of the method; ignored by
#'   \code{predict} and \code{print}.
#'
#' @return An object of class \code{"wafc_competitor"}: a list with the
#'   levels \code{cc}, the functions \code{beta(u)} (the \eqn{n} by
#'   \eqn{p} matrix of functional coefficients) and \code{g(grid)} (the
#'   \eqn{p} by \eqn{q} list of additive components, \code{NULL} when the
#'   method does not decompose), the \eqn{p} by \eqn{q} logical matrix
#'   \code{blocks} of the blocks kept (\code{NULL} when the method selects
#'   nothing), the \code{fitted} values, the \code{intercept} that no
#'   constant covariate carries (when the method has one), the object
#'   \code{fit} of the engine, the \code{design} for the methods fitted on
#'   the design of \code{\link{wafc_design}}, what is particular to the
#'   method in \code{extra}, and the \code{method},
#'   \code{time}, \code{n}, \code{p}, \code{q}, \code{xnames} and
#'   \code{unames}. \code{predict} returns
#'   \code{rowSums(newx * beta(newu))} plus the intercept, and \code{print}
#'   returns \code{x} invisibly.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' fit <- wafc_competitor("linear", d$x, d$u, d$y)
#' head(fit$beta(d$u))
#'
#' @export
wafc_competitor <- function(method = c("gam", "bsgl", "klopp", "aspline",
                                       "vcbart", "linear", "oracle"),
                            x, u, y, active = NULL, ...) {
  method <- match.arg(method)
  x <- wafc_as_matrix(x, "x")
  u <- wafc_as_matrix(u, "u")
  y <- wafc_check_y(y, nrow(x))
  if (nrow(u) != nrow(x)) {
    stop("'x' and 'u' must have the same number of rows.", call. = FALSE)
  }
  fitter <- switch(method, gam = wafc_fit_gam, bsgl = wafc_fit_bsgl,
                   klopp = wafc_fit_klopp, aspline = wafc_fit_aspline,
                   vcbart = wafc_fit_vcbart, linear = wafc_fit_linear,
                   oracle = wafc_fit_oracle)
  ## The pilot calls every method with one list of arguments, so an
  ## argument that belongs to another fitter is dropped here instead of
  ## raising an error. A fitter with a '...' of its own keeps everything
  ## else; only a name that is some other fitter's argument is removed, and
  ## so is a wafc_design() argument unless the '...' of the fitter reaches
  ## wafc_design(). The second rule is what the pilot needed: it passes
  ## 'wavelet.table' to every method (decision D31), gam and bsgl swallowed
  ## it through the '...' of their engines, and VCBART_ind() stopped with
  ## "unused argument", which removed the column from the table.
  args <- list(x = x, u = u, y = y)
  extra <- list(...)
  extra[["active"]] <- NULL
  own <- setdiff(names(formals(fitter)), "...")
  if ("..." %in% names(formals(fitter))) {
    drop <- setdiff(wafc_fitter_args(), own)
    if (!(method %in% wafc_design_fitters)) {
      drop <- union(drop, setdiff(wafc_design_args(), own))
    }
    extra <- extra[!(names(extra) %in% drop)]
  } else {
    extra <- extra[names(extra) %in% own]
  }
  if (method == "oracle") extra["active"] <- list(active)
  t0 <- proc.time()[["elapsed"]]
  out <- do.call(fitter, c(args, extra))
  out[["method"]] <- method
  out[["time"]] <- proc.time()[["elapsed"]] - t0
  out[["n"]] <- nrow(x)
  out[["p"]] <- ncol(x)
  out[["q"]] <- ncol(u)
  if (is.null(out[["xnames"]])) out[["xnames"]] <- wafc_names(x, ncol(x), "x")
  if (is.null(out[["unames"]])) out[["unames"]] <- wafc_names(u, ncol(u), "u")
  class(out) <- "wafc_competitor"
  out
}

#' @rdname wafc_competitor
#' @param object An object of class \code{"wafc_competitor"}.
#' @param newx,newu New covariates. When both are missing, the training
#'   sample is used.
#' @export
predict.wafc_competitor <- function(object, newx, newu, ...) {
  if (missing(newx) && missing(newu)) {
    return(object[["fitted"]])
  }
  if (missing(newx) || missing(newu)) {
    stop("Supply both 'newx' and 'newu', or neither.", call. = FALSE)
  }
  newx <- wafc_as_matrix(newx, "newx")
  newu <- wafc_as_matrix(newu, "newu")
  if (ncol(newx) != object[["p"]] || ncol(newu) != object[["q"]]) {
    stop("'newx' and 'newu' must have ", object[["p"]], " and ",
         object[["q"]], " column(s).", call. = FALSE)
  }
  a0 <- if (is.null(object[["intercept"]])) 0 else object[["intercept"]]
  as.numeric(a0 + rowSums(newx * object[["beta"]](newu)))
}

#' @rdname wafc_competitor
#' @param digits Number of significant digits printed.
#' @export
print.wafc_competitor <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  cat(sprintf("WAFC competitor \"%s\": n = %d, p = %d, q = %d, %.2f s\n",
              x[["method"]], x[["n"]], x[["p"]], x[["q"]], x[["time"]]))
  cat("  levels c:", paste(format(x[["cc"]], digits = digits), collapse = ", "),
      "\n")
  if (!is.null(x[["blocks"]])) {
    cat("  blocks kept:", sum(x[["blocks"]]), "of",
        length(x[["blocks"]]), "\n")
  }
  if (is.null(x[["g"]])) {
    cat("  no additive decomposition (beta_l(u) only)\n")
  }
  invisible(x)
}

#' Additive components of a fit on a grid
#'
#' Evaluates \eqn{\hat g_{\ell m}} of a \code{"wafc"} or a
#' \code{"wafc_competitor"} fit on a grid, in the single form the pilot
#' compares: a \eqn{p} by \eqn{q} list of numeric vectors, each one centred
#' on the grid, since every method fixes the level by its own convention and
#' only the shape is comparable.
#'
#' @param object A \code{"wafc"}, \code{"cv.wafc"} or
#'   \code{"wafc_competitor"} object. A \code{"cv.wafc"} object is read as
#'   \code{\link{coef.cv.wafc}} reads it: at \code{lambda.min} after its
#'   threshold, when it has one.
#' @param grid Matrix with \eqn{q} columns at which the components are
#'   evaluated, on the original scale of the modulating covariates.
#' @param s For a \code{"wafc"} object, the penalty level; for a
#'   \code{"cv.wafc"} object, as in \code{\link{coef.cv.wafc}}, with
#'   \code{NULL} for \code{"lambda.min"}.
#' @param design For a \code{"wafc"} object, an optional design already
#'   built on \code{grid} with \code{spec} equal to the fitted design,
#'   which is how the pilot avoids rebuilding it once per method.
#'
#' @return A \eqn{p} by \eqn{q} list of numeric vectors of length
#'   \code{nrow(grid)}, or \code{NULL} when the method has no additive
#'   decomposition.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' g <- wafc_grid_components(wafc_competitor("gam", d$x, d$u, d$y),
#'                           cbind(seq(0, 1, length.out = 64),
#'                                 seq(0, 1, length.out = 64)))
#' length(g[[1, 1]])
#'
#' @export
wafc_grid_components <- function(object, grid, s = NULL, design = NULL) {
  grid <- wafc_as_matrix(grid, "grid")
  if (inherits(object, "cv.wafc")) {
    f <- wafc_cv_fit(object, if (is.null(s)) "lambda.min" else s)
    object <- f[["fit"]]
    s <- f[["s"]]
  }
  if (inherits(object, "wafc_competitor")) {
    if (is.null(object[["g"]])) return(NULL)
    g <- object[["g"]](grid)
  } else if (inherits(object, "wafc")) {
    des <- object[["design"]]
    p <- des[["p"]]
    q <- des[["q"]]
    if (is.null(design)) {
      design <- wafc_design(matrix(1, nrow(grid), p), grid, spec = des)
    }
    b <- coef.wafc(object, s = wafc_single_s(object, s))[-1L, 1L]
    g <- vector("list", p * q)
    dim(g) <- c(p, q)
    for (l in seq_len(p)) {
      for (m in seq_len(q)) {
        idx <- des[["blocks"]][[wafc_block_name(des, l, m)]]
        g[[l, m]] <- as.numeric(design[["Z"]][, idx, drop = FALSE] %*% b[idx])
      }
    }
  } else {
    stop("'object' must be a \"wafc\", a \"cv.wafc\" or a ",
         "\"wafc_competitor\" fit.", call. = FALSE)
  }
  for (i in seq_along(g)) g[[i]] <- g[[i]] - mean(g[[i]])
  g
}

## ---------------------------------------------------------------------------
## The fitters
## ---------------------------------------------------------------------------

#' Spline basis dimension matched to a WAFC expansion
#'
#' The basis dimension the smooth of each modulating covariate needs for
#' \code{\link{wafc_fit_gam}} to be compared with \code{\link{wafc}} at
#' equal dimension and not at equal convenience: \eqn{2^J} per modulator,
#' the number of wavelet columns a block of resolution \eqn{J} has with
#' \eqn{j_0 = 0} and the constant scaling function discarded, truncated at
#' what the data admit, since a smooth cannot have more basis functions
#' than its covariate has distinct values.
#'
#' Step E6.1a is the reason this exists: on the three real candidates the
#' WAFC seemed to beat \code{mgcv} by 6.6\% to 15.2\% at the default
#' \code{k = 10} and tied with it at matched dimension, so the whole
#' apparent gain was one of dimension. The pilot of step E2.4 ran at the
#' default, and its verdict against this competitor has to be read again
#' at matched dimension.
#'
#' @param u Matrix of modulating covariates.
#' @param J The resolution level, or one per modulating covariate.
#' @param kmin Smallest dimension returned.
#'
#' @return An integer vector of length \code{ncol(u)}.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' wafc_k_matched(d$u, J = 4)
#'
#' @export
wafc_k_matched <- function(u, J, kmin = 3L) {
  u <- wafc_as_matrix(u, "u")
  q <- ncol(u)
  J <- wafc_recycle(J, q, "J")
  if (any(!is.finite(J)) || any(J != round(J)) || any(J < 1)) {
    stop("'J' must contain integer values of at least 1.", call. = FALSE)
  }
  ndist <- vapply(seq_len(q), function(m) length(unique(u[, m])), 0L)
  as.integer(pmax(kmin, pmin(2^as.integer(J), ndist - 1L)))
}

#' Additive varying coefficient model by penalized splines
#'
#' \code{mgcv::gam} on the formula
#' \eqn{y \sim 0 + x_1 + \cdots + x_p + \sum_{\ell m} s(u_m, by = x_\ell)},
#' with the smoothing parameters chosen by REML. The constant linear
#' covariate carries the intercept, exactly as in \code{\link{wafc}}, and
#' \pkg{mgcv} centres each smooth over the sample, which is the same
#' identifiability constraint the discarded scaling function imposes on the
#' wavelet side (Lemma 1 of E1.2), up to the difference between the
#' empirical and the Lebesgue measure (decision D22).
#'
#' Blocks are declared kept when the effective degrees of freedom of the
#' smooth exceed \code{edf.tol}: \pkg{mgcv} shrinks a null smooth towards a
#' straight line and, with \code{select = TRUE}, towards zero, but it never
#' returns an exactly zero fit, so any structure count for this method
#' needs a threshold and is reported with it.
#'
#' Two things step E2.4b changed, both from what step E6.1a measured on
#' real data (\file{docs/aplicacao-candidatas.md}). First, \code{k} may
#' now be given one value per modulating covariate, so that the basis
#' dimension can be matched to the \eqn{2^J} of the WAFC expansion
#' (\code{\link{wafc_k_matched}}): with the single default \code{k = 10}
#' the comparison is one of dimension and not one of basis, and on the
#' three real candidates that difference was the whole apparent gain of the
#' WAFC, 6.6\% to 15.2\% at \code{k = 10} against a tie at matched
#' dimension. Second, \code{engine = "bam"} fits the same model with
#' \code{\link[mgcv]{bam}}, \code{method = "fREML"} and discretized
#' covariates, which at these dimensions is two orders of magnitude faster
#' (6.1 s against 1453 s for eight smooths of \eqn{k = 23}). That is a
#' change of fitting algorithm plus a binning of the covariates, not a
#' change of estimator, and the number it produces is labelled with it
#' wherever it is reported.
#'
#' Step E2.5h added the choice of the basis dimension from a grid, with
#' \code{k.select}. The grid is 5, 10, 20, 40 and 80 (decision D41), one
#' grid for the three criteria, and every candidate is fitted and scored by
#' one criterion, the full search of Ruppert (2002) and option 3 of Pya and
#' Wood (2016, section 2); powers of two are not used, so that the spline
#' is not tuned on the grid of the wavelets. The grid is the one of Ruppert
#' (2002, section 3) without its last value, 120: his additive search
#' (section 6) stops at 40, 80 is kept because step E2.5d saw
#' \code{k = 64} bind in the inhomogeneous scenario, and 120 cost three to
#' four times the rest of the grid. Two things differ from Ruppert
#' and are declared, not corrected: his \eqn{K} is a number of knots of a
#' truncated power basis, while the \code{k} of \pkg{mgcv} is the dimension
#' of a thin plate regression spline basis, so the same number is a slightly
#' smaller space here; and his candidates stop below \eqn{n - p - 1}, while
#' here a candidate is truncated, modulator by modulator, at the number of
#' distinct values of the modulator minus one, as \code{\link{wafc_k_matched}}
#' does, and the candidates that the truncation makes equal are fitted once.
#' The value is common to the smooths, which is what Ruppert (2002, section
#' 6) recommends for additive models, since a search over one value per
#' term multiplies the cost and buys little once each value is large
#' enough.
#'
#' The two criteria are the two smoothing criteria of \pkg{mgcv}, and each
#' search uses the same criterion for the smoothing parameters and for
#' \code{k}. With \code{k.select = "reml"} the smoothing parameters are
#' estimated by REML, as without a search, and \code{k} minimises the
#' restricted negative log-likelihood, profiled over the scale, which is
#' the criterion of Kauermann and Opsomer (2011) with REML in place of
#' maximum likelihood. Their reason for maximum likelihood does not apply:
#' the fixed effects of the restricted likelihood are the \eqn{p} level
#' terms (with \code{select = TRUE} every coefficient of every smooth is
#' penalized), and they do not change with \code{k}, so the scores of two
#' candidates are values of one function of the same data. The score is
#' computed here from the fitted smoothing parameters
#' (\code{wafc_gam_reml()}) and not read from the engine: with
#' \code{engine = "gam"} the two coincide, but with \code{engine = "bam"}
#' the score \code{bam} reports under \code{discrete = TRUE} is off the
#' exact one by an amount that changes with \code{k} (1.09, 0.93 and 0.92
#' at \code{k} = 5, 10 and 20 in the test of this step), which would bias
#' the choice; the one the engine reports is kept beside it. With
#' \code{k.select = "gcv"} the smoothing parameters are chosen by
#' generalized cross-validation (\code{method = "GCV.Cp"}) and \code{k}
#' minimises the GCV score \eqn{n\,RSS/(n - \tau)^2} of the fit,
#' \eqn{\tau} the trace of the influence matrix, as in Ruppert (2002).
#' \code{bam} discretizes only under REML, so \code{engine = "bam"} fits
#' this search without \code{discrete = TRUE}.
#'
#' With \code{k.select = "cv"} the score of a candidate is the squared
#' error of prediction on the folds \code{foldid}, the mean over folds of
#' the mean over the observations a fold leaves out, which is the loss of
#' \code{\link{cv.wafc}}: given its folds, the spline and the WAFC choose
#' their dimension by the same criterion on the same partition. Inside each
#' fold the smoothing parameters are estimated by REML on the engine given,
#' and the fit returned is the one on the whole sample at the \code{k}
#' chosen. It costs one fit per fold and candidate, plus one.
#'
#' @param x,u,y The data, as in \code{\link{wafc_competitor}}.
#' @param k Basis dimension of each smooth: one value for every smooth, or
#'   one value per modulating covariate. With \code{k.select} other than
#'   \code{"none"}, the grid of candidates, each one common to the
#'   smooths. \code{NULL} (the default) is 10 without a search and the grid
#'   of decision D41 with one.
#' @param k.select \code{"none"} (the default) fits at \code{k};
#'   \code{"reml"}, \code{"gcv"} and \code{"cv"} fit every candidate of the
#'   grid and keep the one with the smallest score of that criterion.
#' @param nfolds,foldid Folds of \code{k.select = "cv"}, as in
#'   \code{\link{cv.wafc}}; ignored by the other criteria.
#' @param select Passed to \code{\link[mgcv]{gam}}: \code{TRUE} adds the
#'   extra penalty on the null space, which is what lets a smooth be shrunk
#'   away entirely and is the fair setting when half the blocks are zero.
#' @param edf.tol Threshold on the effective degrees of freedom above which
#'   a block counts as kept.
#' @param engine \code{"gam"} (the default) fits with
#'   \code{\link[mgcv]{gam}} and \code{method = "REML"}, which is what
#'   decision D30 asks for; \code{"bam"} fits with
#'   \code{\link[mgcv]{bam}}, \code{method = "fREML"} and
#'   \code{discrete = TRUE}.
#' @param ... Further arguments to \code{\link[mgcv]{gam}} or to
#'   \code{\link[mgcv]{bam}}.
#'
#' @return An object of class \code{"wafc_competitor"}, whose \code{extra}
#'   carries the \code{edf} matrix, the \code{edf.tol} used, the vector
#'   \code{k} (one value per modulating covariate), the \code{engine}, the
#'   \code{k.select} and the smoothing criterion \code{smooth.method}. With
#'   a search it also carries \code{k.table}, one row per candidate fitted
#'   (the value of the grid, the dimensions used after the truncation, the
#'   \code{score} of the criterion, its standard error \code{cvsd} for
#'   \code{"cv"}, the \code{score.engine} the engine reports, the total
#'   effective degrees of freedom, the number of coefficients and the
#'   seconds; for \code{"cv"} the last three are the means over the folds
#'   and \code{score.engine} is \code{NA}), and \code{k.top}, \code{TRUE}
#'   when the candidate chosen is the largest one fitted; with
#'   \code{"cv"}, also the \code{foldid} used.
#'
#' @references Kauermann, G. and Opsomer, J. D. (2011). Data-driven
#'   selection of the spline dimension in penalized spline regression.
#'   \emph{Biometrika} 98(1), 225-230.
#'
#'   Pya, N. and Wood, S. N. (2016). A note on basis dimension selection in
#'   generalized additive modelling. arXiv:1602.06696.
#'
#'   Ruppert, D. (2002). Selecting the number of knots for penalized
#'   splines. \emph{Journal of Computational and Graphical Statistics}
#'   11(4), 735-757.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' wafc_fit_gam(d$x, d$u, d$y)$blocks
#' wafc_fit_gam(d$x, d$u, d$y, k = c(5, 10), k.select = "reml")$extra$k.table
#'
#' @export
wafc_fit_gam <- function(x, u, y, k = NULL,
                         k.select = c("none", "reml", "gcv", "cv"),
                         select = TRUE, edf.tol = 0.1,
                         engine = c("gam", "bam"), nfolds = 10L,
                         foldid = NULL, ...) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("method = \"gam\" needs the package 'mgcv'.", call. = FALSE)
  }
  engine <- match.arg(engine)
  k.select <- match.arg(k.select)
  p <- ncol(x)
  q <- ncol(u)
  if (k.select == "none") {
    if (is.null(k)) k <- 10L
    k <- wafc_recycle(k, q, "k")
    if (any(!is.finite(k)) || any(k != round(k)) || any(k < 3)) {
      stop("'k' must contain integer values of at least 3.", call. = FALSE)
    }
    cand <- list(as.integer(k))
    grid <- NA_integer_
  } else {
    if (is.null(k)) k <- wafc_k_grid
    cand <- wafc_k_candidates(u, k)
    grid <- attr(cand, "grid")
  }
  if (k.select == "cv") {
    foldid <- wafc_foldid(foldid, nrow(x), nfolds)
    nfolds <- max(foldid)
  }
  smooth.method <- if (k.select == "gcv") "GCV.Cp" else "REML"
  xn <- wafc_names(x, p, "x")
  un <- wafc_names(u, q, "u")
  dat <- wafc_frame(x, u, xn, un)
  dat[["y"]] <- y
  ## One smooth per block, in the lexicographic order of D12, plus the p
  ## parametric terms that carry the levels. The 'by' variable is what
  ## predict(type = "terms") multiplies the smooth by, so evaluating the
  ## terms with every covariate set to one is what turns the term of the
  ## block (l, m) into the component g_{lm} itself.
  fit_at <- function(kk, rows = seq_len(nrow(dat))) {
    terms <- character(0)
    for (l in seq_len(p)) {
      for (m in seq_len(q)) {
        terms <- c(terms, sprintf("s(%s, k = %d, by = %s)", un[m], kk[m],
                                  xn[l]))
      }
    }
    fo <- stats::as.formula(paste("y ~ 0 +",
                                  paste(c(xn, terms), collapse = " + ")))
    d <- dat[rows, , drop = FALSE]
    if (engine == "gam") {
      mgcv::gam(fo, data = d, method = smooth.method, select = select, ...)
    } else if (smooth.method == "REML") {
      mgcv::bam(fo, data = d, method = "fREML", discrete = TRUE,
                select = select, ...)
    } else {
      mgcv::bam(fo, data = d, method = smooth.method, select = select, ...)
    }
  }
  ## The cross-validated error of one candidate: the fold fits on the rows
  ## kept, the squared error on the rows left out, the mean of the fold
  ## means and its standard error, as in wafc_cv_design().
  cv_at <- function(kk) {
    fm <- numeric(nfolds)
    edf <- numeric(nfolds)
    nc <- numeric(nfolds)
    for (f in seq_len(nfolds)) {
      out <- which(foldid == f)
      ff <- fit_at(kk, rows = -out)
      pr <- as.numeric(stats::predict(ff, newdata = dat[out, , drop = FALSE]))
      fm[f] <- mean((y[out] - pr)^2)
      edf[f] <- sum(ff[["edf"]])
      nc[f] <- length(stats::coef(ff))
      rm(ff)
    }
    list(cvm = mean(fm), cvsd = stats::sd(fm) / sqrt(nfolds),
         edf = mean(edf), ncoef = mean(nc))
  }
  ## The search keeps only the best fit so far: a fit of mgcv carries its
  ## model matrix, and at the top of the grid in the cell "mixed" that is
  ## 16 smooths of 80 columns. The cross-validation fits the whole sample
  ## once, at the k it chooses.
  fit <- NULL
  best <- Inf
  ibest <- 1L
  tab <- vector("list", length(cand))
  for (i in seq_along(cand)) {
    t0 <- proc.time()[["elapsed"]]
    if (k.select == "cv") {
      z <- cv_at(cand[[i]])
      score <- z[["cvm"]]
      row <- list(cvsd = z[["cvsd"]], score.engine = NA_real_,
                  edf = z[["edf"]], ncoef = z[["ncoef"]])
    } else {
      fi <- fit_at(cand[[i]])
      score <- switch(k.select, none = NA_real_,
                      reml = as.numeric(wafc_gam_reml(fi, y)),
                      gcv = as.numeric(fi[["gcv.ubre"]]))
      row <- list(cvsd = NA_real_,
                  score.engine = as.numeric(fi[["gcv.ubre"]]),
                  edf = sum(fi[["edf"]]), ncoef = length(stats::coef(fi)))
    }
    tab[[i]] <- data.frame(k = grid[i],
                           k.used = paste(cand[[i]], collapse = ","),
                           score = score, cvsd = row[["cvsd"]],
                           score.engine = row[["score.engine"]],
                           edf = row[["edf"]], ncoef = row[["ncoef"]],
                           time = proc.time()[["elapsed"]] - t0,
                           stringsAsFactors = FALSE)
    if (is.na(best) || i == 1L || (!is.na(score) && score < best)) {
      if (k.select != "cv") fit <- fi
      best <- score
      ibest <- i
    }
    if (k.select != "cv") rm(fi)
  }
  k <- cand[[ibest]]
  if (k.select == "cv") fit <- fit_at(k)
  ## Effective degrees of freedom by smooth, in the order the terms were
  ## written, which is the lexicographic order of D12.
  edf <- vapply(fit[["smooth"]],
                function(sm) sum(fit[["edf"]][sm[["first.para"]]:sm[["last.para"]]]),
                0)
  blocks <- matrix(edf > edf.tol, p, q, byrow = TRUE,
                   dimnames = list(xn, un))
  cc <- as.numeric(stats::coef(fit)[xn])
  names(cc) <- xn

  term_at <- function(newu) {
    nd <- wafc_frame(matrix(1, nrow(newu), p), newu, xn, un)
    ## The generic, and not mgcv::predict.gam: a bam fitted with
    ## discrete = TRUE has its own method.
    stats::predict(fit, newdata = nd, type = "terms")
  }
  beta_fun <- function(newu) {
    newu <- wafc_as_matrix(newu, "newu")
    tm <- term_at(newu)
    out <- matrix(rep(cc, each = nrow(newu)), nrow(newu), p,
                  dimnames = list(NULL, xn))
    b <- 0L
    for (l in seq_len(p)) {
      for (m in seq_len(q)) {
        b <- b + 1L
        out[, l] <- out[, l] + tm[, wafc_gam_term(colnames(tm), b)]
      }
    }
    out
  }
  g_fun <- function(grid) {
    grid <- wafc_as_matrix(grid, "grid")
    tm <- term_at(grid)
    g <- vector("list", p * q)
    dim(g) <- c(p, q)
    b <- 0L
    for (l in seq_len(p)) {
      for (m in seq_len(q)) {
        b <- b + 1L
        g[[l, m]] <- as.numeric(tm[, wafc_gam_term(colnames(tm), b)])
      }
    }
    dimnames(g) <- list(xn, un)
    g
  }
  list(fit = fit, cc = cc, beta = beta_fun, g = g_fun, blocks = blocks,
       fitted = as.numeric(stats::fitted(fit)), xnames = xn, unames = un,
       extra = list(edf = matrix(edf, p, q, byrow = TRUE,
                                 dimnames = list(xn, un)),
                    edf.tol = edf.tol, k = k, engine = engine,
                    k.select = k.select, smooth.method = smooth.method,
                    k.table = if (k.select == "none") NULL else
                      do.call(rbind, tab),
                    k.top = if (k.select == "none") NA else
                      ibest == length(cand),
                    foldid = if (k.select == "cv") foldid else NULL))
}

#' B-spline sieve with a group LASSO by block
#'
#' The same architecture as \code{\link{wafc}} with the wavelet basis
#' replaced by a cubic B-spline one and the LASSO replaced by the group
#' LASSO of \code{grpreg}, one group per block \eqn{(\ell, m)}: it is the
#' comparison that isolates the basis, since the design, the unpenalized
#' level terms and the selection unit are the ones of decision D12.
#'
#' The basis of each modulating covariate is centred on the training sample,
#' which is the spline analogue of discarding \eqn{\phi_{00}}, and the
#' penalty level is chosen by cross-validation. The number of basis
#' functions is chosen from a grid of the form \eqn{2^J}, so the two
#' expansions are compared at comparable dimensions rather than at a number
#' pulled from the air. The default candidates are \eqn{2^J} for
#' \eqn{J = 2, \ldots, 8}, the grid of the WAFC since decision D34 (decision
#' D36). They were 4, 8 and 16 before, the dimensions of \eqn{J \le 4},
#' and that truncated this competitor as the old grid truncated the WAFC:
#' in the inhomogeneous scenario at \eqn{n = 1000} the extended grid
#' selects 64 functions and lowers the cross-validated error from 2.129 to
#' 1.958, at the price of about 60 s against 0.9 s. A candidate larger
#' than the number of distinct values of some modulating covariate has no
#' data between its knots and is dropped.
#'
#' @param x,u,y The data.
#' @param df Number of B-spline basis functions per block, or a vector of
#'   candidates scored by the same cross-validated loss.
#' @param nfolds,foldid Folds of the cross-validation. Supplying
#'   \code{foldid} is how the pilot makes every method use one partition.
#' @param penalty Passed to \code{\link[grpreg]{cv.grpreg}}:
#'   \code{"grLasso"} is the group LASSO.
#' @param ... Further arguments to \code{\link[grpreg]{cv.grpreg}}.
#'
#' @return An object of class \code{"wafc_competitor"}, whose \code{extra}
#'   carries the dimension chosen, \code{df}, and \code{J = log2(df)}, the
#'   level of the WAFC block with as many columns (\code{NA} when \code{df}
#'   is not a power of two).
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' wafc_fit_bsgl(d$x, d$u, d$y, df = 8)$blocks
#'
#' @export
wafc_fit_bsgl <- function(x, u, y, df = 2L^(2:8), nfolds = 10L,
                          foldid = NULL, penalty = "grLasso", ...) {
  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop("method = \"bsgl\" needs the package 'grpreg'.", call. = FALSE)
  }
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(u)
  xn <- wafc_names(x, p, "x")
  un <- wafc_names(u, q, "u")
  foldid <- wafc_foldid(foldid, n, nfolds)
  df <- sort(unique(as.integer(df)))
  if (any(df < 3L)) stop("'df' must be at least 3.", call. = FALSE)
  ## the same truncation as wafc_k_matched(): a block cannot have more basis
  ## functions than its covariate has distinct values
  ndist <- min(vapply(seq_len(q), function(m) length(unique(u[, m])), 0L))
  df <- if (any(df <= ndist - 1L)) df[df <= ndist - 1L] else min(df)

  best <- NULL
  for (d in df) {
    sp <- wafc_bs_spec(u, d, un)
    Z <- wafc_bs_design(x, u, sp)
    grp <- c(rep(0L, p), rep(seq_len(p * q), each = d - 1L))
    cv <- grpreg::cv.grpreg(Z, y, group = grp, penalty = penalty,
                            fold = foldid, ...)
    val <- min(cv[["cve"]])
    if (is.null(best) || val < best[["cve"]]) {
      best <- list(cve = val, df = d, spec = sp, cv = cv, group = grp)
    }
    rm(Z)
  }
  d <- best[["df"]]
  sp <- best[["spec"]]
  b <- as.numeric(stats::coef(best[["cv"]]))
  a0 <- b[1L]
  b <- b[-1L]
  cc <- b[seq_len(p)]
  names(cc) <- xn
  ## grpreg always fits its own intercept; the design has a constant column
  ## when the model has one, so the intercept is folded into its level, the
  ## convention of coef.wafc()
  const <- which(vapply(seq_len(p), function(l) diff(range(x[, l])) == 0, TRUE))
  if (length(const) == 1L) {
    cc[const] <- cc[const] + a0 / x[1L, const]
    a0 <- 0
  }
  idx <- function(l, m) p + (d - 1L) * ((l - 1L) * q + m - 1L) + seq_len(d - 1L)
  nz <- matrix(FALSE, p, q, dimnames = list(xn, un))
  for (l in seq_len(p)) {
    for (m in seq_len(q)) nz[l, m] <- any(b[idx(l, m)] != 0)
  }
  g_of <- function(newu) {
    newu <- wafc_as_matrix(newu, "newu")
    B <- lapply(seq_len(q), function(m) wafc_bs_eval(newu[, m], sp, m))
    g <- vector("list", p * q)
    dim(g) <- c(p, q)
    for (l in seq_len(p)) {
      for (m in seq_len(q)) {
        g[[l, m]] <- as.numeric(B[[m]] %*% b[idx(l, m)])
      }
    }
    dimnames(g) <- list(xn, un)
    g
  }
  ## Without a constant covariate the intercept of grpreg has no level to be
  ## folded into, and it stays out of beta: adding it to beta_1, as this
  ## function did, multiplied it by X_1 in every prediction.
  beta_fun <- function(newu) {
    g <- g_of(newu)
    nr <- length(g[[1L, 1L]])
    out <- matrix(rep(cc, each = nr), nr, p, dimnames = list(NULL, xn))
    for (l in seq_len(p)) {
      for (m in seq_len(q)) out[, l] <- out[, l] + g[[l, m]]
    }
    out
  }
  ## The dimension chosen is also recorded as the level J of the WAFC block
  ## with as many columns (2^J basis functions, 2^J - 1 after centring), so
  ## the pilot can read it in the column it reads the J of the others from;
  ## a candidate that is not a power of two has no such level.
  Jd <- log2(d)
  list(fit = best[["cv"]], cc = cc, beta = beta_fun, g = g_of, blocks = nz,
       fitted = as.numeric(a0 + rowSums(x * beta_fun(u))), intercept = a0,
       xnames = xn, unames = un,
       extra = list(df = d,
                    J = if (Jd == round(Jd)) as.integer(Jd) else NA_integer_,
                    lambda = best[["cv"]][["lambda.min"]],
                    cve = best[["cve"]]))
}

#' Block LASSO of Klopp and Pensky on the WAFC design
#'
#' The estimator whose theory covers the case \eqn{q = 1} with
#' \eqn{X \perp U}, run on the design of \code{\link{wafc_design}} so that
#' the comparison with \code{\link{wafc}} is a comparison of penalties and
#' not of designs. Two things differ from decision D3, and both are theirs:
#' the wavelet coefficients of each block \eqn{(\ell, m)} are grouped in
#' consecutive chunks whose size is about \eqn{\log n} and penalized by the
#' Euclidean norm of the chunk, which is their norm (3.1); and the level
#' terms \eqn{c_\ell} are penalized as a block of their own, where the WAFC
#' leaves them free.
#'
#' By default the chunks respect the order of decision D12, so a chunk
#' never straddles two resolution levels (the coarse chunk of
#' \code{merge.coarse} and \code{balanced}, below, is the exception): the
#' levels \eqn{j} with \eqn{2^j} at most the chunk size are chunks of their
#' own, and a finer level is cut into consecutive pieces of the chunk size. That is the reading of "blocks of
#' size about \eqn{\log n} inside each functional coefficient" that keeps
#' the block a set of neighbouring translates at one scale.
#'
#' The weight of a chunk in the penalty is the one open question 34 of
#' \code{docs/ESTADO.md} asks about (step E2.5c). By default it is the one
#' of grpreg, the square root of the size of the chunk, and then the
#' penalty level the theory needs is set by the chunks of size one of the
#' coarse levels, which brings the bound back to the order of the LASSO;
#' their norm (3.1), and the theory of step E1.11, have weight one on every
#' chunk. \code{merge.coarse} is the third reading: the coarse levels of a
#' block, the ones with \eqn{2^j} below the chunk size, become one chunk,
#' which removes the chunks of size one of level 0 without touching the
#' weight of the full chunks.
#'
#' \code{balanced} is the fourth reading (step E2.5e): the coarse levels
#' are merged as with \code{merge.coarse}, and the short piece at the end
#' of each finer level, the part the chunk size \eqn{b_n} does not divide,
#' is absorbed into the chunk before it in the same level. Every chunk of a
#' finer level then has between \eqn{b_n} and \eqn{2 b_n - 1} columns.
#' The coarse chunk has \eqn{2^{j^* + 1} - 1} columns, with \eqn{j^*} the
#' finest level of the block with \eqn{2^{j^*} < b_n}; that is at most
#' \eqn{2 b_n - 3}, and it falls below \eqn{b_n} exactly when
#' \eqn{2^{j^* + 1} \le b_n}, which happens in two cases: \eqn{b_n} a power
#' of two (the chunk has \eqn{b_n - 1} columns), and a block with no level
#' of \eqn{2^j \ge b_n} and \eqn{2^J \le b_n}, which is then one chunk of
#' \eqn{2^J - 1} columns (with the \eqn{b_n} of 6 and 7 of the pilot, only
#' \eqn{J = 2}, a chunk of 3). With the weights of grpreg the ratio of the
#' largest weight of a wavelet chunk to the smallest is then at most
#' \eqn{\sqrt{(2 b_n - 1)/(b_n - 1)}}, bounded in \eqn{n}, where the other
#' two readings can leave a chunk of one column next to chunks of
#' \eqn{b_n} or more (the default always does, at level 0), a ratio of the
#' order of \eqn{\sqrt{b_n}} that is not bounded in \eqn{n}; whether the
#' theory of step E1.11 accepts weights of bounded ratio is answered in
#' \code{derivations/08a-sondagem-blocos.md}, section 11.
#'
#' \code{free.coarse} is the fifth reading (step E2.5f, open question 37):
#' the coarse levels of each block, the ones with \eqn{2^j} below the chunk
#' size, are left unpenalized, in group 0 of grpreg beside the level terms
#' when those are free, instead of being one penalized chunk. With
#' \code{balanced = TRUE} the finer levels keep the balanced chunks. In the
#' theory those \eqn{2^{j^* + 1} - 1} columns per block join the
#' unpenalized part, as the scaling coefficients of step E1.8 do; see
#' \code{derivations/08a-sondagem-blocos.md}, section 12. The price is that
#' no block is ever zero: a block with a coarse level always has a nonzero
#' coefficient, so \code{blocks} marks every block, and the selection is
#' read from \code{extra$blocks.fine}, the blocks whose penalized part is
#' nonzero, and \code{extra$nzero}, the number of nonzero penalized
#' coefficients. A resolution at which no block has a level with
#' \eqn{2^j} at least the chunk size (with the chunk sizes 6 and 7 of the
#' pilot, \eqn{J \le 3}) has nothing to penalize; grpreg does not fit
#' that, and the candidate is fitted by least squares instead, scored by the
#' same cross-validated squared error grpreg reports, with \code{lambda}
#' recorded as \code{NA}.
#'
#' @param x,u,y The data.
#' @param J Resolution level, or a vector of candidates scored by the same
#'   cross-validated loss. \code{NULL} uses the grid of
#'   \code{\link{cv.wafc}}.
#' @param block.size Size of the chunks. \code{NULL} is
#'   \code{ceiling(log(n))}, the choice of their norm (3.1).
#' @param penalize.levels Whether the level terms are penalized, as they are
#'   in Klopp and Pensky. \code{FALSE} gives the block LASSO with the
#'   unpenalized levels of decision D3, which isolates the effect of the
#'   grouping alone.
#' @param chunk.weights Weight of each chunk. \code{"sqrt"}, the default, is
#'   the default of grpreg, which is the square root of the rank of the
#'   chunk once grpreg has dropped its null columns, that is
#'   \eqn{\sqrt{|G|}} for a chunk of full rank; nothing is passed to
#'   grpreg, so the fit is the one of before the argument existed.
#'   \code{"unit"} gives weight one to every chunk, as their norm (3.1)
#'   does.
#' @param merge.coarse Whether the levels \eqn{j} of a block with
#'   \eqn{2^j} below the chunk size are merged into one chunk; see
#'   \code{wafc_kp_groups()}.
#' @param balanced Whether the chunks are balanced as described above. It
#'   merges the coarse levels whatever \code{merge.coarse} says, and the
#'   object records \code{merge.coarse = TRUE} then.
#' @param free.coarse Whether the coarse levels are left unpenalized, as
#'   described above. It overrides the merge of the coarse levels, which
#'   then are not a chunk, and leaves the finer levels to \code{balanced}.
#' @param nfolds,foldid Folds of the cross-validation.
#' @param ... Passed to \code{\link{wafc_design}}.
#'
#' @return An object of class \code{"wafc_competitor"}.
#'
#' @references Klopp, O. and Pensky, M. (2015). Sparse high-dimensional
#'   varying coefficient model: nonasymptotic minimax study. \emph{The
#'   Annals of Statistics} 43(3), 1273-1299.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' wafc_fit_klopp(d$x, d$u, d$y, J = 3)$extra$block.size
#'
#' @export
wafc_fit_klopp <- function(x, u, y, J = NULL, block.size = NULL,
                           penalize.levels = TRUE,
                           chunk.weights = c("sqrt", "unit"),
                           merge.coarse = FALSE, balanced = FALSE,
                           free.coarse = FALSE,
                           nfolds = 10L, foldid = NULL, ...) {
  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop("method = \"klopp\" needs the package 'grpreg'.", call. = FALSE)
  }
  chunk.weights <- match.arg(chunk.weights)
  merge.coarse <- merge.coarse || balanced
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(u)
  xn <- wafc_names(x, p, "x")
  un <- wafc_names(u, q, "u")
  foldid <- wafc_foldid(foldid, n, nfolds)
  J <- wafc_J_grid(J, n)
  if (is.null(block.size)) block.size <- max(1L, as.integer(ceiling(log(n))))
  block.size <- as.integer(block.size)

  best <- NULL
  conv <- vector("list", length(J))
  for (i in seq_along(J)) {
    Ji <- J[i]
    des <- wafc_design(x, u, J = Ji, ...)
    grp <- wafc_kp_groups(des, block.size, penalize.levels, merge.coarse,
                          balanced, free.coarse)
    Z <- as.matrix(des[["Z"]])
    ## Nothing penalized (free.coarse at a coarse J): least squares, scored
    ## as cv.grpreg scores, in the slot the cv.grpreg object would take.
    if (max(grp) == 0L) {
      cv <- wafc_kp_cv_ols(Z, y, foldid)
      conv[[i]] <- wafc_kp_conv(Ji, NULL)
      if (is.null(best) || cv[["cve"]] < best[["cve"]]) {
        best <- list(cve = cv[["cve"]], J = Ji, design = des, cv = cv,
                     group = grp)
      }
      next
    }
    ## grpreg reads a missing group.multiplier as its default, and there is
    ## no value that means "missing", hence the two calls. The groups are
    ## numbered 1, ..., max(grp) with no gap, which is the order grpreg
    ## expects the multipliers in.
    cv <- if (chunk.weights == "unit") {
      grpreg::cv.grpreg(Z, y, group = grp, penalty = "grLasso",
                        fold = foldid, group.multiplier = rep(1, max(grp)))
    } else {
      grpreg::cv.grpreg(Z, y, group = grp, penalty = "grLasso",
                        fold = foldid)
    }
    conv[[i]] <- wafc_kp_conv(Ji, cv)
    val <- min(cv[["cve"]])
    if (is.null(best) || val < best[["cve"]]) {
      best <- list(cve = val, J = Ji, design = des, cv = cv, group = grp)
    }
  }
  des <- best[["design"]]
  ols <- max(best[["group"]]) == 0L
  b <- if (ols) best[["cv"]][["coef"]] else as.numeric(stats::coef(best[["cv"]]))
  a0 <- b[1L]
  b <- b[-1L]
  cc <- b[des[["unpenalized"]]]
  names(cc) <- xn
  const <- des[["constant"]]
  if (length(const) == 1L) {
    cc[const] <- cc[const] + a0 / as.numeric(des[["Z"]][1L, const])
    a0 <- 0
  }
  nz <- matrix(FALSE, p, q, dimnames = list(xn, un))
  nz_fine <- nz
  pen <- best[["group"]] > 0L
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      idx <- des[["blocks"]][[wafc_block_name(des, l, m)]]
      nz[l, m] <- any(b[idx] != 0)
      nz_fine[l, m] <- any(b[idx[pen[idx]]] != 0)
    }
  }
  g_of <- function(newu) {
    newu <- wafc_as_matrix(newu, "newu")
    d_new <- wafc_design(matrix(1, nrow(newu), p), newu, spec = des)
    g <- vector("list", p * q)
    dim(g) <- c(p, q)
    for (l in seq_len(p)) {
      for (m in seq_len(q)) {
        idx <- des[["blocks"]][[wafc_block_name(des, l, m)]]
        g[[l, m]] <- as.numeric(d_new[["Z"]][, idx, drop = FALSE] %*% b[idx])
      }
    }
    dimnames(g) <- list(xn, un)
    g
  }
  ## The intercept of grpreg stays out of beta when no constant covariate
  ## can carry it; see the same note in wafc_fit_bsgl().
  beta_fun <- function(newu) {
    g <- g_of(newu)
    nr <- length(g[[1L, 1L]])
    out <- matrix(rep(cc, each = nr), nr, p, dimnames = list(NULL, xn))
    for (l in seq_len(p)) {
      for (m in seq_len(q)) out[, l] <- out[, l] + g[[l, m]]
    }
    out
  }
  list(fit = best[["cv"]], cc = cc, beta = beta_fun, g = g_of, blocks = nz,
       fitted = as.numeric(a0 + rowSums(x * beta_fun(u))), intercept = a0,
       xnames = xn, unames = un,
       design = des,
       extra = list(J = best[["J"]], block.size = block.size,
                    ngroups = length(unique(best[["group"]][best[["group"]] > 0L])),
                    penalize.levels = penalize.levels,
                    chunk.weights = chunk.weights,
                    merge.coarse = merge.coarse, balanced = balanced,
                    free.coarse = free.coarse,
                    lambda = best[["cv"]][["lambda.min"]], cve = best[["cve"]],
                    conv = wafc_kp_conv_table(conv, best[["J"]]),
                    blocks.fine = nz_fine,
                    nzero = sum(b[pen & seq_along(b) %in% unlist(des[["blocks"]])] != 0)))
}

## What grpreg returned of the path at one J (step E2.5j, open question
## 41(d) of docs/ESTADO.md), read and not acted upon. grpreg counts
## 'max.iter' over the whole path: when the budget runs out the remaining
## penalty levels are dropped, so the signs of a cut are a total of
## iterations at the budget and a path shorter than asked; a level with
## 'iter' at 'max.iter' is the case its own warning names. cv.grpreg then
## drops every level at which some fold has no finite error, so a cut in a
## fold shows as a cross-validated path shorter than the one of the fit.
## 'nlambda' and 'max.iter' are the ones the path was asked with, the
## defaults of grpreg for wafc_fit_klopp(), which passes neither.
wafc_kp_conv <- function(J, cv,
                         nlambda = wafc_grpreg_default("nlambda"),
                         max.iter = wafc_grpreg_default("max.iter")) {
  if (is.null(cv)) {
    return(data.frame(J = J, nlambda = NA_integer_, nreturned = NA_integer_,
                      ncv = NA_integer_, iter.total = NA_integer_,
                      max.iter = NA_integer_, n.maxiter = NA_integer_,
                      lambda.maxiter = NA_real_, lambda.last = NA_real_,
                      lambda.cv.last = NA_real_, lambda.min = NA_real_))
  }
  fit <- cv[["fit"]]
  mi <- as.integer(max.iter)
  it <- fit[["iter"]]
  at <- it == mi
  data.frame(J = J,
             nlambda = as.integer(nlambda),
             nreturned = length(fit[["lambda"]]),
             ncv = length(cv[["lambda"]]), iter.total = as.integer(sum(it)),
             max.iter = mi, n.maxiter = sum(at),
             lambda.maxiter = if (any(at)) max(fit[["lambda"]][at]) else
               NA_real_,
             lambda.last = min(fit[["lambda"]]),
             lambda.cv.last = min(cv[["lambda"]]),
             lambda.min = cv[["lambda.min"]])
}

## The rows of wafc_kp_conv() over the grid, with the flags read from them:
## the budget of iterations spent ('budget'), the cross-validated path
## shorter than the fit's ('cv.cut'), and whether the selected level is the
## last one of a path that was cut ('cut.at.min'), where the minimum may
## lie past the cut.
wafc_kp_conv_table <- function(rows, J.min) {
  d <- do.call(rbind, rows)
  d[["chosen"]] <- d[["J"]] == J.min
  d[["budget"]] <- d[["iter.total"]] >= d[["max.iter"]]
  d[["cv.cut"]] <- d[["ncv"]] < d[["nreturned"]]
  d[["cut.at.min"]] <- (d[["budget"]] | d[["cv.cut"]] |
                          d[["n.maxiter"]] > 0L) &
    d[["lambda.min"]] == d[["lambda.cv.last"]]
  d
}

## Least squares with an intercept, the fit of the block LASSO when no
## column is penalized, scored as cv.grpreg scores a gaussian fit: the
## squared error of each observation predicted from the folds without it,
## averaged over the n observations. Columns aliased with the intercept
## (the constant covariate) get coefficient zero, which leaves the fitted
## values as they are. 'coef' is (intercept, coefficients of Z), the order
## of coef() on a cv.grpreg object, and lambda.min is NA.
wafc_kp_cv_ols <- function(Z, y, foldid) {
  ls_coef <- function(Zs, ys) {
    b <- stats::lm.fit(cbind(1, Zs), ys)[["coefficients"]]
    b[is.na(b)] <- 0
    unname(b)
  }
  err <- numeric(length(y))
  for (k in unique(foldid)) {
    out <- foldid == k
    b <- ls_coef(Z[!out, , drop = FALSE], y[!out])
    err[out] <- (y[out] - b[1L] - Z[out, , drop = FALSE] %*% b[-1L])^2
  }
  list(cve = mean(err), coef = ls_coef(Z, y), lambda.min = NA_real_)
}

#' Spline with adaptive knots, one knot set per block
#'
#' The competitor that adapts locally without wavelets, in the spirit of
#' Wang, Jiang and Liu (2024): every block \eqn{(\ell, m)} gets its own set
#' of interior knots, chosen from the data rather than fixed in advance, so
#' that a component with a jump can spend its knots where the jump is.
#'
#' What differs from their paper, and has to be said before the pilot
#' reports a number: they search the knot set by an exact dynamic program,
#' and the search here is greedy forward insertion scored by the BIC, run
#' inside a block coordinate descent over the \eqn{pq} terms (their model
#' has a single modulating covariate and no additive structure over
#' \eqn{m}, so their algorithm does not apply verbatim to this design in
#' any case). The greedy search is an upper bound on their criterion, never
#' a lower one, so a win of the WAFC over this column is a win over
#' something no better than their method and a loss is not conclusive. The
#' faithful comparison needs their code and is listed as an open question
#' in the handoff of step E2.4.
#'
#' @param x,u,y The data.
#' @param n.candidate Number of candidate knots per block, placed at
#'   equally spaced quantiles of the modulating covariate.
#' @param max.knots Largest number of interior knots a block may use.
#' @param degree Degree of the spline.
#' @param maxit Number of sweeps of the block coordinate descent.
#' @param tol Relative tolerance on the residual sum of squares between
#'   sweeps.
#'
#' @return An object of class \code{"wafc_competitor"}.
#'
#' @references Wang, X., Jiang, Y. and Liu, Y. (2024). Varying coefficient
#'   model via adaptive spline fitting. \emph{Journal of Computational and
#'   Graphical Statistics} 33(2), 614-624.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "inhomogeneous", seed = 1)
#' wafc_fit_aspline(d$x, d$u, d$y, n.candidate = 16, max.knots = 4)$extra$knots
#'
#' @export
wafc_fit_aspline <- function(x, u, y, n.candidate = 40L, max.knots = 12L,
                             degree = 3L, maxit = 5L, tol = 1e-4) {
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(u)
  xn <- wafc_names(x, p, "x")
  un <- wafc_names(u, q, "u")
  n.candidate <- as.integer(n.candidate)
  max.knots <- as.integer(max.knots)

  rng <- apply(u, 2L, range)
  cand <- lapply(seq_len(q), function(m) {
    pr <- seq_len(n.candidate) / (n.candidate + 1)
    unique(as.numeric(stats::quantile(u[, m], pr, names = FALSE)))
  })

  ## Level terms by least squares, the starting point of the descent.
  qrX <- qr(x)
  cc <- wafc_qr_coef(qrX, y)
  cc[is.na(cc)] <- 0
  cc <- as.numeric(cc)
  names(cc) <- xn
  gfit <- vector("list", p * q)
  dim(gfit) <- c(p, q)
  contrib <- matrix(0, n, p * q)
  fitted <- as.numeric(x %*% cc)
  rss_old <- Inf

  for (it in seq_len(maxit)) {
    b <- 0L
    for (l in seq_len(p)) {
      for (m in seq_len(q)) {
        b <- b + 1L
        r <- y - fitted + contrib[, b]
        sel <- wafc_knot_search(r, x[, l], u[, m], rng[, m], cand[[m]],
                                max.knots, degree, n)
        gfit[[l, m]] <- sel
        new <- x[, l] * sel[["value"]]
        fitted <- fitted - contrib[, b] + new
        contrib[, b] <- new
      }
    }
    ## levels re-fitted on the residual of the components
    cc <- as.numeric(wafc_qr_coef(qrX, y - rowSums(contrib)))
    cc[is.na(cc)] <- 0
    names(cc) <- xn
    fitted <- as.numeric(x %*% cc) + rowSums(contrib)
    rss <- sum((y - fitted)^2)
    if (is.finite(rss_old) && abs(rss_old - rss) <= tol * rss_old) break
    rss_old <- rss
  }

  nkn <- matrix(0L, p, q, dimnames = list(xn, un))
  nz <- matrix(FALSE, p, q, dimnames = list(xn, un))
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      nkn[l, m] <- length(gfit[[l, m]][["knots"]])
      nz[l, m] <- any(gfit[[l, m]][["coef"]] != 0)
    }
  }
  g_of <- function(newu) {
    newu <- wafc_as_matrix(newu, "newu")
    g <- vector("list", p * q)
    dim(g) <- c(p, q)
    for (l in seq_len(p)) {
      for (m in seq_len(q)) {
        s <- gfit[[l, m]]
        g[[l, m]] <- wafc_spline_eval(s, newu[, m], rng[, m], degree)
      }
    }
    dimnames(g) <- list(xn, un)
    g
  }
  beta_fun <- function(newu) {
    g <- g_of(newu)
    nr <- length(g[[1L, 1L]])
    out <- matrix(rep(cc, each = nr), nr, p, dimnames = list(NULL, xn))
    for (l in seq_len(p)) {
      for (m in seq_len(q)) out[, l] <- out[, l] + g[[l, m]]
    }
    out
  }
  list(fit = gfit, cc = cc, beta = beta_fun, g = g_of, blocks = nz,
       fitted = fitted, xnames = xn, unames = un,
       extra = list(knots = nkn, iter = it, rss = rss_old))
}

#' Bayesian trees for varying coefficients (VCBART)
#'
#' The competitor that does not assume the additive structure over the
#' modulating covariates: it estimates each \eqn{\beta_\ell} as a sum of
#' regression trees in the whole vector \eqn{U}, so the scenarios of
#' \code{\link{simulate_wafc}}, which are additive, are inside its model and
#' not inside its inductive bias. Only \eqn{\hat\beta_\ell(u)} is defined
#' for it; there is no \eqn{\hat g_{\ell m}}, and the pilot reads the
#' \code{NULL} in the \code{g} slot and skips the component metric.
#'
#' The engine is \code{VCBART::VCBART_ind}, the version with independent
#' errors, called with one subject per observation. VCBART carries its own
#' intercept function \eqn{\beta_0(U)}, so the constant linear covariate of
#' the design is not passed as a covariate: it is that intercept, divided
#' by the value of the column.
#'
#' @param x,u,y The data.
#' @param burn,nd Length of the burn-in and number of retained draws.
#' @param M Number of trees per coefficient.
#' @param ... Further arguments to \code{\link[VCBART]{VCBART_ind}}.
#'
#' @return An object of class \code{"wafc_competitor"}.
#'
#' @references Deshpande, S. K., Bai, R., Balocchi, C., Starling, J. E. and
#'   Weiss, J. (2026). VCBART: Bayesian trees for varying coefficients.
#'   \emph{Bayesian Analysis} 21(1), 281-308.
#'
#' @examples
#' \dontrun{
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' wafc_fit_vcbart(d$x, d$u, d$y, burn = 100, nd = 100)
#' }
#'
#' @export
wafc_fit_vcbart <- function(x, u, y, burn = 500L, nd = 500L, M = 50L, ...) {
  if (!requireNamespace("VCBART", quietly = TRUE)) {
    stop("method = \"vcbart\" needs the package 'VCBART' ",
         "(install.packages(\"VCBART\")).", call. = FALSE)
  }
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(u)
  xn <- wafc_names(x, p, "x")
  un <- wafc_names(u, q, "u")
  const <- which(vapply(seq_len(p), function(l) diff(range(x[, l])) == 0, TRUE))
  if (length(const) > 1L) {
    stop("More than one constant linear covariate; keep one.", call. = FALSE)
  }
  keep <- setdiff(seq_len(p), const)
  fit <- VCBART::VCBART_ind(Y_train = y, subj_id_train = seq_len(n),
                            ni_train = rep(1L, n),
                            X_train = x[, keep, drop = FALSE],
                            Z_cont_train = u, M = as.integer(M),
                            nd = as.integer(nd), burn = as.integer(burn),
                            verbose = FALSE, ...)
  beta_fun <- function(newu) {
    newu <- wafc_as_matrix(newu, "newu")
    draws <- VCBART::predict_betas(fit, Z_cont = newu, verbose = FALSE)
    bm <- apply(draws, c(2L, 3L), mean)
    out <- matrix(0, nrow(newu), p, dimnames = list(NULL, xn))
    if (length(const) == 1L) out[, const] <- bm[, 1L] / x[1L, const]
    for (i in seq_along(keep)) out[, keep[i]] <- bm[, 1L + i]
    out
  }
  bhat <- beta_fun(u)
  cc <- colMeans(bhat)
  names(cc) <- xn
  list(fit = fit, cc = cc, beta = beta_fun, g = NULL, blocks = NULL,
       fitted = as.numeric(rowSums(x * bhat)), xnames = xn, unames = un,
       extra = list(burn = burn, nd = nd, M = M))
}

#' Linear regression with constant coefficients
#'
#' Ordinary least squares of \eqn{y} on \eqn{x} alone. It is the correctly
#' specified model of the null scenario, where it is the oracle, and the
#' floor of every other scenario: a method that does not beat this column
#' has bought nothing with its basis.
#'
#' @param x,u,y The data. \code{u} is used only for its dimension.
#'
#' @return An object of class \code{"wafc_competitor"}.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "null", seed = 1,
#'                    sigma = 0.6)
#' round(wafc_fit_linear(d$x, d$u, d$y)$cc, 3)
#'
#' @export
wafc_fit_linear <- function(x, u, y) {
  p <- ncol(x)
  q <- ncol(u)
  xn <- wafc_names(x, p, "x")
  un <- wafc_names(u, q, "u")
  cc <- as.numeric(wafc_qr_coef(qr(x), y))
  cc[is.na(cc)] <- 0
  names(cc) <- xn
  beta_fun <- function(newu) {
    newu <- wafc_as_matrix(newu, "newu")
    matrix(rep(cc, each = nrow(newu)), nrow(newu), p,
           dimnames = list(NULL, xn))
  }
  g_of <- function(grid) {
    grid <- wafc_as_matrix(grid, "grid")
    g <- vector("list", p * q)
    dim(g) <- c(p, q)
    for (i in seq_along(g)) g[[i]] <- numeric(nrow(grid))
    dimnames(g) <- list(xn, un)
    g
  }
  list(fit = NULL, cc = cc, beta = beta_fun, g = g_of,
       blocks = matrix(FALSE, p, q, dimnames = list(xn, un)),
       fitted = as.numeric(x %*% cc), xnames = xn, unames = un,
       extra = list())
}

#' The WAFC restricted to the blocks that are really active
#'
#' Not a competitor but the reference the pilot reads "the price of not
#' knowing the structure" against: the same estimator, on the same design,
#' with the columns of the inactive blocks removed before the fit. The
#' penalty level is cross-validated over the same folds, and the resolution
#' over the same grid of \eqn{J}, so the only thing this fit knows and
#' \code{\link{wafc}} does not is which blocks are zero.
#'
#' With \code{penalty = "block"} (step E4.1c, decision D63) the estimator
#' is the block LASSO of \code{\link{wafc}} with the same name: the
#' balanced chunks of size \code{block.size} inside each active block, the
#' weights of \pkg{grpreg} and its tolerance \code{thresh}, with
#' \eqn{(J, \lambda)} chosen as \code{\link{cv.wafc}} chooses them, by the
#' cross-validated error of \code{grpreg::cv.grpreg()} on the folds, and
#' read at \code{lambda.min} without threshold. It is then the oracle of
#' the estimator of decision D44, which the default is not: the default
#' stays \code{"lasso"}, the oracle of the coordinatewise LASSO of decision
#' D3, and fits what it fitted before the option existed. A candidate of
#' the grid of \eqn{J} whose levels, once \code{\link{wafc_design}} has
#' capped them, are those of a candidate already fitted is the same design
#' on the same folds and is not refitted, as in \code{\link{cv.wafc}}; the
#' first of the two is kept, as it was when both were fitted.
#'
#' @param x,u,y The data.
#' @param active Logical \eqn{p} by \eqn{q} matrix of the active blocks,
#'   or a vector of length \eqn{pq} read by column. With no active block the
#'   fit is the least squares of \code{\link{wafc_fit_linear}}.
#' @param J Resolution level or grid of candidates; \code{NULL} uses the
#'   grid of \code{\link{cv.wafc}}.
#' @param penalty \code{"lasso"} (the default), \code{"sglasso"} or
#'   \code{"block"}.
#' @param lambda,nlambda,lambda.min.ratio The penalty path, as in
#'   \code{\link{wafc}}. Giving \code{lambda} explicitly is how the test of
#'   exact recovery reaches the unpenalized end: with a noiseless response
#'   \pkg{glmnet} stops its own path early, because the deviance has stopped
#'   moving.
#' @param nfolds,foldid Folds of the cross-validation.
#' @param block.size The chunk size of \code{penalty = "block"}, as in
#'   \code{\link{wafc}}; \code{NULL} is \code{ceiling(log(n))}. Ignored by
#'   the other penalties.
#' @param thresh Convergence tolerance of the engine, as in
#'   \code{\link{wafc}}; \code{NULL} is the default of \code{\link{wafc}}
#'   for the penalty.
#' @param ... Passed to \code{\link{wafc_design}}.
#'
#' @return An object of class \code{"wafc_competitor"}. Its \code{extra}
#'   has the selected \code{J}, \code{lambda} and \code{cvm} and the
#'   \code{penalty}; with \code{penalty = "block"} also the
#'   \code{block.size}, the number \code{nzero} of nonzero penalized
#'   coefficients, and \code{conv}, the path of \pkg{grpreg} at every
#'   \eqn{J} of the grid as \code{\link{wafc_fit_klopp}} records it, and
#'   the object keeps in \code{intercept} the intercept of \pkg{grpreg}
#'   when no constant covariate carries it.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' wafc_fit_oracle(d$x, d$u, d$y, active = nzchar(d$structure), J = 3)$blocks
#' wafc_fit_oracle(d$x, d$u, d$y, active = nzchar(d$structure), J = 3,
#'                 penalty = "block")$extra$nzero
#'
#' @export
wafc_fit_oracle <- function(x, u, y, active = NULL, J = NULL,
                            penalty = c("lasso", "sglasso", "block"),
                            lambda = NULL, nlambda = 100L,
                            lambda.min.ratio = NULL, nfolds = 10L,
                            foldid = NULL, block.size = NULL, thresh = NULL,
                            ...) {
  penalty <- match.arg(penalty)
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(u)
  xn <- wafc_names(x, p, "x")
  un <- wafc_names(u, q, "u")
  active <- wafc_oracle_active(active, p, q)
  foldid <- wafc_foldid(foldid, n, nfolds)
  J <- wafc_J_grid(J, n)
  if (!any(active)) {
    out <- wafc_fit_linear(x, u, y)
    out[["extra"]] <- list(J = NA_integer_, lambda = NA_real_,
                           note = "no active block: least squares on x")
    return(out)
  }
  block <- penalty == "block"
  ## the default tolerance is the one wafc() has for the penalty, read from
  ## its formals so that the two cannot drift apart
  thr <- if (is.null(thresh)) {
    eval(formals(wafc)[["thresh"]], list(penalty = penalty))
  } else thresh
  fit_args <- list(thresh = thr)
  if (!is.null(block.size)) fit_args[["block.size"]] <- block.size
  best <- NULL
  seen <- character(0)
  conv <- vector("list", length(J))
  for (i in seq_along(J)) {
    Ji <- J[i]
    des <- wafc_design(x, u, J = Ji, ...)
    key <- paste(des[["J"]], collapse = ",")
    prev <- match(key, seen)
    seen[i] <- key
    if (!is.na(prev)) {
      if (block) {
        conv[[i]] <- conv[[prev]]
        conv[[i]][["J"]] <- Ji
      }
      next
    }
    sub <- wafc_design_keep(des, active)
    if (block) {
      zz <- wafc_cv_block(sub, y, foldid, lambda, nlambda, lambda.min.ratio,
                          "mse", fit_args, Ji)
      full <- zz[["fit"]]
      z <- zz[["cv"]]
      conv[[i]] <- z[["conv"]]
    } else {
      full <- wafc(design = sub, y = y, penalty = penalty, lambda = lambda,
                   nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                   thresh = thr)
      z <- wafc_cv_design(sub, y, full, foldid, function(e) e^2, penalty)
    }
    if (is.null(best) || z[["cvm.min"]] < best[["cvm"]]) {
      best <- list(cvm = z[["cvm.min"]], J = Ji, fit = full,
                   lambda = z[["lambda.min"]], design = sub,
                   nzero = z[["nzero.min"]])
    }
  }
  fit <- best[["fit"]]
  des <- best[["design"]]
  lam <- best[["lambda"]]
  cf <- coef.wafc(fit, s = lam)
  b <- cf[-1L, 1L]
  ## the block LASSO always has the intercept of grpreg, which coef.wafc()
  ## has already read into the constant covariate when there is one; the
  ## LASSO fits one only then (wafc()), so it has none to keep
  a0 <- if (block) as.numeric(cf[1L, 1L]) else 0
  cc <- b[des[["unpenalized"]]]
  names(cc) <- xn
  nz <- matrix(FALSE, p, q, dimnames = list(xn, un))
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      nm <- wafc_block_name(des, l, m)
      if (!is.null(des[["blocks"]][[nm]])) {
        nz[l, m] <- any(b[des[["blocks"]][[nm]]] != 0)
      }
    }
  }
  g_of <- function(newu) {
    newu <- wafc_as_matrix(newu, "newu")
    d_new <- wafc_design(matrix(1, nrow(newu), p), newu,
                         spec = des[["parent"]])
    g <- vector("list", p * q)
    dim(g) <- c(p, q)
    for (l in seq_len(p)) {
      for (m in seq_len(q)) {
        nm <- wafc_block_name(des, l, m)
        idx <- des[["blocks"]][[nm]]
        g[[l, m]] <- if (is.null(idx)) numeric(nrow(newu)) else {
          as.numeric(d_new[["Z"]][, des[["parent.blocks"]][[nm]],
                                  drop = FALSE] %*% b[idx])
        }
      }
    }
    dimnames(g) <- list(xn, un)
    g
  }
  beta_fun <- function(newu) {
    g <- g_of(newu)
    nr <- length(g[[1L, 1L]])
    out <- matrix(rep(cc, each = nr), nr, p, dimnames = list(NULL, xn))
    for (l in seq_len(p)) {
      for (m in seq_len(q)) out[, l] <- out[, l] + g[[l, m]]
    }
    out
  }
  extra <- list(J = best[["J"]], lambda = lam, cvm = best[["cvm"]],
                penalty = penalty)
  if (!block) {
    return(list(fit = fit, cc = cc, beta = beta_fun, g = g_of, blocks = nz,
                fitted = as.numeric(rowSums(x * beta_fun(u))), xnames = xn,
                unames = un, extra = extra))
  }
  extra[["block.size"]] <- fit[["group"]][["block.size"]]
  extra[["nzero"]] <- as.integer(best[["nzero"]])
  extra[["conv"]] <- wafc_kp_conv_table(conv, best[["J"]])
  list(fit = fit, cc = cc, beta = beta_fun, g = g_of, blocks = nz,
       fitted = as.numeric(a0 + rowSums(x * beta_fun(u))), intercept = a0,
       xnames = xn, unames = un, extra = extra)
}

## ---------------------------------------------------------------------------
## Internals
## ---------------------------------------------------------------------------

## Every formal argument of every fitter, which is the list of names that
## are known to belong to some method and may therefore be dropped when
## they reach another one.
wafc_fitter_args <- function() {
  fs <- list(wafc_fit_gam, wafc_fit_bsgl, wafc_fit_klopp, wafc_fit_aspline,
             wafc_fit_vcbart, wafc_fit_linear, wafc_fit_oracle)
  setdiff(unique(unlist(lapply(fs, function(f) names(formals(f))))),
          c("...", "x", "u", "y"))
}

## The grid of basis dimensions searched by wafc_fit_gam() when k.select is
## not "none" (step E2.5h, decision D41): the one of Ruppert (2002, section
## 3) without its 120, one grid for the three criteria.
wafc_k_grid <- c(5L, 10L, 20L, 40L, 80L)

## The candidates of a search over k: each value of the grid, common to the
## smooths, truncated modulator by modulator at the number of distinct
## values minus one (the rule of wafc_k_matched()) and at 3 from below, and
## the vectors the truncation makes equal kept once, at the first value of
## the grid that produced them. The grid kept is the attribute "grid".
wafc_k_candidates <- function(u, k) {
  u <- wafc_as_matrix(u, "u")
  if (!is.numeric(k) || length(k) == 0L || any(!is.finite(k)) ||
      any(k != round(k)) || any(k < 3)) {
    stop("'k' must contain integer values of at least 3.", call. = FALSE)
  }
  k <- sort(unique(as.integer(k)))
  ndist <- vapply(seq_len(ncol(u)), function(m) length(unique(u[, m])), 0L)
  cand <- lapply(k, function(kk) as.integer(pmax(3L, pmin(kk, ndist - 1L))))
  keep <- !duplicated(vapply(cand, paste, "", collapse = ","))
  out <- cand[keep]
  attr(out, "grid") <- k[keep]
  out
}

## The restricted negative log-likelihood of a Gaussian fit of mgcv, at its
## smoothing parameters and profiled over the scale (step E2.5h). With the
## penalty S = sum_j sp_j S_j, block diagonal over the smooths, and Mp the
## number of unpenalized coefficients (the parametric ones, plus the null
## spaces left unpenalized when select = FALSE),
##
##   V = (RSS + b'Sb) / (2 s2) + (n - Mp)/2 log(2 pi s2)
##       - log|S|_+ / 2 + log|X'X + S| / 2,
##
## with b = (X'X + S)^{-1} X'y and s2 = (RSS + b'Sb) / (n - Mp), which is the
## REML score of Wood (2011) and of mgcv::gam(method = "REML") at its own
## optimum. It is the Gaussian restricted likelihood of y with the
## parametric terms as fixed effects, so it is comparable between fits that
## share those terms; the test of this step checks it against that
## likelihood written with the n by n covariance. b and s2 are recomputed
## from the model matrix rather than read from the fit, because bam with
## discrete = TRUE reports a score of its own (see wafc_fit_gam()).
wafc_gam_reml <- function(fit, y) {
  X <- stats::predict(fit, type = "lpmatrix")
  n <- nrow(X)
  P <- ncol(X)
  S <- matrix(0, P, P)
  ldS <- 0
  Mp <- fit[["nsdf"]]
  j <- 0L
  for (sm in fit[["smooth"]]) {
    ii <- sm[["first.para"]]:sm[["last.para"]]
    Sb <- matrix(0, length(ii), length(ii))
    for (Sj in sm[["S"]]) {
      j <- j + 1L
      Sb <- Sb + fit[["sp"]][j] * Sj
    }
    S[ii, ii] <- Sb
    r <- length(ii) - sm[["null.space.dim"]]
    ev <- eigen(Sb, symmetric = TRUE, only.values = TRUE)[["values"]]
    if (r > 0L && ev[r] <= 0) {
      stop("the penalty of a smooth has rank below the one mgcv declares.",
           call. = FALSE)
    }
    ldS <- ldS + sum(log(ev[seq_len(r)]))
    Mp <- Mp + sm[["null.space.dim"]]
  }
  R <- chol(crossprod(X) + S)
  b <- backsolve(R, forwardsolve(t(R), crossprod(X, y)))
  e <- as.numeric(y - X %*% b)
  pen <- sum(e^2) + sum(b * (S %*% b))
  s2 <- pen / (n - Mp)
  out <- pen / (2 * s2) + (n - Mp) / 2 * log(2 * pi * s2) - ldS / 2 +
    sum(log(diag(R)))
  attr(out, "sigma2") <- s2
  out
}

## The methods whose '...' is passed to wafc_design(); the '...' of the
## others goes to an engine (mgcv, grpreg, VCBART) that does not know the
## arguments of the basis.
wafc_design_fitters <- c("klopp", "oracle")

## The arguments of wafc_design() that describe the basis and the rescaling.
wafc_design_args <- function() {
  setdiff(names(formals(wafc_design)), c("x", "u", "J", "spec"))
}

## qr.coef on a matrix or a vector, with the NA of a rank-deficient fit
## turned into zero, which is what every caller here wants.
wafc_qr_coef <- function(qrx, v) {
  cf <- qr.coef(qrx, v)
  cf[is.na(cf)] <- 0
  cf
}

## Data frame with one column per covariate, named as the design names, for
## the formula interface of mgcv.
wafc_frame <- function(x, u, xn, un) {
  d <- as.data.frame(x)
  names(d) <- xn
  du <- as.data.frame(u)
  names(du) <- un
  cbind(d, du)
}

## Column of predict(type = "terms") of the b-th smooth. mgcv names the
## columns "s(u1):x2" and keeps the order of the formula, which is the
## lexicographic order of D12, so the position is enough and the name is
## only checked.
wafc_gam_term <- function(nms, b) {
  sm <- grep("^s\\(", nms)
  if (length(sm) < b) {
    stop("mgcv returned ", length(sm), " smooth term(s); expected at least ",
         b, ".", call. = FALSE)
  }
  sm[b]
}

## B-spline specification of every modulating covariate: the knots and the
## boundary, fixed on the training sample, plus the column means that
## centre the basis (the spline analogue of dropping phi_{00}).
wafc_bs_spec <- function(u, df, un) {
  q <- ncol(u)
  sp <- vector("list", q)
  for (m in seq_len(q)) {
    rg <- range(u[, m])
    nk <- df - 4L
    kn <- if (nk > 0L) {
      as.numeric(stats::quantile(u[, m], seq_len(nk) / (nk + 1), names = FALSE))
    } else numeric(0)
    B <- splines::bs(u[, m], knots = kn, degree = 3L, Boundary.knots = rg,
                     intercept = TRUE)
    sp[[m]] <- list(knots = kn, boundary = rg, centre = colMeans(B))
  }
  names(sp) <- un
  sp
}

## The centred basis of one modulating covariate, with the first column
## dropped so that the block has df - 1 columns and does not repeat the
## unpenalized level term.
wafc_bs_eval <- function(v, sp, m) {
  s <- sp[[m]]
  v <- pmin(pmax(v, s[["boundary"]][1L]), s[["boundary"]][2L])
  B <- splines::bs(v, knots = s[["knots"]], degree = 3L,
                   Boundary.knots = s[["boundary"]], intercept = TRUE)
  B <- sweep(as.matrix(B), 2L, s[["centre"]], "-")
  B[, -1L, drop = FALSE]
}

## The design of the B-spline sieve, in the column order of D12.
wafc_bs_design <- function(x, u, sp) {
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(u)
  B <- lapply(seq_len(q), function(m) wafc_bs_eval(u[, m], sp, m))
  cols <- vector("list", 1L + p * q)
  cols[[1L]] <- x
  b <- 0L
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      b <- b + 1L
      cols[[1L + b]] <- B[[m]] * x[, l]
    }
  }
  do.call(cbind, cols)
}

## Group vector of the block LASSO of Klopp and Pensky: zero on a term that
## is not penalized, and otherwise one group per chunk. Inside a block the
## chunks never straddle a resolution level, because the columns are
## ordered by increasing j and then k (D12): a level with at most
## 'block.size' translates is a chunk, a finer level is cut into
## consecutive pieces of that size.
##
## With merge.coarse = TRUE the levels with 2^j below 'block.size' are one
## chunk instead of one chunk each (step E2.5c, open question 34): with the
## default ceiling(log n) of 6 or 7 that is levels 0 to 2, seven columns,
## and it removes the chunk of size one of level 0. The finer levels are cut
## as before, so the full chunks, and the short piece at the end of a level
## that the chunk size does not divide, are unchanged.
##
## With balanced = TRUE (step E2.5e) the coarse levels are merged as above,
## whatever merge.coarse says, and that short piece is absorbed into the
## chunk before it in the same level, so every chunk of a finer level has
## between 'block.size' and 2 'block.size' - 1 columns. A level with 2^j
## at least 'block.size' always has one full chunk to absorb into. The
## sizes of the coarse chunk are in the note on wafc_fit_klopp().
##
## With free.coarse = TRUE (step E2.5f) the coarse levels are not a chunk
## but group 0, unpenalized, whatever merge.coarse says; the finer levels
## are cut as 'balanced' says. A block with no level of 2^j at least
## 'block.size' is then wholly in group 0.
wafc_kp_groups <- function(design, block.size, penalize.levels,
                           merge.coarse = FALSE, balanced = FALSE,
                           free.coarse = FALSE) {
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
    if (merge.coarse || balanced || free.coarse) {
      coarse <- levs[2^levs < block.size]
      if (length(coarse) > 0L) {
        pos <- as.integer(sum(2^coarse))
        if (!free.coarse) {
          g <- g + 1L
          grp[idx[seq_len(pos)]] <- g
        }
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

## The 'active' of wafc_fit_oracle(): a p by q logical matrix, or a vector of
## length pq read by column, with no NA. An empty or absent one is an error
## and not "no active block", which is the matrix of FALSE.
wafc_oracle_active <- function(active, p, q) {
  if (is.null(active) || length(active) == 0L) {
    stop("method = \"oracle\" needs 'active', the p by q matrix of the ",
         "blocks that are really active.", call. = FALSE)
  }
  if (!(is.logical(active) || is.numeric(active))) {
    stop("'active' must be a logical p by q matrix.", call. = FALSE)
  }
  if (is.matrix(active) && !identical(dim(active), c(p, q))) {
    stop("'active' must be ", p, " by ", q, " (p by q); it is ",
         nrow(active), " by ", ncol(active), ".", call. = FALSE)
  }
  if (length(active) != p * q) {
    stop("'active' must have p q = ", p * q, " entries; it has ",
         length(active), ".", call. = FALSE)
  }
  if (anyNA(active)) stop("'active' must not contain NA.", call. = FALSE)
  matrix(as.logical(active), p, q)
}

## A design with the columns of the inactive blocks removed, used by the
## oracle of the structure. The object keeps the parent design, so that a
## fit on it can still be evaluated on new data through the spec route.
wafc_design_keep <- function(design, active) {
  p <- design[["p"]]
  q <- design[["q"]]
  keep <- design[["unpenalized"]]
  blocks <- list()
  parent <- list()
  pos <- length(keep)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      if (!active[l, m]) next
      nm <- wafc_block_name(design, l, m)
      idx <- design[["blocks"]][[nm]]
      keep <- c(keep, idx)
      blocks[[nm]] <- pos + seq_along(idx)
      parent[[nm]] <- idx
      pos <- pos + length(idx)
    }
  }
  out <- design
  out[["Z"]] <- design[["Z"]][, keep, drop = FALSE]
  out[["penalty.factor"]] <- design[["penalty.factor"]][keep]
  out[["blocks"]] <- blocks
  out[["nvars"]] <- length(keep)
  out[["parent"]] <- design
  out[["parent.blocks"]] <- parent
  out[["keep"]] <- keep
  out
}

## Forward knot insertion for one block of the adaptive spline, with the
## number of knots chosen at the end and not on the way. The path inserts,
## at each step, the candidate that most reduces the residual sum of
## squares of the weighted fit r ~ w * s(v), up to 'max.knots'; the BIC
## then picks a point of that path, the empty set of interior knots
## included, which is how a null block ends as a plain polynomial.
##
## Choosing the size at the end and not by a stopping rule is not a detail:
## with a stopping rule the search halts on a component whose first knot
## alone does not pay, which is exactly a spiky component such as 'bumps',
## the one the comparison is about.
##
## The candidates of a step are scored at once, not refitted one by one.
## Adding the knot k to a spline of the given degree adds the truncated
## power (v - k)_+^degree to its span, so the residual sum of squares of the
## fit with that knot is the current one minus (c'e)^2 / c'c, where e is the
## current residual and c is the weighted truncated power projected off the
## current span. That is the same number the refit gives, one projection per
## step instead of one B-spline basis and one least squares fit per
## candidate: the search chose the same knots as the refit in all eighteen
## cases compared (three scenarios, n = 250 and 1000, three seeds), with
## identical fitted values, and was seven to eleven times faster. A candidate
## whose projection is at the level of rounding (no observation between two
## knots, where the refit is rank deficient and the sum of squares does not
## move) is given no gain; the threshold, 1e-24 of the squared norm of the
## column, is far below the legitimate small projections (a knot between
## two close knots adds little, and 1e-10 already discarded such knots in
## the comparison) and far above rounding, which is of order 1e-32.
wafc_knot_search <- function(r, w, v, rng, cand, max.knots, degree, n) {
  vc <- pmin(pmax(v, rng[1L]), rng[2L])
  basis <- function(kn) {
    as.matrix(splines::bs(vc, knots = kn, degree = degree,
                          Boundary.knots = rng, intercept = TRUE))
  }
  score <- function(kn) {
    fitk <- stats::.lm.fit(basis(kn) * w, r)
    rss <- sum(fitk[["residuals"]]^2)
    df <- fitk[["rank"]]
    list(bic = n * log(max(rss, .Machine[["double.eps"]]) / n) + df * log(n),
         coef = wafc_zero_na(fitk[["coefficients"]]), rss = rss, knots = kn)
  }
  tpow <- w * pmax(outer(vc, cand, "-"), 0)^degree
  tnorm <- colSums(tpow^2)
  kn <- numeric(0)
  path <- list(score(kn))
  while (length(kn) < max.knots) {
    pool <- which(!(cand %in% kn))
    if (length(pool) == 0L) break
    Q <- qr.Q(qr(basis(kn) * w))
    e <- r - as.numeric(Q %*% crossprod(Q, r))
    Cp <- tpow[, pool, drop = FALSE]
    Cp <- Cp - Q %*% crossprod(Q, Cp)
    nrm <- colSums(Cp^2)
    gain <- ifelse(nrm > 1e-24 * tnorm[pool],
                   as.numeric(crossprod(Cp, e))^2 / nrm, 0)
    ## which.max takes the first of tied candidates, as the refit loop did
    best <- score(sort(c(kn, cand[pool[which.max(gain)]])))
    kn <- best[["knots"]]
    path[[length(path) + 1L]] <- best
  }
  cur <- path[[which.min(vapply(path, `[[`, 0, "bic"))]]
  value <- as.numeric(basis(cur[["knots"]]) %*% cur[["coef"]])
  list(knots = cur[["knots"]], coef = cur[["coef"]], value = value,
       bic = cur[["bic"]])
}

wafc_zero_na <- function(v) {
  v[is.na(v)] <- 0
  v
}

## Evaluation of a block of the adaptive spline at new points.
wafc_spline_eval <- function(s, v, rng, degree) {
  B <- splines::bs(pmin(pmax(v, rng[1L]), rng[2L]), knots = s[["knots"]],
                   degree = degree, Boundary.knots = rng, intercept = TRUE)
  as.numeric(as.matrix(B) %*% s[["coef"]])
}
