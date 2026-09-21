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
##
## Every fitter returns an object of class "wafc_competitor" with the same
## four accessors, so the pilot loops over methods and not over special
## cases: beta(u) gives the n by p matrix of functional coefficients,
## g(grid) gives the p by q list of additive components (NULL when the
## method does not decompose, which is the case of vcbart), blocks gives
## the p by q matrix of selected blocks (NULL when the method selects
## nothing), and predict() gives rowSums(newx * beta(newu)).
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
#' @param u Matrix of modulating covariates, \eqn{n} by \eqn{q}.
#' @param y Numeric response of length \eqn{n}.
#' @param active Logical \eqn{p} by \eqn{q} matrix of the blocks that are
#'   really active. Required by \code{method = "oracle"} and ignored by
#'   every other method.
#' @param ... Passed to the fitter of the method.
#'
#' @return An object of class \code{"wafc_competitor"}.
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
  ## The pilot loops over the methods with one list of arguments, so an
  ## argument meant for another method is dropped here rather than raising
  ## an error. A fitter with a '...' of its own keeps receiving everything.
  args <- list(x = x, u = u, y = y)
  extra <- list(...)
  extra[["active"]] <- NULL
  ## The pilot calls every method with one list of arguments, so an
  ## argument that belongs to another fitter is dropped here instead of
  ## raising an error. A fitter with a '...' of its own keeps everything
  ## else, because that is how wafc_design() arguments such as 'rescale'
  ## reach it; only a name that is some other fitter's argument is removed.
  own <- setdiff(names(formals(fitter)), "...")
  if ("..." %in% names(formals(fitter))) {
    extra <- extra[!(names(extra) %in% setdiff(wafc_fitter_args(), own))]
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
  as.numeric(rowSums(newx * object[["beta"]](newu)))
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
#' @param object A \code{"wafc"} or \code{"wafc_competitor"} object.
#' @param grid Matrix with \eqn{q} columns at which the components are
#'   evaluated, on the original scale of the modulating covariates.
#' @param s For a \code{"wafc"} object, the penalty level.
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
    stop("'object' must be a \"wafc\" or a \"wafc_competitor\" fit.",
         call. = FALSE)
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
#' @param x,u,y The data, as in \code{\link{wafc_competitor}}.
#' @param k Basis dimension of each smooth: one value for every smooth, or
#'   one value per modulating covariate.
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
#'   \code{k} and the \code{engine}.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' wafc_fit_gam(d$x, d$u, d$y)$blocks
#'
#' @export
wafc_fit_gam <- function(x, u, y, k = 10L, select = TRUE, edf.tol = 0.1,
                         engine = c("gam", "bam"), ...) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("method = \"gam\" needs the package 'mgcv'.", call. = FALSE)
  }
  engine <- match.arg(engine)
  p <- ncol(x)
  q <- ncol(u)
  k <- wafc_recycle(k, q, "k")
  if (any(!is.finite(k)) || any(k != round(k)) || any(k < 3)) {
    stop("'k' must contain integer values of at least 3.", call. = FALSE)
  }
  k <- as.integer(k)
  xn <- wafc_names(x, p, "x")
  un <- wafc_names(u, q, "u")
  dat <- wafc_frame(x, u, xn, un)
  dat[["y"]] <- y
  ## One smooth per block, in the lexicographic order of D12, plus the p
  ## parametric terms that carry the levels. The 'by' variable is what
  ## predict(type = "terms") multiplies the smooth by, so evaluating the
  ## terms with every covariate set to one is what turns the term of the
  ## block (l, m) into the component g_{lm} itself.
  terms <- character(0)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      terms <- c(terms, sprintf("s(%s, k = %d, by = %s)", un[m], k[m], xn[l]))
    }
  }
  fo <- stats::as.formula(paste("y ~ 0 +", paste(c(xn, terms), collapse = " + ")))
  fit <- if (engine == "gam") {
    mgcv::gam(fo, data = dat, method = "REML", select = select, ...)
  } else {
    mgcv::bam(fo, data = dat, method = "fREML", discrete = TRUE,
              select = select, ...)
  }
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
                    edf.tol = edf.tol, k = k, engine = engine))
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
#' functions is chosen from the same grid the WAFC uses for \eqn{2^J}, so
#' the two sieves are compared at comparable dimensions rather than at a
#' number pulled from the air.
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
#' @return An object of class \code{"wafc_competitor"}.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' wafc_fit_bsgl(d$x, d$u, d$y, df = 8)$blocks
#'
#' @export
wafc_fit_bsgl <- function(x, u, y, df = c(4L, 8L, 16L), nfolds = 10L,
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

  best <- NULL
  for (d in df) {
    sp <- wafc_bs_spec(u, d, un)
    Z <- wafc_bs_design(x, u, sp)
    grp <- c(rep(0L, p), rep(seq_len(p * q), each = d - 1L))
    cv <- grpreg::cv.grpreg(Z, y, group = grp, penalty = penalty,
                            fold = foldid, ...)
    val <- min(cv[["cve"]])
    if (is.null(best) || val < best[["cve"]]) {
      best <- list(cve = val, df = d, spec = sp, cv = cv, group = grp,
                   Z = Z)
    }
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
  beta_fun <- function(newu) {
    g <- g_of(newu)
    nr <- length(g[[1L, 1L]])
    out <- matrix(rep(cc, each = nr), nr, p, dimnames = list(NULL, xn))
    out[, 1L] <- out[, 1L] + a0
    for (l in seq_len(p)) {
      for (m in seq_len(q)) out[, l] <- out[, l] + g[[l, m]]
    }
    out
  }
  list(fit = best[["cv"]], cc = cc, beta = beta_fun, g = g_of, blocks = nz,
       fitted = as.numeric(rowSums(x * beta_fun(u))), xnames = xn, unames = un,
       extra = list(df = d, lambda = best[["cv"]][["lambda.min"]],
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
#' The chunks respect the order of decision D12, so a chunk never straddles
#' two resolution levels: the levels \eqn{j} with \eqn{2^j} at most the
#' chunk size are chunks of their own, and a finer level is cut into
#' consecutive pieces of the chunk size. That is the reading of "blocks of
#' size about \eqn{\log n} inside each functional coefficient" that keeps
#' the block a set of neighbouring translates at one scale.
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
                           penalize.levels = TRUE, nfolds = 10L,
                           foldid = NULL, ...) {
  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop("method = \"klopp\" needs the package 'grpreg'.", call. = FALSE)
  }
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
  for (Ji in J) {
    des <- wafc_design(x, u, J = Ji, ...)
    grp <- wafc_kp_groups(des, block.size, penalize.levels)
    Z <- as.matrix(des[["Z"]])
    cv <- grpreg::cv.grpreg(Z, y, group = grp, penalty = "grLasso",
                            fold = foldid)
    val <- min(cv[["cve"]])
    if (is.null(best) || val < best[["cve"]]) {
      best <- list(cve = val, J = Ji, design = des, cv = cv, group = grp)
    }
  }
  des <- best[["design"]]
  b <- as.numeric(stats::coef(best[["cv"]]))
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
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      nz[l, m] <- any(b[des[["blocks"]][[wafc_block_name(des, l, m)]]] != 0)
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
  beta_fun <- function(newu) {
    g <- g_of(newu)
    nr <- length(g[[1L, 1L]])
    out <- matrix(rep(cc, each = nr), nr, p, dimnames = list(NULL, xn))
    out[, 1L] <- out[, 1L] + a0
    for (l in seq_len(p)) {
      for (m in seq_len(q)) out[, l] <- out[, l] + g[[l, m]]
    }
    out
  }
  list(fit = best[["cv"]], cc = cc, beta = beta_fun, g = g_of, blocks = nz,
       fitted = as.numeric(rowSums(x * beta_fun(u))), xnames = xn, unames = un,
       design = des,
       extra = list(J = best[["J"]], block.size = block.size,
                    ngroups = length(unique(best[["group"]][best[["group"]] > 0L])),
                    penalize.levels = penalize.levels,
                    lambda = best[["cv"]][["lambda.min"]], cve = best[["cve"]]))
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
#' @param x,u,y The data.
#' @param active Logical \eqn{p} by \eqn{q} matrix of the active blocks.
#' @param J Resolution level or grid of candidates; \code{NULL} uses the
#'   grid of \code{\link{cv.wafc}}.
#' @param penalty \code{"lasso"} or \code{"sglasso"}.
#' @param lambda,nlambda,lambda.min.ratio The penalty path, as in
#'   \code{\link{wafc}}. Giving \code{lambda} explicitly is how the test of
#'   exact recovery reaches the unpenalized end: with a noiseless response
#'   \pkg{glmnet} stops its own path early, because the deviance has stopped
#'   moving.
#' @param nfolds,foldid Folds of the cross-validation.
#' @param ... Passed to \code{\link{wafc_design}}.
#'
#' @return An object of class \code{"wafc_competitor"}.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' wafc_fit_oracle(d$x, d$u, d$y, active = nzchar(d$structure), J = 3)$blocks
#'
#' @export
wafc_fit_oracle <- function(x, u, y, active = NULL, J = NULL,
                            penalty = c("lasso", "sglasso"), lambda = NULL,
                            nlambda = 100L, lambda.min.ratio = NULL,
                            nfolds = 10L, foldid = NULL, ...) {
  penalty <- match.arg(penalty)
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(u)
  xn <- wafc_names(x, p, "x")
  un <- wafc_names(u, q, "u")
  if (is.null(active)) {
    stop("method = \"oracle\" needs 'active', the p by q matrix of the ",
         "blocks that are really active.", call. = FALSE)
  }
  active <- matrix(as.logical(active), p, q)
  foldid <- wafc_foldid(foldid, n, nfolds)
  J <- wafc_J_grid(J, n)
  if (!any(active)) {
    out <- wafc_fit_linear(x, u, y)
    out[["extra"]] <- list(J = NA_integer_, lambda = NA_real_,
                           note = "no active block: least squares on x")
    return(out)
  }
  best <- NULL
  for (Ji in J) {
    des <- wafc_design(x, u, J = Ji, ...)
    sub <- wafc_design_keep(des, active)
    full <- wafc(design = sub, y = y, penalty = penalty, lambda = lambda,
                 nlambda = nlambda, lambda.min.ratio = lambda.min.ratio)
    z <- wafc_cv_design(sub, y, full, foldid, function(e) e^2, penalty)
    if (is.null(best) || z[["cvm.min"]] < best[["cvm"]]) {
      best <- list(cvm = z[["cvm.min"]], J = Ji, fit = full,
                   lambda = z[["lambda.min"]], design = sub)
    }
  }
  fit <- best[["fit"]]
  des <- best[["design"]]
  lam <- best[["lambda"]]
  b <- coef.wafc(fit, s = lam)[-1L, 1L]
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
  list(fit = fit, cc = cc, beta = beta_fun, g = g_of, blocks = nz,
       fitted = as.numeric(rowSums(x * beta_fun(u))), xnames = xn, unames = un,
       extra = list(J = best[["J"]], lambda = lam, cvm = best[["cvm"]],
                    penalty = penalty))
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
wafc_kp_groups <- function(design, block.size, penalize.levels) {
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
    for (j in j0:(Jm - 1L)) {
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
wafc_knot_search <- function(r, w, v, rng, cand, max.knots, degree, n) {
  basis <- function(kn, at = v) {
    B <- splines::bs(pmin(pmax(at, rng[1L]), rng[2L]), knots = kn,
                     degree = degree, Boundary.knots = rng, intercept = TRUE)
    as.matrix(B)
  }
  score <- function(kn) {
    fitk <- stats::.lm.fit(basis(kn) * w, r)
    rss <- sum(fitk[["residuals"]]^2)
    df <- fitk[["rank"]]
    list(bic = n * log(max(rss, .Machine[["double.eps"]]) / n) + df * log(n),
         coef = wafc_zero_na(fitk[["coefficients"]]), rss = rss, knots = kn)
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
