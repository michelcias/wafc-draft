## wafc/R/threshold.R -- estimation followed by a threshold on the blocks
## (step E2.5g): the estimator of Corollary 8 of
## derivations/06-selecao-limiar.tex (decision D32), and the calibration of
## its threshold t_n (open question 11 of docs/ESTADO.md).
##
## A fit, of the WAFC or of the block LASSO of Klopp and Pensky in any of
## its forms, is read block by block: the statistic of a block (l, m) is the
## Euclidean norm of its wavelet coefficients, the one wafc_blocks() reports,
## which is the L2 norm of g_{lm} exactly when eps = 0 and rescale = FALSE,
## and that norm times a factor sqrt(scale_m) common to the blocks of U_m
## when rescale = TRUE (step E1.7c). The thresholded estimator keeps the
## blocks whose norm is strictly above t and zeroes the others, which is the
## set S-hat(t) = {(l, m) : N-hat_{lm} > t} of Lemma 13 with its strict
## inequality; nothing else changes, so the levels and the coefficients of
## the blocks kept are the ones of the fit. That is the estimator of
## Corollary 8. A refit by least squares of what is kept is offered as a
## secondary option, in the spirit of the thresholded and refitted LASSO of
## van de Geer, Buhlmann and Zhou (2011): on the nonzero coefficients of the
## blocks kept ("support", their reading, carried to blocks), or on every
## column of those blocks ("block").
##
## Three rules choose t, and only the first two are estimators:
##
##   "max"     t = c max_{lm} N-hat_{lm}, invariant to the scale of the
##             response, with c = 0.15, the value step E1.7c started from;
##   "cv"      t by cross-validation on the folds of the fit: each fold is
##             refitted at the (J, lambda) already chosen, with no new
##             search, thresholded at every candidate t and scored on the
##             observations it left out;
##   "oracle"  the t with the smallest prediction error against the true
##             regression function on a test sample, which no estimator
##             can see; it is the reference the other two are read
##             against, and is labelled as such.
##
## The candidates of "cv" and "oracle" are exhaustive and not a grid. A
## thresholded fit is constant in t between two consecutive norms of its
## blocks, so the cross-validated error is constant between two consecutive
## norms of the fold fits, and the error of the oracle between two
## consecutive norms of the fit; each interval is scored once, and the t
## returned is its midpoint, or the largest norm when the interval is the
## one in which every block is zero.
##
## The object returned is a "wafc_competitor" (with the class
## "wafc_threshold" in front), so the pilot reads it with the code it reads
## every competitor with: predict(), beta(u), g(grid) and blocks.

#' Estimation followed by a threshold on the blocks
#'
#' Zeroes the blocks \eqn{(\ell, m)} of a fit whose coefficient norm is not
#' above a threshold \eqn{t}, with \eqn{t} given or chosen by one of three
#' rules, and returns the result in the form of a
#' \code{\link{wafc_competitor}}. Without refit this is the estimator of
#' Corollary 8 of \file{derivations/06-selecao-limiar.tex}.
#'
#' @param object A \code{"cv.wafc"} object (the fit at \code{lambda.min} of
#'   the selected \eqn{J}, with its folds), a \code{"wafc"} object (then
#'   \code{s} is needed, and \code{foldid} for \code{rule = "cv"}), or a
#'   \code{"wafc_competitor"} of method \code{"klopp"}, in any of its forms.
#' @param t The threshold. When given, \code{rule} is not used.
#' @param rule \code{"max"}, \code{"cv"} or \code{"oracle"}; see the header
#'   of \file{wafc/R/threshold.R}.
#' @param c The fraction of the largest norm used by \code{rule = "max"}.
#' @param refit \code{"none"} (the default, the estimator of Corollary 8),
#'   \code{"support"} (least squares on the level terms and the nonzero
#'   coefficients of the blocks kept) or \code{"block"} (least squares on
#'   the level terms and every column of the blocks kept).
#' @param s For a \code{"wafc"} object, the penalty level of the fit.
#' @param y The response, when the object does not carry it (a
#'   \code{"klopp"} fit at a resolution with nothing to penalize).
#' @param foldid The folds, when the object does not carry them.
#' @param truth For \code{rule = "oracle"}: a list with the test sample
#'   \code{x}, \code{u} and the true regression function \code{f} there.
#' @param fold.fits The fits of the folds, as returned by
#'   \code{wafc_threshold_folds()}, so that several calls on the same
#'   object fit the folds once.
#' @param ... Passed to \code{\link{wafc}} when a fold of a WAFC fit is
#'   refitted (arguments of the engine that the object does not record).
#'
#' @return An object of class \code{c("wafc_threshold",
#'   "wafc_competitor")}. Its \code{extra} has the base \code{J} and
#'   \code{lambda}, the threshold \code{t}, the \code{rule}, \code{c},
#'   \code{refit}, the \eqn{p \times q} matrices \code{norm} (of the fit
#'   before the threshold) and \code{kept}, \code{nzero} (nonzero wavelet
#'   coefficients after the threshold), and, for the rules that score
#'   candidates, \code{candidates} (a data frame of the intervals and their
#'   error).
#'
#' @references van de Geer, S., Buhlmann, P. and Zhou, S. (2011). The
#'   adaptive and the thresholded Lasso for potentially misspecified models
#'   (and a lower bound for the Lasso). \emph{Electronic Journal of
#'   Statistics} 5, 688-749.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' cv <- cv.wafc(d$x, d$u, d$y, J = 3:4, nfolds = 5)
#' th <- wafc_threshold(cv, rule = "max")
#' th$blocks
#'
#' @export
wafc_threshold <- function(object, t = NULL, rule = c("max", "cv", "oracle"),
                           c = 0.15, refit = c("none", "support", "block"),
                           s = NULL, y = NULL, foldid = NULL, truth = NULL,
                           fold.fits = NULL, ...) {
  t0 <- proc.time()[["elapsed"]]
  rule <- match.arg(rule)
  refit <- match.arg(refit)
  base <- wafc_thr_base(object, s = s, y = y, foldid = foldid)
  des <- base[["design"]]
  nrm <- wafc_thr_norms(base[["b"]], des)
  cand <- NULL

  if (!is.null(t)) {
    if (!is.numeric(t) || length(t) != 1L || !is.finite(t) || t < 0) {
      stop("'t' must be a single non-negative value.", call. = FALSE)
    }
    rule <- "fixed"
  } else if (rule == "max") {
    if (!is.numeric(c) || length(c) != 1L || !is.finite(c) || c < 0) {
      stop("'c' must be a single non-negative value.", call. = FALSE)
    }
    t <- c * max(nrm)
  } else if (rule == "cv") {
    if (is.null(fold.fits)) {
      fold.fits <- wafc_threshold_folds(object, s = s, y = y,
                                        foldid = foldid, ...)
    }
    cand <- wafc_thr_cv(base, fold.fits, refit)
    t <- wafc_thr_pick(cand, max(c(nrm, cand[["lower"]])))
  } else {
    if (is.null(truth) || is.null(truth[["x"]]) || is.null(truth[["u"]]) ||
        is.null(truth[["f"]])) {
      stop("rule = \"oracle\" needs 'truth', a list with the test sample ",
           "'x', 'u' and the true regression function 'f' there.",
           call. = FALSE)
    }
    cand <- wafc_thr_oracle(base, nrm, refit, truth)
    t <- wafc_thr_pick(cand, max(nrm))
  }

  kept <- nrm > t
  est <- wafc_thr_apply(base[["a0"]], base[["b"]], des, kept)
  if (refit != "none") {
    est <- wafc_thr_refit(est[["b"]], des, kept, base[["y"]],
                          seq_len(des[["n"]]), refit, base[["intercept"]])
  }
  out <- wafc_thr_object(est[["a0"]], est[["b"]], des, kept)
  out[["method"]] <- paste0(base[["method"]], "+", rule,
                            if (refit != "none") paste0("+", refit) else "")
  out[["extra"]] <- c(out[["extra"]],
                      list(J = base[["J"]], lambda = base[["s"]], t = t,
                           rule = rule, c = if (rule == "max") c else
                             if (max(nrm) > 0) t / max(nrm) else NA_real_,
                           refit = refit, norm = nrm, candidates = cand))
  out[["time"]] <- proc.time()[["elapsed"]] - t0
  out
}

#' Fits of the folds of a fit, at the resolution and penalty it chose
#'
#' Refits each fold of the cross-validation of a fit at the pair
#' \eqn{(J, \lambda)} the fit selected, with no new search, and returns the
#' coefficients. It is what \code{\link{wafc_threshold}} does for
#' \code{rule = "cv"}, exposed so that several rules or refits on the same
#' object pay for the folds once.
#'
#' For a WAFC fit the fold is refitted as \code{\link{cv.wafc}} refits it:
#' the columns of the design of the whole sample, the penalty path of the
#' fit on the whole sample (down to the selected level), and the
#' coefficients read at that level, so the fold fits are the ones the
#' cross-validation scored. For the block LASSO, grpreg standardizes each
#' fold with the fold's own sample, where \code{cv.grpreg} uses the
#' standardization of the whole sample; the penalty path and the level read
#' are the same.
#'
#' @inheritParams wafc_threshold
#'
#' @return A list with \code{foldid} and \code{fits}, one list
#'   \code{(a0, b)} per fold, in the coordinates of the design (intercept
#'   and one coefficient per column).
#'
#' @export
wafc_threshold_folds <- function(object, s = NULL, y = NULL, foldid = NULL,
                                 ...) {
  base <- wafc_thr_base(object, s = s, y = y, foldid = foldid)
  if (is.null(base[["foldid"]])) {
    stop("rule = \"cv\" needs the folds: pass 'foldid', or a \"cv.wafc\" ",
         "object, which carries them.", call. = FALSE)
  }
  fid <- base[["foldid"]]
  fits <- lapply(seq_len(max(fid)), function(k) {
    base[["fold_fit"]](which(fid != k), ...)
  })
  list(foldid = fid, fits = fits)
}

## ---------------------------------------------------------------------------
## Internals
## ---------------------------------------------------------------------------

## What the threshold needs from a fit, whatever its class: the design, the
## response, the coefficients (a0, b) in the coordinates of the design, the
## penalty level, the folds, whether the fit has an intercept, and a function
## that refits a set of rows at the same (J, lambda).
wafc_thr_base <- function(object, s = NULL, y = NULL, foldid = NULL) {
  if (inherits(object, "cv.wafc")) {
    if (is.null(foldid)) foldid <- object[["foldid"]]
    if (is.null(s)) s <- object[["lambda.min"]]
    object <- object[["wafc.fit"]]
  }
  if (inherits(object, "wafc")) {
    return(wafc_thr_base_wafc(object, s, y, foldid))
  }
  if (inherits(object, "wafc_competitor") &&
      identical(object[["method"]], "klopp")) {
    return(wafc_thr_base_klopp(object, y, foldid))
  }
  stop("'object' must be a \"cv.wafc\" or \"wafc\" fit, or a ",
       "\"wafc_competitor\" of method \"klopp\".", call. = FALSE)
}

wafc_thr_base_wafc <- function(object, s, y, foldid) {
  if (is.null(s)) {
    stop("A \"wafc\" object needs 's', the penalty level of the fit.",
         call. = FALSE)
  }
  s <- wafc_single_s(object, s)
  des <- object[["design"]]
  if (is.null(y)) y <- object[["y"]]
  if (!is.null(foldid)) foldid <- wafc_foldid(foldid, des[["n"]], 10L)
  cf <- wafc_raw_coef(object, s = s)[, 1L]
  lam <- object[["lambda"]]
  ## the path down to s and one point past it, which is all the
  ## interpolation at s reads; the points before are computed exactly as on
  ## the whole path, since the engines walk it in order
  path <- lam[seq_len(min(length(lam), sum(lam >= s) + 1L))]
  pen <- object[["penalty"]]
  fold_fit <- function(rows, ...) {
    a <- list(design = wafc_subset_design(des, rows), y = y[rows],
              penalty = pen, lambda = path, intercept = object[["intercept"]])
    if (pen == "sglasso") a[["asparse"]] <- object[["asparse"]]
    fi <- do.call(wafc, c(a, list(...)))
    v <- wafc_raw_coef(fi, s = s)[, 1L]
    list(a0 = v[[1L]], b = unname(v[-1L]))
  }
  list(design = des, y = y, a0 = cf[[1L]], b = unname(cf[-1L]), s = s,
       J = des[["J"]][1L], foldid = foldid,
       intercept = isTRUE(object[["intercept"]]),
       method = paste0("wafc.", pen), fold_fit = fold_fit)
}

wafc_thr_base_klopp <- function(object, y, foldid) {
  des <- object[["design"]]
  ex <- object[["extra"]]
  cvo <- object[["fit"]]
  ols <- is.null(cvo[["fit"]])
  if (is.null(y)) {
    if (ols) {
      stop("This \"klopp\" fit is a least squares fit and does not carry ",
           "the response; pass 'y'.", call. = FALSE)
    }
    ## grpreg keeps the response as given, with its mean as an attribute
    y <- as.numeric(cvo[["fit"]][["y"]])
  }
  y <- wafc_check_y(y, des[["n"]])
  if (is.null(foldid) && !ols) foldid <- cvo[["fold"]]
  if (!is.null(foldid)) foldid <- wafc_foldid(foldid, des[["n"]], 10L)
  cf <- if (ols) cvo[["coef"]] else as.numeric(stats::coef(cvo))
  grp <- wafc_kp_groups(des, ex[["block.size"]], ex[["penalize.levels"]],
                        ex[["merge.coarse"]], ex[["balanced"]],
                        ex[["free.coarse"]])
  Z <- as.matrix(des[["Z"]])
  s <- ex[["lambda"]]
  path <- if (ols) NULL else cvo[["lambda"]][seq_len(cvo[["min"]])]
  fold_fit <- function(rows, ...) {
    if (ols) {
      b <- stats::lm.fit(cbind(1, Z[rows, , drop = FALSE]),
                         y[rows])[["coefficients"]]
      b[is.na(b)] <- 0
      b <- unname(b)
      return(list(a0 = b[1L], b = b[-1L]))
    }
    a <- list(Z[rows, , drop = FALSE], y[rows], group = grp,
              penalty = "grLasso", lambda = path)
    if (identical(ex[["chunk.weights"]], "unit")) {
      a[["group.multiplier"]] <- rep(1, max(grp))
    }
    fi <- do.call(grpreg::grpreg, a)
    k <- which.min(abs(fi[["lambda"]] - s))
    v <- as.numeric(fi[["beta"]][, k])
    list(a0 = v[1L], b = v[-1L])
  }
  list(design = des, y = y, a0 = cf[1L], b = cf[-1L], s = s,
       J = ex[["J"]], foldid = foldid, intercept = TRUE,
       method = "klopp", fold_fit = fold_fit)
}

## Norm of the wavelet coefficients of each block, the statistic of
## wafc_blocks(), as a p by q matrix.
wafc_thr_norms <- function(b, design) {
  p <- design[["p"]]
  q <- design[["q"]]
  out <- matrix(0, p, q, dimnames = list(design[["xnames"]],
                                         design[["unames"]]))
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      idx <- design[["blocks"]][[wafc_block_name(design, l, m)]]
      out[l, m] <- sqrt(sum(b[idx]^2))
    }
  }
  out
}

## The fit with the blocks not kept set to zero.
wafc_thr_apply <- function(a0, b, design, kept) {
  for (l in seq_len(design[["p"]])) {
    for (m in seq_len(design[["q"]])) {
      if (kept[l, m]) next
      b[design[["blocks"]][[wafc_block_name(design, l, m)]]] <- 0
    }
  }
  list(a0 = a0, b = b)
}

## Least squares on the rows 'rows', over the level terms and, in each block
## kept, the nonzero coefficients of 'b' ("support") or every column
## ("block"). The constant covariate, when the fit has an intercept, is left
## out of the columns, since the intercept carries its level (the
## convention of both engines); without an intercept every level term is
## in. Aliased columns get coefficient zero.
wafc_thr_refit <- function(b, design, kept, y, rows, refit, intercept) {
  cols <- design[["unpenalized"]]
  if (intercept) cols <- setdiff(cols, design[["constant"]])
  for (l in seq_len(design[["p"]])) {
    for (m in seq_len(design[["q"]])) {
      if (!kept[l, m]) next
      idx <- design[["blocks"]][[wafc_block_name(design, l, m)]]
      if (refit == "support") idx <- idx[b[idx] != 0]
      cols <- c(cols, idx)
    }
  }
  W <- as.matrix(design[["Z"]][rows, cols, drop = FALSE])
  if (intercept) W <- cbind(1, W)
  v <- if (ncol(W) == 0L) numeric(0) else {
    z <- qr.coef(qr(W), y[rows])
    z[is.na(z)] <- 0
    unname(z)
  }
  out <- numeric(length(b))
  if (intercept) {
    a0 <- v[1L]
    v <- v[-1L]
  } else {
    a0 <- 0
  }
  out[cols] <- v
  list(a0 = a0, b = out)
}

## The thresholded fit as a "wafc_competitor": the levels with the intercept
## folded into the constant covariate when there is one, and the components
## through the spec route, as wafc_fit_klopp() builds its object.
wafc_thr_object <- function(a0, b, design, kept) {
  p <- design[["p"]]
  q <- design[["q"]]
  xn <- design[["xnames"]]
  un <- design[["unames"]]
  Z <- design[["Z"]]
  fitted <- as.numeric(a0 + Z %*% b)
  cc <- b[design[["unpenalized"]]]
  names(cc) <- xn
  const <- design[["constant"]]
  if (length(const) == 1L) {
    cc[const] <- cc[const] + a0 / as.numeric(Z[1L, const])
    a0 <- 0
  }
  wav <- unlist(design[["blocks"]], use.names = FALSE)
  nz <- matrix(FALSE, p, q, dimnames = list(xn, un))
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      idx <- design[["blocks"]][[wafc_block_name(design, l, m)]]
      nz[l, m] <- any(b[idx] != 0)
    }
  }
  g_of <- function(newu) {
    newu <- wafc_as_matrix(newu, "newu")
    d_new <- wafc_design(matrix(1, nrow(newu), p), newu, spec = design)
    g <- vector("list", p * q)
    dim(g) <- c(p, q)
    for (l in seq_len(p)) {
      for (m in seq_len(q)) {
        idx <- design[["blocks"]][[wafc_block_name(design, l, m)]]
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
    for (l in seq_len(p)) {
      for (m in seq_len(q)) out[, l] <- out[, l] + g[[l, m]]
    }
    out
  }
  out <- list(cc = cc, beta = beta_fun, g = g_of, blocks = nz,
              fitted = fitted, intercept = a0, coef = c(a0, b),
              xnames = xn, unames = un, design = design,
              n = design[["n"]], p = p, q = q,
              extra = list(kept = kept, nzero = sum(b[wav] != 0)))
  class(out) <- c("wafc_threshold", "wafc_competitor")
  out
}

## The cross-validated error of the thresholded fit at every candidate. The
## breakpoints are the norms of the blocks of every fold fit; between two of
## them no fold fit changes, so the error of the interval [lower, upper) is
## the error at its lower end. Each fold scores only its distinct kept sets
## (at most pq + 1), and the error is the mean of the fold means, the
## measure of cv.wafc().
wafc_thr_cv <- function(base, folds, refit) {
  des <- base[["design"]]
  fid <- folds[["foldid"]]
  K <- length(folds[["fits"]])
  y <- base[["y"]]
  Z <- des[["Z"]]
  fnorm <- lapply(folds[["fits"]], function(f) wafc_thr_norms(f[["b"]], des))
  brk <- sort(unique(c(0, unlist(fnorm))))
  fold_err <- matrix(NA_real_, K, length(brk))
  for (k in seq_len(K)) {
    out <- which(fid == k)
    rows <- which(fid != k)
    f <- folds[["fits"]][[k]]
    nk <- fnorm[[k]]
    lev <- sort(unique(c(0, as.numeric(nk))))
    ## the kept set at each breakpoint is the one at the largest norm of
    ## this fold not above it
    which_lev <- findInterval(brk, lev)
    e <- vapply(seq_along(lev), function(i) {
      kept <- nk > lev[i]
      est <- wafc_thr_apply(f[["a0"]], f[["b"]], des, kept)
      if (refit != "none") {
        est <- wafc_thr_refit(est[["b"]], des, kept, y, rows, refit,
                              base[["intercept"]])
      }
      eta <- as.numeric(Z[out, , drop = FALSE] %*% est[["b"]]) + est[["a0"]]
      mean((y[out] - eta)^2)
    }, 0)
    fold_err[k, ] <- e[which_lev]
  }
  data.frame(lower = brk, upper = c(brk[-1L], Inf),
             error = colMeans(fold_err))
}

## The error of the thresholded fit against the truth at every candidate:
## the intervals between consecutive norms of the fit, the error the root
## mean squared one on the test sample.
wafc_thr_oracle <- function(base, nrm, refit, truth) {
  des <- base[["design"]]
  Zt <- wafc_design(truth[["x"]], truth[["u"]], spec = des)[["Z"]]
  f <- as.numeric(truth[["f"]])
  if (length(f) != nrow(Zt)) {
    stop("'truth$f' must have one value per row of 'truth$x'.", call. = FALSE)
  }
  brk <- sort(unique(c(0, as.numeric(nrm))))
  e <- vapply(brk, function(tt) {
    kept <- nrm > tt
    est <- wafc_thr_apply(base[["a0"]], base[["b"]], des, kept)
    if (refit != "none") {
      est <- wafc_thr_refit(est[["b"]], des, kept, base[["y"]],
                            seq_len(des[["n"]]), refit, base[["intercept"]])
    }
    sqrt(mean((as.numeric(Zt %*% est[["b"]]) + est[["a0"]] - f)^2))
  }, 0)
  data.frame(lower = brk, upper = c(brk[-1L], Inf), error = e)
}

## The threshold of the best interval: its midpoint, or 'top' when it is the
## last one, in which every block is zero. Ties go to the first interval,
## the smallest t.
wafc_thr_pick <- function(cand, top) {
  i <- which.min(cand[["error"]])
  if (is.infinite(cand[["upper"]][i])) return(top)
  (cand[["lower"]][i] + cand[["upper"]][i]) / 2
}
