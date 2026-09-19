## wafc/R/design.R -- rescaling of the modulating covariates and the design
## matrix of the WAFC model (step E2.1).
##
## The design is the one of wall() in WaveBased with every block multiplied
## by a linear covariate: wall() is read, not copied, and the functions here
## are written from scratch with their own names (docs/instrucoes.md,
## section 6). The blocks by covariate, the rescaling to [eps, 1 - eps] and
## the discarding of the constant scaling function follow .wall_design(),
## .wall_rescale_pars() and .wall_penalty() of WaveBased/R/wall.R.
##
## Column order is frozen by decision D12 (docs/ESTADO.md): first the p
## unpenalized columns X_1, ..., X_p, then the blocks (l, m) in
## lexicographic order, and inside each block the wavelets by increasing j
## and, within a level, by increasing k. The constant scaling function
## phi_{00} of each block is discarded, which is what imposes the
## identifiability constraint at the level of the basis (Lemma 1 of
## derivations/01-identificabilidade.md).

#' Rescale modulating covariates to the unit interval
#'
#' Maps each column of \code{u} linearly onto \eqn{[\epsilon, 1-\epsilon]},
#' the transformation used before evaluating the wavelet basis. The
#' periodized basis carries an artifact on the boundary strip, and
#' \eqn{\epsilon > 0} keeps the data away from it. The margin has a
#' declared role (decision D23): it periodizes the basis on an interval
#' strictly larger than the support of the data, which is what allows the
#' functional coefficient to be replaced by an extension that emends into
#' \eqn{0 \equiv 1}, so \eqn{\epsilon} is a fixed feature of the target
#' and not a function of the resolution level \eqn{J}.
#'
#' @param u Matrix (or data frame, or vector) of modulating covariates, with
#'   \eqn{n} rows and \eqn{q} columns.
#' @param eps Numeric in \eqn{[0, 0.5)}, of length 1 or \eqn{q}: the target
#'   image is \eqn{[\epsilon, 1-\epsilon]}.
#' @param location,scale Optional numeric vectors of length \eqn{q} giving
#'   the transformation \eqn{(u_m - location_m)/scale_m} directly. When
#'   \code{NULL} (the default) they are computed from the column ranges of
#'   \code{u}, so that the image is exactly \eqn{[\epsilon, 1-\epsilon]}.
#' @param clip Logical. If \code{TRUE}, the rescaled values are truncated to
#'   \eqn{[\epsilon, 1-\epsilon]}. Used when \code{location} and \code{scale}
#'   come from a previous sample and new observations may fall outside its
#'   range.
#'
#' @return A list with the rescaled matrix \code{u} and the components
#'   \code{location}, \code{scale} and \code{eps} that define the
#'   transformation.
#'
#' @examples
#' r <- wafc_rescale(matrix(rnorm(20), 10, 2), eps = 0.05)
#' range(r$u)
#'
#' @export
wafc_rescale <- function(u, eps = 0, location = NULL, scale = NULL,
                         clip = FALSE) {
  u <- wafc_as_matrix(u, "u")
  q <- ncol(u)
  eps <- wafc_recycle(eps, q, "eps")
  if (any(!is.finite(eps)) || any(eps < 0) || any(eps >= 0.5)) {
    stop("'eps' must belong to [0, 0.5).", call. = FALSE)
  }
  if (is.null(location) || is.null(scale)) {
    location <- numeric(q)
    scale <- numeric(q)
    for (m in seq_len(q)) {
      rg <- range(u[, m])
      if (!all(is.finite(rg))) {
        stop("Column ", m, " of 'u' has non-finite values.", call. = FALSE)
      }
      if (rg[1L] == rg[2L]) {
        stop("Column ", m, " of 'u' is constant and cannot be rescaled.",
             call. = FALSE)
      }
      ## image [eps, 1 - eps]: widen the observed range by a on each side,
      ## with a/(2a + range) = eps.
      a <- eps[m] * diff(rg) / (1 - 2 * eps[m])
      location[m] <- rg[1L] - a
      scale[m] <- 2 * a + diff(rg)
    }
  } else {
    location <- wafc_recycle(location, q, "location")
    scale <- wafc_recycle(scale, q, "scale")
    if (any(!is.finite(scale)) || any(scale <= 0)) {
      stop("'scale' must be positive and finite.", call. = FALSE)
    }
  }
  for (m in seq_len(q)) {
    u[, m] <- (u[, m] - location[m]) / scale[m]
  }
  if (clip) {
    for (m in seq_len(q)) {
      u[, m] <- pmin(pmax(u[, m], eps[m]), 1 - eps[m])
    }
  }
  list(u = u, location = location, scale = scale, eps = eps)
}

#' Design matrix of the WAFC model
#'
#' Builds the matrix \eqn{Z} whose rows are
#' \eqn{(X_{i1}, \ldots, X_{ip}, X_{i\ell}\psi_{jk}(U_{im}))}, that is, the
#' additive wavelet expansion of every functional coefficient multiplied by
#' its linear covariate, together with the penalty factors that leave the
#' \eqn{p} level terms \eqn{c_\ell} unpenalized.
#'
#' The wavelet basis of each modulating covariate is evaluated once by
#' \code{WaveBased::wbasis} and reused by the \eqn{p} blocks that share it.
#' With the default periodized basis and \code{j0 = 0} the single scaling
#' function \eqn{\phi_{00}} is constant equal to one, so its column would
#' repeat the unpenalized column \eqn{X_\ell} in every block and make the
#' Gram matrix singular with nullity \eqn{pq}: it is discarded, which
#' imposes \eqn{\int_0^1 g_{\ell m} = 0} at the level of the basis. See
#' Lemma 1 of \file{derivations/01-identificabilidade.md}.
#'
#' @param x Matrix (or data frame, or vector) of linear covariates, with
#'   \eqn{n} rows and \eqn{p} columns. A constant column (the usual
#'   \eqn{X_1 \equiv 1}) is allowed and is the way to obtain an additive
#'   intercept.
#' @param u Matrix (or data frame, or vector) of modulating covariates, with
#'   \eqn{n} rows and \eqn{q} columns.
#' @param J Resolution level of the sieve: a single integer used for every
#'   modulating covariate, or one entry per covariate. It must satisfy
#'   \code{J > j0}.
#' @param j0 Coarsest resolution level of the basis. Default \code{0}, the
#'   case in which the periodized basis discards a single constant column
#'   per block.
#' @param family,filter.size,prec.wavelet,wavelet.filter The wavelet basis,
#'   as in \code{\link[WaveBased]{wbasis}}. The default is the Daublet with
#'   filter size 8 (four vanishing moments), which is the basis used in the
#'   numerical checks of \file{derivations/check/}.
#' @param boundary Boundary treatment of the basis. Only \code{"periodic"}
#'   is implemented; \code{"interval"} is left for a later step (open
#'   question 6 of \file{docs/ESTADO.md}), and asking for it raises an
#'   informative error.
#' @param rescale Logical. If \code{TRUE} (default), each modulating
#'   covariate is mapped to \eqn{[\epsilon, 1-\epsilon]} by
#'   \code{\link{wafc_rescale}} before the basis is evaluated, and the
#'   transformation is stored in the returned object. If \code{FALSE}, the
#'   modulating covariates are assumed to lie in \eqn{[0,1]}.
#' @param eps Value in \eqn{[0, 0.5)} used by the rescaling, of length 1 or
#'   \eqn{q}. The default does not depend on \eqn{J} (decision D23): it is
#'   \eqn{0} when \code{boundary = "interval"}, where there is no
#'   periodization to keep the data away from, and a fixed constant,
#'   currently \eqn{0.05}, in the periodized case. That constant is
#'   provisional and is the quantity step E2.4 measures, by reading
#'   \eqn{\lambda_{\min}(G_\epsilon)} and the integrated squared error
#'   against \eqn{\epsilon}. The rule \eqn{1.9^{-J}} of
#'   \code{WaveBased::wall}, used until now, tied the margin to the finest
#'   scale: it made every candidate of \code{\link{cv.wafc}} estimate a
#'   slightly different target, and it is not admissible at \eqn{J = 1}.
#' @param use.table Whether the basis is evaluated by table lookup, which is
#'   faster than the exact Daubechies-Lagarias algorithm on large samples:
#'   one of \code{"auto"} (default), \code{"always"} or \code{"never"}. See
#'   \code{\link[WaveBased]{wtable}}.
#' @param wavelet.table An optional table built by
#'   \code{\link[WaveBased]{wtable}}, which overrides \code{use.table}.
#' @param sparse Whether the design is stored in the sparse format of the
#'   \pkg{Matrix} package: one of \code{"auto"} (default), \code{"always"}
#'   or \code{"never"}. The unpenalized columns are dense by nature and the
#'   rule looks only at the wavelet part.
#' @param spec An object returned by a previous call, typically the design
#'   of the training sample. When supplied, every argument describing the
#'   basis and the rescaling is taken from it, and the rescaled modulating
#'   covariates are truncated to the interval used there, so that the two
#'   designs have the same columns in the same order.
#'
#' @return An object of class \code{"wafc_design"}: a list with the design
#'   matrix \code{Z} (dense or sparse, \eqn{n \times (p + pqN_J)}), the
#'   vector \code{penalty.factor} (zero on the \eqn{p} first columns, one on
#'   the wavelet columns), the column indices \code{unpenalized} and
#'   \code{blocks} (one entry per pair \eqn{(\ell, m)}, in the order of
#'   D12), the indices \code{constant} of the linear covariates that do not
#'   vary (whose level is carried by the intercept of \pkg{glmnet}, which
#'   drops constant columns from the fit), the basis and rescaling
#'   specification, and the dimensions
#'   \code{n}, \code{p}, \code{q}, \code{NJ} and \code{nvars}.
#'
#' @examples
#' n <- 200
#' x <- cbind(1, rnorm(n))
#' u <- matrix(runif(2 * n), n, 2)
#' d <- wafc_design(x, u, J = 3)
#' dim(d$Z)
#' head(colnames(d$Z))
#'
#' @export
wafc_design <- function(x, u, J, j0 = 0L, family = "Daublets",
                        filter.size = 8L, prec.wavelet = 30L,
                        wavelet.filter = NULL,
                        boundary = c("periodic", "interval"),
                        rescale = TRUE, eps = NULL,
                        use.table = c("auto", "always", "never"),
                        wavelet.table = NULL,
                        sparse = c("auto", "always", "never"),
                        spec = NULL) {

  this_call <- match.call()
  x <- wafc_as_matrix(x, "x")
  u <- wafc_as_matrix(u, "u")
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(u)
  if (nrow(u) != n) {
    stop("'x' and 'u' must have the same number of rows (", n, " and ",
         nrow(u), ").", call. = FALSE)
  }
  if (n == 0L) stop("'x' and 'u' must have at least one row.", call. = FALSE)

  if (is.null(spec)) {
    boundary <- match.arg(tolower(boundary[1L]), c("periodic", "interval"))
    use.table <- match.arg(use.table)
    sparse <- match.arg(sparse)
    if (boundary == "interval") {
      stop("boundary = \"interval\" is not implemented in wafc_design() yet. ",
           "The reparametrisation of the scaling block that keeps the ",
           "identifiability constraint is described in section 5 of ",
           "derivations/01-identificabilidade.md and is open question 6 of ",
           "docs/ESTADO.md.", call. = FALSE)
    }
    if (missing(J)) {
      stop("The resolution level 'J' must be provided.", call. = FALSE)
    }
    j0 <- wafc_check_j0(j0)
    J <- wafc_check_J(J, j0, q)
    eps <- wafc_eps(eps, q, rescale, boundary)
    if (rescale) {
      rs <- wafc_rescale(u, eps = eps, clip = FALSE)
    } else {
      rs <- list(u = u, location = rep_len(0, q), scale = rep_len(1, q),
                 eps = eps)
      if (min(u) < 0 || max(u) > 1) {
        warning("Some modulating covariates lie outside [0,1] and rescale = FALSE. ",
                "The wavelet basis is only defined on the unit interval.",
                call. = FALSE)
      }
    }
    wtab <- wafc_table(use.table, wavelet.table, workload = n * q,
                       family = family, filter.size = filter.size,
                       prec.wavelet = prec.wavelet,
                       wavelet.filter = wavelet.filter)
    obj <- list(J = J, j0 = j0, boundary = boundary, family = family,
                filter.size = filter.size, prec.wavelet = prec.wavelet,
                wavelet.filter = wavelet.filter, wavelet.table = wtab,
                rescale = rescale, location = rs[["location"]],
                scale = rs[["scale"]], eps = rs[["eps"]],
                xnames = wafc_names(x, p, "x"), unames = wafc_names(u, q, "u"))
    obj[["sparse"]] <- wafc_sparse(sparse, obj[["J"]], j0,
                                   wafc_filter_length(family, filter.size,
                                                      wavelet.filter), p)
    u_eval <- rs[["u"]]
  } else {
    if (!inherits(spec, "wafc_design")) {
      stop("'spec' must be an object returned by wafc_design().", call. = FALSE)
    }
    if (ncol(u) != length(spec[["J"]])) {
      stop("'u' must have ", length(spec[["J"]]),
           " column(s), as the data used to build 'spec'.", call. = FALSE)
    }
    if (p != spec[["p"]]) {
      stop("'x' must have ", spec[["p"]],
           " column(s), as the data used to build 'spec'.", call. = FALSE)
    }
    obj <- spec[c("J", "j0", "boundary", "family", "filter.size",
                  "prec.wavelet", "wavelet.filter", "wavelet.table",
                  "rescale", "location", "scale", "eps", "xnames", "unames",
                  "sparse")]
    u_eval <- if (obj[["rescale"]]) {
      wafc_rescale(u, eps = obj[["eps"]], location = obj[["location"]],
                   scale = obj[["scale"]], clip = TRUE)[["u"]]
    } else {
      u
    }
  }

  NJ <- as.integer(2^obj[["J"]] - 2^obj[["j0"]])

  ## One basis evaluation per modulating covariate, shared by the p blocks
  ## that use it. The constant column phi_{00} is dropped by wafc_wbasis().
  basis <- vector("list", q)
  for (m in seq_len(q)) {
    basis[[m]] <- wafc_wbasis(u_eval[, m], obj, obj[["J"]][m])
  }

  ## Columns in the order of D12: the p unpenalized terms, then the blocks
  ## (l, m) lexicographically, each one X_l times the basis of U_m.
  cols <- vector("list", 1L + p * q)
  cols[[1L]] <- if (obj[["sparse"]]) Matrix::Matrix(x, sparse = TRUE) else x
  colnames(cols[[1L]]) <- obj[["xnames"]]
  blocks <- vector("list", p * q)
  block_names <- character(p * q)
  pos <- p
  b <- 0L
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      b <- b + 1L
      B <- wafc_scale_rows(basis[[m]], x[, l], obj[["sparse"]])
      colnames(B) <- wafc_colnames(obj[["xnames"]][l], obj[["unames"]][m],
                                   obj[["j0"]], obj[["J"]][m])
      cols[[1L + b]] <- B
      blocks[[b]] <- pos + seq_len(NJ[m])
      block_names[b] <- paste0(obj[["xnames"]][l], ":", obj[["unames"]][m])
      pos <- pos + NJ[m]
    }
  }
  names(blocks) <- block_names
  Z <- if (obj[["sparse"]]) Reduce(Matrix::cbind2, cols) else do.call(cbind, cols)

  obj[["Z"]] <- Z
  obj[["penalty.factor"]] <- c(rep.int(0, p), rep.int(1, sum(NJ) * p))
  obj[["unpenalized"]] <- seq_len(p)
  ## Constant linear covariates, typically X_1 = 1. glmnet drops columns
  ## with zero variance from the fit and lets its own intercept carry them,
  ## so the level c_l of a constant covariate is read from the intercept
  ## and not from the coefficient of the column (the same phenomenon that
  ## E1.4 recorded in wall()). Recorded here for the fitting step.
  obj[["constant"]] <- which(vapply(seq_len(p),
                                    function(l) diff(range(x[, l])) == 0,
                                    TRUE))
  obj[["blocks"]] <- blocks
  obj[["n"]] <- n
  obj[["p"]] <- p
  obj[["q"]] <- q
  obj[["NJ"]] <- NJ
  obj[["nvars"]] <- ncol(Z)
  obj[["call"]] <- this_call
  class(obj) <- "wafc_design"
  obj
}

## ---------------------------------------------------------------------------
## Internals
## ---------------------------------------------------------------------------

## Coerces to a numeric matrix, keeping the column names when there are any.
wafc_as_matrix <- function(z, what) {
  if (is.data.frame(z)) z <- as.matrix(z)
  if (is.vector(z) && !is.list(z)) z <- matrix(z, ncol = 1L)
  if (!is.matrix(z) || !is.numeric(z)) {
    stop("'", what, "' must be a numeric matrix, data frame or vector.",
         call. = FALSE)
  }
  if (anyNA(z)) stop("'", what, "' must not contain NA.", call. = FALSE)
  z
}

wafc_recycle <- function(v, len, what) {
  if (!is.numeric(v) || length(v) == 0L ||
      (length(v) != 1L && length(v) != len)) {
    stop("'", what, "' must be numeric of length 1 or ", len, ".",
         call. = FALSE)
  }
  rep_len(v, len)
}

wafc_names <- function(z, k, prefix) {
  nm <- colnames(z)
  if (is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
    nm <- paste0(prefix, seq_len(k))
  }
  nm
}

wafc_check_j0 <- function(j0) {
  if (length(j0) != 1L || !is.finite(j0) || j0 < 0 || j0 != round(j0)) {
    stop("'j0' must be a single non-negative integer.", call. = FALSE)
  }
  if (j0 != 0L) {
    stop("Only j0 = 0 is implemented in wafc_design(): with the periodized ",
         "basis it is the case in which a single constant scaling function ",
         "is discarded per block (decision D2, Lemma 1 of ",
         "derivations/01-identificabilidade.md).", call. = FALSE)
  }
  as.integer(j0)
}

## Validates the resolution level(s) and recycles them to one per
## modulating covariate.
wafc_check_J <- function(J, j0, q) {
  if (!is.numeric(J) || any(!is.finite(J)) || any(J != round(J))) {
    stop("'J' must contain only integer values.", call. = FALSE)
  }
  if (length(J) != 1L && length(J) != q) {
    stop("'J' must have length 1 or one entry per modulating covariate (",
         q, ").", call. = FALSE)
  }
  if (any(J <= j0)) {
    stop("'J' must be larger than 'j0' (= ", j0, ").", call. = FALSE)
  }
  rep_len(as.integer(J), q)
}

## Provisional default margin of the periodized case. The value is the
## smallest fixed margin the numerical check of E1.3b measured as buying
## the full approximation rate (part C of derivations/check/
## 02-aproximacao-besov.R: at eps = 0.05 the restricted projection error
## falls 4 to 13 bits per level, against 1/2 bit at eps = 0 and 0.96 bit
## under the inherited rule 1.9^(-J)). It is PROVISIONAL, and provisional
## in a direction that is already known: the same check reads
## lambda_min(G_eps) collapsing to 10^(-12) at J = 6 under this margin,
## because wavelets supported inside the excluded strip become invisible.
## The margin that buys the rate is the one that degenerates the restricted
## Gram matrix, and the arbitration between the two is what step E2.4
## measures, over eps in {0, 2^(-(J+1)), 1.9^(-J), 0.02, 0.05, 0.10}
## (decision D23).
wafc_eps_periodic <- 0.05

## Margin of the rescaling, one entry per modulating covariate. The value
## does not depend on J (decision D23): the margin is what buys the
## extension of g_{lm} to [0,1] that emends into 0 = 1 (E1.3b), a
## populational device, and not a numerical convenience of the finest
## scale. The rule 1.9^(-J) of WaveBased::wall, inherited by E2.1, made
## every candidate of cv.wafc() estimate a slightly different target,
## because the rescaled support moves with J; put the margin at the order
## of one cell of the finest scale, 2^(-J), which is exactly where the
## boundary wavelets of that level lose observations; and was not even
## admissible at J = 1, where it gives 0.526, outside the [0, 0.5) that
## wafc_rescale() requires. With boundary = "interval" there is no
## periodization to keep the data away from and the default margin is zero.
wafc_eps <- function(eps, q, rescale, boundary = "periodic") {
  if (!rescale) return(rep_len(0, q))
  if (is.null(eps)) {
    return(rep_len(if (boundary == "interval") 0 else wafc_eps_periodic, q))
  }
  eps <- wafc_recycle(eps, q, "eps")
  if (any(!is.finite(eps)) || any(eps < 0) || any(eps >= 0.5)) {
    stop("'eps' must belong to [0, 0.5).", call. = FALSE)
  }
  eps
}

wafc_filter_length <- function(family, filter.size, wavelet.filter) {
  if (is.null(wavelet.filter)) as.integer(filter.size)
  else length(wavelet.filter)
}

## Lookup-table policy, with the rule of .wall_table(): the workload here is
## n * q, because the basis is evaluated once per modulating covariate and
## reused by the p blocks that share it.
wafc_table <- function(use.table, wavelet.table, workload, family,
                       filter.size, prec.wavelet, wavelet.filter) {
  if (!is.null(wavelet.table)) return(wavelet.table)
  L <- wafc_filter_length(family, filter.size, wavelet.filter)
  build <- switch(use.table,
                  always = TRUE,
                  never = FALSE,
                  auto = workload >= 2000 * L)
  if (!build) return(NULL)
  if (is.null(wavelet.filter)) {
    WaveBased::wtable(family = family, filter.size = filter.size,
                      prec.wavelet = prec.wavelet, check = FALSE)
  } else {
    WaveBased::wtable(family = "Own", prec.wavelet = prec.wavelet,
                      wavelet.filter = wavelet.filter, check = FALSE)
  }
}

## Sparse storage of the design. At level j at most min(2^j, L - 1)
## translates are non-zero at a given point, which bounds the density of the
## wavelet part; the p unpenalized columns are dense and are counted as
## such.
wafc_sparse <- function(sparse, J, j0, L, p) {
  if (sparse == "always") return(TRUE)
  if (sparse == "never") return(FALSE)
  ncol_w <- 0
  nnz_w <- 0
  for (Jm in J) {
    ncol_w <- ncol_w + 2^Jm - 2^j0
    nnz_w <- nnz_w + sum(pmin(2^(j0:(Jm - 1L)), L - 1))
  }
  ncols <- p + p * ncol_w
  nnz <- p + p * nnz_w
  ncols >= 128 && nnz / ncols < 0.4
}

## Evaluates the wavelet basis of one (already rescaled) modulating
## covariate and drops the constant scaling function phi_{00}.
wafc_wbasis <- function(u, obj, J) {
  B <- if (is.null(obj[["wavelet.filter"]])) {
    WaveBased::wbasis(u, j0 = obj[["j0"]], J = J, family = obj[["family"]],
                      filter.size = obj[["filter.size"]],
                      prec.wavelet = obj[["prec.wavelet"]],
                      boundary = obj[["boundary"]],
                      wavelet.table = obj[["wavelet.table"]])
  } else {
    WaveBased::wbasis(u, j0 = obj[["j0"]], J = J, family = "Own",
                      filter.size = obj[["filter.size"]],
                      prec.wavelet = obj[["prec.wavelet"]],
                      wavelet.filter = obj[["wavelet.filter"]],
                      boundary = obj[["boundary"]],
                      wavelet.table = obj[["wavelet.table"]])
  }
  B <- B[, -1L, drop = FALSE]
  if (obj[["sparse"]]) Matrix::Matrix(B, sparse = TRUE) else B
}

## Multiplies every row of the basis block by the corresponding entry of the
## linear covariate.
wafc_scale_rows <- function(B, v, sparse) {
  if (sparse) Matrix::Diagonal(x = v) %*% B else B * v
}

## Labels of the columns of one block, e.g. "x2:u1:psi3.5" for
## X_2 psi_{3,5}(U_1).
wafc_colnames <- function(xname, uname, j0, J) {
  nm <- character(0L)
  for (j in j0:(J - 1L)) {
    nm <- c(nm, paste0(xname, ":", uname, ":psi", j, ".", seq_len(2^j) - 1L))
  }
  nm
}
