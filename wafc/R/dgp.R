## wafc/R/dgp.R -- data generating processes for the WAFC model (step E2.1).
##
## Model D1 (docs/notacao.md, section 2):
##
##   Y = sum_l beta_l(U) X_l + eps,   beta_l(u) = c_l + sum_m g_{lm}(u_m),
##
## with the components centred in Lebesgue measure on [0,1], which is the
## centring the basis imposes when the constant scaling function is
## discarded (Lemma 1 of derivations/01-identificabilidade.md, and open
## question 4 of docs/ESTADO.md).
##
## Three scenarios, in the terms of docs/plano-projeto.md, E2.1:
##
##   "smooth"          sine, cosine and a cubic: functions in every Besov
##                     space, where a linear sieve is already efficient;
##   "inhomogeneous"   bumps, blocks and heavisine of Donoho and Johnstone
##                     (1994), rescaled to [0,1]: spatially inhomogeneous
##                     regularity, which is where the wavelet LASSO is
##                     supposed to win;
##   "null"            every g_{lm} is zero, so the model is a linear
##                     regression and the whole penalized part must be
##                     shrunk away.
##
## The active structure is the same in the first two: beta_1 depends on two
## modulating covariates (the additive structure that is the point of the
## model), beta_2 on one, and beta_l is constant for l >= 3.

#' Additive components of the data generating processes
#'
#' Returns one of the univariate functions used by
#' \code{\link{simulate_wafc}}, centred and normalised on \eqn{[0,1]}: the
#' returned function \eqn{g} satisfies \eqn{\int_0^1 g = 0} and
#' \eqn{\int_0^1 g^2 = 1}, both up to the accuracy of the grid used to
#' compute the two constants (\eqn{2^{16}} midpoints).
#'
#' @param name One of \code{"sine"}, \code{"cosine"}, \code{"cubic"},
#'   \code{"bumps"}, \code{"blocks"}, \code{"heavisine"} or \code{"zero"}.
#'   The last three are the test functions of Donoho and Johnstone (1994),
#'   with their usual constants, seen as functions on the unit interval.
#'
#' @return A function of a numeric vector in \eqn{[0,1]}.
#'
#' @references Donoho, D. L. and Johnstone, I. M. (1994). Ideal spatial
#'   adaptation by wavelet shrinkage. \emph{Biometrika} 81(3), 425-455.
#'
#' @examples
#' g <- wafc_component("bumps")
#' mean(g((seq_len(1024) - 0.5)/1024))
#'
#' @export
wafc_component <- function(name = c("sine", "cosine", "cubic", "bumps",
                                    "blocks", "heavisine", "zero")) {
  name <- match.arg(name)
  if (name == "zero") return(function(u) rep_len(0, length(u)))
  raw <- switch(
    name,
    sine = function(u) sin(2 * pi * u),
    cosine = function(u) cos(4 * pi * u),
    cubic = function(u) u^3 - 1.4 * u^2 + 0.4 * u,
    bumps = function(u) {
      t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
      h <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
      w <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005,
             0.008, 0.005)
      out <- numeric(length(u))
      for (i in seq_along(t)) {
        out <- out + h[i] * (1 + abs(u - t[i]) / w[i])^(-4)
      }
      out
    },
    blocks = function(u) {
      t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
      h <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
      out <- numeric(length(u))
      for (i in seq_along(t)) {
        out <- out + h[i] * (1 + sign(u - t[i])) / 2
      }
      out
    },
    heavisine = function(u) {
      4 * sin(4 * pi * u) - sign(u - 0.3) - sign(0.72 - u)
    }
  )
  ## Centring and normalisation constants on a fixed midpoint grid: the
  ## spikes of bumps are narrow (w = 0.005), so the grid has to be fine.
  grid <- (seq_len(2^16) - 0.5) / 2^16
  v <- raw(grid)
  mu <- mean(v)
  sdv <- sqrt(mean((v - mu)^2))
  function(u) (raw(u) - mu) / sdv
}

#' Active structure of a scenario
#'
#' The matrix of component names of a scenario, with one row per linear
#' covariate and one column per modulating covariate. An empty string marks
#' a block \eqn{(\ell, m)} whose component is zero.
#'
#' @param scenario One of \code{"smooth"}, \code{"inhomogeneous"} or
#'   \code{"null"}.
#' @param p,q Numbers of linear and of modulating covariates.
#'
#' @return A character matrix of dimension \eqn{p \times q}.
#'
#' @examples
#' wafc_scenario("smooth", 3, 2)
#'
#' @export
wafc_scenario <- function(scenario = c("smooth", "inhomogeneous", "null"),
                          p, q) {
  scenario <- match.arg(scenario)
  if (p < 1L || q < 1L) stop("'p' and 'q' must be at least 1.", call. = FALSE)
  out <- matrix("", p, q)
  if (scenario == "null") return(out)
  nm <- if (scenario == "smooth") c("sine", "cubic", "cosine")
        else c("bumps", "blocks", "heavisine")
  ## beta_1 additive in two modulating covariates, beta_2 in one, the
  ## remaining coefficients constant.
  out[1L, 1L] <- nm[1L]
  if (q >= 2L) out[1L, 2L] <- nm[2L]
  if (p >= 2L) out[2L, 1L] <- nm[3L]
  out
}

#' Simulate from the WAFC model
#'
#' Draws a sample from \eqn{Y = \sum_\ell \beta_\ell(U) X_\ell + \varepsilon}
#' with \eqn{\beta_\ell(u) = c_\ell + \sum_m g_{\ell m}(u_m)} and Gaussian
#' errors, in one of the scenarios of \code{\link{wafc_scenario}}.
#'
#' @param n Sample size.
#' @param p Number of linear covariates. The first one is constant equal to
#'   one when \code{intercept = TRUE}, which makes the additive intercept
#'   \eqn{c_1 + \sum_m g_{1m}(U_m)} part of the model.
#' @param q Number of modulating covariates.
#' @param scenario One of \code{"smooth"}, \code{"inhomogeneous"} or
#'   \code{"null"}; see \code{\link{wafc_scenario}}.
#' @param seed Optional integer passed to \code{\link{set.seed}} before the
#'   draw.
#' @param snr Signal to noise ratio: the error standard deviation is set to
#'   \eqn{\mathrm{sd}(f)/snr}, where \eqn{f} is the realised regression
#'   function of the sample. Ignored when \code{sigma} is supplied, and
#'   unusable when the regression function is constant (\code{p = 1} with
#'   an intercept in the null scenario), where \code{sigma} must be given.
#' @param sigma Error standard deviation, overriding \code{snr}.
#' @param x_dist Distribution of the non-constant linear covariates:
#'   \code{"gaussian"} (standard normal) or \code{"uniform"} on
#'   \eqn{[-1,1]}.
#' @param intercept Logical: should \eqn{X_1 \equiv 1}?
#' @param u_dist Distribution of the modulating covariates:
#'   \code{"uniform"} on \eqn{[0,1]} or \code{"beta"}, in which case the
#'   even coordinates are Beta(2,3), a density bounded away from zero and
#'   infinity that is far from uniform.
#' @param u_rho Correlation of the Gaussian copula that couples the
#'   modulating covariates. Zero (the default) makes them independent; a
#'   value in \eqn{(-1,1)} keeps the marginals and makes them dependent,
#'   which is allowed by hypothesis (A2) and by D14.
#' @param cc Vector of levels \eqn{c_\ell}, of length \eqn{p}. The default
#'   recycles \code{c(1, 2, -1.5, 0.5)}.
#' @param amplitude Multiplier of every non-zero component, of length 1 or
#'   \eqn{p \times q}. Since the components are normalised to unit
#'   \eqn{L_2[0,1]} norm, it is the standard deviation each one contributes
#'   to its functional coefficient.
#'
#' @return An object of class \code{"wafc_dgp"}: a list with the response
#'   \code{y}, the covariates \code{x} and \code{u}, the regression function
#'   \code{f} and the matrix \code{beta} of the \eqn{\beta_\ell(U_i)}
#'   evaluated at the sample, the true \code{cc}, the \eqn{p \times q} list
#'   \code{g} of components (\code{NULL} where the block is zero), the
#'   \code{structure} matrix of names, and \code{sigma}, \code{scenario},
#'   \code{seed} and \code{call}.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' str(d$y)
#' colMeans(d$beta)
#'
#' @export
simulate_wafc <- function(n, p = 3L, q = 2L,
                          scenario = c("smooth", "inhomogeneous", "null"),
                          seed = NULL, snr = 4, sigma = NULL,
                          x_dist = c("gaussian", "uniform"),
                          intercept = TRUE,
                          u_dist = c("uniform", "beta"), u_rho = 0,
                          cc = NULL, amplitude = 1) {

  this_call <- match.call()
  scenario <- match.arg(scenario)
  x_dist <- match.arg(x_dist)
  u_dist <- match.arg(u_dist)
  if (length(n) != 1L || !is.finite(n) || n < 1 || n != round(n)) {
    stop("'n' must be a single positive integer.", call. = FALSE)
  }
  p <- as.integer(p)
  q <- as.integer(q)
  if (p < 1L || q < 1L) stop("'p' and 'q' must be at least 1.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(cc)) cc <- rep_len(c(1, 2, -1.5, 0.5), p)
  cc <- wafc_recycle(cc, p, "cc")
  amp <- matrix(wafc_recycle(amplitude, p * q, "amplitude"), p, q)

  ## Modulating covariates: marginals through a Gaussian copula, so that
  ## dependence between them can be switched on without changing the
  ## marginal laws (D14).
  if (length(u_rho) != 1L || !is.finite(u_rho) || abs(u_rho) >= 1) {
    stop("'u_rho' must be a single value in (-1, 1).", call. = FALSE)
  }
  if (u_rho == 0 || q == 1L) {
    w <- matrix(stats::runif(n * q), n, q)
  } else {
    R <- matrix(u_rho, q, q)
    diag(R) <- 1
    L <- chol(R)
    w <- stats::pnorm(matrix(stats::rnorm(n * q), n, q) %*% L)
  }
  u <- w
  if (u_dist == "beta" && q >= 2L) {
    even <- seq(2L, q, by = 2L)
    u[, even] <- stats::qbeta(w[, even], 2, 3)
  }
  colnames(u) <- paste0("u", seq_len(q))

  ## Linear covariates.
  x <- switch(x_dist,
              gaussian = matrix(stats::rnorm(n * p), n, p),
              uniform = matrix(stats::runif(n * p, -1, 1), n, p))
  if (intercept) x[, 1L] <- 1
  colnames(x) <- paste0("x", seq_len(p))

  ## Functional coefficients.
  struct <- wafc_scenario(scenario, p, q)
  g <- vector("list", p * q)
  dim(g) <- c(p, q)
  beta <- matrix(rep(cc, each = n), n, p)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      if (!nzchar(struct[l, m])) next
      gl <- wafc_component(struct[l, m])
      a <- amp[l, m]
      g[[l, m]] <- local({
        gg <- gl
        aa <- a
        function(v) aa * gg(v)
      })
      beta[, l] <- beta[, l] + g[[l, m]](u[, m])
    }
  }
  colnames(beta) <- colnames(x)

  f <- as.numeric(rowSums(x * beta))
  if (is.null(sigma)) {
    sf <- stats::sd(f)
    if (!is.finite(sf) || sf <= 0) {
      stop("The regression function is constant in this sample; supply 'sigma' ",
           "instead of 'snr'.", call. = FALSE)
    }
    if (length(snr) != 1L || !is.finite(snr) || snr <= 0) {
      stop("'snr' must be a single positive value.", call. = FALSE)
    }
    sigma <- sf / snr
  }
  if (length(sigma) != 1L || !is.finite(sigma) || sigma < 0) {
    stop("'sigma' must be a single non-negative value.", call. = FALSE)
  }
  y <- f + stats::rnorm(n, sd = sigma)

  out <- list(y = y, x = x, u = u, f = f, beta = beta, cc = cc, g = g,
              structure = struct, sigma = sigma, scenario = scenario,
              seed = seed, n = as.integer(n), p = p, q = q, call = this_call)
  class(out) <- "wafc_dgp"
  out
}

#' Functional coefficients of a simulated model
#'
#' Evaluates \eqn{\beta_\ell(u) = c_\ell + \sum_m g_{\ell m}(u_m)} of an
#' object returned by \code{\link{simulate_wafc}} at a new matrix of
#' modulating covariates.
#'
#' @param object An object of class \code{"wafc_dgp"}.
#' @param u Matrix (or vector) of modulating covariates with \eqn{q}
#'   columns, the number used in the simulation.
#'
#' @return A matrix with one row per row of \code{u} and one column per
#'   linear covariate.
#'
#' @examples
#' d <- simulate_wafc(100, seed = 1)
#' wafc_beta(d, matrix(c(0.2, 0.8), 1, 2))
#'
#' @export
wafc_beta <- function(object, u) {
  if (!inherits(object, "wafc_dgp")) {
    stop("'object' must be an object returned by simulate_wafc().",
         call. = FALSE)
  }
  u <- wafc_as_matrix(u, "u")
  if (ncol(u) != object[["q"]]) {
    stop("'u' must have ", object[["q"]], " column(s).", call. = FALSE)
  }
  p <- object[["p"]]
  beta <- matrix(rep(object[["cc"]], each = nrow(u)), nrow(u), p)
  for (l in seq_len(p)) {
    for (m in seq_len(object[["q"]])) {
      gm <- object[["g"]][[l, m]]
      if (is.null(gm)) next
      beta[, l] <- beta[, l] + gm(u[, m])
    }
  }
  colnames(beta) <- colnames(object[["x"]])
  beta
}
