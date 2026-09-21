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
## Four scenarios, in the terms of docs/plano-projeto.md, E2.1, plus
## the one decision D30 adds:
##
##   "smooth"          sine, cosine and a cubic: functions in every Besov
##                     space, where an expansion of growing dimension is
##                     already efficient;
##   "uneven"          a mixture of a wide and a narrow Gaussian, a chirp
##                     and a cosine: every component is C^infinity, and so
##                     inside the hypothesis of Xue and Yang (2006), but
##                     the scale of two of them varies along the domain.
##                     It is the scenario decision D30 asks for: the smooth
##                     case where a single smoothing parameter has to
##                     compromise and an adaptive basis need not;
##   "inhomogeneous"   bumps, blocks and heavisine of Donoho and Johnstone
##                     (1994), rescaled to [0,1]: spatially inhomogeneous
##                     regularity, which is where the wavelet LASSO is
##                     supposed to win;
##   "null"            every g_{lm} is zero, so the model is a linear
##                     regression and the whole penalized part must be
##                     shrunk away.
##
## The active structure is the same in the first three: beta_1 depends on
## two modulating covariates (the additive structure that is the point of
## the model), beta_2 on one, and beta_l is constant for l >= 3.
##
## Each scenario also declares the effective regularity s' of decision D27,
## through wafc_sprime(); see the note on that function for why the number
## has to carry the regime it is read in.

#' A window that vanishes to every order at the endpoints
#'
#' The \eqn{C^\infty} partition-of-unity profile of Step 2 of Lemma 10
#' (\file{derivations/02-aproximacao-besov.tex}), used here to make a
#' component and all of its derivatives vanish at \eqn{0} and at \eqn{1}.
#' A component built this way has a \eqn{C^\infty} periodic extension, so
#' the effective regularity the wavelet basis reads is the one of the
#' function itself and not the one of a seam.
#'
#' @param u A numeric vector in \eqn{[0,1]}.
#' @param d Width of each transition. The window is one on
#'   \eqn{[d, 1-d]}.
#'
#' @return A numeric vector of the same length as \code{u}, with values in
#'   \eqn{[0,1]}.
#'
#' @examples
#' wafc_window(c(0, 0.04, 0.5, 1))
#'
#' @export
wafc_window <- function(u, d = 0.08) {
  if (length(d) != 1L || !is.finite(d) || d <= 0 || d >= 0.5) {
    stop("'d' must be a single value in (0, 0.5).", call. = FALSE)
  }
  step <- function(t) {
    a <- exp(-1 / pmax(t, 1e-300))
    b <- exp(-1 / pmax(1 - t, 1e-300))
    ifelse(t <= 0, 0, ifelse(t >= 1, 1, a / (a + b)))
  }
  step(u / d) * step((1 - u) / d)
}

#' Additive components of the data generating processes
#'
#' Returns one of the univariate functions used by
#' \code{\link{simulate_wafc}}, centred and normalised on \eqn{[0,1]}: the
#' returned function \eqn{g} satisfies \eqn{\int_0^1 g = 0} and
#' \eqn{\int_0^1 g^2 = 1}, both up to the accuracy of the grid used to
#' compute the two constants (\eqn{2^{16}} midpoints).
#'
#' \code{"gaussians"} and \code{"chirp"} are the components of uneven
#' curvature decision D30 asks for. Both are \eqn{C^\infty} on the whole
#' line, so both are inside the hypothesis of Xue and Yang (2006) that a
#' penalized spline is efficient under, and both have a scale that changes
#' along the domain: \code{"gaussians"} superposes a bump of standard
#' deviation \eqn{0.10} on one of \eqn{0.018}, and \code{"chirp"} sweeps
#' the frequency from one to eight cycles over the interval. Both are
#' multiplied by the \eqn{C^\infty} window \code{\link{wafc_window}},
#' which vanishes to every order at the two endpoints, so that the periodic
#' extension is \eqn{C^\infty} too and the scenario measures uneven
#' curvature and not a seam: the corner of the cubic is the reason the
#' smooth scenario reads \eqn{s' = 3/2} (see \code{\link{wafc_sprime}}),
#' and repeating it here would confound the two effects.
#'
#' @param name One of \code{"sine"}, \code{"cosine"}, \code{"cubic"},
#'   \code{"gaussians"}, \code{"chirp"}, \code{"bumps"},
#'   \code{"blocks"}, \code{"heavisine"} or \code{"zero"}. The last three
#'   are the test functions of Donoho and Johnstone (1994), with their usual
#'   constants, seen as functions on the unit interval.
#'
#' @return A function of a numeric vector in \eqn{[0,1]}.
#'
#' @references Donoho, D. L. and Johnstone, I. M. (1994). Ideal spatial
#'   adaptation by wavelet shrinkage. \emph{Biometrika} 81(3), 425-455.
#'
#'   Xue, L. and Yang, L. (2006). Additive coefficient modeling via
#'   polynomial spline. \emph{Statistica Sinica} 16(4), 1423-1446.
#'
#' @examples
#' g <- wafc_component("bumps")
#' mean(g((seq_len(1024) - 0.5)/1024))
#'
#' @export
wafc_component <- function(name = c("sine", "cosine", "cubic", "gaussians",
                                    "chirp", "bumps", "blocks", "heavisine",
                                    "zero")) {
  name <- match.arg(name)
  if (name == "zero") return(function(u) rep_len(0, length(u)))
  raw <- switch(
    name,
    sine = function(u) sin(2 * pi * u),
    cosine = function(u) cos(4 * pi * u),
    cubic = function(u) u^3 - 1.4 * u^2 + 0.4 * u,
    ## Uneven curvature, decision D30. The wide bump is five and a half
    ## times the narrow one, which is the ratio a single smoothing
    ## parameter has to split the difference between.
    gaussians = function(u) {
      wafc_window(u) * (exp(-(u - 0.40)^2 / (2 * 0.10^2)) -
                          0.8 * exp(-(u - 0.72)^2 / (2 * 0.018^2)))
    },
    ## Linear chirp: the instantaneous frequency goes from f0 = 1 to
    ## f1 = 8 cycles per unit, so the local scale varies by a factor of
    ## eight while the function stays analytic.
    chirp = function(u) {
      wafc_window(u) * sin(2 * pi * (1 * u + (8 - 1) * u^2 / 2))
    },
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
#' The returned matrix carries the effective regularity of the scenario in
#' the attribute \code{"sprime"}, with the regime it is read in in
#' \code{"regime"}; see \code{\link{wafc_sprime}}.
#'
#' @param scenario One of \code{"smooth"}, \code{"uneven"},
#'   \code{"inhomogeneous"} or \code{"null"}.
#' @param p,q Numbers of linear and of modulating covariates.
#' @param regime The regime the attribute \code{"sprime"} is read in, as
#'   in \code{\link{wafc_sprime}}.
#'
#' @return A character matrix of dimension \eqn{p \times q}, with the
#'   attributes \code{"sprime"} and \code{"regime"}.
#'
#' @examples
#' wafc_scenario("smooth", 3, 2)
#' attr(wafc_scenario("uneven", 3, 2), "sprime")
#'
#' @export
wafc_scenario <- function(scenario = c("smooth", "uneven", "inhomogeneous",
                                       "null"),
                          p, q, regime = c("periodic", "margin", "interval")) {
  scenario <- match.arg(scenario)
  regime <- match.arg(regime)
  if (p < 1L || q < 1L) stop("'p' and 'q' must be at least 1.", call. = FALSE)
  out <- matrix("", p, q)
  if (scenario != "null") {
    nm <- switch(scenario,
                 smooth = c("sine", "cubic", "cosine"),
                 uneven = c("gaussians", "chirp", "cosine"),
                 inhomogeneous = c("bumps", "blocks", "heavisine"))
    ## beta_1 additive in two modulating covariates, beta_2 in one, the
    ## remaining coefficients constant.
    out[1L, 1L] <- nm[1L]
    if (q >= 2L) out[1L, 2L] <- nm[2L]
    if (p >= 2L) out[2L, 1L] <- nm[3L]
  }
  attr(out, "sprime") <- wafc_sprime_value(scenario, regime)
  attr(out, "regime") <- regime
  out
}

## The table of declared regularities, and the only place it is written.
## wafc_scenario() reads it through this internal name and not through the
## exported wafc_sprime(), so that a script that binds an object of its own
## to that name, as wafc/scripts/04-pilot.R still does, shadows the
## accessor and not the data generating process.
wafc_sprime_value <- function(scenario, regime) {
  tab <- rbind(smooth        = c(periodic = 3/2, margin = 4, interval = 4),
               uneven        = c(4,              4,          4),
               inhomogeneous = c(1/2,            1/2,        1/2),
               null          = c(3/2,            3/2,        3/2))
  out <- tab[scenario, regime]
  attr(out, "regime") <- regime
  out
}

#' Effective regularity declared by a scenario
#'
#' The \eqn{s'} of decision D27: the regularity the wavelet expansion
#' actually reads off the components of a scenario, which is what the
#' penalty level of \code{\link{wafc_lambda_theory}} and the resolution
#' rule of \code{\link{wafc_J_theory}} consume. It is the smallest value
#' over the components the scenario activates, since the slowest one sets
#' the approximation error of the sum.
#'
#' The number cannot be stated without the regime it is read in, which is
#' the point decision D27 and Lemma 12 (\file{derivations/
#' 07-rota-intervalo.tex}) make: a component may be smooth on the interval
#' and still have a corner, or a jump, where the periodic extension joins
#' \eqn{1} to \eqn{0}. Three regimes are distinguished, and they are the
#' three of that lemma:
#'
#' \describe{
#'   \item{\code{"periodic"}}{the periodized basis with no margin, which
#'     is where the theory of the manuscript is stated (decisions D23 and
#'     D26) and the regime the numbers of step E2.3 were read in;}
#'   \item{\code{"margin"}}{the periodized basis with the fixed margin of
#'     \code{\link{wafc_eps}}, where the component is extended off its
#'     support by the \eqn{C^\infty} cut of Lemma 10, so a corner at the
#'     seam is bought away;}
#'   \item{\code{"interval"}}{the boundary-corrected basis of
#'     \code{boundary = "interval"}, where there is no seam to begin with.}
#' }
#'
#' The declared values are the ones measured in part (C2) of
#' \file{derivations/check/07-rota-intervalo.R}, rounded to the value of
#' the table of that file: \code{"smooth"} reads \eqn{3/2} under
#' \code{"periodic"}, because the cubic has a corner at the seam
#' (\eqn{g'(0) = 0.4} against \eqn{g'(1) = 0.6}), and \eqn{4} under the
#' other two, where that corner is bought away; \code{"inhomogeneous"}
#' reads \eqn{1/2} under all three, because \code{blocks} and
#' \code{heavisine} jump inside the interval, where no basis helps; and
#' \code{"uneven"} reads \eqn{4} under all three, because its components
#' are windowed by \code{\link{wafc_window}} and have no seam. The cap of
#' \eqn{4} is the one of the filter in use, not a property of the
#' functions: an analytic component cannot be read as decaying faster than
#' the number of vanishing moments allows.
#'
#' The null scenario has no component, so no regularity is defined for it.
#' The value returned, \eqn{3/2}, is the convention step E2.4 used when it
#' had to give \code{rule = "theory"} some number in that cell, and it is
#' recorded here so that the two are not chosen independently.
#'
#' @param scenario One of \code{"smooth"}, \code{"uneven"},
#'   \code{"inhomogeneous"} or \code{"null"}.
#' @param regime \code{"periodic"} (the default), \code{"margin"} or
#'   \code{"interval"}.
#'
#' @return A single value, with the attribute \code{"regime"}.
#'
#' @seealso \code{\link{wafc_J_theory}},
#'   \code{\link{wafc_lambda_theory}}.
#'
#' @examples
#' wafc_sprime("smooth")
#' wafc_sprime("smooth", "margin")
#'
#' @export
wafc_sprime <- function(scenario = c("smooth", "uneven", "inhomogeneous",
                                     "null"),
                        regime = c("periodic", "margin", "interval")) {
  scenario <- match.arg(scenario)
  regime <- match.arg(regime)
  wafc_sprime_value(scenario, regime)
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
#' @param scenario One of \code{"smooth"}, \code{"uneven"},
#'   \code{"inhomogeneous"} or \code{"null"}; see
#'   \code{\link{wafc_scenario}}.
#' @param regime The regime the declared \code{sprime} is read in; see
#'   \code{\link{wafc_sprime}}. It labels a number, and changes nothing
#'   about the sample.
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
#'   \code{structure} matrix of names, the declared effective regularity
#'   \code{sprime} of \code{\link{wafc_sprime}} with its \code{regime},
#'   and \code{sigma}, \code{scenario}, \code{seed} and \code{call}.
#'
#' @examples
#' d <- simulate_wafc(200, p = 3, q = 2, scenario = "smooth", seed = 1)
#' str(d$y)
#' colMeans(d$beta)
#'
#' @export
simulate_wafc <- function(n, p = 3L, q = 2L,
                          scenario = c("smooth", "uneven", "inhomogeneous",
                                       "null"),
                          seed = NULL, snr = 4, sigma = NULL,
                          x_dist = c("gaussian", "uniform"),
                          intercept = TRUE,
                          u_dist = c("uniform", "beta"), u_rho = 0,
                          cc = NULL, amplitude = 1,
                          regime = c("periodic", "margin", "interval")) {

  this_call <- match.call()
  scenario <- match.arg(scenario)
  regime <- match.arg(regime)
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
  struct <- wafc_scenario(scenario, p, q, regime = regime)
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
              structure = struct, sprime = as.numeric(attr(struct, "sprime")),
              regime = regime, sigma = sigma, scenario = scenario,
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
