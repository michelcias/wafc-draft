## wafc/R/plot.R -- the plots of a WAFC fit (step E3.2), in base graphics and
## with no package beyond grDevices, graphics and utils.
##
## Four panels, each drawn by an internal function from numbers the methods
## prepare, so that plot.wafc() and plot.cv.wafc() draw the same panel from
## the same quantities:
##
##   "components"  one panel per block (l, m) with the reconstructed
##                 component g-hat_{lm} of wafc_functions(), the thresholded
##                 one when the object has a threshold, the blocks that are
##                 zero shaded, and the true component when it is supplied;
##   "path"        the norm nu-hat_{lm} of every block, the statistic of
##                 wafc_blocks() and of the threshold, against log(lambda)
##                 along the path;
##   "cv"          the cross-validated loss of cv.wafc() against log(lambda),
##                 one curve per candidate J, with lambda.min and lambda.1se
##                 of the selected J;
##   "threshold"   the cross-validated error of the thresholded fit against
##                 t, with the t chosen and the band of one standard error
##                 that the rule "cv1se" reads.
##
## The true component is shifted by a constant before it is drawn, so that
## its mean on the grid is the mean of the estimate there. The level of a
## component is a convention and not part of its shape: the estimate
## integrates to zero over the rescaled range of U_m and the components of
## simulate_wafc() over [0, 1], and the difference belongs to c_l (header of
## wafc/R/reconstruct.R). The estimate itself is drawn as wafc_functions()
## returns it. wafc/scripts/01-smoke.R centres both, which is the same
## comparison.
##
## Every method returns invisibly what it drew, a list with one entry per
## panel, so that a figure can be redrawn or checked without the device.

#' Plot a WAFC fit
#'
#' Two plots of a fit of \code{\link{wafc}}: the reconstructed components
#' \eqn{\hat g_{\ell m}}, one panel per block \eqn{(\ell, m)}, at one penalty
#' level, and the path of the block norms \eqn{\hat\nu_{\ell m}} against
#' \eqn{\lambda}.
#'
#' With \code{which = "components"} the device is divided into \eqn{p} rows
#' and \eqn{q} columns, one per block, and restored afterwards. Each panel
#' draws \eqn{\hat g_{\ell m}} on the grid of \code{\link{wafc_functions}},
#' with its norm \eqn{\hat\nu_{\ell m}} in the upper right corner; a block
#' whose coefficients are all zero is shaded and drawn in grey. When
#' \code{truth} is supplied the true component is drawn dashed, shifted by
#' a constant so that its mean on the grid is the one of the estimate: the
#' level of a component is a convention, and only its shape is compared
#' (the header of \file{wafc/R/plot.R}).
#'
#' With \code{which = "path"} a single panel draws
#' \eqn{\hat\nu_{\ell m}}, the Euclidean norm of the wavelet coefficients of
#' the block, the statistic of \code{\link{wafc_blocks}} and of
#' \code{\link{wafc_threshold}}, against \eqn{\log\lambda}, one line per
#' block, with \code{s} marked when it is given.
#'
#' @param x An object of class \code{"wafc"}.
#' @param which \code{"components"} (the default), \code{"path"}, or both.
#' @param s A single penalty level, on the scale of the objective of
#'   \code{\link{wafc}}. The components are drawn there; \code{NULL} takes
#'   the smallest level of the path, as \code{\link{wafc_functions}} does.
#'   The path marks it with a vertical line when it is given.
#' @param truth The true components, for a simulation: an object returned
#'   by \code{\link{simulate_wafc}}, or a \eqn{p} by \eqn{q} list of
#'   functions of one variable, with \code{NULL} where the component is
#'   zero. \code{NULL} (the default) draws the estimate alone.
#' @param grid,n_grid The grid of the components, as in
#'   \code{\link{wafc_functions}}.
#' @param ask Whether to ask before each new page; \code{NULL} asks when
#'   there are more panels than the layout of the device holds and the
#'   device is interactive, as \code{\link[stats]{plot.lm}} does.
#' @param ... Graphical parameters passed to \code{\link[graphics]{plot}}
#'   for every panel, which override the defaults (\code{main},
#'   \code{xlab}, \code{ylim}, ...).
#'
#' @return Invisibly, a list with one entry per panel drawn.
#'   \code{components}: the penalty level \code{s}, the \code{grid}, the
#'   \eqn{p} by \eqn{q} lists \code{g} of the estimate and \code{truth} of
#'   the true components as drawn (\code{NULL} without \code{truth}), the
#'   \eqn{p} by \eqn{q} matrices \code{zero} (the blocks drawn as zero) and
#'   \code{norm} (\eqn{\hat\nu_{\ell m}}), and the threshold \code{t}
#'   (\code{NULL} here). \code{path}: the path \code{lambda}, the matrix
#'   \code{norm} of the block norms (one row per penalty level, one column
#'   per block), the levels \code{marked} with a vertical line and the
#'   threshold \code{t} (\code{NULL} here).
#'
#' @seealso \code{\link{plot.cv.wafc}}, \code{\link{wafc_functions}},
#'   \code{\link{wafc_blocks}}.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' fit <- wafc(d$x, d$u, d$y, J = 4)
#' plot(fit, which = "path", s = fit$lambda[30])
#' plot(fit, s = fit$lambda[30], truth = d)
#'
#' @export
plot.wafc <- function(x, which = "components", s = NULL, truth = NULL,
                      grid = NULL, n_grid = 512L, ask = NULL, ...) {
  which <- wafc_plot_which(which, c("components", "path"))
  dots <- list(...)
  ask <- wafc_plot_ask(ask, length(which))
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }
  out <- list()
  for (w in which) {
    if (w == "components") {
      fn <- wafc_functions(x, s = s, grid = grid, n_grid = n_grid)
      th <- x[["threshold"]]
      if (is.null(th)) {
        zero <- fn[["nonzero"]] == 0L
        note <- "shaded: zero block"
      } else {
        ## the thresholded fit of a "cv.wafc" read as a "wafc" object
        ## (wafc_cv_thr_view() of wafc/R/tune.R)
        zero <- !th[["kept"]]
        note <- sprintf("shaded: set to zero by the threshold \"%s\"",
                        th[["rule"]])
      }
      out[["components"]] <- wafc_plot_components(
        fn, zero = zero, norm = fn[["norm"]],
        truth = wafc_plot_truth(truth, fn),
        title = wafc_plot_title(fn[["s"]], th, note, !is.null(truth),
                                J = wafc_J_text(x[["design"]], unique = TRUE)),
        t = if (is.null(th)) NULL else th[["t"]], dots = dots)
    } else {
      marked <- if (is.null(s)) NULL else c(s = wafc_single_s(x, s))
      out[["path"]] <- wafc_plot_path(
        x, marked = marked, t = NULL, dots = dots,
        main = sprintf("Block norms along the path, J = %s",
                       wafc_J_text(x[["design"]], unique = TRUE)))
    }
  }
  invisible(out)
}

#' Plot a cross-validated WAFC fit
#'
#' Four plots of a fit of \code{\link{cv.wafc}}: the cross-validation over
#' \eqn{(J, \lambda)}, the cross-validation of the threshold, the
#' components of the fit, and the path of the block norms at the selected
#' \eqn{J}.
#'
#' \code{"cv"} draws the cross-validated loss against \eqn{\log\lambda},
#' one curve per candidate \eqn{J}, the minimum of each curve as a point,
#' the band of one standard error (\code{cvlo} to \code{cvup}) of the
#' selected \eqn{J}, and \code{lambda.min} and \code{lambda.1se} of that
#' \eqn{J} as vertical lines.
#'
#' \code{"threshold"} draws the cross-validated squared error of the
#' thresholded fit against the threshold \eqn{t}, a step function that is
#' constant between two consecutive norms of the fold fits (the candidates
#' of \code{\link{wafc_threshold}}), with the \eqn{t} chosen as a vertical
#' line, the norms \eqn{\hat\nu_{\ell m}} of the blocks of the fit as ticks
#' on the horizontal axis (black for the blocks kept, grey for the others),
#' and the band from the smallest error to one standard error above it,
#' whose largest \eqn{t} is the one the rule \code{"cv1se"} takes. It needs
#' a rule that scores candidates, \code{"cv1se"} or \code{"cv"}; the rule
#' \code{"max"} scores none, and \code{threshold = "none"} has no threshold.
#'
#' \code{"components"} draws the components at \code{s}, as
#' \code{\link{plot.wafc}} does: at \code{lambda.min} the thresholded ones
#' when the object has a threshold (the WAFC estimator), with the blocks the
#' threshold set to zero shaded and the norm of each block before the
#' threshold in its corner. \code{"path"} draws the path of the block norms
#' of \code{wafc.fit}, the fit at the selected \eqn{J}, with
#' \code{lambda.min} and \code{lambda.1se} as vertical lines and the
#' threshold \eqn{t} as a horizontal one: the blocks whose norm at
#' \code{lambda.min} is above that line are the ones kept.
#'
#' @param x An object of class \code{"cv.wafc"}.
#' @param which Any of \code{"cv"}, \code{"threshold"}, \code{"components"}
#'   and \code{"path"}. \code{NULL} (the default) is \code{"cv"} and, when
#'   the threshold scored candidates, \code{"threshold"}.
#' @param s For \code{"components"}, \code{"lambda.min"} (the default),
#'   \code{"lambda.1se"} or a numeric penalty level, as in
#'   \code{\link{coef.cv.wafc}}.
#' @param thresholded For \code{"components"}, as in
#'   \code{\link{coef.cv.wafc}}.
#' @param truth,grid,n_grid,ask,... As in \code{\link{plot.wafc}}.
#'
#' @return Invisibly, a list with one entry per panel drawn.
#'   \code{cv}: a data frame \code{table} with the candidate \code{J}, the
#'   \code{lambda}, \code{cvm}, \code{cvsd} and \code{nzero} of every point
#'   of every curve and the flag \code{selected}, and the \code{J.min},
#'   \code{lambda.min} and \code{lambda.1se} marked. \code{threshold}: the
#'   \code{candidates} of the rule (\code{lower}, \code{upper},
#'   \code{error}, \code{se}), the \code{t} chosen, the \code{rule}, the
#'   matrices \code{norm} and \code{kept}, and the \code{band}, the smallest
#'   error and that error plus its standard error. \code{components} and
#'   \code{path}: as in \code{\link{plot.wafc}}, with the threshold
#'   \code{t} when the components drawn are the thresholded ones (for
#'   \code{components}) or when the object has a threshold (for
#'   \code{path}).
#'
#' @seealso \code{\link{plot.wafc}}, \code{\link{cv.wafc}},
#'   \code{\link{wafc_threshold}}.
#'
#' @examples
#' d <- simulate_wafc(300, p = 3, q = 2, scenario = "smooth", seed = 1)
#' cvfit <- cv.wafc(d$x, d$u, d$y, J = 2:4, nfolds = 5)
#' op <- par(mfrow = c(1, 2))
#' plot(cvfit)
#' par(op)
#' plot(cvfit, which = "components", truth = d)
#' plot(cvfit, which = "path")
#'
#' @export
plot.cv.wafc <- function(x, which = NULL, s = "lambda.min",
                         thresholded = NULL, truth = NULL, grid = NULL,
                         n_grid = 512L, ask = NULL, ...) {
  th <- x[["threshold"]]
  scored <- !is.null(th) && !is.null(th[["candidates"]])
  if (is.null(which)) which <- c("cv", if (scored) "threshold")
  which <- wafc_plot_which(which, c("cv", "threshold", "components", "path"))
  if ("threshold" %in% which && !scored) {
    stop(if (is.null(th)) {
      "This \"cv.wafc\" object has no threshold (threshold = \"none\"); "
    } else {
      sprintf("The threshold rule \"%s\" scores no candidate t; ",
              th[["rule"]])
    }, "there is no threshold curve to draw.", call. = FALSE)
  }
  dots <- list(...)
  ask <- wafc_plot_ask(ask, length(which))
  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }
  out <- list()
  for (w in which) {
    if (w == "cv") {
      out[["cv"]] <- wafc_plot_cv(x, dots)
    } else if (w == "threshold") {
      out[["threshold"]] <- wafc_plot_thr(th, dots)
    } else if (w == "components") {
      f <- wafc_cv_fit(x, s, thresholded)
      fn <- wafc_functions(f[["fit"]], s = f[["s"]], grid = grid,
                           n_grid = n_grid)
      thr_used <- !is.null(f[["fit"]][["threshold"]])
      if (thr_used) {
        zero <- !th[["kept"]]
        norm <- th[["norm"]]
        note <- sprintf("shaded: set to zero by the threshold \"%s\"",
                        th[["rule"]])
      } else {
        zero <- fn[["nonzero"]] == 0L
        norm <- fn[["norm"]]
        note <- "shaded: zero block"
      }
      out[["components"]] <- wafc_plot_components(
        fn, zero = zero, norm = norm, truth = wafc_plot_truth(truth, fn),
        title = wafc_plot_title(fn[["s"]], if (thr_used) th, note,
                                !is.null(truth),
                                J = wafc_J_text(x[["wafc.fit"]][["design"]],
                                                unique = TRUE)),
        t = if (thr_used) th[["t"]] else NULL, dots = dots)
    } else {
      out[["path"]] <- wafc_plot_path(
        x[["wafc.fit"]],
        marked = c(lambda.min = x[["lambda.min"]],
                   lambda.1se = x[["lambda.1se"]]),
        t = if (is.null(th)) NULL else th[["t"]], dots = dots,
        main = sprintf("Block norms along the path, J = %s",
                       wafc_J_text(x[["wafc.fit"]][["design"]],
                                   unique = TRUE)))
    }
  }
  invisible(out)
}

## ---------------------------------------------------------------------------
## Internals
## ---------------------------------------------------------------------------

## The colours: the estimate in black, the truth and the threshold in the
## vermilion of Okabe and Ito, the curves of several J or blocks in the
## qualitative palette "Dark 3" of grDevices::hcl.colors().
wafc_plot_accent <- "#D55E00"
wafc_plot_palette <- function(k) grDevices::hcl.colors(max(k, 1L), "Dark 3")[seq_len(k)]

wafc_plot_which <- function(which, choices) {
  if (!is.character(which) || length(which) == 0L || anyNA(which)) {
    stop("'which' must be a character vector of ",
         paste0("\"", choices, "\"", collapse = ", "), ".", call. = FALSE)
  }
  bad <- setdiff(which, choices)
  if (length(bad) > 0L) {
    stop("unknown panel(s) in 'which': ", paste(bad, collapse = ", "),
         ". Choose among ", paste0("\"", choices, "\"", collapse = ", "),
         ".", call. = FALSE)
  }
  unique(which)
}

wafc_plot_ask <- function(ask, npanels) {
  if (is.null(ask)) {
    return(npanels > 1L && prod(graphics::par("mfcol")) < npanels &&
             grDevices::dev.interactive())
  }
  if (!is.logical(ask) || length(ask) != 1L || is.na(ask)) {
    stop("'ask' must be NULL, TRUE or FALSE.", call. = FALSE)
  }
  ask
}

## The range of the finite values: a penalty level of zero, which a path
## passed by the caller may hold, has no logarithm and is left off the axis.
wafc_plot_range <- function(v) range(v[is.finite(v)])

## Vertical lines at the penalty levels 'marked', labelled at the top of the
## panel: the first label to the left of its line and the second to the
## right, since lambda.min is at or below lambda.1se and the two are often
## close.
wafc_plot_marks <- function(marked) {
  if (length(marked) == 0L) return(invisible(NULL))
  usr <- graphics::par("usr")
  graphics::abline(v = log(marked), lty = c(2, 3)[seq_along(marked)],
                   col = "grey30")
  graphics::text(log(marked), usr[4L] - 0.04 * diff(usr[3:4]), names(marked),
                 pos = c(2, 4)[seq_along(marked)], offset = 0.2, cex = 0.7,
                 col = "grey30")
  invisible(NULL)
}

## The arguments of a call to graphics::plot(), with the ones the caller
## passed through '...' taking precedence over the defaults.
wafc_plot_args <- function(defaults, dots) utils::modifyList(defaults, dots)

## The true components on the grid of the estimate, each shifted to the
## mean of the estimate there (header of this file). 'truth' is a
## "wafc_dgp" object or a p by q list of functions with NULL for a zero
## component.
wafc_plot_truth <- function(truth, fn) {
  if (is.null(truth)) return(NULL)
  p <- length(fn[["xnames"]])
  q <- length(fn[["unames"]])
  if (inherits(truth, "wafc_dgp")) truth <- truth[["g"]]
  if (!is.list(truth) || !identical(as.integer(dim(truth)), c(p, q))) {
    stop("'truth' must be an object returned by simulate_wafc() or a ",
         p, " by ", q, " list of functions, with NULL for a zero component.",
         call. = FALSE)
  }
  out <- vector("list", p * q)
  dim(out) <- c(p, q)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      gm <- truth[[l, m]]
      v <- fn[["grid"]][, m]
      gt <- if (is.null(gm)) numeric(length(v)) else {
        if (!is.function(gm)) {
          stop("'truth[[", l, ", ", m, "]]' must be a function or NULL.",
               call. = FALSE)
        }
        as.numeric(gm(v))
      }
      if (length(gt) != length(v)) {
        stop("'truth[[", l, ", ", m, "]]' must return one value per point ",
             "of the grid.", call. = FALSE)
      }
      out[[l, m]] <- gt - mean(gt) + mean(fn[["g"]][[l, m]])
    }
  }
  dimnames(out) <- list(fn[["xnames"]], fn[["unames"]])
  out
}

## The line of text above the grid of components.
wafc_plot_title <- function(s, th, note, has_truth, J = NULL) {
  first <- sprintf("WAFC components at lambda = %s%s%s",
                   format(s, digits = 3),
                   if (is.null(J)) "" else
                     sprintf(", J = %s", paste(unique(J), collapse = ", ")),
                   if (is.null(th)) "" else
                     sprintf(", threshold t = %s", format(th[["t"]], digits = 3)))
  second <- paste0("solid: estimate", if (has_truth) "; dashed: truth" else "",
                   "; ", note)
  c(first, second)
}

## One panel per block (l, m), on a p by q layout that is restored on exit.
wafc_plot_components <- function(fn, zero, norm, truth, title, t, dots) {
  p <- length(fn[["xnames"]])
  q <- length(fn[["unames"]])
  op <- graphics::par(mfrow = c(p, q), mar = c(3.6, 3.6, 2.2, 0.8),
                      mgp = c(2.2, 0.7, 0), oma = c(0, 0, 3, 0))
  on.exit(graphics::par(op))
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      v <- fn[["grid"]][, m]
      gh <- fn[["g"]][[l, m]]
      gt <- if (is.null(truth)) NULL else truth[[l, m]]
      yl <- range(c(gh, gt))
      if (diff(yl) < 1e-8) yl <- yl + c(-0.5, 0.5)
      do.call(graphics::plot,
              wafc_plot_args(list(x = range(v), y = yl, type = "n",
                                  xlab = fn[["unames"]][m], ylab = "",
                                  main = sprintf("%s : %s", fn[["xnames"]][l],
                                                 fn[["unames"]][m])),
                             dots))
      if (zero[l, m]) {
        usr <- graphics::par("usr")
        graphics::rect(usr[1L], usr[3L], usr[2L], usr[4L], col = "grey93",
                       border = NA)
        graphics::box()
      }
      graphics::abline(h = 0, col = "grey80")
      if (!is.null(gt)) {
        graphics::lines(v, gt, col = wafc_plot_accent, lty = 2, lwd = 1.5)
      }
      graphics::lines(v, gh, col = if (zero[l, m]) "grey50" else "black",
                      lwd = 1.5)
      graphics::mtext(bquote(hat(nu) == .(format(norm[l, m], digits = 3))),
                      side = 3, line = 0.1, adj = 1, cex = 0.7,
                      col = if (zero[l, m]) "grey40" else "black")
    }
  }
  graphics::mtext(title[1L], side = 3, line = 1.5, outer = TRUE, cex = 0.9)
  graphics::mtext(title[2L], side = 3, line = 0.3, outer = TRUE, cex = 0.75)
  list(s = fn[["s"]], grid = fn[["grid"]], g = fn[["g"]], truth = truth,
       zero = zero, norm = norm, t = t)
}

## The block norms along the path of a "wafc" object, against log(lambda).
wafc_plot_path <- function(object, marked, t, dots, main) {
  des <- object[["design"]]
  cf <- wafc_raw_coef(object)[-1L, , drop = FALSE]
  nb <- names(des[["blocks"]])
  nrm <- vapply(nb, function(b) {
    idx <- des[["blocks"]][[b]]
    sqrt(colSums(cf[idx, , drop = FALSE]^2))
  }, numeric(ncol(cf)))
  nrm <- matrix(nrm, ncol = length(nb), dimnames = list(NULL, nb))
  lam <- object[["lambda"]]
  xl <- log(lam)
  cols <- wafc_plot_palette(length(nb))
  do.call(graphics::plot,
          wafc_plot_args(list(x = wafc_plot_range(c(xl, if (length(marked)) log(marked))),
                              y = range(0, nrm, t), type = "n",
                              xlab = "log(lambda)",
                              ylab = expression("block norm " * hat(nu)[l * m]),
                              main = main),
                         dots))
  graphics::matlines(xl, nrm, lty = 1, lwd = 1.5, col = cols,
                     type = if (length(lam) == 1L) "p" else "l", pch = 19)
  wafc_plot_marks(marked)
  if (!is.null(t)) {
    graphics::abline(h = t, lty = 2, col = wafc_plot_accent)
  }
  if (length(nb) <= 12L) {
    graphics::legend("topright", legend = nb, col = cols, lty = 1, lwd = 1.5,
                     bg = grDevices::adjustcolor("white", 0.85),
                     box.col = "grey80", cex = 0.7, inset = c(0.01, 0.08))
  }
  list(lambda = lam, norm = nrm, marked = marked, t = t)
}

## The cross-validated loss of every candidate J against log(lambda).
wafc_plot_cv <- function(object, dots) {
  cvl <- object[["cv"]]
  J <- object[["J"]]
  sel <- which(J == object[["J.min"]])
  tab <- do.call(rbind, lapply(seq_along(cvl), function(i) {
    z <- cvl[[i]]
    data.frame(J = J[i], lambda = z[["lambda"]], cvm = z[["cvm"]],
               cvsd = z[["cvsd"]], nzero = as.integer(z[["nzero"]]),
               selected = i == sel)
  }))
  zs <- cvl[[sel]]
  cols <- wafc_plot_palette(length(J))
  ## the vertical axis runs from the lowest point drawn to the largest loss
  ## at the start of a path, the fit with every block at zero: a curve that
  ## overfits at small lambda leaves the panel there instead of flattening
  ## the minima, which are what the panel is for
  top <- max(vapply(cvl, function(z) z[["cvm"]][1L], 0), zs[["cvup"]][
    which(zs[["lambda"]] == object[["lambda.1se"]])[1L]], na.rm = TRUE)
  do.call(graphics::plot,
          wafc_plot_args(list(x = wafc_plot_range(log(tab[["lambda"]])),
                              y = c(min(tab[["cvm"]], zs[["cvlo"]]), top),
                              type = "n", xlab = "log(lambda)",
                              ylab = sprintf("cross-validated %s",
                                             object[["type.measure"]]),
                              main = "Cross-validation over (J, lambda)"),
                         dots))
  xs <- log(zs[["lambda"]])
  graphics::polygon(c(xs, rev(xs)), c(zs[["cvlo"]], rev(zs[["cvup"]])),
                    col = grDevices::adjustcolor(cols[sel], 0.2), border = NA)
  for (i in seq_along(cvl)) {
    z <- cvl[[i]]
    graphics::lines(log(z[["lambda"]]), z[["cvm"]], col = cols[i],
                    lwd = if (i == sel) 2.2 else 1)
    graphics::points(log(z[["lambda.min"]]), z[["cvm.min"]], pch = 19,
                     cex = 0.7, col = cols[i])
  }
  marked <- c(lambda.min = object[["lambda.min"]],
              lambda.1se = object[["lambda.1se"]])
  wafc_plot_marks(marked)
  graphics::legend("topleft", legend = sprintf("J = %d%s", J,
                                               ifelse(seq_along(J) == sel,
                                                      " (selected)", "")),
                   col = cols, lwd = ifelse(seq_along(J) == sel, 2.2, 1),
                   bg = grDevices::adjustcolor("white", 0.85),
                   box.col = "grey80", cex = 0.7, inset = c(0.01, 0.08))
  list(table = tab, J.min = object[["J.min"]],
       lambda.min = object[["lambda.min"]],
       lambda.1se = object[["lambda.1se"]])
}

## The cross-validated error of the threshold against t, a step function on
## the intervals the rule scored, with the band of one standard error.
wafc_plot_thr <- function(th, dots) {
  cand <- th[["candidates"]]
  nrm <- th[["norm"]]
  kept <- th[["kept"]]
  top <- max(c(nrm, cand[["lower"]], th[["t"]])) * 1.05
  if (!(top > 0)) top <- 1
  imin <- which.min(cand[["error"]])
  band <- c(cand[["error"]][imin],
            cand[["error"]][imin] + cand[["se"]][imin])
  xs <- c(cand[["lower"]], top)
  ys <- c(cand[["error"]], cand[["error"]][nrow(cand)])
  do.call(graphics::plot,
          wafc_plot_args(list(x = c(0, top), y = range(ys, band), type = "n",
                              xlab = "threshold t",
                              ylab = "cross-validated mse",
                              main = sprintf("Threshold \"%s\": t = %s",
                                             th[["rule"]],
                                             format(th[["t"]], digits = 3))),
                         dots))
  usr <- graphics::par("usr")
  graphics::rect(usr[1L], band[1L], usr[2L], band[2L],
                 col = grDevices::adjustcolor("grey60", 0.3), border = NA)
  graphics::abline(h = band, lty = 3, col = "grey40")
  graphics::lines(xs, ys, type = "s", lwd = 1.5)
  graphics::abline(v = th[["t"]], lty = 2, col = wafc_plot_accent)
  ## two calls, because rug() takes a single colour
  graphics::rug(as.numeric(nrm)[!as.logical(kept)], ticksize = 0.04, lwd = 2,
                col = "grey60")
  graphics::rug(as.numeric(nrm)[as.logical(kept)], ticksize = 0.04, lwd = 2,
                col = "black")
  list(candidates = cand, t = th[["t"]], rule = th[["rule"]], norm = nrm,
       kept = kept, band = band)
}
