## R/application_report.R -- from the units of the applications to the
## tables and figures.
##
## For each application, under <out>/tables/ and <out>/figures/:
##
##   prediction-<app>.csv  the root mean squared error on the held-out
##                         blocks, mean and standard error over the
##                         partitions; the ratio to the better of the two
##                         splines in each partition; the paired difference
##                         to the spline with the smaller mean error, with
##                         its naive standard error and the one corrected by
##                         Nadeau and Bengio (2003); the fraction of
##                         partitions won against that spline; the R^2 on
##                         the test sample;
##   structure-<app>.csv   for every block (l, m): the frequency with which
##                         each method keeps it over the partitions and
##                         whether the fit on the whole sample keeps it (a
##                         smooth of mgcv is kept when its effective degrees
##                         of freedom exceed 0.1, since select = TRUE shrinks
##                         a smooth towards zero but never to zero exactly);
##                         the median norm of the block in the WAFC; the
##                         range of the component in the fit on the whole
##                         sample; the effective degrees of freedom of the
##                         splines; where the largest increment of the
##                         component is;
##   choices-<app>.csv     J, lambda and the threshold of the WAFC, k of the
##                         splines, split by split, with the k.check() of
##                         the splines;
##   readings-<app>.csv    the quantities the text of the article reads off
##                         the figure (see application_readings below);
##   partitions-<app>.csv  one row per method and split;
##   <app>.pdf, <app>.png  beta_l along each modulator: the fits on the
##                         whole sample and the range of the WAFC over the
##                         partitions, broken over the long gaps in the data
##                         of a modulator, with a rug of the observed values.

#' The units of one application, as a list of rows
#'
#' @return A list with `rows` (one entry per split and label: the split,
#'   the label, the row of R/application.R and the sizes of the split) and
#'   the grid.
collect_application <- function(cfg, out, base) {
  files <- list.files(file.path(out, "units", base), pattern = "^split[0-9]+\\.rds$",
                      recursive = TRUE, full.names = TRUE)
  rows <- list()
  grid <- NULL
  for (f in files) {
    u <- readRDS(f)
    if (!identical(u[["settings"]], application_settings(
      cfg, base, u[["settings"]][["split"]], u[["settings"]][["method"]]))) {
      stop("the unit ", f, " was run with other settings; rerun the fit ",
           "part with --refresh.", call. = FALSE)
    }
    grid <- u[["grid"]]
    for (lab in names(u[["rows"]])) {
      rows[[length(rows) + 1L]] <- list(split = u[["settings"]][["split"]],
                                        label = lab, row = u[["rows"]][[lab]],
                                        sizes = u[["sizes"]])
    }
  }
  list(rows = rows, grid = grid)
}

## The value of one field of the row of a label at a split.
pick_row <- function(col, label, split) {
  for (r in col[["rows"]]) {
    if (r[["label"]] == label && r[["split"]] == split) return(r[["row"]])
  }
  NULL
}

## beta_l along modulator m with the other components at their mean over
## the training sample: c_l + g_lm(u_m) + sum over m' != m of mean g_lm',
## with g_lm the component in the convention of the fit (the centred
## component plus its offset).
partial_curve <- function(row, l, m) {
  if (is.null(row) || is.null(row[["components"]])) return(NULL)
  row[["cc"]][[l]] + row[["components"]][[l, m]] + row[["offset"]][l, m] +
    sum(row[["gmean"]][l, -m])
}

## The gaps in the observed values `u` of a modulator longer than `width`:
## one row per gap, its first and last end.
data_gaps <- function(u, width) {
  s <- sort(unique(u))
  i <- which(diff(s) > width)
  cbind(start = s[i], end = s[i + 1L])
}

## Where the largest increment of a curve on the grid is.
largest_increment <- function(g, v) {
  if (is.null(g) || all(g == g[1L])) return(NA_real_)
  k <- which.max(abs(diff(g)))
  (v[k] + v[k + 1L]) / 2
}

#' The prediction table
application_prediction <- function(col, labels, splines) {
  K <- sort(unique(vapply(col[["rows"]], `[[`, 0, "split")))
  K <- K[K > 0]
  if (length(K) == 0L) return(NULL)
  R <- vapply(labels, function(lab) vapply(K, function(s) {
    r <- pick_row(col, lab, s)
    if (is.null(r) || is.null(r[["rmse"]])) NA_real_ else r[["rmse"]]
  }, 0), numeric(length(K)))
  R <- matrix(R, nrow = length(K), dimnames = list(K, labels))
  sizes <- lapply(K, function(s) {
    for (r in col[["rows"]]) if (r[["split"]] == s) return(r[["sizes"]])
  })
  rho <- mean(vapply(sizes, function(z) z[["ntest"]] / z[["ntrain"]], 0))
  sdt <- mean(vapply(sizes, `[[`, 0, "sd_test"))
  sp <- intersect(splines, labels)
  best <- apply(R[, sp, drop = FALSE], 1L, min)
  ref <- sp[which.min(colMeans(R[, sp, drop = FALSE]))]
  n <- colSums(!is.na(R))
  ## The partitions share most of their training observations, so the
  ## naive standard error of a mean over them is too small; Nadeau and
  ## Bengio (2003) correct the variance by (1/K + n_test/n_train) / (1/K).
  nb <- sqrt((1 / nrow(R) + rho) / (1 / nrow(R)))
  d <- R - R[, ref]
  data.frame(label = labels, partitions = n,
             rmse = colMeans(R, na.rm = TRUE),
             se = apply(R, 2L, stats::sd, na.rm = TRUE) / sqrt(n),
             ratio_best_spline = colMeans(R / best, na.rm = TRUE),
             versus = ref,
             diff = colMeans(d, na.rm = TRUE),
             diff_se = apply(d, 2L, stats::sd, na.rm = TRUE) / sqrt(n),
             diff_se_nb = nb * apply(d, 2L, stats::sd, na.rm = TRUE) / sqrt(n),
             wins = colMeans(R < R[, ref], na.rm = TRUE),
             r2 = 1 - (colMeans(R, na.rm = TRUE) / sdt)^2,
             nb_factor = nb, row.names = NULL, stringsAsFactors = FALSE)
}

#' The structure table
application_structure <- function(col, labels, data) {
  xn <- colnames(data[["x"]])
  un <- colnames(data[["u"]])
  splits <- sort(unique(vapply(col[["rows"]], `[[`, 0, "split")))
  parts <- splits[splits > 0]
  tab <- expand.grid(x = xn, u = un, stringsAsFactors = FALSE)
  tab <- tab[order(match(tab[["x"]], xn), match(tab[["u"]], un)), ]
  rownames(tab) <- NULL
  at <- cbind(match(tab[["x"]], xn), match(tab[["u"]], un))
  grid <- col[["grid"]]
  for (lab in labels) {
    kept <- vapply(seq_len(nrow(at)), function(i) {
      v <- vapply(parts, function(s) {
        r <- pick_row(col, lab, s)
        if (is.null(r[["blocks"]])) NA else r[["blocks"]][at[i, 1L], at[i, 2L]]
      }, NA)
      if (length(v) == 0L) NA_real_ else mean(v, na.rm = TRUE)
    }, 0)
    full <- pick_row(col, lab, 0L)
    tab[[paste0("kept_", lab)]] <- kept
    tab[[paste0("full_", lab)]] <- if (is.null(full[["blocks"]])) NA else
      full[["blocks"]][at]
    tab[[paste0("range_", lab)]] <- vapply(seq_len(nrow(at)), function(i) {
      g <- full[["components"]]
      if (is.null(g)) NA_real_ else diff(range(g[[at[i, 1L], at[i, 2L]]]))
    }, 0)
    if (!is.null(full[["edf"]])) {
      tab[[paste0("edf_", lab)]] <- full[["edf"]][at]
    }
    if (!is.null(full[["norm"]])) {
      tab[[paste0("norm_", lab)]] <- vapply(seq_len(nrow(at)), function(i) {
        stats::median(vapply(parts, function(s) {
          r <- pick_row(col, lab, s)
          if (is.null(r[["norm"]])) NA_real_ else r[["norm"]][at[i, 1L], at[i, 2L]]
        }, 0))
      }, 0)
    }
    tab[[paste0("jump_", lab)]] <- vapply(seq_len(nrow(at)), function(i) {
      g <- full[["components"]]
      if (is.null(g)) NA_real_ else
        largest_increment(g[[at[i, 1L], at[i, 2L]]], grid[, at[i, 2L]])
    }, 0)
  }
  tab
}

#' The choices of every method at every split
application_choices <- function(col) {
  one <- function(r) {
    row <- r[["row"]]
    g <- function(v) if (is.null(row[[v]]) || length(row[[v]]) != 1L)
      NA else row[[v]]
    data.frame(split = r[["split"]], label = r[["label"]],
               rmse = g("rmse"), rmse_train = g("rmse_train"),
               time = g("time"), J = g("J"), top = g("top"),
               lambda = g("lambda"), t = g("t"), t_frac = g("t_frac"),
               nzero = g("nzero"),
               k = if (is.null(row[["k"]])) NA_character_ else
                 paste(row[["k"]], collapse = ","),
               k_top = g("k_top"), kcheck_min_p = g("kcheck_min_p"),
               kcheck_n_low = g("kcheck_n_low"),
               blocks_kept = if (is.null(row[["blocks"]])) NA else
                 sum(row[["blocks"]]),
               ntrain = r[["sizes"]][["ntrain"]],
               ntest = r[["sizes"]][["ntest"]],
               error = if (is.null(row[["error"]])) NA_character_ else
                 row[["error"]],
               stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, lapply(col[["rows"]], one))
  out[order(out[["label"]], out[["split"]]), , drop = FALSE]
}

#' The readings of the article
#'
#' Marylebone Road ("step"): the step of beta_nox along the date. The level
#' before is the mean of the curve on the grid from `margin` days after its
#' start up to `before`; the level after, from `after` to `margin` days
#' before the end of the grid (the ends are left out because there the
#' periodized basis joins the end of the series to its start); the rise is
#' their difference, read as percentage points of the primary fraction. For every date of `at`, the fraction of the rise the
#' curve has reached there; the first date after `from` at which it reaches
#' 90 percent; the peak between `after` minus six months and `after` plus
#' six months. Beside it, the range of the wind component of the slope,
#' absolute and as a fraction of the rise. Each quantity for the fit on the
#' whole sample and as median, minimum and maximum over the partitions.
#'
#' Beijing ("season"): for every component along the day of the year, the
#' mean of the curve inside the heating season minus the mean outside it.
application_readings <- function(col, labels, data, spec) {
  if (is.null(spec)) return(NULL)
  splits <- sort(unique(vapply(col[["rows"]], `[[`, 0, "split")))
  grid <- col[["grid"]]
  summarise <- function(what, lab, f) {
    v <- vapply(splits, function(s) {
      r <- pick_row(col, lab, s)
      if (is.null(r) || !is.null(r[["error"]])) NA_real_ else f(r)
    }, 0)
    pt <- v[splits > 0]
    data.frame(quantity = what, label = lab,
               full = if (any(splits == 0)) v[splits == 0] else NA_real_,
               median = if (length(pt)) stats::median(pt, na.rm = TRUE) else NA,
               min = if (length(pt)) min(pt, na.rm = TRUE) else NA,
               max = if (length(pt)) max(pt, na.rm = TRUE) else NA,
               partitions = sum(!is.na(pt)), stringsAsFactors = FALSE)
  }
  out <- list()
  st <- spec[["step"]]
  if (!is.null(st)) {
    l <- match(st[["block"]][[1L]], colnames(data[["x"]]))
    m <- match(st[["block"]][[2L]], colnames(data[["u"]]))
    o <- data[["def"]][["origin"]][[colnames(data[["u"]])[m]]]
    v <- grid[, m]
    day <- function(s) as.numeric(as.Date(s) - o)
    pre <- v >= min(v) + st[["margin"]] & v < day(st[["before"]])
    post <- v >= day(st[["after"]]) & v <= max(v) - st[["margin"]]
    pk <- v >= day(st[["after"]]) - 183 & v <= day(st[["after"]]) + 183
    lev <- function(r) {
      g <- partial_curve(r, l, m)
      c(pre = mean(g[pre]), post = mean(g[post]))
    }
    wl <- match(st[["wind"]][[1L]], colnames(data[["x"]]))
    wm <- match(st[["wind"]][[2L]], colnames(data[["u"]]))
    for (lab in labels) {
      out[[length(out) + 1L]] <- summarise("level_before", lab,
                                           function(r) lev(r)[["pre"]])
      out[[length(out) + 1L]] <- summarise("level_after", lab,
                                           function(r) lev(r)[["post"]])
      out[[length(out) + 1L]] <- summarise("rise", lab, function(r)
        diff(lev(r)))
      for (a in unlist(st[["at"]])) {
        out[[length(out) + 1L]] <- summarise(
          paste0("fraction_at_", a), lab, function(r) {
            g <- partial_curve(r, l, m)
            lv <- lev(r)
            (stats::approx(v, g, xout = day(a))[["y"]] - lv[["pre"]]) /
              diff(lv)
          })
      }
      out[[length(out) + 1L]] <- summarise("date_90", lab, function(r) {
        g <- partial_curve(r, l, m)
        lv <- lev(r)
        i <- which(v >= day(st[["from"]]) & g >= lv[["pre"]] + 0.9 * diff(lv))[1L]
        if (is.na(i)) NA_real_ else v[i]
      })
      out[[length(out) + 1L]] <- summarise("peak_above_before", lab,
                                           function(r) {
        g <- partial_curve(r, l, m)
        max(g[pk]) - lev(r)[["pre"]]
      })
      out[[length(out) + 1L]] <- summarise("date_peak", lab, function(r) {
        g <- partial_curve(r, l, m)
        v[pk][which.max(g[pk])]
      })
      out[[length(out) + 1L]] <- summarise("wind_range", lab, function(r)
        diff(range(r[["components"]][[wl, wm]])))
      out[[length(out) + 1L]] <- summarise("wind_range_over_rise", lab,
                                           function(r)
        diff(range(r[["components"]][[wl, wm]])) / diff(lev(r)))
    }
  }
  se <- spec[["season"]]
  if (!is.null(se)) {
    m <- match(se[["modulator"]], colnames(data[["u"]]))
    v <- grid[, m]
    on <- unlist(se[["inside"]])
    inside <- if (on[1L] > on[2L]) v >= on[1L] | v < on[2L] else
      v >= on[1L] & v < on[2L]
    for (lab in labels) {
      for (l in seq_len(ncol(data[["x"]]))) {
        out[[length(out) + 1L]] <- summarise(
          paste0("season_contrast_", colnames(data[["x"]])[l]), lab,
          function(r) {
            g <- r[["components"]][[l, m]]
            mean(g[inside]) - mean(g[!inside])
          })
        out[[length(out) + 1L]] <- summarise(
          paste0("largest_increment_", colnames(data[["x"]])[l]), lab,
          function(r) largest_increment(r[["components"]][[l, m]], v))
      }
    }
  }
  do.call(rbind, out)
}

#' The figure of one application
#'
#' One panel per block (l, m): beta_l along u_m with the other components
#' at their mean, for the labels of `figure$labels` fitted on the whole
#' sample, over the range of the label `figure$band` over the partitions
#' (grey). The ends of the grid of a modulator listed in `figure$trim` are
#' not drawn (in the units of the modulator, from the start and from the
#' end): there the periodized basis joins the end of a non-periodic series
#' to its start. Along a modulator listed in `figure$gaps`, the curves and
#' the band are broken over the gaps in its observed values longer than the
#' spacing of the finest level of the WAFC (the label `figure$band`) fitted
#' on the whole sample, 2^-J of the range of the modulator, which is the
#' period of the basis: the wavelets of that level centred in such a gap
#' are set by the penalty and the boundaries, not by the data. Along a
#' modulator listed in `figure$rug`, the panel carries a rug of its observed
#' values, rounded to the resolution given there (in the units of the
#' modulator).
application_figure <- function(col, data, fig, path) {
  xn <- colnames(data[["x"]])
  un <- colnames(data[["u"]])
  def <- data[["def"]]
  grid <- col[["grid"]]
  splits <- sort(unique(vapply(col[["rows"]], `[[`, 0, "split")))
  labs <- unlist(fig[["labels"]])
  lty <- c(1, 2, 4, 3)[seq_along(labs)]
  colr <- unlist(fig[["colours"]])[seq_along(labs)]
  draw <- function() {
    op <- graphics::par(mfrow = c(length(xn), length(un)),
                        mar = c(3.2, 3.6, 1.8, 0.8), mgp = c(2, 0.6, 0),
                        oma = c(0, 0, 1.8, 0), cex = 0.8)
    on.exit(graphics::par(op))
    for (l in seq_along(xn)) {
      for (m in seq_along(un)) {
        v <- grid[, m]
        tr <- unlist(fig[["trim"]][[un[m]]])
        show <- if (is.null(tr)) rep(TRUE, length(v)) else
          v >= min(v) + tr[1L] & v <= max(v) - tr[2L]
        um <- data[["u"]][, m]
        hole <- rep(FALSE, length(v))
        J <- pick_row(col, fig[["band"]], 0L)[["J"]]
        if (un[m] %in% unlist(fig[["gaps"]]) && length(J) > 0L) {
          gp <- data_gaps(um, diff(range(um)) / 2^J[min(m, length(J))])
          for (i in seq_len(nrow(gp))) {
            hole <- hole | (v > gp[i, "start"] & v < gp[i, "end"])
          }
        }
        band <- do.call(cbind, lapply(splits[splits > 0], function(s)
          partial_curve(pick_row(col, fig[["band"]], s), l, m)))
        if (!is.null(band)) band[hole, ] <- NA
        curves <- lapply(labs, function(lab) {
          g <- partial_curve(pick_row(col, lab, 0L), l, m)
          if (!is.null(g)) g[hole] <- NA
          g
        })
        yl <- range(c(unlist(lapply(curves, function(g) g[show])),
                      if (!is.null(band)) band[show, ]), na.rm = TRUE)
        o <- def[["origin"]][[un[m]]]
        xv <- if (is.null(o)) v else o + v
        graphics::plot(xv[show], v[show], type = "n", ylim = yl,
                       xlab = def[["ulab"]][[un[m]]],
                       ylab = def[["betalab"]][[xn[l]]],
                       main = sprintf("%s along %s", def[["xlab"]][[xn[l]]],
                                      tolower(sub(" \\(.*\\)$", "",
                                                  def[["ulab"]][[un[m]]]))))
        if (!is.null(band)) {
          ## one polygon per run of the grid drawn outside the gaps
          run <- rle(show & !hole)
          end <- cumsum(run[["lengths"]])
          for (k in which(run[["values"]])) {
            ii <- (end[k] - run[["lengths"]][k] + 1L):end[k]
            graphics::polygon(c(xv[ii], rev(xv[ii])),
                              c(apply(band[ii, , drop = FALSE], 1L, min),
                                rev(apply(band[ii, , drop = FALSE], 1L, max))),
                              col = "grey85", border = NA)
          }
        }
        for (i in seq_along(labs)) {
          if (!is.null(curves[[i]])) {
            graphics::lines(xv[show], curves[[i]][show], lty = lty[i],
                            col = colr[i], lwd = 1.6)
          }
        }
        mk <- def[["marks"]][[un[m]]]
        if (!is.null(mk)) {
          graphics::abline(v = if (is.null(o)) mk else o + mk, lty = 3,
                           col = "grey40")
        }
        res <- fig[["rug"]][[un[m]]]
        if (!is.null(res)) {
          ru <- unique(round(um / res) * res)
          ru <- ru[ru >= min(v[show]) & ru <= max(v[show])]
          graphics::rug(if (is.null(o)) ru else o + ru, ticksize = 0.025,
                        lwd = 0.4, col = "grey30")
        }
      }
    }
    graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                  mar = c(0, 0, 0, 0), new = TRUE)
    graphics::plot.new()
    graphics::legend("top", legend = c(unlist(fig[["names"]])[labs],
                                       "WAFC over the partitions"),
                     lty = c(lty, NA), col = c(colr, "grey80"), lwd = 1.6,
                     pch = c(rep(NA, length(labs)), 15), pt.cex = 1.8,
                     bty = "n", horiz = TRUE, cex = 0.85, inset = 0.005,
                     xpd = NA)
  }
  w <- 3.6 * length(un)
  h <- 2.8 * length(xn)
  grDevices::pdf(paste0(path, ".pdf"), width = w, height = h)
  draw()
  grDevices::dev.off()
  grDevices::png(paste0(path, ".png"), width = w, height = h, units = "in",
                 res = 150)
  draw()
  grDevices::dev.off()
  invisible(path)
}

#' Write the tables and figures of the applications
application_report <- function(cfg, out, bases = NULL) {
  rp <- cfg[["report"]]
  labels <- unlist(lapply(cfg[["methods"]], function(m)
    application_methods[[m]][["rows"]]))
  tdir <- file.path(out, "tables")
  fdir <- file.path(out, "figures")
  dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  write <- function(d, name) {
    if (!is.null(d)) utils::write.csv(d, file.path(tdir, name),
                                      row.names = FALSE)
  }
  for (b in if (is.null(bases)) cfg[["bases"]] else bases) {
    col <- collect_application(cfg, out, b)
    if (length(col[["rows"]]) == 0L) {
      message("[report] ", b, ": no unit yet")
      next
    }
    data <- application_data(cfg, b)
    pred <- application_prediction(col, labels, unlist(rp[["splines"]]))
    write(pred, sprintf("prediction-%s.csv", b))
    write(application_structure(col, labels, data),
          sprintf("structure-%s.csv", b))
    ch <- application_choices(col)
    write(ch, sprintf("choices-%s.csv", b))
    ## the linear model has constant coefficients: nothing to read along a
    ## modulator
    rd <- application_readings(col, setdiff(labels, "linear"), data,
                               rp[["readings"]][[b]])
    write(rd, sprintf("readings-%s.csv", b))
    if (any(ch[["split"]] == 0L)) {
      application_figure(col, data, rp[["figure"]], file.path(fdir, b))
    }
    message(sprintf("[report] %s: %d row(s) from %d split(s)%s", b, nrow(ch),
                    length(unique(ch[["split"]])),
                    if (is.null(pred)) "" else
                      sprintf("; %d partition(s) in the prediction table",
                              max(pred[["partitions"]]))))
    if (!is.null(pred)) {
      print(format(pred[c("label", "rmse", "se", "ratio_best_spline", "diff",
                          "diff_se_nb", "wins")], digits = 4),
            row.names = FALSE)
    }
  }
  invisible(TRUE)
}
