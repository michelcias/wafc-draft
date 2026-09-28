## Question 31 of docs/ESTADO.md -- why the cross-validation of the sparse
## group LASSO, on the null scenario and the deep grid of decision D34,
## selects fits at the end of the lambda path.
##
##     Rscript wafc/scripts/08-sgl-null.R [parts] [ncores]
##
## Three parts, each a hypothesis test, with the seeds of the investigation
## of 2026-09-28 so that the numbers recorded there are reproduced:
##
##   curve    One replicate, the one wafc/scripts/07-accel.R flagged (part
##            "sgl", cell "null", n = 500, replicate 1). Along the path, at
##            J = 2 and J = 8, the cross-validated error, the error on an
##            independent test sample of 5000 and the training error; and
##            cv.sparsegl() of the package itself on the same design, path
##            and folds, which is an independent implementation of the loop
##            of folds: if it gives the same curve, the loop of wafc/ is not
##            the cause.
##   cross    Leakage or chance. Writing etahat for the held-out prediction
##            and f for the truth, the excess of the cross-validated error
##            of a point of the path over the null point is
##            mean((etahat - f)^2) - 2 mean(eps (etahat - f)) (differences
##            taken against the null point). Without leakage the held-out
##            prediction is independent of the held-out noise and the cross
##            term has mean zero; a cross term that is positive on average
##            would mean the folds see the held-out responses. Thirty null
##            replicates at n = 500, J = 8, for both penalties.
##   cost     What the selection costs. cv.wafc() over the whole grid, the
##            excess of the true error of the selected fit over the truth
##            at lambda.min and at lambda.1se, for both penalties: 20 null
##            replicates at each of n = 250, 500, 1000, and 10 of the
##            inhomogeneous scenario at n = 500 for contrast.
##
## Default: every part, all cores but two. The part 'cost' takes about 20
## minutes on 14 cores and 'cross' about 5.

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

args <- commandArgs(trailingOnly = TRUE)
parts <- if (length(args) >= 1L && nzchar(args[1L])) {
  strsplit(args[1L], ",", fixed = TRUE)[[1L]]
} else c("curve", "cross", "cost")
ncores <- if (length(args) >= 2L && nzchar(args[2L])) as.integer(args[2L]) else {
  max(1L, parallel::detectCores() - 2L)
}
out_dir <- Sys.getenv("WAFC_OUT", ".")

## Decision D31: the basis is evaluated by a table built once.
tab <- WaveBased::wtable(family = "Daublets", filter.size = 8L,
                         prec.wavelet = 30L, check = FALSE)

## Fitted values of a path on a design.
path_eta <- function(Z, cf) {
  sweep(as.matrix(Z %*% cf[-1L, , drop = FALSE]), 2L, cf[1L, ], "+")
}

## ---------------------------------------------------------------------------
## Part "curve"
## ---------------------------------------------------------------------------

if ("curve" %in% parts) {
  ## the seed of wafc/scripts/07-accel.R, part "sgl", cell "null" (the fifth),
  ## n = 500 (the second), replicate 1
  seed <- 20260928L + 200000L + 10000L * 5L + 1000L * 2L + 1L
  d <- simulate_wafc(500, p = 3, q = 2, scenario = "null", seed = seed,
                     sigma = 0.62)
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:10, 500))
  te <- simulate_wafc(5000, p = 3, q = 2, scenario = "null", seed = seed + 999L,
                      sigma = 0.62)
  cat("\n== curve: one null replicate, n = 500 ==\n")
  for (pen in c("sglasso", "lasso")) for (J in c(2L, 8L)) {
    des <- wafc_design(d[["x"]], d[["u"]], J = J, wavelet.table = tab)
    full <- wafc(design = des, y = d[["y"]], penalty = pen)
    z <- wafc_cv_design(des, d[["y"]], full, foldid, function(e) e^2, pen)
    cf <- wafc_raw_coef(full)
    dt <- wafc_design(te[["x"]], te[["u"]], spec = des)
    test_mse <- colMeans((te[["y"]] - path_eta(dt[["Z"]], cf))^2)
    train_mse <- colMeans((d[["y"]] - path_eta(des[["Z"]], cf))^2)
    k <- unique(c(1, seq(10, length(full[["lambda"]]), by = 10),
                  length(full[["lambda"]])))
    cat(sprintf("\n%s J = %d: lambda end / lambda max = %.4g; cv argmin %d, test argmin %d\n",
                pen, J, utils::tail(full[["lambda"]], 1) / full[["lambda"]][1],
                which.min(z[["cvm"]]), which.min(test_mse)))
    print(data.frame(k = k, nzero = full[["nzero"]][k],
                     cvm = round(z[["cvm"]][k], 4), test = round(test_mse[k], 4),
                     train = round(train_mse[k], 4)), row.names = FALSE)
    if (pen == "sglasso") {
      g <- wafc_groups(des)
      ref <- sparsegl::cv.sparsegl(as.matrix(des[["Z"]]), d[["y"]],
                                   group = g[["group"]], family = "gaussian",
                                   lambda = full[["lambda"]], foldid = foldid,
                                   pf_group = g[["pf_group"]],
                                   pf_sparse = g[["pf_sparse"]], asparse = 0.05,
                                   intercept = TRUE, standardize = FALSE,
                                   eps = 1e-9)
      ## the first point is left out: wafc() writes the least squares fit on
      ## the levels there, and sparsegl leaves the zero vector (note in
      ## wafc/R/fit.R)
      cat(sprintf("   cv.sparsegl on the same design, path and folds: argmin %d; largest |difference| of cvm from the second point on: %.3g\n",
                  which.min(ref[["cvm"]]),
                  max(abs(z[["cvm"]][-1L] - ref[["cvm"]][-1L]))))
    }
  }
}

## ---------------------------------------------------------------------------
## Part "cross"
## ---------------------------------------------------------------------------

## Held-out predictions at every level of the path of the fit on the whole
## sample, with the folds of wafc_cv_design().
heldout <- function(des, y, full, foldid, pen) {
  lam <- full[["lambda"]]
  eta <- matrix(NA_real_, length(y), length(lam))
  for (k in seq_len(max(foldid))) {
    o <- foldid == k
    fit <- wafc(design = wafc_subset_design(des, !o), y = y[!o], penalty = pen,
                lambda = lam)
    eta[o, ] <- path_eta(des[["Z"]][o, , drop = FALSE],
                         wafc_raw_coef(fit, s = lam))
  }
  eta
}

run_cross <- function(r) {
  seed <- 20260928L + 700000L + r
  d <- simulate_wafc(500, p = 3, q = 2, scenario = "null", seed = seed,
                     sigma = 0.62)
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:10, 500))
  e <- d[["y"]] - d[["f"]]
  out <- NULL
  for (pen in c("sglasso", "lasso")) {
    des <- wafc_design(d[["x"]], d[["u"]], J = 8L, wavelet.table = tab)
    full <- wafc(design = des, y = d[["y"]], penalty = pen)
    eta <- heldout(des, d[["y"]], full, foldid, pen)
    cvm <- colMeans((d[["y"]] - eta)^2)
    tr <- colMeans((d[["y"]] - path_eta(des[["Z"]], wafc_raw_coef(full)))^2)
    ## the end of the path for the sparse group LASSO, and the LASSO point
    ## with the same training error as that end
    k <- if (pen == "sglasso") length(cvm) else which.min(abs(tr - out[["train"]][1L]))
    g <- eta[, k] - d[["f"]]
    g0 <- eta[, 1L] - d[["f"]]
    out <- rbind(out, data.frame(
      rep = r, pen = pen, k = k, nzero = full[["nzero"]][k], train = tr[k],
      cv_excess = cvm[k] - cvm[1L], mean_g2 = mean(g^2) - mean(g0^2),
      cross = mean(e * g) - mean(e * g0), cv_argmin = which.min(cvm),
      nl = length(cvm)))
  }
  out
}

if ("cross" %in% parts) {
  cat("\n== cross: 30 null replicates, n = 500, J = 8 ==\n")
  res <- do.call(rbind, parallel::mclapply(1:30, run_cross, mc.cores = ncores,
                                           mc.preschedule = FALSE))
  saveRDS(res, file.path(out_dir, "e25-sgl-cross.rds"))
  for (pen in c("sglasso", "lasso")) {
    a <- res[res[["pen"]] == pen, ]
    cat(sprintf("\n%s, point: %s\n", pen,
                if (pen == "sglasso") "end of the path" else
                  "same training error as the sparse group LASSO end"))
    cat(sprintf("  median nzero %d, median training mse %.3f\n",
                as.integer(stats::median(a[["nzero"]])),
                stats::median(a[["train"]])))
    cat(sprintf("  cv excess over the null point: mean %+.4f (sd %.4f), negative in %d of %d\n",
                mean(a[["cv_excess"]]), stats::sd(a[["cv_excess"]]),
                sum(a[["cv_excess"]] < 0), nrow(a)))
    cat(sprintf("  = excess of mean(ghat^2) %+.4f minus twice the cross term %+.4f\n",
                mean(a[["mean_g2"]]), mean(a[["cross"]])))
    cat(sprintf("  cross term: mean %+.4f, standard error %.4f, t = %.1f\n",
                mean(a[["cross"]]), stats::sd(a[["cross"]]) / sqrt(nrow(a)),
                mean(a[["cross"]]) / (stats::sd(a[["cross"]]) / sqrt(nrow(a)))))
    cat(sprintf("  cv argmin at the last point of the path in %d of %d\n",
                sum(a[["cv_argmin"]] == a[["nl"]]), nrow(a)))
  }
}

## ---------------------------------------------------------------------------
## Part "cost"
## ---------------------------------------------------------------------------

run_cost <- function(j) {
  seed <- 20260928L + 800000L + 1000L * match(j[["n"]], c(250L, 500L, 1000L)) +
    100L * (j[["sc"]] == "null") + j[["r"]]
  a <- list(n = j[["n"]], p = 3, q = 2, scenario = j[["sc"]], seed = seed)
  if (j[["sc"]] == "null") a[["sigma"]] <- 0.62 else a[["snr"]] <- 3
  d <- do.call(simulate_wafc, a)
  a[["n"]] <- 5000L
  a[["seed"]] <- seed + 5000L
  te <- do.call(simulate_wafc, a)
  te[["y"]] <- te[["f"]] + stats::rnorm(5000, sd = d[["sigma"]])
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:10, j[["n"]]))
  out <- NULL
  for (pen in c("lasso", "sglasso")) {
    cv <- cv.wafc(d[["x"]], d[["u"]], d[["y"]], penalty = pen, foldid = foldid,
                  wavelet.table = tab)
    fit <- cv[["wafc.fit"]]
    dt <- wafc_design(te[["x"]], te[["u"]], spec = fit[["design"]])
    excess <- function(s) {
      cf <- wafc_raw_coef(fit, s = s)
      mean((as.numeric(dt[["Z"]] %*% cf[-1L, 1L]) + cf[1L, 1L] - te[["f"]])^2)
    }
    z <- cv[["cv"]][[which(cv[["J"]] == cv[["J.min"]])]]
    im <- which(z[["lambda"]] == cv[["lambda.min"]])
    out <- rbind(out, data.frame(
      sc = j[["sc"]], n = j[["n"]], rep = j[["r"]], pen = pen,
      J = cv[["J.min"]], imin = im, nl = length(z[["lambda"]]),
      nzero = fit[["nzero"]][im],
      dip = cv[["cvm.min"]] - cv[["cv"]][[1L]][["cvm"]][1L],
      excess_min = excess(cv[["lambda.min"]]),
      excess_1se = excess(cv[["lambda.1se"]]),
      nzero_1se = fit[["nzero"]][which(z[["lambda"]] == cv[["lambda.1se"]])],
      sigma2 = d[["sigma"]]^2))
  }
  out
}

if ("cost" %in% parts) {
  cat("\n== cost: the selected fit against the truth ==\n")
  jobs <- c(lapply(1:20, function(r) list(sc = "null", n = 250L, r = r)),
            lapply(1:20, function(r) list(sc = "null", n = 500L, r = r)),
            lapply(1:20, function(r) list(sc = "null", n = 1000L, r = r)),
            lapply(1:10, function(r) list(sc = "inhomogeneous", n = 500L, r = r)))
  res <- do.call(rbind, parallel::mclapply(jobs, run_cost, mc.cores = ncores,
                                           mc.preschedule = FALSE))
  saveRDS(res, file.path(out_dir, "e25-sgl-cost.rds"))
  res[["rel_min"]] <- res[["excess_min"]] / res[["sigma2"]]
  res[["rel_1se"]] <- res[["excess_1se"]] / res[["sigma2"]]
  for (sc in unique(res[["sc"]])) for (n in unique(res[["n"]][res[["sc"]] == sc])) {
    cat(sprintf("\n%s, n = %d: excess of the true error over the truth, as a fraction of sigma^2 (median [90th percentile])\n",
                sc, n))
    for (pen in c("lasso", "sglasso")) {
      a <- res[res[["sc"]] == sc & res[["n"]] == n & res[["pen"]] == pen, ]
      cat(sprintf("  %-8s lambda.min %.4f [%.4f]  lambda.1se %.4f [%.4f] | J: %s | lambda.min in the last 20 points: %d of %d | nzero median %d (lambda.1se %d)\n",
                  pen, stats::median(a[["rel_min"]]),
                  stats::quantile(a[["rel_min"]], 0.9),
                  stats::median(a[["rel_1se"]]),
                  stats::quantile(a[["rel_1se"]], 0.9),
                  paste(names(table(a[["J"]])), table(a[["J"]]), sep = ":",
                        collapse = " "),
                  sum(a[["imin"]] > a[["nl"]] - 20), nrow(a),
                  as.integer(stats::median(a[["nzero"]])),
                  as.integer(stats::median(a[["nzero_1se"]]))))
    }
  }
}

cat("\nOK\n")
