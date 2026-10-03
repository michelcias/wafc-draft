## Question 29 of docs/ESTADO.md -- two accelerations of the cross-validation
## of the WAFC, measured before either is adopted.
##
##     Rscript wafc/scripts/07-accel.R [n_rep] [parts] [ncores] [ns] [cells]
##
## Since decision D34 the default grid of J is 2:8, and the cross-validation
## costs 10 to 18 times what it cost on the grid of cv.wall(). Two changes
## would cut that, and each can move the selected pair, so neither enters
## wafc/R/ unless it leaves J.min, lambda.min and lambda.1se unchanged in
## every replicate measured here.
##
##   staged  The fold fits of each candidate J are computed on a prefix of
##           the lambda path, which is extended only when the minimum of the
##           cross-validated error lies near the end of the prefix. glmnet
##           and sparsegl solve the path in sequence, warm starting from the
##           previous level, so the solutions on a prefix are the ones the
##           whole path gives at those levels: the cross-validated error is
##           the same number there, and what can change is only a minimum
##           lying beyond the point where the search stopped. lambda.1se is
##           the first level with cvm <= cvm[imin] + cvsd[imin], which lies
##           before imin, so it moves only if lambda.min moves. The prefix
##           starts at 'm0' levels, grows by 'step', and stops once the
##           minimum is at least 'guard' levels before its end.
##   sgl     The sparse group LASSO at thresh = 1e-8, the default of
##           sparsegl, against the 1e-9 it inherited from the measurement
##           step E2.4b made for glmnet; and the staged path above at 1e-9.
##           The path of sparsegl does not depend on the tolerance, so the
##           two runs choose among the same penalty levels.
##
## Default: 20 replicates, both parts, all cores but two. The part 'sgl'
## caps the replicates at 10, and the cell "mixed" (p q = 16 blocks, up to
## 4084 columns at J = 8) runs 5 replicates in 'staged' and 3 in 'sgl',
## because one cross-validation of it costs minutes. Every arm of a
## replicate sees the same data and the same folds.
##
## The cells are the four of wafc/scripts/04-pilot.R plus the scenario
## "uneven" of decision D30, at the signal to noise ratio of "smooth";
## the seeds are this script's own.

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1L && nzchar(args[1L])) as.integer(args[1L]) else 20L
parts <- if (length(args) >= 2L && nzchar(args[2L])) {
  strsplit(args[2L], ",", fixed = TRUE)[[1L]]
} else c("staged", "sgl")
ncores <- if (length(args) >= 3L && nzchar(args[3L])) as.integer(args[3L]) else {
  max(1L, parallel::detectCores() - 2L)
}
ns <- if (length(args) >= 4L && nzchar(args[4L])) {
  as.integer(strsplit(args[4L], ",", fixed = TRUE)[[1L]])
} else c(250L, 500L, 1000L)
out_dir <- Sys.getenv("WAFC_OUT", ".")
seed0 <- 20260928L

## Decision D31: the basis is evaluated by a table built once.
accel_table <- WaveBased::wtable(family = "Daublets", filter.size = 8L,
                                 prec.wavelet = 30L, check = FALSE)

cells <- list(
  list(name = "smooth", scenario = "smooth", p = 3L, q = 2L, snr = 4),
  list(name = "uneven", scenario = "uneven", p = 3L, q = 2L, snr = 4),
  list(name = "inhomogeneous", scenario = "inhomogeneous", p = 3L, q = 2L,
       snr = 3),
  list(name = "mixed", scenario = "inhomogeneous", p = 4L, q = 4L, snr = 3,
       reps = c(staged = 5L, sgl = 3L)),
  list(name = "null", scenario = "null", p = 3L, q = 2L, sigma = 0.62)
)
if (length(args) >= 5L && nzchar(args[5L])) {
  want <- strsplit(args[5L], ",", fixed = TRUE)[[1L]]
  cells <- cells[vapply(cells, `[[`, "", "name") %in% want]
}

## Parameters of the staged path.
m0 <- 40L
step <- 20L
guard <- 10L

## ---------------------------------------------------------------------------
## The staged cross-validation
## ---------------------------------------------------------------------------

## The cross-validation of one design on a prefix of the path of the fit on
## the whole sample. It goes through wafc_cv_design() itself, handed the fit
## with its path truncated, so the fold fits and the loss are the ones of
## the production code and only the length of the path differs.
cv_design_staged <- function(design, y, full, foldid, penalty, ...) {
  lam <- full[["lambda"]]
  nl <- length(lam)
  m <- min(nl, m0)
  repeat {
    part <- full
    part[["lambda"]] <- lam[seq_len(m)]
    part[["nzero"]] <- full[["nzero"]][seq_len(m)]
    z <- wafc_cv_design(design, y, part, foldid, function(e) e^2, penalty,
                        ...)
    imin <- which.min(z[["cvm"]])
    if (m == nl || imin <= m - guard) break
    m <- min(nl, m + step)
  }
  z[["m"]] <- m
  z[["nl"]] <- nl
  z
}

## The selection of cv.wafc() over the grid of J, with the staged path in
## place of the whole one: one design and one fit on the whole sample per
## candidate, ties to the smallest J.
cv_wafc_staged <- function(x, u, y, foldid, penalty, thresh) {
  J <- wafc_J_grid(NULL, nrow(x))
  zs <- vector("list", length(J))
  for (i in seq_along(J)) {
    des <- wafc_design(x, u, J = J[i], wavelet.table = accel_table)
    full <- wafc(design = des, y = y, penalty = penalty, thresh = thresh)
    zs[[i]] <- cv_design_staged(des, y, full, foldid, penalty, thresh = thresh)
  }
  best <- which.min(vapply(zs, `[[`, 0, "cvm.min"))
  list(J = J, cv = zs, J.min = J[best], lambda.min = zs[[best]][["lambda.min"]],
       lambda.1se = zs[[best]][["lambda.1se"]])
}

## ---------------------------------------------------------------------------
## Replicates
## ---------------------------------------------------------------------------

draw <- function(cell, n, r, part) {
  seed <- seed0 + 100000L * match(part, c("staged", "sgl")) +
    10000L * match(cell[["name"]], vapply(cells, `[[`, "", "name")) +
    1000L * match(n, ns) + r
  a <- list(n = n, p = cell[["p"]], q = cell[["q"]],
            scenario = cell[["scenario"]], seed = seed)
  if (!is.null(cell[["sigma"]])) a[["sigma"]] <- cell[["sigma"]]
  else a[["snr"]] <- cell[["snr"]]
  dgp <- do.call(simulate_wafc, a)
  set.seed(seed + 77L)
  list(dgp = dgp, foldid = sample(rep_len(1:10, n)))
}

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  v <- expr
  list(value = v, secs = proc.time()[["elapsed"]] - t0)
}

## One comparison of the whole path (cv.wafc, the production code) with the
## staged one, at one penalty and one tolerance. Returns one row per
## candidate J and one summary row.
compare_staged <- function(cell, n, r, dat, penalty, thresh, part) {
  x <- dat[["dgp"]][["x"]]
  u <- dat[["dgp"]][["u"]]
  y <- dat[["dgp"]][["y"]]
  full <- timed(cv.wafc(x, u, y, penalty = penalty, foldid = dat[["foldid"]],
                        wavelet.table = accel_table, thresh = thresh,
                        threshold = "none"))
  stg <- timed(cv_wafc_staged(x, u, y, dat[["foldid"]], penalty, thresh))
  a <- full[["value"]]
  b <- stg[["value"]]
  byJ <- do.call(rbind, lapply(seq_along(a[["J"]]), function(i) {
    za <- a[["cv"]][[i]]
    zb <- b[["cv"]][[i]]
    m <- zb[["m"]]
    data.frame(part = part, cell = cell[["name"]], n = n, rep = r,
               penalty = penalty, thresh = thresh, J = a[["J"]][i],
               nl = length(za[["lambda"]]), m = m,
               imin_full = which.min(za[["cvm"]]),
               imin_staged = which.min(zb[["cvm"]]),
               same_path = isTRUE(all.equal(za[["lambda"]][seq_len(m)],
                                            zb[["lambda"]], tolerance = 0)),
               cvm_maxdiff = max(abs(za[["cvm"]][seq_len(m)] - zb[["cvm"]])),
               stringsAsFactors = FALSE)
  }))
  summ <- data.frame(part = part, cell = cell[["name"]], n = n, rep = r,
                     penalty = penalty, thresh = thresh, arm = "staged",
                     J_full = a[["J.min"]], J_arm = b[["J.min"]],
                     same_J = a[["J.min"]] == b[["J.min"]],
                     same_min = a[["lambda.min"]] == b[["lambda.min"]],
                     same_1se = a[["lambda.1se"]] == b[["lambda.1se"]],
                     secs_full = full[["secs"]], secs_arm = stg[["secs"]],
                     stringsAsFactors = FALSE)
  list(byJ = byJ, summ = summ, full = a)
}

run_staged <- function(cell, n, r) {
  dat <- draw(cell, n, r, "staged")
  z <- compare_staged(cell, n, r, dat, "lasso", 1e-9, "staged")
  list(byJ = z[["byJ"]], summ = z[["summ"]])
}

run_sgl <- function(cell, n, r) {
  dat <- draw(cell, n, r, "sgl")
  x <- dat[["dgp"]][["x"]]
  u <- dat[["dgp"]][["u"]]
  y <- dat[["dgp"]][["y"]]
  ## the staged path at the current tolerance, against the whole one
  z <- compare_staged(cell, n, r, dat, "sglasso", 1e-9, "sgl")
  a <- z[["full"]]
  ## the looser tolerance, on the whole path
  loose <- timed(cv.wafc(x, u, y, penalty = "sglasso",
                         foldid = dat[["foldid"]],
                         wavelet.table = accel_table, thresh = 1e-8,
                         threshold = "none"))
  b <- loose[["value"]]
  same_paths <- all(vapply(seq_along(a[["J"]]), function(i) {
    isTRUE(all.equal(a[["cv"]][[i]][["lambda"]], b[["cv"]][[i]][["lambda"]],
                     tolerance = 0))
  }, TRUE))
  kkt_bad <- function(fit, lam) {
    k <- wafc_kkt(fit, s = lam)
    sum(!k[["ok"]])
  }
  summ2 <- data.frame(part = "sgl", cell = cell[["name"]], n = n, rep = r,
                      penalty = "sglasso", thresh = 1e-8, arm = "thresh",
                      J_full = a[["J.min"]], J_arm = b[["J.min"]],
                      same_J = a[["J.min"]] == b[["J.min"]],
                      same_min = a[["lambda.min"]] == b[["lambda.min"]],
                      same_1se = a[["lambda.1se"]] == b[["lambda.1se"]],
                      secs_full = z[["summ"]][["secs_full"]],
                      secs_arm = loose[["secs"]],
                      stringsAsFactors = FALSE)
  summ2[["same_paths"]] <- same_paths
  summ2[["kkt_bad_full"]] <- kkt_bad(a[["wafc.fit"]], a[["lambda.min"]])
  summ2[["kkt_bad_arm"]] <- kkt_bad(b[["wafc.fit"]], b[["lambda.min"]])
  s1 <- z[["summ"]]
  s1[["same_paths"]] <- NA
  s1[["kkt_bad_full"]] <- NA_integer_
  s1[["kkt_bad_arm"]] <- NA_integer_
  list(byJ = z[["byJ"]], summ = rbind(s1, summ2))
}

## ---------------------------------------------------------------------------
## Driver
## ---------------------------------------------------------------------------

run_part <- function(fun, part, reps) {
  jobs <- list()
  for (cell in cells) {
    nrep <- if (is.null(cell[["reps"]])) reps else min(reps, cell[["reps"]][[part]])
    for (n in ns) for (r in seq_len(nrep)) {
      jobs[[length(jobs) + 1L]] <- list(cell = cell, n = n, r = r)
    }
  }
  cat(sprintf("\n== %s: %d jobs on %d core(s) ==\n", part, length(jobs), ncores))
  t0 <- proc.time()[["elapsed"]]
  res <- parallel::mclapply(jobs, function(j) {
    tryCatch(fun(j[["cell"]], j[["n"]], j[["r"]]),
             error = function(e) {
               structure(list(cell = j[["cell"]][["name"]], n = j[["n"]],
                              r = j[["r"]], msg = conditionMessage(e)),
                         class = "accel_failure")
             })
  }, mc.cores = ncores, mc.preschedule = FALSE)
  cat(sprintf("   %.1f s\n", proc.time()[["elapsed"]] - t0))
  ok <- vapply(res, function(v) is.list(v) && !is.null(v[["summ"]]), TRUE)
  if (any(!ok)) {
    cat(sprintf("   %d of %d job(s) failed:\n", sum(!ok), length(jobs)))
    for (v in res[!ok]) {
      cat("     ", if (inherits(v, "accel_failure")) {
        sprintf("%s n = %d rep %d: %s", v[["cell"]], v[["n"]], v[["r"]],
                substr(v[["msg"]], 1L, 100L))
      } else "worker returned no result", "\n")
    }
  }
  list(byJ = do.call(rbind, lapply(res[ok], `[[`, "byJ")),
       summ = do.call(rbind, lapply(res[ok], `[[`, "summ")),
       failed = sum(!ok))
}

report <- function(z, label) {
  s <- z[["summ"]]
  b <- z[["byJ"]]
  cat(sprintf("\n-- %s --\n", label))
  cat(sprintf("prefix equal to the whole path in every candidate: %s; largest |cvm difference| on the prefix: %.3g\n",
              all(b[["same_path"]]), max(b[["cvm_maxdiff"]])))
  agg <- stats::aggregate(
    cbind(same_J, same_min, same_1se, secs_full, secs_arm) ~ arm + cell + n,
    data = s, FUN = function(v) c(sum = sum(v), med = stats::median(v)))
  for (i in seq_len(nrow(agg))) {
    k <- sum(s[["arm"]] == agg[["arm"]][i] & s[["cell"]] == agg[["cell"]][i] &
               s[["n"]] == agg[["n"]][i])
    cat(sprintf("%-7s %-14s n = %4d: J %2d/%d, lambda.min %2d/%d, lambda.1se %2d/%d; median %.1f s -> %.1f s\n",
                agg[["arm"]][i], agg[["cell"]][i], agg[["n"]][i],
                as.integer(agg[["same_J"]][i, "sum"]), k,
                as.integer(agg[["same_min"]][i, "sum"]), k,
                as.integer(agg[["same_1se"]][i, "sum"]), k,
                agg[["secs_full"]][i, "med"], agg[["secs_arm"]][i, "med"]))
  }
  for (a in unique(s[["arm"]])) {
    sa <- s[s[["arm"]] == a, ]
    cat(sprintf("%s: identical selection in %d of %d replicates; total %.0f s -> %.0f s (%.1fx)\n",
                a, sum(sa[["same_J"]] & sa[["same_min"]] & sa[["same_1se"]]),
                nrow(sa), sum(sa[["secs_full"]]), sum(sa[["secs_arm"]]),
                sum(sa[["secs_full"]]) / sum(sa[["secs_arm"]])))
  }
  if (!is.null(s[["kkt_bad_arm"]]) && any(!is.na(s[["kkt_bad_arm"]]))) {
    cat(sprintf("KKT at lambda.min, failures: thresh 1e-9 %d, thresh 1e-8 %d\n",
                sum(s[["kkt_bad_full"]], na.rm = TRUE),
                sum(s[["kkt_bad_arm"]], na.rm = TRUE)))
    cat(sprintf("same lambda paths at the two tolerances: %s\n",
                all(s[["same_paths"]], na.rm = TRUE)))
  }
  cat("where the staged path stopped (median m of nl, and imin of the whole path), by J:\n")
  st <- stats::aggregate(cbind(m, nl, imin_full) ~ J, data = b,
                         FUN = stats::median)
  print(st, row.names = FALSE)
}

if ("staged" %in% parts) {
  z <- run_part(run_staged, "staged", R)
  saveRDS(z, file.path(out_dir, "e25-accel-staged.rds"))
  report(z, "staged path, LASSO, thresh 1e-9")
}
if ("sgl" %in% parts) {
  z <- run_part(run_sgl, "sgl", min(R, 10L))
  saveRDS(z, file.path(out_dir, "e25-accel-sgl.rds"))
  report(z, "sparse group LASSO: staged path at 1e-9, and thresh 1e-8")
}

cat("\nOK\n")
