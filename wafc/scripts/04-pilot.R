## E2.4 -- the pilot: the WAFC against the competitors of proposal L2d, the
## six penalty rules against the oracle of the grid, and the margin eps
## against the two things it has to arbitrate.
##
##     Rscript wafc/scripts/04-pilot.R [n_rep] [parts] [ncores] [ns] [cells]
##
## 'parts' is a comma separated subset of
##
##   competitors  the table of plano-projeto.md E2.4: scenarios x n x
##                methods x replicates, with the integrated squared error of
##                each component, the out-of-sample prediction error, the
##                structure recovered block by block, and the seconds;
##   lambda       the five rules of E2.3 plus the QUT of Giacobino et al.
##                (2017), against the oracle of the grid, which is the
##                column that gives the word "cost" a meaning;
##   margin       lambda_min(G_eps) and the error against eps, on the grid
##                of open question 6 of docs/ESTADO.md plus the family
##                eps = a (L-1) 2^(-J); it is what fixes wafc_eps_periodic,
##                provisional at 0.05;
##   j1           whether the grid of cv.wafc() should start at J = 1, now
##                that the margin no longer depends on J (step E2.1b).
##
## 'ns' and 'cells' are comma separated and restrict the sweep, which is
## how a part is rerun on one cell without rerunning the rest.
##
## Default: 50 replicates, every part, as many cores as the machine has
## minus two. The two parts that answer a question about the code rather
## than about the method, 'margin' and 'j1', are capped at 20 replicates:
## what they measure is a ranking of margins and a frequency of selection,
## and 20 paired replicates settle both well inside the resolution of the
## decision they feed. The parts that resample write one RDS each into
## WAFC_OUT (the working directory by default), so a part can be rerun
## alone and the tables recomputed from the saved draws.
##
## Every method of a replicate sees the same data, the same folds and the
## same test sample: the comparison is paired, and a difference of medians
## is read across replicates and not across fits.

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1L && nzchar(args[1L])) as.integer(args[1L]) else 50L
parts <- if (length(args) >= 2L && nzchar(args[2L])) {
  strsplit(args[2L], ",", fixed = TRUE)[[1L]]
} else c("competitors", "lambda", "margin", "j1")
if (identical(parts, "all")) parts <- c("competitors", "lambda", "margin", "j1")
ncores <- if (length(args) >= 3L && nzchar(args[3L])) as.integer(args[3L]) else {
  max(1L, parallel::detectCores() - 2L)
}
out_dir <- Sys.getenv("WAFC_OUT", ".")
seed0 <- 20260919L
n_test <- 2000L
n_grid <- 256L
ns <- if (length(args) >= 4L && nzchar(args[4L])) {
  as.integer(strsplit(args[4L], ",", fixed = TRUE)[[1L]])
} else c(250L, 500L, 1000L)
R_small <- min(R, 20L)

## Scenarios of dgp.R, with the effective regularity s' each one declares.
## Decision D27 fixed the two values: the cubic of "smooth" has a corner in
## its periodic extension and the theory is stated at eps = 0 (D26), so the
## smooth scenario reads 3/2 and not the 4 it would read with a margin and
## an extension; blocks and heavisine jump inside the interval, where no
## basis helps, so the inhomogeneous one reads 1/2.
##
## The number belongs in dgp.R as an attribute of the scenario, which is
## what D27 asks for and what step E2.4 may not edit (docs/TAREFA.md,
## section 3); the line to add is in the handoff. Until then it lives here,
## and here only.
wafc_sprime <- c(smooth = 3/2, inhomogeneous = 1/2, null = 3/2)

## The four cells of the pilot. The third is the "mixture with half the
## components null" of plano-projeto.md E2.4, as far as dgp.R reaches:
## wafc_scenario() activates three blocks whatever p and q are, so p = 4 and
## q = 4 give 3 active blocks of 16, which is sparser than a half and not
## looser. Making the fraction a parameter is a change to dgp.R and is in
## the handoff.
cells <- list(
  list(name = "smooth", scenario = "smooth", p = 3L, q = 2L, snr = 4),
  list(name = "inhomogeneous", scenario = "inhomogeneous", p = 3L, q = 2L,
       snr = 3),
  list(name = "mixed", scenario = "inhomogeneous", p = 4L, q = 4L, snr = 3),
  list(name = "null", scenario = "null", p = 3L, q = 2L, sigma = 0.62)
)

if (length(args) >= 5L && nzchar(args[5L])) {
  want <- strsplit(args[5L], ",", fixed = TRUE)[[1L]]
  cells <- cells[vapply(cells, `[[`, "", "name") %in% want]
}

methods <- c("wafc.lasso", "wafc.sglasso", "gam", "bsgl", "klopp", "aspline",
             "vcbart", "linear", "oracle")

## ---------------------------------------------------------------------------
## Shared machinery
## ---------------------------------------------------------------------------

draw_cell <- function(cell, n, seed) {
  a <- list(n = n, p = cell[["p"]], q = cell[["q"]],
            scenario = cell[["scenario"]], seed = seed)
  if (!is.null(cell[["sigma"]])) a[["sigma"]] <- cell[["sigma"]]
  else a[["snr"]] <- cell[["snr"]]
  do.call(simulate_wafc, a)
}

## The test sample of a replicate, with the noise drawn at the error scale
## of the training sample and not at the one the test sample would set
## itself through the signal to noise ratio: the two differ by a little,
## and the difference would end up inside the out-of-sample error.
test_for <- function(cell, dgp, seed) {
  test <- draw_cell(cell, n_test, seed + 500000L)
  test[["y"]] <- test[["f"]] + stats::rnorm(n_test, sd = dgp[["sigma"]])
  test[["sigma"]] <- dgp[["sigma"]]
  test
}

## The grid on which the components are compared: the observed range of
## each modulating covariate, which is where every method is identified and
## where the WAFC integrates to zero (the caveat of reconstruct.R).
grid_of <- function(dgp, n_grid) {
  q <- dgp[["q"]]
  g <- matrix(0, n_grid, q)
  for (m in seq_len(q)) {
    rg <- range(dgp[["u"]][, m])
    g[, m] <- seq(rg[1L], rg[2L], length.out = n_grid)
  }
  g
}

## Integrated squared error of each component on the grid, with both the
## estimate and the truth centred there: every method fixes the level by its
## own convention, so only the shape is comparable (decision D22 and the
## caveat of wafc/R/reconstruct.R).
ise_components <- function(ghat, dgp, grid) {
  p <- dgp[["p"]]
  q <- dgp[["q"]]
  ise <- matrix(NA_real_, p, q)
  if (is.null(ghat)) return(ise)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      gm <- dgp[["g"]][[l, m]]
      gt <- if (is.null(gm)) numeric(nrow(grid)) else gm(grid[, m])
      gh <- ghat[[l, m]] - mean(ghat[[l, m]])
      gt <- gt - mean(gt)
      ise[l, m] <- mean((gh - gt)^2) * diff(range(grid[, m]))
    }
  }
  ise
}

## The metric every method shares, including the one that does not
## decompose: the mean squared error of the functional coefficients at the
## test sample, summed over the p coefficients. Levels are compared here,
## unlike in the component metric, because beta_l(u) is identified.
mse_beta <- function(bhat, test) {
  mean(rowSums((bhat - test[["beta"]])^2))
}

one_row <- function(cell, n, r, method, dgp, test, grid, fit, active,
                    blocks, ghat, extra = list()) {
  bhat <- fit[["beta_test"]]
  ise <- ise_components(ghat, dgp, grid)
  sel <- blocks
  data.frame(
    cell = cell[["name"]], scenario = cell[["scenario"]], n = n, rep = r,
    method = method,
    rmse_f = sqrt(mean((fit[["f_test"]] - test[["f"]])^2)),
    rmse_y = sqrt(mean((fit[["f_test"]] - test[["y"]])^2)),
    mse_beta = mse_beta(bhat, test),
    ise = sum(ise), ise_active = sum(ise[active]),
    ise_null = sum(ise[!active]),
    n_true = if (is.null(sel)) NA_integer_ else sum(sel[active]),
    n_false = if (is.null(sel)) NA_integer_ else sum(sel[!active]),
    n_active = sum(active), n_block = length(active),
    J = if (is.null(extra[["J"]])) NA_integer_ else extra[["J"]],
    lambda = if (is.null(extra[["lambda"]])) NA_real_ else extra[["lambda"]],
    nzero = if (is.null(extra[["nzero"]])) NA_integer_ else extra[["nzero"]],
    sigma = dgp[["sigma"]], time = fit[["time"]],
    stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## Part "competitors"
## ---------------------------------------------------------------------------

run_competitors <- function(cell, n, r) {
  seed <- seed0 + 100000L * match(cell[["name"]], vapply(cells, `[[`, "", "name")) +
    1000L * match(n, ns) + r
  dgp <- draw_cell(cell, n, seed)
  test <- test_for(cell, dgp, seed)
  p <- cell[["p"]]
  grid <- grid_of(dgp, n_grid)
  active <- nzchar(dgp[["structure"]])
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:10, n))
  rows <- list()

  ## The two WAFC variants, tuned by cross-validation over (J, lambda),
  ## which decision D20 made the default rule.
  for (pen in c("lasso", "sglasso")) {
    t0 <- proc.time()[["elapsed"]]
    cv <- try(cv.wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], penalty = pen,
                      foldid = foldid), silent = TRUE)
    if (inherits(cv, "try-error")) next
    el <- proc.time()[["elapsed"]] - t0
    f <- cv[["wafc.fit"]]
    lam <- cv[["lambda.min"]]
    d_test <- wafc_design(test[["x"]], test[["u"]], spec = f[["design"]])
    d_grid <- wafc_design(matrix(1, n_grid, p), grid, spec = f[["design"]])
    cf <- wafc_raw_coef(f, s = lam)
    fh <- as.numeric(d_test[["Z"]] %*% cf[-1L, 1L]) + cf[1L, 1L]
    bh <- predict(f, newu = test[["u"]], s = lam, type = "beta")
    gh <- wafc_grid_components(f, grid, s = lam, design = d_grid)
    blk <- wafc_blocks(f, s = lam)[["nonzero"]] > 0L
    rows[[length(rows) + 1L]] <- one_row(
      cell, n, r, paste0("wafc.", pen), dgp, test, grid,
      list(f_test = fh, beta_test = bh, time = el), active, blk, gh,
      list(J = cv[["J.min"]], lambda = lam,
           nzero = sum(cf[-1L, 1L][-f[["design"]][["unpenalized"]]] != 0)))
  }

  for (mth in c("gam", "bsgl", "klopp", "aspline", "vcbart", "linear",
                "oracle")) {
    f <- try(wafc_competitor(mth, dgp[["x"]], dgp[["u"]], dgp[["y"]],
                             active = active, foldid = foldid), silent = TRUE)
    if (inherits(f, "try-error")) next
    gh <- wafc_grid_components(f, grid)
    rows[[length(rows) + 1L]] <- one_row(
      cell, n, r, mth, dgp, test, grid,
      list(f_test = predict(f, test[["x"]], test[["u"]]),
           beta_test = f[["beta"]](test[["u"]]), time = f[["time"]]),
      active, f[["blocks"]], gh,
      list(J = f[["extra"]][["J"]], lambda = f[["extra"]][["lambda"]]))
  }
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------------
## Part "lambda"
## ---------------------------------------------------------------------------

## The six rules and the oracle of the grid, all on the LASSO and on the
## same replicate. The oracle is the pair (J, lambda) of the whole grid that
## minimises the realised out-of-sample error, which no rule can see: the
## cost of a rule is its error divided by that one.
run_lambda <- function(cell, n, r) {
  seed <- seed0 + 200000L +
    10000L * match(cell[["name"]], vapply(cells, `[[`, "", "name")) +
    1000L * match(n, ns) + r
  dgp <- draw_cell(cell, n, seed)
  test <- test_for(cell, dgp, seed)
  p <- cell[["p"]]
  grid <- grid_of(dgp, n_grid)
  active <- nzchar(dgp[["structure"]])
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:10, n))
  sp <- wafc_sprime[[cell[["scenario"]]]]
  rows <- list()
  dcache <- list()

  designs <- function(fit) {
    key <- paste(fit[["design"]][["J"]], collapse = "-")
    if (is.null(dcache[[key]])) {
      dcache[[key]] <<- list(
        test = wafc_design(test[["x"]], test[["u"]], spec = fit[["design"]]),
        grid = wafc_design(matrix(1, n_grid, p), grid, spec = fit[["design"]]))
    }
    dcache[[key]]
  }
  record <- function(rule, fit, lam, el, J) {
    dd <- designs(fit)
    cf <- wafc_raw_coef(fit, s = lam)
    fh <- as.numeric(dd[["test"]][["Z"]] %*% cf[-1L, 1L]) + cf[1L, 1L]
    bh <- predict(fit, newu = test[["u"]], s = lam, type = "beta")
    gh <- wafc_grid_components(fit, grid, s = lam, design = dd[["grid"]])
    blk <- wafc_blocks(fit, s = lam)[["nonzero"]] > 0L
    one_row(cell, n, r, rule, dgp, test, grid,
            list(f_test = fh, beta_test = bh, time = el), active, blk, gh,
            list(J = J, lambda = lam,
                 nzero = sum(cf[-1L, 1L][-fit[["design"]][["unpenalized"]]] != 0)))
  }

  for (rule in c("cv.min", "cv.1se", "bic", "ebic", "theory")) {
    t0 <- proc.time()[["elapsed"]]
    tn <- try(wafc_tune(dgp[["x"]], dgp[["u"]], dgp[["y"]], rule = rule,
                        foldid = foldid, s = sp, sigma = dgp[["sigma"]]),
              silent = TRUE)
    if (inherits(tn, "try-error")) next
    el <- proc.time()[["elapsed"]] - t0
    rows[[length(rows) + 1L]] <- record(rule, tn[["fit"]], tn[["lambda"]], el,
                                        tn[["J"]])
  }

  ## QUT: the resolution is not part of the rule, so it is taken from the
  ## same cross-validation the other rules use, and only lambda changes.
  t0 <- proc.time()[["elapsed"]]
  cvq <- try(cv.wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], foldid = foldid),
             silent = TRUE)
  if (!inherits(cvq, "try-error")) {
    fq <- cvq[["wafc.fit"]]
    lq <- wafc_lambda_qut(fq[["design"]], dgp[["y"]], nsim = 200L,
                          seed = seed + 5L)
    fq2 <- wafc(design = fq[["design"]], y = dgp[["y"]],
                lambda = wafc_path_to(lq, fq[["design"]], dgp[["y"]]))
    el <- proc.time()[["elapsed"]] - t0
    rows[[length(rows) + 1L]] <- record("qut", fq2, lq, el, cvq[["J.min"]])
  }

  ## The oracle of the grid.
  t0 <- proc.time()[["elapsed"]]
  best <- NULL
  for (Ji in wafc_J_grid(NULL, n)) {
    fj <- wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], J = Ji)
    dd <- designs(fj)
    cfj <- wafc_raw_coef(fj)
    fh <- sweep(as.matrix(dd[["test"]][["Z"]] %*% cfj[-1L, , drop = FALSE]),
                2L, cfj[1L, ], "+")
    rm_j <- sqrt(colMeans((fh - test[["f"]])^2))
    k <- which.min(rm_j)
    if (is.null(best) || rm_j[k] < best[["rmse"]]) {
      best <- list(rmse = rm_j[k], J = Ji, lambda = fj[["lambda"]][k], fit = fj)
    }
  }
  el <- proc.time()[["elapsed"]] - t0
  rows[[length(rows) + 1L]] <- record("oracle.grid", best[["fit"]],
                                      best[["lambda"]], el, best[["J"]])
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------------
## Part "margin"
## ---------------------------------------------------------------------------

## lambda_min of the Gram of the basis restricted to the support, the
## constant that replaces c_U in the design condition of E1.4 once the
## margin is positive (decision D23). Deterministic: it depends on eps, J
## and the filter, not on the data. The construction is the one of part (7)
## of derivations/check/02-aproximacao-besov.R, rewritten here because a
## check script is not a library.
gram_restricted <- function(J, eps, filter.size = 8L, ngrid = 2^14) {
  uG <- (seq_len(ngrid) - 0.5) / ngrid
  W <- WaveBased::wbasis(uG, j0 = 0L, J = J, family = "Daublets",
                         filter.size = filter.size)
  idx <- which(uG >= eps & uG <= 1 - eps)
  ev <- sort(eigen(crossprod(W[idx, , drop = FALSE]) / ngrid,
                   symmetric = TRUE, only.values = TRUE)[["values"]])
  c(lmin = ev[1L], nzero = sum(ev < 1e-8))
}

## The grid of margins of open question 6: the three of decision D23, the
## two fixed values, and the family a (L-1) 2^(-J), which is the width of
## the strip the periodization contaminates. 'a' around 1 is the crossing
## that D26 showed to be exact: above it the margin buys the extension, below
## it the restricted Gram survives, and the two demands cannot both hold.
margin_grid <- function(J, L = 8L) {
  c("0" = 0, "2^-(J+1)" = 2^(-J - 1), "1.9^-J" = 1.9^(-J),
    "0.02" = 0.02, "0.05" = 0.05, "0.10" = 0.10,
    "a=0.5" = 0.5 * (L - 1) * 2^(-J), "a=0.75" = 0.75 * (L - 1) * 2^(-J),
    "a=1" = (L - 1) * 2^(-J), "a=1.5" = 1.5 * (L - 1) * 2^(-J),
    "a=2" = 2 * (L - 1) * 2^(-J))
}

run_margin_fit <- function(cell, n, r) {
  seed <- seed0 + 300000L +
    10000L * match(cell[["name"]], vapply(cells, `[[`, "", "name")) +
    1000L * match(n, ns) + r
  dgp <- draw_cell(cell, n, seed)
  test <- test_for(cell, dgp, seed)
  p <- cell[["p"]]
  grid <- grid_of(dgp, n_grid)
  active <- nzchar(dgp[["structure"]])
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:10, n))
  rows <- list()
  for (J in 2:5) {
    eg <- margin_grid(J)
    for (en in names(eg)) {
      e <- eg[[en]]
      if (e >= 0.5) next
      fit <- try(wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], J = J, eps = e),
                 silent = TRUE)
      if (inherits(fit, "try-error")) next
      z <- wafc_cv_design(fit[["design"]], dgp[["y"]], fit, foldid,
                          function(v) v^2, "lasso")
      lam <- z[["lambda.min"]]
      d_test <- wafc_design(test[["x"]], test[["u"]], spec = fit[["design"]])
      d_grid <- wafc_design(matrix(1, n_grid, p), grid, spec = fit[["design"]])
      cf <- wafc_raw_coef(fit, s = lam)
      fh <- as.numeric(d_test[["Z"]] %*% cf[-1L, 1L]) + cf[1L, 1L]
      gh <- wafc_grid_components(fit, grid, s = lam, design = d_grid)
      ise <- ise_components(gh, dgp, grid)
      rows[[length(rows) + 1L]] <- data.frame(
        cell = cell[["name"]], n = n, rep = r, J = J, eps.rule = en, eps = e,
        lambda = lam, cvm = z[["cvm.min"]],
        rmse_f = sqrt(mean((fh - test[["f"]])^2)),
        ise = sum(ise), ise_active = sum(ise[active]),
        nzero = sum(cf[-1L, 1L][-fit[["design"]][["unpenalized"]]] != 0),
        nullcol = sum(Matrix::colSums(abs(fit[["design"]][["Z"]])) == 0),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------------
## Part "j1"
## ---------------------------------------------------------------------------

## Whether the grid of cv.wafc() should start at J = 1. The margin no
## longer depends on J (step E2.1b), so the design exists there; the
## question is whether the candidate is ever selected and what it costs to
## carry it. Measured by cross-validating twice on the same folds, once on
## 1:Jmax and once on 2:Jmax.
run_j1 <- function(cell, n, r) {
  seed <- seed0 + 400000L +
    10000L * match(cell[["name"]], vapply(cells, `[[`, "", "name")) +
    1000L * match(n, ns) + r
  dgp <- draw_cell(cell, n, seed)
  test <- test_for(cell, dgp, seed)
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:10, n))
  Jmax <- max(wafc_J_grid(NULL, n))
  rows <- list()
  for (lo in 1:2) {
    t0 <- proc.time()[["elapsed"]]
    cv <- try(cv.wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], J = lo:Jmax,
                      foldid = foldid), silent = TRUE)
    if (inherits(cv, "try-error")) next
    el <- proc.time()[["elapsed"]] - t0
    f <- cv[["wafc.fit"]]
    d_test <- wafc_design(test[["x"]], test[["u"]], spec = f[["design"]])
    cf <- wafc_raw_coef(f, s = cv[["lambda.min"]])
    fh <- as.numeric(d_test[["Z"]] %*% cf[-1L, 1L]) + cf[1L, 1L]
    rows[[length(rows) + 1L]] <- data.frame(
      cell = cell[["name"]], n = n, rep = r, from = lo, J = cv[["J.min"]],
      lambda = cv[["lambda.min"]], cvm = cv[["cvm.min"]],
      rmse_f = sqrt(mean((fh - test[["f"]])^2)), time = el,
      cvm_J1 = if (lo == 1L) cv[["cvtab"]][["mse"]][1L] else NA_real_,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------------
## Driver
## ---------------------------------------------------------------------------

sweep_part <- function(fun, label, cells_used = cells, ns_used = ns,
                       reps = R) {
  jobs <- list()
  for (cell in cells_used) {
    for (n in ns_used) {
      for (r in seq_len(reps)) {
        jobs[[length(jobs) + 1L]] <- list(cell = cell, n = n, r = r)
      }
    }
  }
  cat(sprintf("\n== %s: %d jobs on %d core(s) ==\n", label, length(jobs),
              ncores))
  t0 <- proc.time()[["elapsed"]]
  res <- parallel::mclapply(jobs, function(j) {
    tryCatch(fun(j[["cell"]], j[["n"]], j[["r"]]),
             error = function(e) {
               message("job failed: ", conditionMessage(e))
               NULL
             })
  }, mc.cores = ncores, mc.preschedule = FALSE)
  cat(sprintf("   %.1f s\n", proc.time()[["elapsed"]] - t0))
  do.call(rbind, res[!vapply(res, is.null, TRUE)])
}

## Median of each variable by the grouping columns. The ordering goes
## through unname(): one of the grouping columns is called "method", which
## order() would take for its own 'method' argument.
med <- function(d, vars, by) {
  a <- stats::aggregate(d[vars], d[by], FUN = stats::median, na.rm = TRUE)
  a[do.call(order, unname(as.list(a[by]))), ]
}

if ("competitors" %in% parts) {
  res <- sweep_part(run_competitors, "competitors")
  saveRDS(res, file.path(out_dir, "e24-competitors.rds"))
  res[["method"]] <- factor(res[["method"]], levels = methods)
  cat("\nmedians by cell, n and method (rmse_f out of sample, ISE of the",
      "components, blocks kept of the active and of the null ones, s):\n")
  print(med(res, c("rmse_f", "mse_beta", "ise", "ise_active", "ise_null",
                   "n_true", "n_false", "time"),
            c("cell", "n", "method")),
        row.names = FALSE, digits = 3)
  cat("\nrelative to the WAFC with the LASSO (rmse_f, median of the ratios",
      "within replicate):\n")
  base <- res[res[["method"]] == "wafc.lasso",
              c("cell", "n", "rep", "rmse_f", "ise")]
  names(base)[4:5] <- c("rmse0", "ise0")
  cmp <- merge(res, base, by = c("cell", "n", "rep"))
  cmp[["r_rmse"]] <- cmp[["rmse_f"]] / cmp[["rmse0"]]
  ## a replicate of the null scenario in which the WAFC zeroes everything
  ## has ise0 = 0, and the ratio there is not a number: it is dropped, not
  ## turned into an infinity that would then be the median
  cmp[["r_ise"]] <- ifelse(cmp[["ise0"]] > 0, cmp[["ise"]] / cmp[["ise0"]],
                           NA_real_)
  print(med(cmp, c("r_rmse", "r_ise"), c("cell", "n", "method")),
        row.names = FALSE, digits = 3)
}

if ("lambda" %in% parts) {
  res <- sweep_part(run_lambda, "lambda rules")
  saveRDS(res, file.path(out_dir, "e24-lambda.rds"))
  cat("\nmedians by cell, n and rule:\n")
  print(med(res, c("J", "lambda", "nzero", "rmse_f", "ise", "n_true",
                   "n_false", "time"), c("cell", "n", "method")),
        row.names = FALSE, digits = 3)
  base <- res[res[["method"]] == "oracle.grid",
              c("cell", "n", "rep", "rmse_f", "lambda")]
  names(base)[4:5] <- c("rmse0", "lam0")
  cmp <- merge(res, base, by = c("cell", "n", "rep"))
  cmp[["cost"]] <- cmp[["rmse_f"]] / cmp[["rmse0"]]
  cmp[["lam_ratio"]] <- cmp[["lambda"]] / cmp[["lam0"]]
  cat("\ncost over the oracle of the grid, and lambda over its lambda:\n")
  print(med(cmp, c("cost", "lam_ratio"), c("cell", "n", "method")),
        row.names = FALSE, digits = 3)
}

if ("margin" %in% parts) {
  cat("\n== margin: lambda_min(G_eps), deterministic ==\n")
  gm <- list()
  for (J in 2:8) {
    eg <- margin_grid(J)
    for (en in names(eg)) {
      e <- eg[[en]]
      if (e >= 0.5) next
      z <- gram_restricted(J, e)
      gm[[length(gm) + 1L]] <- data.frame(J = J, eps.rule = en, eps = e,
                                          lmin = z[["lmin"]],
                                          nzero = z[["nzero"]],
                                          stringsAsFactors = FALSE)
    }
  }
  gm <- do.call(rbind, gm)
  saveRDS(gm, file.path(out_dir, "e24-gram.rds"))
  cat("\nlambda_min(G_eps):\n")
  print(stats::reshape(gm[c("J", "eps.rule", "lmin")], idvar = "J",
                       timevar = "eps.rule", direction = "wide"),
        row.names = FALSE, digits = 3)
  cat("\nnull directions of G_eps:\n")
  print(stats::reshape(gm[c("J", "eps.rule", "nzero")], idvar = "J",
                       timevar = "eps.rule", direction = "wide"),
        row.names = FALSE)

  ## The margin is a question about the basis on [0,1], not about p and q,
  ## and the sweep over J and eps is the most expensive of the four: it runs
  ## on the two cells with q = 2 that have components to resolve, and on the
  ## ends of the range of n.
  res <- sweep_part(run_margin_fit, "margin: error against eps",
                    cells_used = cells[vapply(cells, `[[`, "", "name") %in%
                                         c("smooth", "inhomogeneous")],
                    ns_used = range(ns), reps = R_small)
  saveRDS(res, file.path(out_dir, "e24-margin.rds"))
  cat("\nmedians by cell, n, J and margin:\n")
  print(med(res, c("eps", "rmse_f", "ise", "ise_active", "nzero", "nullcol"),
            c("cell", "n", "J", "eps.rule")),
        row.names = FALSE, digits = 3)
  ## The ranking of the margins is read at the J the cross-validation would
  ## pick for that margin, and not averaged over J: eps and J are not
  ## separable questions, because the rescaling squeezes the sample into
  ## [eps, 1 - eps] and a wide margin is a finer effective resolution at the
  ## same J.
  key <- paste(res[["cell"]], res[["n"]], res[["rep"]], res[["eps.rule"]])
  sel <- unlist(lapply(split(seq_len(nrow(res)), key),
                       function(i) i[which.min(res[["cvm"]][i])]))
  at_cv <- res[sort(sel), ]
  saveRDS(at_cv, file.path(out_dir, "e24-margin-atcv.rds"))
  cat("\nat the J the cross-validation picks for each margin:\n")
  print(med(at_cv, c("J", "eps", "rmse_f", "ise", "ise_active", "nzero",
                     "nullcol"),
            c("cell", "n", "eps.rule")),
        row.names = FALSE, digits = 3)
  cat("\nranking of the margins, by median rmse_f at that J:\n")
  a <- med(at_cv, c("rmse_f", "ise"), c("cell", "n", "eps.rule"))
  for (cn in unique(a[["cell"]])) {
    for (nn in unique(a[["n"]])) {
      b <- a[a[["cell"]] == cn & a[["n"]] == nn, ]
      b <- b[order(b[["rmse_f"]]), ]
      cat(sprintf("  %-14s n = %4d: %s\n", cn, nn,
                  paste(sprintf("%s (%.4f)", b[["eps.rule"]], b[["rmse_f"]]),
                        collapse = "  ")))
    }
  }
}

if ("j1" %in% parts) {
  res <- sweep_part(run_j1, "J = 1 in the grid", reps = R_small)
  saveRDS(res, file.path(out_dir, "e24-j1.rds"))
  cat("\nmedians by cell, n and lower end of the grid:\n")
  print(med(res, c("J", "cvm", "rmse_f", "time"), c("cell", "n", "from")),
        row.names = FALSE, digits = 3)
  cat("\nhow often J = 1 is selected when it is in the grid, and how often",
      "the two grids disagree:\n")
  w <- stats::reshape(res[c("cell", "n", "rep", "from", "J", "rmse_f")],
                      idvar = c("cell", "n", "rep"), timevar = "from",
                      direction = "wide")
  tab <- stats::aggregate(
    cbind(sel1 = w[["J.1"]] == 1L, differ = w[["J.1"]] != w[["J.2"]],
          worse = w[["rmse_f.1"]] > w[["rmse_f.2"]]) ~ cell + n,
    data = w, FUN = mean)
  print(tab[order(tab[["cell"]], tab[["n"]]), ], row.names = FALSE, digits = 3)
}

cat("\nOK\n")
