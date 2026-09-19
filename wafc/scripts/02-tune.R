## E2.3 -- the four selection rules compared on the scenarios of E2.1.
##
##     Rscript 02-tune.R [n_rep]
##
## Five columns: the two readings of the cross-validation (lambda.min and
## lambda.1se), the BIC, the EBIC, and the rule of the theory (J_n of
## equation (Jn) of E1.6 with lambda_n of Corollary 2 of E1.5, sigma known).
## For each scenario and each n it reports the selected (J, lambda), the
## out-of-sample prediction error, the integrated squared error by block and
## the time.

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1L) as.integer(args[1L]) else 30L
p <- 3L
q <- 2L
n_test <- 2000L
n_grid <- 256L
seed0 <- 20260919L
rules <- c("cv.min", "cv.1se", "bic", "ebic", "theory", "theory.sigma")
## effective regularity s' = s - (1/pi - 1/2)_+ of the components of each
## scenario: the cubic of "smooth" has a corner in its periodic extension
## (s' = 3/2, as the bumps of E1.3), and blocks and heavisine jump inside
## the interval (s = 1, pi = 1, s' = 1/2). The null scenario has no
## component; it is run at the value of "smooth" and the column is there to
## show what the rule costs when there is nothing to resolve.
sprime <- c(smooth = 3/2, inhomogeneous = 1/2, null = 3/2)
scenarios <- names(sprime)
ns <- c(250L, 1000L)

ise_of <- function(fit, lambda, dgp, grid, d_grid) {
  design <- fit[["design"]]
  cf <- coef(fit, s = lambda)[, 1L]
  b <- cf[-1L]
  ise <- matrix(NA_real_, p, q)
  nz <- matrix(0L, p, q)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      idx <- design[["blocks"]][[paste0(design[["xnames"]][l], ":",
                                        design[["unames"]][m])]]
      gh <- as.numeric(d_grid[["Z"]][, idx, drop = FALSE] %*% b[idx])
      gm <- dgp[["g"]][[l, m]]
      gt <- if (is.null(gm)) numeric(nrow(grid)) else gm(grid[, m])
      gh <- gh - mean(gh)
      gt <- gt - mean(gt)
      ise[l, m] <- mean((gh - gt)^2) * diff(range(grid[, m]))
      nz[l, m] <- sum(b[idx] != 0)
    }
  }
  list(ise = ise, nz = nz)
}

out <- list()
for (scen in scenarios) {
  for (n in ns) {
    cat(sprintf("\n== scenario %s, n = %d ==\n", scen, n))
    for (r in seq_len(R)) {
      seed <- seed0 + 1000L * match(scen, scenarios) + 10L * match(n, ns) + r
      dgp <- if (scen == "null") {
        simulate_wafc(n, p = p, q = q, scenario = scen, seed = seed,
                      sigma = 0.62)
      } else {
        simulate_wafc(n, p = p, q = q, scenario = scen, seed = seed, snr = 4)
      }
      test <- simulate_wafc(n_test, p = p, q = q, scenario = scen,
                            seed = seed + 500000L, sigma = dgp[["sigma"]])
      foldid <- sample(rep_len(1:10, n))
      active <- nzchar(dgp[["structure"]])
      grid <- matrix(0, n_grid, q)
      for (m in seq_len(q)) {
        rg <- range(dgp[["u"]][, m])
        grid[, m] <- seq(rg[1L], rg[2L], length.out = n_grid)
      }
      ## The design of the test sample and of the grid depend only on J, so
      ## they are built once per J and shared by the rules that land there.
      dcache <- list()
      record <- function(rule, J, lambda, nzero, fit, el) {
        key <- as.character(J)
        if (is.null(dcache[[key]])) {
          dcache[[key]] <<- list(
            grid = wafc_design(matrix(1, n_grid, p), grid,
                               spec = fit[["design"]]),
            test = wafc_design(test[["x"]], test[["u"]],
                               spec = fit[["design"]]))
        }
        d_grid <- dcache[[key]][["grid"]]
        z <- ise_of(fit, lambda, dgp, grid, d_grid)
        cf <- wafc_raw_coef(fit, s = lambda)
        fhat <- as.numeric(dcache[[key]][["test"]][["Z"]] %*% cf[-1L, 1L]) +
          cf[1L, 1L]
        data.frame(scenario = scen, n = n, rep = r, rule = rule,
                   J = J, lambda = lambda, nzero = nzero,
                   nblocks = sum(z[["nz"]] > 0L),
                   nblocks_false = sum(z[["nz"]][!active] > 0L),
                   nblocks_true = sum(z[["nz"]][active] > 0L),
                   rmse_f = sqrt(mean((fhat - test[["f"]])^2)),
                   rmse_y = sqrt(mean((fhat - test[["y"]])^2)),
                   ise = sum(z[["ise"]]),
                   ise_active = sum(z[["ise"]][active]),
                   ise_null = sum(z[["ise"]][!active]),
                   sigma = dgp[["sigma"]], time = el,
                   stringsAsFactors = FALSE)
      }
      for (rule in rules) {
        t0 <- proc.time()[["elapsed"]]
        tn <- wafc_tune(dgp[["x"]], dgp[["u"]], dgp[["y"]],
                        rule = sub("\\.sigma$", "", rule),
                        foldid = foldid, s = sprime[[scen]],
                        sigma = if (rule == "theory.sigma") NULL
                                else dgp[["sigma"]])
        el <- proc.time()[["elapsed"]] - t0
        out[[length(out) + 1L]] <- record(rule, tn[["J"]], tn[["lambda"]],
                                          tn[["nzero"]], tn[["fit"]], el)
      }
      ## The oracle of the grid: the pair (J, lambda) that minimises the
      ## out-of-sample prediction error itself, which no rule can see. It is
      ## the reference against which the cost of every rule is read, and the
      ## J it selects is the "J of minimum realised error" of E1.6.
      t0 <- proc.time()[["elapsed"]]
      best <- NULL
      for (Ji in wafc_J_grid(NULL, n)) {
        fj <- wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], J = Ji)
        key <- as.character(Ji)
        if (is.null(dcache[[key]])) {
          dcache[[key]] <- list(
            grid = wafc_design(matrix(1, n_grid, p), grid,
                               spec = fj[["design"]]),
            test = wafc_design(test[["x"]], test[["u"]],
                               spec = fj[["design"]]))
        }
        cfj <- wafc_raw_coef(fj)
        fh <- sweep(as.matrix(dcache[[key]][["test"]][["Z"]] %*%
                                cfj[-1L, , drop = FALSE]), 2L, cfj[1L, ], "+")
        rm_j <- sqrt(colMeans((fh - test[["f"]])^2))
        k <- which.min(rm_j)
        if (is.null(best) || rm_j[k] < best[["rmse"]]) {
          best <- list(rmse = rm_j[k], J = Ji, lambda = fj[["lambda"]][k],
                       nzero = fj[["nzero"]][k], fit = fj)
        }
      }
      el <- proc.time()[["elapsed"]] - t0
      out[[length(out) + 1L]] <- record("oracle", best[["J"]],
                                        best[["lambda"]], best[["nzero"]],
                                        best[["fit"]], el)
      if (r %% 5L == 0L) cat(sprintf("  %d/%d\n", r, R))
    }
  }
}
res <- do.call(rbind, out)
saveRDS(res, file.path(Sys.getenv("WAFC_OUT", "."), "e23-rules.rds"))

agg <- aggregate(cbind(J, lambda, nzero, nblocks, nblocks_false, rmse_f,
                       rmse_y, ise, ise_active, ise_null, time) ~
                   scenario + n + rule, data = res, FUN = median)
agg <- agg[order(agg[["scenario"]], agg[["n"]],
                 match(agg[["rule"]], c(rules, "oracle"))), ]
print(agg, row.names = FALSE, digits = 3)
cat("\nOK\n")
