## wafc/scripts/01-smoke.R -- smoke test of step E2.1: one scenario, n = 500,
## the design of wafc_design() fed to glmnet, and the estimated components
## against the truth. From the root of the repository:
##
##     Rscript wafc/scripts/01-smoke.R [n] [J] [scenario]
##
## Defaults: n = 500, J = 4, scenario = "smooth". It prints a table of
## integrated squared errors, the prediction error on a fresh sample and the
## selected structure, and writes wafc/cache/01-smoke.png (not versioned).
##
## There is no wafc() yet: the fit here is a direct call to glmnet, which is
## what E2.2 will wrap. Two things this script pins down for that step:
## glmnet drops the constant column X_1 and carries c_1 in its own
## intercept, so the fit uses intercept = TRUE; and the estimated component
## is centred over the rescaled range of U_m, not over [0,1], so the
## comparison with the truth recentres both on the comparison grid.

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 1L) as.integer(args[1L]) else 500L
J <- if (length(args) >= 2L) as.integer(args[2L]) else 4L
scenario <- if (length(args) >= 3L) args[3L] else "smooth"
p <- 3L
q <- 2L
seed <- 20260919L

cat(sprintf("WAFC smoke test: scenario %s, n = %d, p = %d, q = %d, J = %d\n",
            scenario, n, p, q, J))

## ---- data ------------------------------------------------------------------

dgp <- simulate_wafc(n, p = p, q = q, scenario = scenario, seed = seed,
                     snr = 4)
test <- simulate_wafc(2000L, p = p, q = q, scenario = scenario,
                      seed = seed + 1L, sigma = dgp[["sigma"]])
cat(sprintf("  sigma = %.4f (snr 4), sd(f) = %.4f\n",
            dgp[["sigma"]], stats::sd(dgp[["f"]])))

## ---- design ----------------------------------------------------------------

t0 <- proc.time()[["elapsed"]]
d <- wafc_design(dgp[["x"]], dgp[["u"]], J = J)
t_design <- proc.time()[["elapsed"]] - t0
cat(sprintf("  design: %d x %d (%s), %d unpenalized column(s), eps = %.4f, %.2f s\n",
            d[["n"]], d[["nvars"]], if (d[["sparse"]]) "sparse" else "dense",
            length(d[["unpenalized"]]), d[["eps"]][1L], t_design))
cat(sprintf("  constant linear covariate(s): %s\n",
            if (length(d[["constant"]]) == 0L) "none"
            else paste(d[["xnames"]][d[["constant"]]], collapse = ", ")))

## ---- fit -------------------------------------------------------------------

## intercept = TRUE because glmnet drops the constant column X_1 and lets
## its intercept carry c_1; standardize = FALSE because the wavelet basis is
## already orthonormal, as in wall().
t0 <- proc.time()[["elapsed"]]
cvfit <- glmnet::cv.glmnet(d[["Z"]], dgp[["y"]], family = "gaussian",
                           intercept = TRUE, standardize = FALSE,
                           penalty.factor = d[["penalty.factor"]], nfolds = 10L)
t_fit <- proc.time()[["elapsed"]] - t0
cat(sprintf("  fit: %.2f s, lambda.min = %.5f, lambda.1se = %.5f\n",
            t_fit, cvfit[["lambda.min"]], cvfit[["lambda.1se"]]))

## ---- evaluation at one value of lambda --------------------------------------

d_test <- wafc_design(test[["x"]], test[["u"]], spec = d)

## g_hat_{lm}(u) is the block (l, m) of the design evaluated with X_l = 1 at
## a grid of values of U_m, times the estimated coefficients of the block.
## One design call on a matrix whose column m is the grid of U_m gives every
## block at once.
G <- 512L
grid <- matrix(0, G, q)
for (m in seq_len(q)) {
  rg <- range(dgp[["u"]][, m])
  grid[, m] <- seq(rg[1L], rg[2L], length.out = G)
}
d_grid <- wafc_design(matrix(1, G, p), grid, spec = d)

evaluate <- function(s) {
  cf <- as.numeric(stats::coef(cvfit, s = s))
  a0 <- cf[1L]
  b_hat <- cf[-1L]
  c_hat <- b_hat[d[["unpenalized"]]]
  c_hat[d[["constant"]]] <- c_hat[d[["constant"]]] + a0
  f_hat <- as.numeric(a0 + d_test[["Z"]] %*% b_hat)
  ise <- matrix(NA_real_, p, q, dimnames = list(d[["xnames"]], d[["unames"]]))
  nz <- matrix(0L, p, q, dimnames = dimnames(ise))
  g_hat <- vector("list", p * q); dim(g_hat) <- c(p, q)
  g_true <- g_hat
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      idx <- d[["blocks"]][[paste0(d[["xnames"]][l], ":", d[["unames"]][m])]]
      gh <- as.numeric(d_grid[["Z"]][, idx, drop = FALSE] %*% b_hat[idx])
      gm <- dgp[["g"]][[l, m]]
      gt <- if (is.null(gm)) numeric(G) else gm(grid[, m])
      ## both recentred on the comparison grid: the estimate integrates to
      ## zero over the rescaled range of U_m, the truth over [0,1]
      g_hat[[l, m]] <- gh - mean(gh)
      g_true[[l, m]] <- gt - mean(gt)
      ise[l, m] <- mean((g_hat[[l, m]] - g_true[[l, m]])^2) *
        diff(range(grid[, m]))
      nz[l, m] <- sum(b_hat[idx] != 0)
    }
  }
  list(a0 = a0, c_hat = c_hat, nz = nz, ise = ise, g_hat = g_hat,
       g_true = g_true, nzero = sum(b_hat[-d[["unpenalized"]]] != 0),
       rmse_f = sqrt(mean((f_hat - test[["f"]])^2)),
       rmse_y = sqrt(mean((f_hat - test[["y"]])^2)))
}

res <- lapply(c("lambda.min", "lambda.1se"), evaluate)
names(res) <- c("lambda.min", "lambda.1se")
npen <- d[["nvars"]] - length(d[["unpenalized"]])
for (s in names(res)) {
  r <- res[[s]]
  cat(sprintf("\n  %s: %d of %d wavelet coefficients, sum of ISE = %.4f, RMSE(f_hat, f) = %.4f on %d fresh observations (sigma = %.4f)\n",
              s, r[["nzero"]], npen, sum(r[["ise"]]), r[["rmse_f"]],
              test[["n"]], dgp[["sigma"]]))
  cat(sprintf("    c_hat = (%s) against c = (%s)\n",
              paste(sprintf("%.3f", r[["c_hat"]]), collapse = ", "),
              paste(sprintf("%.3f", dgp[["cc"]]), collapse = ", ")))
  cat("    ISE of g_hat_{lm}, with the non-zero coefficients of the block:\n")
  tab <- matrix(sprintf("%.4f (%d)", r[["ise"]], r[["nz"]]), p, q,
                dimnames = dimnames(r[["ise"]]))
  print(tab, quote = FALSE)
}
cat("\n  True structure (empty = zero component):\n")
struct <- dgp[["structure"]]
dimnames(struct) <- list(d[["xnames"]], d[["unames"]])
print(struct, quote = FALSE)

## ---- plot ------------------------------------------------------------------

outdir <- if (dir.exists("wafc")) file.path("wafc", "cache") else "cache"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, "01-smoke.png")
best <- res[["lambda.min"]]
g_hat <- best[["g_hat"]]
g_true <- best[["g_true"]]
ise <- best[["ise"]]
grDevices::png(outfile, width = 300 * q, height = 250 * p, res = 110)
op <- graphics::par(mfrow = c(p, q), mar = c(4, 4, 2.5, 1))
for (l in seq_len(p)) {
  for (m in seq_len(q)) {
    yl <- range(c(g_hat[[l, m]], g_true[[l, m]]))
    if (diff(yl) < 1e-8) yl <- yl + c(-1, 1) * 0.5
    graphics::plot(grid[, m], g_true[[l, m]], type = "l", ylim = yl,
                   xlab = d[["unames"]][m], ylab = "",
                   main = sprintf("g[%s,%s]  ISE %.4f", d[["xnames"]][l],
                                  d[["unames"]][m], ise[l, m]))
    graphics::lines(grid[, m], g_hat[[l, m]], col = "red", lty = 2)
  }
}
graphics::par(op)
grDevices::dev.off()
cat(sprintf("\n  plot written to %s\n", outfile))
cat("OK\n")
