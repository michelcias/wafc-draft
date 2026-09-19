## Decomposition of the cost of the theoretical rule: how much comes from
## the resolution J_n of equation (Jn) and how much from the penalty level
## lambda_n of Corollary 2. Same seeds as 02-tune.R.
local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  source(hit[1L])
})
R <- 20L; p <- 3L; q <- 2L; n_test <- 2000L; seed0 <- 20260919L
sprime <- c(smooth = 3/2, inhomogeneous = 1/2, null = 3/2)
scenarios <- names(sprime); ns <- c(250L, 1000L)
out <- list()
for (scen in scenarios) for (n in ns) {
  for (r in seq_len(R)) {
    seed <- seed0 + 1000L * match(scen, scenarios) + 10L * match(n, ns) + r
    dgp <- if (scen == "null")
      simulate_wafc(n, p = p, q = q, scenario = scen, seed = seed, sigma = 0.62)
      else simulate_wafc(n, p = p, q = q, scenario = scen, seed = seed, snr = 4)
    test <- simulate_wafc(n_test, p = p, q = q, scenario = scen,
                          seed = seed + 500000L, sigma = dgp[["sigma"]])
    Jn <- wafc_J_theory(n, s = sprime[[scen]])
    grid <- wafc_J_grid(NULL, n)
    fits <- lapply(grid, function(Ji) wafc(dgp$x, dgp$u, dgp$y, J = Ji))
    names(fits) <- as.character(grid)
    rmse_path <- lapply(fits, function(fj) {
      fh <- predict(fj, test$x, test$u)
      sqrt(colMeans((fh - test$f)^2))
    })
    best_each <- vapply(rmse_path, min, 0)
    Jo <- grid[which.min(best_each)]
    lam_n <- function(J) wafc_lambda_theory(fits[[as.character(J)]]$design,
                                            sigma = dgp$sigma)
    at <- function(J, lam) {
      fj <- wafc(design = fits[[as.character(J)]]$design, y = dgp$y,
                 lambda = wafc_path_to(lam, fits[[as.character(J)]]$design, dgp$y))
      cf <- wafc_raw_coef(fj, s = lam)
      Zt <- wafc_design(test$x, test$u, spec = fj$design)$Z
      sqrt(mean((as.numeric(Zt %*% cf[-1L, 1L]) + cf[1L, 1L] - test$f)^2))
    }
    out[[length(out) + 1L]] <- data.frame(
      scenario = scen, n = n, rep = r, Jn = Jn, Jo = Jo,
      rmse_oracle = min(best_each),                      # (J_o, lambda_o)
      rmse_Jn_best = min(rmse_path[[as.character(Jn)]]), # (J_n, lambda_o at J_n)
      rmse_Jo_lamn = at(Jo, lam_n(Jo)),                  # (J_o, lambda_n)
      rmse_Jn_lamn = at(Jn, lam_n(Jn)),                  # the rule itself
      lam_n_Jn = lam_n(Jn), lam_n_Jo = lam_n(Jo),
      lam_o = fits[[as.character(Jo)]]$lambda[which.min(rmse_path[[as.character(Jo)]])],
      stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, out)
saveRDS(res, file.path(Sys.getenv("WAFC_OUT", "."), "e23-decomp.rds"))
res$c_J <- res$rmse_Jn_best / res$rmse_oracle
res$c_lam <- res$rmse_Jo_lamn / res$rmse_oracle
res$c_both <- res$rmse_Jn_lamn / res$rmse_oracle
agg <- aggregate(cbind(Jn, Jo, c_J, c_lam, c_both, lam_n_Jn, lam_o) ~ scenario + n,
                 data = res, FUN = median)
print(agg[order(agg$scenario, agg$n), ], row.names = FALSE, digits = 3)
cat("\nOK\n")
