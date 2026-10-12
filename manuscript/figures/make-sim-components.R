## Figure 1 of ms_6.tex (figures/sim-components.pdf) and Figure S1 of supp_6.tex
## (figures/sim-bands.pdf), E5b, 2026-10-12. Run from the root of wafc-draft on the
## machine that holds the outputs of E4.4: it reads
## wafc-studies/outputs/core/tables/curves.csv (not versioned; its SHA-256 is in
## results/README.md), keeps the inhomogeneous cell at n = 1000 and the components
## g_11 (bumps) and g_12 (blocks) of the WAFC and of GAM (REML), interpolates each
## replicate onto a common grid of 400 points of [0.005, 0.995] (each replicate has
## its own grid of 256 points in the observed range of the modulator) and takes the
## pointwise median and the 5% and 95% quantiles over the 100 replicates.
src <- "wafc-studies/outputs/core/tables/curves.csv"
cv <- read.csv(pipe(paste("awk -F, 'NR==1 || ($1==\"\\\"inhomogeneous\\\"\" && $2==1000)'", src)))
cv <- cv[cv$l == 1 & cv$m %in% 1:2, ]
stopifnot(setequal(unique(cv$method), c("wafc.cv1se", "gam.reml")),
          length(unique(cv$rep)) == 100)
ug <- seq(0.005, 0.995, length.out = 400)
on_grid <- function(meth, m, v) sapply(1:100, function(r) {
  x <- cv[cv$method == meth & cv$m == m & cv$rep == r, ]
  approx(x$u, x[[v]], ug, rule = 2)$y
})
q <- function(meth, m) t(apply(on_grid(meth, m, "ghat"), 1, quantile, c(0.05, 0.5, 0.95)))
truth <- lapply(1:2, function(m) on_grid("wafc.cv1se", m, "truth")[, 1])
B <- list()
for (meth in c("wafc.cv1se", "gam.reml")) for (m in 1:2) B[[paste(meth, m)]] <- q(meth, m)
nm <- c("Bumps", "Blocks")

pdf("manuscript/figures/sim-components.pdf", width = 5.5, height = 2.4, pointsize = 8)
par(mfrow = c(1, 2), mar = c(2.8, 2.6, 1.4, 0.4), mgp = c(1.6, 0.45, 0))
for (m in 1:2) {
  a <- B[[paste("wafc.cv1se", m)]][, 2]; g <- B[[paste("gam.reml", m)]][, 2]
  plot(ug, truth[[m]], type = "l", col = "grey70", lwd = 3, xlab = "u", ylab = "",
       main = nm[m], ylim = range(c(truth[[m]], a, g)))
  lines(ug, a, lwd = 1); lines(ug, g, lty = 2, lwd = 1)
  if (m == 1) legend("topright", c("truth", "WAFC", "GAM (REML)"), lwd = c(3, 1, 1),
                     lty = c(1, 1, 2), col = c("grey70", 1, 1), bty = "n")
}
invisible(dev.off())

pdf("manuscript/figures/sim-bands.pdf", width = 5.5, height = 4.4, pointsize = 8)
par(mfrow = c(2, 2), mar = c(2.8, 2.6, 1.4, 0.4), mgp = c(1.6, 0.45, 0))
for (meth in c("wafc.cv1se", "gam.reml")) for (m in 1:2) {
  b <- B[[paste(meth, m)]]
  plot(ug, b[, 2], type = "n", ylim = range(c(b, truth[[m]])), xlab = "u", ylab = "",
       main = paste(nm[m], if (meth == "wafc.cv1se") "WAFC" else "GAM (REML)", sep = ", "))
  polygon(c(ug, rev(ug)), c(b[, 1], rev(b[, 3])), col = "grey80", border = NA)
  lines(ug, truth[[m]], lty = 2); lines(ug, b[, 2])
}
invisible(dev.off())
