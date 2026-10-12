## Generates the tables of Section S9 of supp_6.tex from results/tables/ (E5b).
## Run from the root of wafc-draft; writes manuscript/tables/sim-*_6.tex.
core <- read.csv("results/tables/core/summary.csv")
pair <- read.csv("results/tables/core/paired.csv")
arms <- rbind(read.csv("results/tables/arms/summary.csv"),
              read.csv("results/tables/scale/summary.csv"))
slopes <- read.csv("results/tables/core/slopes.csv")

lab <- c(wafc.cv1se = "WAFC", wafc.cv = "WAFC (cv)", lasso = "Lasso",
         lasso.cv1se = "Lasso + threshold", gam.reml = "GAM (REML)",
         gam.gcv = "GAM (GCV)", bsgl = "B-spline group lasso",
         aspline = "Adaptive spline", klopp = "Block lasso (KP)",
         oracle = "Oracle support")
cells <- c(smooth = "smooth", uneven = "uneven", inhomogeneous = "inhomogeneous",
           mixed = "mixed", null = "null")
ns <- c(250, 500, 1000, 2000)
val <- function(d, cell, n, m, v) {
  x <- d[d$cell == cell & d$n == n & d$method == m, v]
  if (length(x) == 0 || is.na(x)) NA_real_ else x
}
fmt <- function(x, k) ifelse(is.na(x), "--", formatC(x, format = "f", digits = k))
row <- function(name, xs) paste0(name, " & ", paste(xs, collapse = " & "), " \\\\")
head_float <- function(caption, label, size = "\\footnotesize", sep = "3pt") c(
  "\\begin{table}[tp]", "\\color{colR1}", paste0("\\caption{", caption, "}"),
  paste0("\\label{", label, "}"), "\\centering", size,
  paste0("\\setlength{\\tabcolsep}{", sep, "}"))
tail_float <- c("\\hline", "\\end{tabular}", "\\end{table}", "")
hdr <- "% Gerado por tables/make-sim-supp.R (E5b, 2026-10-12) de results/tables/ (o compendio e a verdade).\n% \\color{colR1} dentro do float (D75(d))."

## Tables S1 to S5: one per cell, RMSE, ISE x 100 and P(exact) at every n.
out <- hdr
for (ci in seq_along(cells)) {
  cell <- cells[[ci]]
  name <- c(smooth = "smooth", uneven = "uneven",
            inhomogeneous = "inhomogeneous", mixed = "mixed",
            null = "null")[[cell]]
  cap <- sprintf(paste0("Simulation, the %s cell, 100 replicates per sample size: ",
    "the root mean squared error of prediction, the integrated squared error of ",
    "the components times 100, and the frequency of the exact structure."), name)
  if (ci == 1) cap <- paste(cap, paste0("The constant coefficients have a ",
    "prediction error of 1.73 in the cells with varying coefficients, and ",
    "coincide with the oracle support in the null cell; GAM (GCV) is not run ",
    "in the mixed cell, and the lasso rows are the coordinatewise lasso of ",
    "Section 4 of the paper."))
  lb <- if (ci == 1) "tab:sim-sizes" else paste0("tab:sim-sizes-", cell)
  out <- c(out, head_float(cap, lb), "\\begin{tabular}{l cccc cccc cccc}", "\\hline",
           " & \\multicolumn{4}{c}{Prediction error} & \\multicolumn{4}{c}{Components ($\\times 100$)} & \\multicolumn{4}{c}{Exact structure} \\\\",
           paste0("$n$ & ", paste(rep(ns, 3), collapse = " & "), " \\\\"), "\\hline")
  for (m in names(lab)) {
    if (m == "gam.gcv" && cell == "mixed") next
    r <- c(fmt(sapply(ns, function(n) val(core, cell, n, m, "rmse_f_mean")), 3),
           fmt(100 * sapply(ns, function(n) val(core, cell, n, m, "ise_mean")), 2),
           fmt(sapply(ns, function(n) val(core, cell, n, m, "exact_support")), 2))
    out <- c(out, row(lab[[m]], r))
  }
  out <- c(out, tail_float)
}
writeLines(out, "manuscript/tables/sim-sizes_6.tex")

## Table S6: structure at n = 1000, null blocks kept and active blocks dropped;
## Table S7: paired comparisons at n = 1000.
out <- c(hdr, head_float(paste0("Simulation at $n = 1000$, the structure: the ",
  "mean number of null blocks kept per replicate; no method drops an active block ",
  "in any replicate at this sample size. ",
  "The smooth, uneven and inhomogeneous cells have 3 active and 3 null blocks, ",
  "the mixed cell 3 active and 13 null, and the null cell 6 null blocks. A ",
  "smooth of \\texttt{mgcv} counts as kept when its effective degrees of ",
  "freedom exceed 0.1."), "tab:sim-structure", "\\small", "6pt"),
  "\\begin{tabular}{lccccc}", "\\hline",
  " & Smooth & Uneven & Inhomogeneous & Mixed & Null \\\\", "\\hline")
stopifnot(all(core$false_neg_mean[core$n == 1000 & core$method %in% setdiff(names(lab), "linear")] == 0, na.rm = TRUE))
for (m in setdiff(names(lab), "oracle")) {
  r <- fmt(sapply(cells, function(cl) val(core, cl, 1000, m, "false_pos_mean")), 2)
  out <- c(out, row(lab[[m]], r))
}
out <- c(out, tail_float)
pv <- function(cl, m) {
  x <- pair[pair$cell == cl & pair$n == 1000 & pair$method == m, ]
  if (!nrow(x)) "--" else sprintf("%.2f (%.2f)", x$rmse_f_ratio_median, x$rmse_f_below_ref)
}
out <- c(out, head_float(paste0("Simulation at $n = 1000$, paired comparisons: ",
  "the median over the replicates of the ratio of the prediction error of each ",
  "method to that of the WAFC in the same replicate, and in parentheses the ",
  "fraction of replicates in which the method has the smaller error."),
  "tab:sim-paired", "\\small", "6pt"),
  "\\begin{tabular}{lccccc}", "\\hline",
  " & Smooth & Uneven & Inhomogeneous & Mixed & Null \\\\", "\\hline")
for (m in setdiff(names(lab), "wafc.cv1se")) out <- c(out, row(lab[[m]], sapply(cells, pv, m = m)))
out <- c(out, tail_float)
writeLines(out, "manuscript/tables/sim-structure_6.tex")

## Table S8: the arms and the larger model at n = 1000.
ac <- c("inhomogeneous.snr1", "inhomogeneous.urho", "inhomogeneous.xu",
        "mixed.snr1", "mixed.urho", "mixed.xu", "scale")
out <- c(hdr, head_float(paste0("Simulation, the arms at $n = 1000$ and the larger ",
  "model, 100 replicates each: (a) the root mean squared error of prediction; (b) ",
  "the frequency of the exact structure. SNR 1: signal to noise ratio 1; copula: ",
  "modulators coupled by a Gaussian copula; $X$ on $U$: linear covariates ",
  "dependent on the modulators; larger: $(p,q) = (6,4)$ with 6 of 24 blocks ",
  "active. GAM (GCV) is run only in the arms of the inhomogeneous cell."),
  "tab:sim-arms", "\\footnotesize", "4pt"),
  "\\begin{tabular}{lccccccc}", "\\hline",
  " & \\multicolumn{3}{c}{Inhomogeneous} & \\multicolumn{3}{c}{Mixed} & \\\\",
  " & SNR 1 & Copula & $X$ on $U$ & SNR 1 & Copula & $X$ on $U$ & Larger \\\\", "\\hline")
for (v in c("rmse_f_mean", "exact_support")) {
  out <- c(out, sprintf("\\multicolumn{8}{l}{%s} \\\\",
    if (v == "rmse_f_mean") "(a) Prediction error" else "(b) Exact structure"))
  for (m in names(lab)) {
    if (v == "exact_support" && m == "oracle") next
    r <- fmt(sapply(ac, function(cl) val(arms, cl, 1000, m, v)), if (v == "rmse_f_mean") 3 else 2)
    out <- c(out, row(lab[[m]], r))
  }
  if (v == "rmse_f_mean") out <- c(out, "\\hline")
}
out <- c(out, tail_float)
writeLines(out, "manuscript/tables/sim-arms_6.tex")

## Tables S9 and S10: tuning and time; Table S11: slopes.
tun <- function(cl_set, lb, first) {
  cap <- paste0("Simulation, tuning and time: the mean resolution $J$ chosen by ",
    "the WAFC and the mean dimension $k$ chosen by the two splines, each followed ",
    "by the fraction of replicates at the top of its grid, and the median time ",
    "of a fit with its tuning, in seconds.")
  if (!first) cap <- "Simulation, tuning and time (continued from Table~\\ref{tab:sim-tuning})."
  o <- c(head_float(cap, lb, "\\footnotesize", "4pt"),
         "\\begin{tabular}{ll cc cc cc ccc}", "\\hline",
         " & & \\multicolumn{2}{c}{WAFC} & \\multicolumn{2}{c}{GAM (REML)} & \\multicolumn{2}{c}{GAM (GCV)} & \\multicolumn{3}{c}{Time (s)} \\\\",
         "Cell & $n$ & $J$ & Top & $k$ & Top & $k$ & Top & WAFC & REML & GCV \\\\", "\\hline")
  for (cl in cl_set) for (n in ns) {
    g <- function(m, v, k) fmt(val(core, cl, n, m, v), k)
    o <- c(o, paste0(if (n == 250) paste0(toupper(substr(cl, 1, 1)), substring(cl, 2)) else "",
      " & ", n, " & ", paste(c(g("wafc.cv1se", "chosen_mean", 2), g("wafc.cv1se", "top_fraction", 2),
      g("gam.reml", "chosen_mean", 0), g("gam.reml", "top_fraction", 2),
      g("gam.gcv", "chosen_mean", 0), g("gam.gcv", "top_fraction", 2),
      g("wafc.cv1se", "time_median", 1), g("gam.reml", "time_median", 1),
      g("gam.gcv", "time_median", 1)), collapse = " & "), " \\\\"))
  }
  c(o, tail_float)
}
out <- c(hdr, tun(c("smooth", "uneven", "null"), "tab:sim-tuning", TRUE),
         tun(c("inhomogeneous", "mixed"), "tab:sim-tuning-b", FALSE))
out <- c(out, head_float(paste0("Simulation, the slope of the least squares line of ",
  "the logarithm of the mean squared error of prediction against $\\log n$, over ",
  "$n = 250$, 500, 1000 and 2000."), "tab:sim-slopes", "\\small", "6pt"),
  "\\begin{tabular}{lccccc}", "\\hline",
  " & Smooth & Uneven & Inhomogeneous & Mixed & Null \\\\", "\\hline")
for (m in names(lab)) {
  r <- sapply(cells, function(cl) { x <- slopes$slope_mse_f[slopes$cell == cl & slopes$method == m]; if (length(x)) sprintf("$-$%.2f", -x) else "--" })
  out <- c(out, row(lab[[m]], r))
}
out <- c(out, tail_float)
writeLines(out, "manuscript/tables/sim-tuning_6.tex")
