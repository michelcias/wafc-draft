## Generates the rows of Table 1 of ms_6.tex (tables/sim-main_6.tex) from
## results/tables/core/summary.csv (E5b). Run from the root of wafc-draft; prints the
## rows of the two panels, which go between the header and the \hline of the float
## (panel (a) without the standard errors, whose range is in the caption; panel (b)
## without the oracle support and the constant coefficients).
d <- read.csv("results/tables/core/summary.csv"); d <- d[d$n == 1000, ]
m <- c(wafc.cv1se = "WAFC", wafc.cv = "WAFC (cv)", gam.reml = "GAM (REML)",
       gam.gcv = "GAM (GCV)", bsgl = "B-spline group lasso",
       aspline = "Adaptive spline", klopp = "Block lasso (KP)",
       oracle = "Oracle support", linear = "Constant coefficients")
cells <- c("smooth", "uneven", "inhomogeneous", "mixed", "null")
get <- function(mm, v) { x <- d[d$method == mm, ]; x[match(cells, x$cell), v] }
f <- function(a, s, k, mult = 1) ifelse(is.na(a), "--",
  sprintf(paste0("%.", k, "f (%.", k, "f)"), a * mult, s * mult))
out <- character(0)
panel <- function(title, rows) c(sprintf("\\multicolumn{6}{l}{%s} \\\\", title), rows)
ra <- vapply(names(m), function(mm) paste(m[[mm]], "&", paste(f(get(mm, "rmse_f_mean"), get(mm, "rmse_f_se"), 3), collapse = " & "), "\\\\"), "")
rb <- vapply(names(m), function(mm) paste(m[[mm]], "&", paste(f(get(mm, "ise_mean"), get(mm, "ise_se"), 2, 100), collapse = " & "), "\\\\"), "")
rc <- vapply(names(m), function(mm) { v <- get(mm, "exact_support"); paste(m[[mm]], "&", paste(ifelse(is.na(v), "--", sprintf("%.2f", v)), collapse = " & "), "\\\\") }, "")
out <- c(panel("(a) Prediction error", ra), "\\hline",
         panel("(b) Error of the components, $\\times 100$", rb), "\\hline",
         panel("(c) Frequency of the exact structure, $\\Prob(\\widehat{\\mathcal{S}} = \\mathcal{S})$", rc))
writeLines(out)
