## E2.4b -- what each method costs, with the tuning separated from the final
## fit.
##
##     Rscript wafc/scripts/06-timing.R [n_rep] [parts] [ns]
##
## 'parts' is a comma separated subset of
##
##   methods    one scenario, three sample sizes, every competitor of the
##              pilot, with the seconds split into the search over the
##              tuning parameters and the fit that search selects;
##   tolerance  what the convergence threshold of wafc() costs and what it
##              buys, which is the second defect step E2.4b was opened to
##              fix.
##
## Why the split. The table of step E2.4 reports one number per method, and
## the numbers are not comparable: the seconds of the WAFC are a whole
## cross-validation over the pair (J, lambda), while the seconds of the gam
## are a single call in which mgcv chooses the smoothing parameters by REML
## on the inside. Reading them side by side says the WAFC is expensive when
## what it says is that one of the two was asked to show its search and the
## other was not. This script reports three columns, secs.tune, secs.fit
## and secs.total, and a fourth that says whether the split is explicit or
## whether the method hides its search inside one call.
##
## Why the gam appears three times. Step E6.1a measured that comparing the
## WAFC with mgcv at the default k = 10, against an expansion with 2^J - 1
## columns per block, is a comparison of dimension and not of basis: on the
## three real candidates the apparent gain of 6.6% to 15.2% became a tie at
## matched dimension. The pilot ran at k = 10. So the cost of the honest
## comparison has to be measured too, and it is not the same number: at
## matched dimension, gam() is the call that stops being affordable, which
## is why mgcv::bam with fREML and discretized covariates is measured beside
## it (step E2.4b added both to wafc_fit_gam()).
##
## This is a timing script. It does not report an error, and no verdict
## about any method belongs to it; that is step E2.5.

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1L && nzchar(args[1L])) as.integer(args[1L]) else 5L
parts <- if (length(args) >= 2L && nzchar(args[2L])) {
  strsplit(args[2L], ",", fixed = TRUE)[[1L]]
} else c("methods", "tolerance")
ns <- if (length(args) >= 3L && nzchar(args[3L])) {
  as.integer(strsplit(args[3L], ",", fixed = TRUE)[[1L]])
} else c(250L, 500L, 1000L)

seed0 <- 20260921L
n_test <- 1000L

## Decision D31, as in step E2.4: the basis is evaluated by table, built
## once. A timing that included the construction of the basis in every fit
## would measure the loader.
wafc_timing_table <- WaveBased::wtable(family = "Daublets", filter.size = 8L,
                                       prec.wavelet = 30L, check = FALSE)

## The cell is the one that carries the verdict of step E2.5: the
## inhomogeneous scenario at p = 3, q = 2, where the wavelet expansion is
## supposed to be worth its cost. One scenario, because the question here
## is the cost of an algorithm and not the difficulty of a function.
cell <- list(scenario = "inhomogeneous", p = 3L, q = 2L, snr = 3)

draw <- function(n, seed) {
  simulate_wafc(n, p = cell[["p"]], q = cell[["q"]],
                scenario = cell[["scenario"]], seed = seed, snr = cell[["snr"]])
}

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, secs = proc.time()[["elapsed"]] - t0)
}

row_of <- function(n, r, method, tune, fit, split, note = "") {
  data.frame(n = n, rep = r, method = method,
             secs.tune = tune, secs.fit = fit,
             secs.total = sum(c(tune, fit), na.rm = TRUE),
             split = split, note = note, stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## Part "methods"
## ---------------------------------------------------------------------------

run_methods <- function(n, r) {
  seed <- seed0 + 1000L * match(n, ns) + r
  dgp <- draw(n, seed)
  x <- dgp[["x"]]; u <- dgp[["u"]]; y <- dgp[["y"]]
  active <- nzchar(dgp[["structure"]])
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:10, n))
  rows <- list()

  ## The two WAFC variants. The search is the cross-validation over the
  ## pair (J, lambda); the fit it selects is one path down to the selected
  ## penalty level at the selected resolution, which is what a user who
  ## already knew the pair would pay.
  Jmin <- NA_integer_
  for (pen in c("lasso", "sglasso")) {
    tn <- try(timed(cv.wafc(x, u, y, penalty = pen, foldid = foldid,
                            wavelet.table = wafc_timing_table)), silent = TRUE)
    if (inherits(tn, "try-error")) next
    cv <- tn[["value"]]
    des <- cv[["wafc.fit"]][["design"]]
    lam <- cv[["lambda.min"]]
    ft <- timed(wafc(design = des, y = y, penalty = pen,
                     lambda = wafc_path_to(lam, des, y)))
    if (pen == "lasso") Jmin <- cv[["J.min"]]
    rows[[length(rows) + 1L]] <- row_of(
      n, r, paste0("wafc.", pen), tn[["secs"]], ft[["secs"]], "explicit",
      sprintf("J = %d of %s", cv[["J.min"]],
              paste(range(cv[["J"]]), collapse = ":")))
  }

  ## The gam, three ways. The smoothing parameters are chosen by REML
  ## inside the one call, so there is no split to report and the whole
  ## number goes to secs.fit, which is exactly the asymmetry this script
  ## exists to make visible.
  if (is.na(Jmin)) Jmin <- 4L
  kmat <- wafc_k_matched(u, Jmin)
  gam_specs <- list(
    list(name = "gam.k10", k = 10L, engine = "gam",
         note = "default basis size of the pilot"),
    list(name = "gam.matched", k = kmat, engine = "gam",
         note = sprintf("k = (%s), matched to 2^%d",
                        paste(kmat, collapse = ", "), Jmin)),
    list(name = "gam.matched.bam", k = kmat, engine = "bam",
         note = sprintf("k = (%s), by bam/fREML/discrete",
                        paste(kmat, collapse = ", "))))
  for (sp in gam_specs) {
    ft <- try(timed(wafc_competitor("gam", x, u, y, k = sp[["k"]],
                                    engine = sp[["engine"]])), silent = TRUE)
    if (inherits(ft, "try-error")) {
      rows[[length(rows) + 1L]] <- row_of(n, r, sp[["name"]], NA_real_,
                                          NA_real_, "failed",
                                          conditionMessage(attr(ft, "condition")))
      next
    }
    rows[[length(rows) + 1L]] <- row_of(n, r, sp[["name"]], NA_real_,
                                        ft[["secs"]], "internal",
                                        sp[["note"]])
  }

  ## The remaining competitors. Every one of them searches on the inside,
  ## by cross-validation (bsgl, klopp), by a BIC scored knot insertion
  ## (aspline) or by an MCMC that has no tuning parameter to search over
  ## (vcbart); the oracle is not timed, because knowing the structure is
  ## not a cost anyone pays.
  ## The basis table goes only to the two methods that build a WAFC
  ## design. A fitter with a '...' of its own forwards what it does not
  ## know to its engine, and an engine that is not the WAFC stops with
  ## "unused argument": that is how the pilot lost its vcbart column
  ## without saying so (see docs/handoff-E2.4b.md).
  for (mth in c("bsgl", "klopp", "aspline", "vcbart", "linear")) {
    a <- list(mth, x, u, y, active = active, foldid = foldid)
    if (mth %in% c("klopp", "oracle")) {
      a[["wavelet.table"]] <- wafc_timing_table
    }
    ft <- try(timed(do.call(wafc_competitor, a)), silent = TRUE)
    if (inherits(ft, "try-error")) {
      rows[[length(rows) + 1L]] <- row_of(
        n, r, mth, NA_real_, NA_real_, "failed",
        substr(conditionMessage(attr(ft, "condition")), 1L, 60L))
      next
    }
    rows[[length(rows) + 1L]] <- row_of(n, r, mth, NA_real_, ft[["secs"]],
                                        if (mth == "linear") "none" else "internal")
  }
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------------
## Part "tolerance"
## ---------------------------------------------------------------------------

## The convergence threshold of wafc(), against the seconds it costs, the
## penalty level it selects and the number of points of the path at which
## wafc_kkt() rejects the optimality conditions. Step E2.4 ran at 1e-10,
## three orders below the glmnet default, and the chat of 2026-09-21
## measured that dropping to 1e-7 is 28 times faster and leaves the tuning
## where it was. What that measurement did not read is the third column
## here: with the design unstandardized (decision D17), the relative
## criterion of the engine is read on a scale the penalty does not share,
## and at 1e-7 the fit stops satisfying the conditions decision D17 made
## the project's own check.
run_tolerance <- function(n, r) {
  seed <- seed0 + 50000L + 1000L * match(n, ns) + r
  dgp <- draw(n, seed)
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:5, n))
  rows <- list()
  for (J in c(3L, 4L, 5L)) {
    des <- wafc_design(dgp[["x"]], dgp[["u"]], J = J,
                       wavelet.table = wafc_timing_table)
    for (th in c(1e-7, 1e-8, 1e-9, 1e-10)) {
      tm <- timed({
        f <- wafc(design = des, y = dgp[["y"]], thresh = th)
        z <- wafc_cv_design(des, dgp[["y"]], f, foldid, function(e) e^2,
                            "lasso", thresh = th)
        list(f = f, z = z)
      })
      k <- wafc_kkt(tm[["value"]][["f"]])
      rows[[length(rows) + 1L]] <- data.frame(
        n = n, rep = r, J = J, thresh = th, secs = tm[["secs"]],
        lambda.min = tm[["value"]][["z"]][["lambda.min"]],
        cvm.min = tm[["value"]][["z"]][["cvm.min"]],
        kkt.bad = sum(!k[["ok"]]), path = nrow(k))
    }
  }
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------------
## The sweep
## ---------------------------------------------------------------------------

sweep_part <- function(fun, label) {
  out <- list()
  for (n in ns) {
    for (r in seq_len(R)) {
      cat(sprintf("  %s: n = %d, replicate %d of %d\n", label, n, r, R))
      utils::flush.console()
      out[[length(out) + 1L]] <- fun(n, r)
    }
  }
  do.call(rbind, out)
}

med <- function(v) if (all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE)

if ("methods" %in% parts) {
  cat("\n== methods ==\n")
  tab <- sweep_part(run_methods, "methods")
  agg <- do.call(rbind, lapply(split(tab, list(tab[["n"]], tab[["method"]]),
                                     drop = TRUE), function(d) {
    data.frame(n = d[["n"]][1L], method = d[["method"]][1L],
               secs.tune = med(d[["secs.tune"]]),
               secs.fit = med(d[["secs.fit"]]),
               secs.total = med(d[["secs.total"]]),
               split = d[["split"]][1L], note = d[["note"]][1L],
               stringsAsFactors = FALSE)
  }))
  agg <- agg[order(agg[["n"]], agg[["secs.total"]]), ]
  cat("\nmedian seconds over", R, "replicates:\n\n")
  print(format(agg[setdiff(names(agg), "note")], digits = 3),
        row.names = FALSE)
  ## The notes go in a legend of their own: printed as a column they wrap
  ## the table and make it unreadable.
  leg <- unique(agg[nzchar(agg[["note"]]), c("n", "method", "note")])
  if (nrow(leg) > 0L) {
    cat("\n")
    print(leg[order(leg[["n"]], leg[["method"]]), ], row.names = FALSE)
  }
}

if ("tolerance" %in% parts) {
  cat("\n== tolerance ==\n")
  tab <- sweep_part(run_tolerance, "tolerance")
  agg <- do.call(rbind, lapply(split(tab, list(tab[["n"]], tab[["J"]],
                                               tab[["thresh"]]), drop = TRUE),
                               function(d) {
    data.frame(n = d[["n"]][1L], J = d[["J"]][1L], thresh = d[["thresh"]][1L],
               secs = med(d[["secs"]]),
               lambda.min = med(d[["lambda.min"]]),
               cvm.min = med(d[["cvm.min"]]),
               kkt.bad = med(d[["kkt.bad"]]), path = med(d[["path"]]))
  }))
  agg <- agg[order(agg[["n"]], agg[["J"]], -agg[["thresh"]]), ]
  cat("\nmedian over", R, "replicates; kkt.bad counts the points of the path\n",
      "at which wafc_kkt() rejects the optimality conditions:\n\n")
  print(format(agg, digits = 4), row.names = FALSE)
}

cat("\ndone.\n")
