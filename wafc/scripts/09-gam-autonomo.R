## E2.5d -- the gam on its own: mgcv with a generous basis dimension fixed
## in advance, against the 'gam.matched' of step E2.5a and against the WAFC
## with the LASSO (question 33(d) of docs/ESTADO.md).
##
##     Rscript wafc/scripts/09-gam-autonomo.R [n_rep] [parts] [ncores] [ns] [cells] [ks]
##
## Why this exists. The 'gam.matched' of the pilot is mgcv at the basis
## dimension 2^J of the J that the cross-validation of 'wafc.lasso' selected
## on the same replicate, so it inherits the search of the WAFC and is not a
## method anyone could run without running the WAFC first. It answered the
## question of step E6.1a (was the gain one of dimension or one of basis?)
## and is the wrong reference for the exit criterion of step E2.5, which
## asks how the WAFC does against the spline a practitioner would fit. The
## usual practice is REML with a basis dimension large enough not to bind,
## and that is what this script fits: the model of wafc_fit_gam(), with
## select = TRUE, on the engine bam (fast REML, discretized covariates, the
## only one that makes these dimensions affordable, see the note on
## wafc_fit_gam()), at k = 64 and k = 128 per modulating covariate, each
## truncated at the number of distinct values minus one, the truncation of
## wafc_k_matched(). The two values are fixed here and not read from any
## fit of the WAFC; running both is what says whether 64 already binds.
##
## 'parts' is a comma separated subset of
##
##   fit     the sweep: for every cell, n and replicate of step E2.5a, the
##           'gam.matched' refitted and one 'gam.k<k>' per value of 'ks';
##           writes <WAFC_TAG>-fits.rds (one row per fit, the columns of the
##           pilot plus the basis dimension, the effective degrees of
##           freedom and the memory) and <WAFC_TAG>-edf.rds (one row per
##           fit and block);
##   report  the check that the refitted 'gam.matched' reproduces the
##           files of step E2.5a, and only if it passes, the tables. Reads
##           what 'fit' wrote, so it can be rerun alone.
##
## Default: 50 replicates, both parts, as many cores as the machine has
## minus two, the three sample sizes, the five cells, ks = 64,128. The cell
## "mixed" runs the 15 replicates it has in step E2.5a: a replicate is
## fitted only if E2.5a has a J for it, since 'gam.matched' cannot be
## refitted without one and the check is the point of refitting it.
##
## Three environment variables place the input and the output: WAFC_E25A is
## the directory of the files of step E2.5a (wafc/cache/e25a by default;
## every *competitors.rds in it is read), WAFC_OUT the directory written to
## (the working directory by default, created if missing) and WAFC_TAG the
## prefix of every file name ("e25d" by default).
##
## The draw. This script reproduces, without editing it, the draw of the
## part 'competitors' of wafc/scripts/04-pilot.R as step E2.5a ran it
## (commit 861582f): the same list of cells in the same order, the same
## sample sizes, the seed seed0 + 100000 cell + 1000 n + r from the
## positions in the full lists (step E2.4c), the same training sample, test
## sample, folds and grid. The functions below that carry the name of one
## of 04-pilot.R are copies of it and must not drift from it; the refitted
## 'gam.matched' is what checks that they have not, because it depends on
## every one of them and on nothing else that could change.
##
## Memory. Every fit records two numbers: 'mem_r', the peak of the R heap
## during the fit (gc(reset = TRUE) before, the "max used" of gc() after),
## which misses what mgcv allocates in C; and 'mem_rss', the peak resident
## set of the process during the fit, read from VmHWM after resetting it
## through /proc/self/clear_refs, which counts everything, the baseline of
## the process included. The second is the one comparable with the peak per
## process that step E2.5a reported; it is NA where the reset is not
## available.

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})
if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("This script needs the package 'mgcv'.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
arg_at <- function(i) if (length(args) >= i && nzchar(args[i])) args[i] else NULL
split_arg <- function(s) strsplit(s, ",", fixed = TRUE)[[1L]]
R <- if (is.null(arg_at(1L))) 50L else as.integer(arg_at(1L))
parts <- if (is.null(arg_at(2L))) c("fit", "report") else split_arg(arg_at(2L))
ncores <- if (is.null(arg_at(3L))) max(1L, parallel::detectCores() - 2L) else {
  as.integer(arg_at(3L))
}
in_dir <- Sys.getenv("WAFC_E25A", file.path("wafc", "cache", "e25a"))
out_dir <- Sys.getenv("WAFC_OUT", ".")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_tag <- Sys.getenv("WAFC_TAG", "e25d")
out_file <- function(what) file.path(out_dir, paste0(out_tag, "-", what, ".rds"))

## ---------------------------------------------------------------------------
## The draw of 04-pilot.R (commit 861582f), copied
## ---------------------------------------------------------------------------

seed0 <- 20260919L
n_test <- 2000L
n_grid <- 256L
ns_default <- c(250L, 500L, 1000L)
ns <- if (is.null(arg_at(4L))) ns_default else as.integer(split_arg(arg_at(4L)))
ns_all <- c(ns_default, setdiff(ns, ns_default))
n_index <- function(n) match(n, ns_all)

cells_all <- list(
  list(name = "smooth", scenario = "smooth", p = 3L, q = 2L, snr = 4),
  list(name = "inhomogeneous", scenario = "inhomogeneous", p = 3L, q = 2L,
       snr = 3),
  list(name = "mixed", scenario = "inhomogeneous", p = 4L, q = 4L, snr = 3,
       reps = 15L),
  list(name = "null", scenario = "null", p = 3L, q = 2L, sigma = 0.62),
  list(name = "uneven", scenario = "uneven", p = 3L, q = 2L, snr = 4)
)
cell_index <- function(cell) {
  match(cell[["name"]], vapply(cells_all, `[[`, "", "name"))
}
cells <- cells_all
if (!is.null(arg_at(5L))) {
  cells <- cells[vapply(cells, `[[`, "", "name") %in% split_arg(arg_at(5L))]
}

draw_cell <- function(cell, n, seed) {
  a <- list(n = n, p = cell[["p"]], q = cell[["q"]],
            scenario = cell[["scenario"]], seed = seed)
  if (!is.null(cell[["sigma"]])) a[["sigma"]] <- cell[["sigma"]]
  else a[["snr"]] <- cell[["snr"]]
  do.call(simulate_wafc, a)
}

test_for <- function(cell, dgp, seed) {
  test <- draw_cell(cell, n_test, seed + 500000L)
  test[["y"]] <- test[["f"]] + stats::rnorm(n_test, sd = dgp[["sigma"]])
  test[["sigma"]] <- dgp[["sigma"]]
  test
}

grid_of <- function(dgp, n_grid) {
  q <- dgp[["q"]]
  g <- matrix(0, n_grid, q)
  for (m in seq_len(q)) {
    rg <- range(dgp[["u"]][, m])
    g[, m] <- seq(rg[1L], rg[2L], length.out = n_grid)
  }
  g
}

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
    sigma = dgp[["sigma"]], time = fit[["time"]], error = NA_character_,
    stringsAsFactors = FALSE)
}

fail_row <- function(cell, n, r, method, dgp, active, err) {
  data.frame(
    cell = cell[["name"]], scenario = cell[["scenario"]], n = n, rep = r,
    method = method, rmse_f = NA_real_, rmse_y = NA_real_,
    mse_beta = NA_real_, ise = NA_real_, ise_active = NA_real_,
    ise_null = NA_real_, n_true = NA_integer_, n_false = NA_integer_,
    n_active = sum(active), n_block = length(active), J = NA_integer_,
    lambda = NA_real_, nzero = NA_integer_, sigma = dgp[["sigma"]],
    time = NA_real_,
    error = gsub("\\s+", " ", conditionMessage(attr(err, "condition"))),
    stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## What this script adds
## ---------------------------------------------------------------------------

ks <- if (is.null(arg_at(6L))) c(64L, 128L) else as.integer(split_arg(arg_at(6L)))
auto_methods <- paste0("gam.k", ks)
methods <- c("gam.matched", auto_methods)

## The rows of step E2.5a, and from them the J of 'wafc.lasso' on each
## replicate, which is the J 'gam.matched' carries in its own row.
read_e25a <- function(dir) {
  fs <- list.files(dir, pattern = "competitors\\.rds$", full.names = TRUE)
  if (length(fs) == 0L) {
    stop("No *competitors.rds in ", dir, "; set WAFC_E25A.", call. = FALSE)
  }
  do.call(rbind, lapply(fs, readRDS))
}
## The block LASSO with free levels of step E2.5b, read only for the
## factor of the smooth case, and only if its files are there (WAFC_E25B,
## wafc/cache/e25b by default): the table without it is complete.
read_e25b <- function(dir) {
  fs <- list.files(dir, pattern = "competitors\\.rds$", full.names = TRUE)
  if (length(fs) == 0L) return(NULL)
  d <- do.call(rbind, lapply(fs, readRDS))
  d[d[["method"]] == "klopp.free", ]
}
key_of <-function(d) paste(d[["cell"]], d[["n"]], d[["rep"]])

vmhwm_reset <- function() {
  isTRUE(tryCatch({
    writeLines("5", "/proc/self/clear_refs")
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE))
}
vmhwm_mb <- function() {
  s <- tryCatch(readLines("/proc/self/status"), error = function(e) character(0))
  h <- grep("^VmHWM:", s, value = TRUE)
  if (length(h) == 0L) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", h)) / 1024
}
gc_max_mb <- function(g) sum(g[, which(colnames(g) == "max used") + 1L])

## The basis dimension of the gam on its own: k per modulating covariate,
## truncated as wafc_k_matched() truncates 2^J, since a smooth cannot have
## more basis functions than its covariate has distinct values.
k_auto <- function(u, k) {
  nd <- vapply(seq_len(ncol(u)), function(m) length(unique(u[, m])), 0L)
  as.integer(pmax(3L, pmin(k, nd - 1L)))
}

run_job <- function(cell, n, r, J_lasso) {
  seed <- seed0 + 100000L * cell_index(cell) + 1000L * n_index(n) + r
  dgp <- draw_cell(cell, n, seed)
  test <- test_for(cell, dgp, seed)
  grid <- grid_of(dgp, n_grid)
  active <- nzchar(dgp[["structure"]])
  ## the folds of the pilot, drawn so that the stream of random numbers
  ## is the one of the pilot at this point; no gam uses them
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:10, n))
  rows <- list()
  edfs <- list()
  for (lab in methods) {
    k <- if (lab == "gam.matched") {
      wafc_k_matched(dgp[["u"]], J_lasso)
    } else {
      k_auto(dgp[["u"]], as.integer(sub("^gam\\.k", "", lab)))
    }
    ## the fit of the previous method is dropped before the baseline is
    ## taken, or its memory would be charged to this one
    f <- gh <- NULL
    gc(reset = TRUE)
    reset_ok <- vmhwm_reset()
    f <- try(wafc_competitor("gam", dgp[["x"]], dgp[["u"]], dgp[["y"]],
                             k = k, engine = "bam"), silent = TRUE)
    mem_r <- gc_max_mb(gc())
    mem_rss <- if (reset_ok) vmhwm_mb() else NA_real_
    if (inherits(f, "try-error")) {
      rows[[length(rows) + 1L]] <- fail_row(cell, n, r, lab, dgp, active, f)
      next
    }
    gh <- wafc_grid_components(f, grid)
    rw <- one_row(
      cell, n, r, lab, dgp, test, grid,
      list(f_test = predict(f, test[["x"]], test[["u"]]),
           beta_test = f[["beta"]](test[["u"]]), time = f[["time"]]),
      active, f[["blocks"]], gh,
      list(J = if (lab == "gam.matched") J_lasso else NA_integer_))
    ## The effective degrees of freedom of each smooth against the number
    ## of coefficients it has, which is what says whether k binds: a smooth
    ## whose edf approaches its coefficient count wanted a larger basis.
    edf <- f[["extra"]][["edf"]]
    ncoef <- vapply(f[["fit"]][["smooth"]],
                    function(sm) sm[["last.para"]] - sm[["first.para"]] + 1, 0)
    ncoef <- matrix(ncoef, nrow(edf), ncol(edf), byrow = TRUE)
    rw[["k"]] <- max(k)
    rw[["edf_total"]] <- sum(edf)
    rw[["edf_active"]] <- sum(edf[active])
    rw[["edf_max"]] <- max(edf)
    rw[["edf_frac_max"]] <- max(edf / ncoef)
    rw[["mem_r"]] <- mem_r
    rw[["mem_rss"]] <- mem_rss
    rows[[length(rows) + 1L]] <- rw
    edfs[[length(edfs) + 1L]] <- data.frame(
      cell = cell[["name"]], n = n, rep = r, method = lab,
      l = as.vector(row(edf)), m = as.vector(col(edf)),
      active = as.vector(active), k = k[as.vector(col(edf))],
      ncoef = as.vector(ncoef), edf = as.vector(edf),
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, lapply(rows, function(d) {
    for (v in c("k", "edf_total", "edf_active", "edf_max", "edf_frac_max",
                "mem_r", "mem_rss")) {
      if (is.null(d[[v]])) d[[v]] <- NA_real_
    }
    d
  }))
  list(rows = out, edf = do.call(rbind, edfs))
}

## ---------------------------------------------------------------------------
## Part "fit"
## ---------------------------------------------------------------------------

if ("fit" %in% parts) {
  e25a <- read_e25a(in_dir)
  Jtab <- e25a[e25a[["method"]] == "gam.matched" & !is.na(e25a[["J"]]),
               c("cell", "n", "rep", "J")]
  Jkey <- stats::setNames(Jtab[["J"]], key_of(Jtab))
  jobs <- list()
  for (cell in cells) {
    nrep <- if (is.null(cell[["reps"]])) R else min(R, cell[["reps"]])
    for (n in ns) {
      for (r in seq_len(nrep)) {
        J <- Jkey[paste(cell[["name"]], n, r)]
        if (is.na(J)) next
        jobs[[length(jobs) + 1L]] <- list(cell = cell, n = n, r = r,
                                          J = as.integer(J))
      }
    }
  }
  ## The longest jobs first, so that the tail of the sweep is not one fit
  ## of 'gam.matched' at J = 8 started last: mclapply without prescheduling
  ## hands the jobs out in this order.
  cost <- vapply(jobs, function(j) {
    j[["cell"]][["q"]] * j[["n"]] * 2^max(j[["J"]], log2(max(ks)))
  }, 0)
  jobs <- jobs[order(-cost)]
  cat(sprintf("\n== fit: %d jobs, methods %s, on %d core(s) ==\n",
              length(jobs), paste(methods, collapse = ", "), ncores))
  t0 <- proc.time()[["elapsed"]]
  res <- parallel::mclapply(jobs, function(j) {
    tryCatch(run_job(j[["cell"]], j[["n"]], j[["r"]], j[["J"]]),
             error = function(e) {
               message("job failed: ", conditionMessage(e))
               NULL
             })
  }, mc.cores = ncores, mc.preschedule = FALSE)
  cat(sprintf("   %.1f s\n", proc.time()[["elapsed"]] - t0))
  ok <- vapply(res, function(z) is.list(z) && is.data.frame(z[["rows"]]), TRUE)
  if (any(!ok)) {
    cat(sprintf("   %d of %d job(s) failed as a whole and left no row\n",
                sum(!ok), length(jobs)))
  }
  fits <- do.call(rbind, lapply(res[ok], `[[`, "rows"))
  edf <- do.call(rbind, lapply(res[ok], `[[`, "edf"))
  saveRDS(fits, out_file("fits"))
  saveRDS(edf, out_file("edf"))
  if (any(!is.na(fits[["error"]]))) {
    bad <- fits[!is.na(fits[["error"]]), c("method", "error")]
    cat("   failed fits, kept as rows with no numbers:\n")
    for (m in unique(bad[["method"]])) {
      cat(sprintf("     %-12s %3d: %s\n", m, sum(bad[["method"]] == m),
                  substr(bad[["error"]][bad[["method"]] == m][1L], 1L, 100L)))
    }
  }
}

## ---------------------------------------------------------------------------
## Part "report"
## ---------------------------------------------------------------------------

med <- function(d, vars, by) {
  a <- stats::aggregate(d[vars], d[by], FUN = stats::median, na.rm = TRUE)
  a[do.call(order, unname(as.list(a[by]))), ]
}

if ("report" %in% parts) {
  e25a <- read_e25a(in_dir)
  fits <- readRDS(out_file("fits"))
  edf <- readRDS(out_file("edf"))

  ## The check, before any table. Every number of the row of 'gam.matched'
  ## except the time has to come back: the fit is deterministic (bam with
  ## one thread, and the reference BLAS), so the tolerance is a rounding
  ## one, and a replicate that does not come back means the draw has
  ## drifted from the pilot.
  num <- c("rmse_f", "rmse_y", "mse_beta", "ise", "ise_active", "ise_null",
           "sigma")
  int <- c("n_true", "n_false", "n_active", "n_block", "J")
  ref <- e25a[e25a[["method"]] == "gam.matched", c("cell", "n", "rep", num, int)]
  new <- fits[fits[["method"]] == "gam.matched", c("cell", "n", "rep", num, int)]
  both <- merge(new, ref, by = c("cell", "n", "rep"), suffixes = c("", ".a"))
  cat(sprintf("\n== check: gam.matched refitted against step E2.5a, %d of %d rows ==\n",
              nrow(both), nrow(new)))
  worst <- vapply(num, function(v) {
    max(abs(both[[v]] - both[[paste0(v, ".a")]]) /
          pmax(abs(both[[paste0(v, ".a")]]), 1e-12), na.rm = TRUE)
  }, 0)
  print(signif(worst, 3))
  same_int <- vapply(int, function(v) {
    all(both[[v]] == both[[paste0(v, ".a")]], na.rm = TRUE)
  }, TRUE)
  cat("integer columns identical:", paste(int, same_int, collapse = ", "), "\n")
  if (nrow(both) != nrow(new) || nrow(new) == 0L || any(worst > 1e-8) ||
      !all(same_int) || any(!is.na(fits[["error"]][fits[["method"]] == "gam.matched"]))) {
    stop("gam.matched does not reproduce step E2.5a: the draw has drifted, ",
         "and no table is printed.", call. = FALSE)
  }
  cat("OK: gam.matched reproduces step E2.5a row by row.\n")

  ## The tables read the WAFC with the LASSO, the gam at k = 10 and the
  ## block LASSO from step E2.5a, and the block LASSO with free levels from
  ## step E2.5b when its files are there, on the replicates this run fitted.
  keys <- unique(key_of(fits))
  old <- e25a[e25a[["method"]] %in% c("wafc.lasso", "gam", "klopp") &
                key_of(e25a) %in% keys, ]
  e25b <- read_e25b(Sys.getenv("WAFC_E25B", file.path("wafc", "cache", "e25b")))
  if (!is.null(e25b)) {
    cat(sprintf("klopp.free read from step E2.5b: %d rows on these replicates\n",
                sum(key_of(e25b) %in% keys)))
    old <- rbind(old, e25b[key_of(e25b) %in% keys, names(old)])
  }
  cols <- intersect(names(old), names(fits))
  all_rows <- rbind(old[cols],
                    fits[fits[["method"]] != "gam.matched", cols],
                    fits[fits[["method"]] == "gam.matched", cols])
  refs <- intersect(c("wafc.lasso", "klopp", "klopp.free"), old[["method"]])
  lev <- c(refs, "gam", "gam.matched", auto_methods)
  all_rows[["method"]] <- factor(all_rows[["method"]], levels = lev)
  cells_order <- c("smooth", "uneven", "inhomogeneous", "mixed", "null")
  all_rows[["cell"]] <- factor(all_rows[["cell"]], levels = cells_order)
  fits[["cell"]] <- factor(fits[["cell"]], levels = cells_order)
  saveRDS(all_rows, out_file("joined"))

  cat("\nmedians by cell, n and method (rmse_f out of sample, ISE, blocks",
      "kept of the active and of the null ones, s):\n")
  print(med(all_rows, c("rmse_f", "ise", "ise_active", "ise_null", "n_true",
                        "n_false", "time"), c("cell", "n", "method")),
        row.names = FALSE, digits = 3)

  ## Ratios within replicate: the median of the ratios, the fraction of
  ## replicates in which the gam on its own is the better of the two, and
  ## the ratio the other way round for the WAFC, which is the "factor" of
  ## the exit criterion in the smooth case (question 33(c)).
  paired <- function(against) {
    b <- all_rows[all_rows[["method"]] == against,
                  c("cell", "n", "rep", "rmse_f", "ise")]
    names(b)[4:5] <- c("rmse0", "ise0")
    d <- merge(all_rows[all_rows[["method"]] %in% auto_methods, ], b,
               by = c("cell", "n", "rep"))
    d[["r_rmse"]] <- d[["rmse_f"]] / d[["rmse0"]]
    d[["r_ise"]] <- ifelse(d[["ise0"]] > 0, d[["ise"]] / d[["ise0"]], NA_real_)
    d[["win_rmse"]] <- d[["rmse_f"]] < d[["rmse0"]]
    d[["win_ise"]] <- d[["ise"]] < d[["ise0"]]
    d[["method"]] <- droplevels(d[["method"]])
    a <- stats::aggregate(cbind(r_rmse, r_ise) ~ cell + n + method, d,
                          FUN = stats::median, na.action = stats::na.pass,
                          na.rm = TRUE)
    w <- stats::aggregate(cbind(win_rmse, win_ise) ~ cell + n + method, d,
                          FUN = mean)
    out <- merge(a, w, by = c("cell", "n", "method"))
    out[order(out[["cell"]], out[["n"]], out[["method"]]), ]
  }
  for (against in c("gam.matched", refs, "gam")) {
    cat(sprintf(paste("\nthe gam on its own against %s (median ratio within",
                      "replicate, < 1 is the gam on its own better; fraction",
                      "of replicates it wins):\n"), against))
    print(paired(against), row.names = FALSE, digits = 3)
  }

  ## The factor of the exit criterion in the smooth case (question 33(c)),
  ## ISE(method) / ISE(gam), for each method that could be the WAFC and
  ## each version of the gam, and the fraction of replicates the gam wins.
  ## A replicate in which the gam has zero ISE has no factor and is dropped.
  gams <- c("gam.matched", auto_methods)
  for (rf in refs) {
    b <- all_rows[all_rows[["method"]] == rf, c("cell", "n", "rep", "ise")]
    names(b)[4L] <- "ise_w"
    d <- merge(all_rows[all_rows[["method"]] %in% gams, ], b,
               by = c("cell", "n", "rep"))
    d <- d[which(d[["ise"]] > 0), ]
    d[["factor"]] <- d[["ise_w"]] / d[["ise"]]
    d[["gam_wins"]] <- d[["ise"]] < d[["ise_w"]]
    d[["method"]] <- droplevels(d[["method"]])
    cat(sprintf(paste("\nthe factor of the smooth case, ISE(%s) / ISE(gam),",
                      "median within replicate, and the fraction of replicates",
                      "the gam wins:\n"), rf))
    a <- merge(med(d, "factor", c("cell", "n", "method")),
               stats::aggregate(gam_wins ~ cell + n + method, d, FUN = mean),
               by = c("cell", "n", "method"))
    a <- a[order(a[["method"]]), ]
    w <- stats::reshape(a, idvar = c("cell", "n"), timevar = "method",
                        direction = "wide")
    print(w[order(w[["cell"]], w[["n"]]), ], row.names = FALSE, digits = 3)
  }

  cat("\nthe two dimensions against each other, and how often the gam on its",
      "own is the same fit as gam.matched (J = log2 k):\n")
  if (length(auto_methods) >= 2L) {
    a1 <- fits[fits[["method"]] == auto_methods[1L], c("cell", "n", "rep", "rmse_f", "ise")]
    a2 <- fits[fits[["method"]] == auto_methods[2L], c("cell", "n", "rep", "rmse_f", "ise")]
    d <- merge(a2, a1, by = c("cell", "n", "rep"), suffixes = c("", ".1"))
    d[["r_rmse"]] <- d[["rmse_f"]] / d[["rmse_f.1"]]
    d[["r_ise"]] <- ifelse(d[["ise.1"]] > 0, d[["ise"]] / d[["ise.1"]], NA_real_)
    cat(sprintf("  %s / %s:\n", auto_methods[2L], auto_methods[1L]))
    print(stats::aggregate(cbind(r_rmse, r_ise) ~ cell + n, d,
                           FUN = function(v) stats::median(v, na.rm = TRUE),
                           na.action = stats::na.pass),
          row.names = FALSE, digits = 4)
  }
  jm <- fits[fits[["method"]] == "gam.matched", c("cell", "n", "rep", "J")]
  for (k in ks) {
    cat(sprintf("  fraction of replicates with J = %d, where gam.k%d is gam.matched:\n",
                as.integer(log2(k)), k))
    print(stats::aggregate(same ~ cell + n,
                           transform(jm, same = J == log2(k)), FUN = mean),
          row.names = FALSE, digits = 3)
  }

  cat("\neffective degrees of freedom: total, of the active blocks, the",
      "largest smooth, and the largest edf / number of coefficients of a",
      "smooth (near 1 is a basis that binds):\n")
  print(med(fits, c("k", "edf_total", "edf_active", "edf_max", "edf_frac_max"),
            c("cell", "n", "method")),
        row.names = FALSE, digits = 3)
  cat("\nlargest edf / coefficients over every smooth and replicate:\n")
  edf[["frac"]] <- edf[["edf"]] / edf[["ncoef"]]
  print(stats::aggregate(frac ~ cell + n + method, edf, FUN = max),
        row.names = FALSE, digits = 3)

  cat("\ntime (s) and memory (MB: peak of the R heap, peak resident set of",
      "the process) per fit, median and maximum:\n")
  tm <- stats::aggregate(cbind(time, mem_r, mem_rss) ~ cell + n + method, fits,
                         FUN = function(v) c(med = stats::median(v), max = max(v)),
                         na.action = stats::na.pass)
  print(tm[order(tm[["cell"]], tm[["n"]], tm[["method"]]), ], row.names = FALSE,
        digits = 3)
  cat("\nprocessor hours by method:\n")
  print(round(tapply(fits[["time"]], fits[["method"]], sum, na.rm = TRUE) / 3600, 2))
}

cat("\nOK\n")
