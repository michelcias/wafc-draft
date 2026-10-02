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
##                eps = a (L-1) 2^(-J); it is what fixed wafc_eps_periodic
##                at 0, decision D35 (it was provisional at 0.05);
##   j1           whether the grid of cv.wafc() should start at J = 1, now
##                that the margin no longer depends on J (step E2.1b);
##   jgrid        what the grid of cv.wafc() costs at the other end. The
##                part 'lambda' found the oracle of the grid sitting on the
##                top of it in every replicate of every cell that has
##                components to resolve, and E1.3 measured that the rate
##                2^{-Js'} of bumps only appears from J = 9. This part
##                cross-validates on 2:Jmax and on a grid twice as deep,
##                and reads what the extra levels buy and cost.
##
## 'ns' and 'cells' are comma separated and restrict the sweep, which is
## how a part is rerun on one cell without rerunning the rest. The seed of a
## replicate is fixed by the position of its cell and of its n in the full
## lists below, and not in the restricted ones, so a restricted run draws
## the same data as the full run does for those cells.
##
## Two environment variables place the output: WAFC_OUT is the directory
## (the working directory by default; created if missing, so that a long
## run does not end on a failed saveRDS), and WAFC_TAG the prefix of every
## file name ("e24" by default, the name of the files of step E2.4). The
## rerun of step E2.5a sets both, and so leaves the files of E2.4 alone.
##
## Two more restrict the part 'competitors' without touching what it draws
## (step E2.5b). WAFC_METHODS is a comma separated subset of the methods
## below: the others are not run, and the data, the folds, the test sample
## and the random stream each method starts from are the ones of the full
## run, so a method added later can be run alone and joined to an earlier
## table, and the rows of a method run twice coincide. 'gam.matched' needs
## the J of 'wafc.lasso' and cannot be run without it. WAFC_REPS_MIXED
## replaces the 15 replicates of the cell "mixed"; the seeds are numbered
## by replicate, so its first 15 are the ones of a run with 15.
##
## Default: 50 replicates, every part, as many cores as the machine has
## minus two, and 15 replicates in the cell "mixed". That budget was fixed
## in step E2.4, when the grid of J stopped at ceiling(log2(n)/2) and one
## replicate of that cell cost 27 minutes of processor over the three
## sample sizes, 23 processor-hours for 50 against 3 for the other three
## cells together. With the grid 2:8 (decision D34) that is no longer so.
## Measured in step E2.5a, part 'competitors': the four cells with q = 2
## took 49.1 processor-hours for 600 replicate-by-n jobs (4 h 08 min on 12
## cores), and the cell "mixed" 3.75 for 45 (1 h on 4 cores), about five
## minutes of processor per job in both, because the cells with q = 2 now
## go deep in J and the WAFC does not choose J = 8 in the mixed one. Of the
## 52 hours, the sparse group LASSO took 22.2 (42%) and the B-spline group
## LASSO 10.9; the WAFC with the LASSO takes 7 to 26 s per fit, search
## included. What the mixed cell still has is memory: 'gam.matched' at
## J = 8 there with n = 1000 ran past 15 minutes and 3 to 4 GB per process
## (step E2.4c), against a peak of 1.35 GB in E2.5a, where the WAFC chose
## J <= 7. The 15 stay the default so that a rerun reproduces E2.5a. With
## WAFC_REPS_MIXED = 50 and every method, step E2.5b took 15.7 processor-
## hours (2 h 02 min on 8 cores) and peaked at 1.49 GB per process; the
## WAFC chose J = 8 in none of the 150 jobs. Step E2.5c ran 'klopp.unit'
## and 'klopp.merged' alone, five cells with 50 replicates each, in 9.7
## processor-hours (1 h 14 min on 8 cores, alone on the machine), with a
## peak of 0.94 GB per process. Step E2.5e ran 'klopp.balanced' alone, the
## same 750 jobs, in 4.8 processor-hours (36 min on 8 cores), with a peak of
## 0.67 GB per process. Step E2.5f ran 'klopp.freecoarse' alone, the same
## 750 jobs, in 6.0 processor-hours (46 min on 8 cores), with a peak of
## 0.67 GB per process.
##
## Step E2.5g adds, in the part 'competitors', the estimator of Corollary 8
## (estimation followed by a threshold on the blocks, wafc_threshold() of
## wafc/R/threshold.R) on two of the methods, 'wafc.lasso' and
## 'klopp.balanced': from the one fit of the replicate, one row per rule of
## the threshold ("max", "cv", "oracle"), without refit and with the two
## least squares refits ("+ls" on the support of the blocks kept, "+lsb" on
## their every column), labelled "<method>+<rule>[+ls|+lsb]". The rows of the
## methods themselves are not touched, so they still coincide with the ones
## of the earlier runs. The time of a thresholded row is the time of the
## fit plus the time of the threshold (and of the fold fits, for "cv"). Two
## side tables go beside the rows: '<tag>-thr-norms.rds', the norm of every
## block of the fit before the threshold with its truth, from which the
## curve of recovery against t is read, and '<tag>-thr-t.rds', the t and
## the fraction c = t / max norm each thresholded row used.
##
## Step E2.5h adds four methods to the part 'competitors'. 'gam.reml',
## 'gam.gcv' and 'gam.cv' are mgcv with the basis dimension chosen on one
## grid, 5, 10, 20, 40 and 80, common to the smooths (decision D41), by the
## REML of the fit (engine bam with discrete = TRUE, as 'gam.matched'), by
## its GCV (smoothing parameters by GCV as well; bam discretizes only under
## REML, so this one runs on the engine named in gam_gcv_engine below), or
## by the squared error on the folds of the replicate, the ones cv.wafc()
## uses (REML and bam inside each fold), with the time of the whole search
## in the row. D41 leaves 'gam.gcv' out of the cell "mixed": the run does
## not ask for it there (WAFC_METHODS), and a fit at k = 80 took more than
## an hour there. 'wafc.gcv' is the WAFC with
## (J, lambda) by generalized cross-validation over the grid 2:8
## (wafc_tune(rule = "gcv"), with its guard: points with df >= n/2 left
## out), and 'wafc.gcv+max' the threshold of step E2.5g on it, with
## c = 0.15 and no refit. Two side tables go beside the rows:
## '<tag>-gam-k.rds', one row per candidate of k fitted, with its score,
## the score the engine reports, its edf, its number of coefficients and its
## seconds, and the flags of the k chosen and of the top of the grid; and
## '<tag>-gcv-guard.rds', how many points of the paths of 'wafc.gcv' the
## guard left out and whether it decided the pair.
##
## The two parts that answer a question about the code rather than about
## the method, 'margin' and 'j1', are capped at 20 replicates:
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
if (identical(parts, "all")) {
  parts <- c("competitors", "lambda", "margin", "j1", "jgrid")
}
ncores <- if (length(args) >= 3L && nzchar(args[3L])) as.integer(args[3L]) else {
  max(1L, parallel::detectCores() - 2L)
}
out_dir <- Sys.getenv("WAFC_OUT", ".")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_tag <- Sys.getenv("WAFC_TAG", "e24")
out_file <- function(what) file.path(out_dir, paste0(out_tag, "-", what, ".rds"))
seed0 <- 20260919L

## Decision D31: in repeated numerical work the basis is fixed and evaluated
## by table, built once, and not by the Daubechies-Lagarias algorithm at
## every fit. The rule use.table = "auto" of wafc_design() only fires at
## n q >= 2000 L, that is n q >= 16000 with the filter of size 8, so it
## never fires at the sample sizes of this pilot: leaving it to decide
## means paying the exact evaluation over the whole sweep. The table is
## built here and passed to every design the script builds from scratch;
## a design built with spec = takes it from the spec. The interpolation
## error is 3.1e-06 (step E2.1), which is why the exception of D31, the
## checks of derivations/check/, does not apply here.
wafc_pilot_table <- WaveBased::wtable(family = "Daublets", filter.size = 8L,
                                      prec.wavelet = 30L, check = FALSE)
n_test <- 2000L
n_grid <- 256L
ns_default <- c(250L, 500L, 1000L)
ns <- if (length(args) >= 4L && nzchar(args[4L])) {
  as.integer(strsplit(args[4L], ",", fixed = TRUE)[[1L]])
} else ns_default
## the position of n in the seeds: the default sizes keep theirs, and a size
## outside them is numbered after them
ns_all <- c(ns_default, setdiff(ns, ns_default))
n_index <- function(n) match(n, ns_all)
R_small <- min(R, 20L)

## The effective regularity s' each scenario declares (decision D27) is read
## from wafc_sprime() of dgp.R, in the regime "periodic" of the theory
## (D26). This script used to carry its own vector under that name, from
## before step E2.4b moved the number into dgp.R; the vector shadowed the
## function and had no entry for "uneven", so a cell of that scenario would
## have failed in the part 'lambda' and been dropped.

## The five cells of the pilot. The third is the "mixture with half the
## components null" of plano-projeto.md E2.4, as far as dgp.R reaches:
## wafc_scenario() activates three blocks whatever p and q are, so p = 4 and
## q = 4 give 3 active blocks of 16, which is sparser than a half and not
## looser. Making the fraction a parameter is a change to dgp.R and is in
## the handoff.
##
## The fifth, "uneven", is the smooth scenario of uneven curvature decision
## D30 asks for (step E2.4b): C-infinity components whose scale varies along
## the domain, inside the hypothesis of the spline competitors, and the cell
## on which the factor of the exit criterion of E2.5 in the smooth case is
## to be fixed. It is appended last so that the four cells of E2.4 keep
## their position, and with it their seeds.
cells_all <- list(
  list(name = "smooth", scenario = "smooth", p = 3L, q = 2L, snr = 4),
  list(name = "inhomogeneous", scenario = "inhomogeneous", p = 3L, q = 2L,
       snr = 3),
  list(name = "mixed", scenario = "inhomogeneous", p = 4L, q = 4L, snr = 3,
       reps = 15L),
  list(name = "null", scenario = "null", p = 3L, q = 2L, sigma = 0.62),
  list(name = "uneven", scenario = "uneven", p = 3L, q = 2L, snr = 4)
)
reps_mixed <- Sys.getenv("WAFC_REPS_MIXED", "")
if (nzchar(reps_mixed)) {
  reps_mixed <- as.integer(reps_mixed)
  if (is.na(reps_mixed) || reps_mixed < 1L) {
    stop("WAFC_REPS_MIXED must be a positive integer.", call. = FALSE)
  }
  im <- match("mixed", vapply(cells_all, `[[`, "", "name"))
  cells_all[[im]][["reps"]] <- reps_mixed
}
cell_index <- function(cell) {
  match(cell[["name"]], vapply(cells_all, `[[`, "", "name"))
}

cells <- cells_all
if (length(args) >= 5L && nzchar(args[5L])) {
  want <- strsplit(args[5L], ",", fixed = TRUE)[[1L]]
  cells <- cells[vapply(cells, `[[`, "", "name") %in% want]
}

## 'gam' is mgcv at the default k = 10, the column of the table of E2.4,
## kept to read the rerun against it; 'gam.matched' is the same model at the
## dimension of the WAFC (step E6.1a showed that the difference between the
## two can be the whole verdict). 'klopp.free' is the block LASSO of Klopp
## and Pensky with the level terms left unpenalized, as decision D3 leaves
## them in the WAFC (step E2.5b): step E2.5a found 'klopp' ahead of the
## WAFC wherever there are components, and behind it in the null scenario,
## where penalizing the levels shrinks them. 'klopp.unit' and
## 'klopp.merged' are 'klopp.free' in the two other forms of open
## question 34 of docs/ESTADO.md (step E2.5c): weight one on every chunk,
## which is their norm (3.1) and the one the theory of step E1.11 covers,
## and the coarse levels of each block merged into one chunk with the
## sqrt(|G|) of grpreg kept. 'klopp.balanced' is the fourth form, the one
## step E2.5c proposed and step E2.5e measures: the coarse levels merged
## and the short piece at the end of each finer level absorbed into the
## chunk before it, so that every chunk of a finer level has between b_n
## and 2 b_n - 1 columns and the sqrt(|G|) of grpreg stay within a bounded
## ratio of one another. 'klopp.freecoarse' is the fifth (step E2.5f, open
## question 37): the balanced chunks on the finer levels, and the coarse
## levels of each block, the ones the balanced form penalizes together,
## left unpenalized with the level terms. No block of it is ever zero, so
## its n_true and n_false count the blocks whose penalized part is nonzero
## (extra$blocks.fine), and its nzero is the number of nonzero penalized
## coefficients, as for the WAFC; at J <= 3 nothing is penalized, the fit is
## least squares and its lambda is NA. The order is the order of the
## tables; a method added later goes after the one it varies, which changes
## no seed, since the random stream of every method is the one it would
## start from alone (run_competitors()).
methods <- c("wafc.lasso", "wafc.sglasso", "wafc.gcv", "gam", "gam.matched",
             "gam.reml", "gam.gcv", "gam.cv", "bsgl",
             "klopp", "klopp.free", "klopp.unit", "klopp.merged",
             "klopp.balanced", "klopp.freecoarse", "aspline", "vcbart",
             "linear", "oracle")
## The engine of 'gam.gcv' (step E2.5h): bam discretizes the covariates only
## under REML, and with method = "GCV.Cp" it warns and fits without the
## discretization, so the choice is between bam and gam without it. bam by
## default: measured on one replicate, the two cost the same with q = 2 at
## n = 250 (419 s and 394 s over the grid up to 120), and bam is eight
## times faster in the cell "mixed" at n = 1000 and k = 20 (55 s against
## 467 s). WAFC_GAM_GCV_ENGINE = "gam" switches.
gam_gcv_engine <- Sys.getenv("WAFC_GAM_GCV_ENGINE", "bam")
## Step E2.5g: the thresholded rows of 'wafc.lasso' and 'klopp.balanced'.
thr_bases <- c("wafc.lasso", "klopp.balanced")
thr_rules <- c("max", "cv", "oracle")
thr_refits <- c(none = "", support = "+ls", block = "+lsb")
thr_labels <- function(base) {
  unlist(lapply(names(thr_refits), function(rf) {
    paste0(base, "+", thr_rules, thr_refits[[rf]])
  }))
}
method_levels <- unlist(lapply(methods, function(m) {
  c(m, if (m %in% thr_bases) thr_labels(m),
    if (m == "wafc.gcv") "wafc.gcv+max")
}))
run_methods <- methods
if (nzchar(Sys.getenv("WAFC_METHODS", ""))) {
  run_methods <- strsplit(Sys.getenv("WAFC_METHODS"), ",", fixed = TRUE)[[1L]]
  bad <- setdiff(run_methods, methods)
  if (length(bad) > 0L) {
    stop("Unknown method(s) in WAFC_METHODS: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  if ("gam.matched" %in% run_methods && !("wafc.lasso" %in% run_methods)) {
    stop("'gam.matched' needs the J of 'wafc.lasso'; add it to WAFC_METHODS.",
         call. = FALSE)
  }
  run_methods <- methods[methods %in% run_methods]
}

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
    sigma = dgp[["sigma"]], time = fit[["time"]], error = NA_character_,
    stringsAsFactors = FALSE)
}

## A method that fails leaves a row with its error message and no numbers,
## instead of no row at all: a column that disappears from the table is
## read as a method that was never run, which is how the failure of vcbart
## went unnoticed after step E2.4b.
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

## The thresholded rows of one fit (step E2.5g): every rule and refit of
## thr_rules and thr_refits, the folds refitted once for the three "cv"
## rows. Besides the rows, two side tables: the norm of every block of the
## fit, with its truth, and the t of every row.
thr_rows <- function(cell, n, r, base, obj, base_time, dgp, test, grid,
                     active, foldid) {
  rows <- list()
  tt <- list()
  truth <- list(x = test[["x"]], u = test[["u"]], f = test[["f"]])
  t0 <- proc.time()[["elapsed"]]
  ff <- try(wafc_threshold_folds(obj, y = dgp[["y"]], foldid = foldid),
            silent = TRUE)
  ff_time <- proc.time()[["elapsed"]] - t0
  if (inherits(ff, "try-error")) ff <- NULL
  nrm <- NULL
  for (rf in names(thr_refits)) {
    for (rule in thr_rules) {
      lab <- paste0(base, "+", rule, thr_refits[[rf]])
      th <- try(wafc_threshold(obj, rule = rule, refit = rf, y = dgp[["y"]],
                               foldid = foldid, truth = truth,
                               fold.fits = ff), silent = TRUE)
      if (inherits(th, "try-error")) {
        rows[[length(rows) + 1L]] <- fail_row(cell, n, r, lab, dgp, active,
                                              th)
        next
      }
      ex <- th[["extra"]]
      if (is.null(nrm)) nrm <- ex[["norm"]]
      rows[[length(rows) + 1L]] <- one_row(
        cell, n, r, lab, dgp, test, grid,
        list(f_test = predict(th, test[["x"]], test[["u"]]),
             beta_test = th[["beta"]](test[["u"]]),
             time = base_time + th[["time"]] +
               if (rule == "cv") ff_time else 0),
        active, th[["blocks"]], wafc_grid_components(th, grid),
        list(J = ex[["J"]], lambda = ex[["lambda"]], nzero = ex[["nzero"]]))
      tt[[length(tt) + 1L]] <- data.frame(
        cell = cell[["name"]], n = n, rep = r, method = lab, t = ex[["t"]],
        c = ex[["c"]], stringsAsFactors = FALSE)
    }
  }
  side <- list(
    "thr-t" = if (length(tt)) do.call(rbind, tt) else NULL,
    "thr-norms" = if (is.null(nrm)) NULL else data.frame(
      cell = cell[["name"]], n = n, rep = r, method = base,
      l = as.vector(row(nrm)), m = as.vector(col(nrm)),
      active = active, norm = as.vector(nrm),
      stringsAsFactors = FALSE))
  list(rows = rows, side = side)
}

## Appends the side tables of a job to the ones it already has.
side_add <- function(side, more) {
  for (nm in names(more)) side[[nm]] <- rbind(side[[nm]], more[[nm]])
  side
}

## ---------------------------------------------------------------------------
## Part "competitors"
## ---------------------------------------------------------------------------

run_competitors <- function(cell, n, r) {
  seed <- seed0 + 100000L * cell_index(cell) + 1000L * n_index(n) + r
  dgp <- draw_cell(cell, n, seed)
  test <- test_for(cell, dgp, seed)
  p <- cell[["p"]]
  grid <- grid_of(dgp, n_grid)
  active <- nzchar(dgp[["structure"]])
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:10, n))
  rows <- list()
  side <- list()
  ## Every method starts from the random stream left by the folds, so that
  ## what one method draws (vcbart is the one that draws) does not depend on
  ## which methods ran before it: that is what lets WAFC_METHODS restrict the
  ## run and a method be inserted in the list without moving anything.
  ## Nothing before vcbart draws, so its rows are the ones of E2.5a.
  rng <- get(".Random.seed", envir = globalenv())
  from_start <- function() assign(".Random.seed", rng, envir = globalenv())

  ## The two WAFC variants, tuned by cross-validation over (J, lambda),
  ## which decision D20 made the default rule. The resolution the LASSO
  ## selects is kept for 'gam.matched' below.
  J_lasso <- NULL
  for (pen in c("lasso", "sglasso")) {
    if (!(paste0("wafc.", pen) %in% run_methods)) next
    from_start()
    t0 <- proc.time()[["elapsed"]]
    cv <- try(cv.wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], penalty = pen,
                      foldid = foldid, wavelet.table = wafc_pilot_table),
              silent = TRUE)
    if (inherits(cv, "try-error")) {
      rows[[length(rows) + 1L]] <- fail_row(cell, n, r, paste0("wafc.", pen),
                                            dgp, active, cv)
      next
    }
    el <- proc.time()[["elapsed"]] - t0
    if (pen == "lasso") J_lasso <- cv[["J.min"]]
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
    if (paste0("wafc.", pen) %in% thr_bases) {
      tr <- thr_rows(cell, n, r, paste0("wafc.", pen), cv, el, dgp, test,
                     grid, active, foldid)
      rows <- c(rows, tr[["rows"]])
      side <- side_add(side, tr[["side"]])
    }
  }

  ## 'wafc.gcv' (step E2.5h): the pair (J, lambda) by generalized
  ## cross-validation over the grid 2:8, with the guard of wafc_gcv(), and
  ## its threshold by the rule "max" of step E2.5g at c = 0.15, no refit. The
  ## time of the row is the time of the search; the one of the thresholded
  ## row adds the threshold.
  if ("wafc.gcv" %in% run_methods) {
    from_start()
    t0 <- proc.time()[["elapsed"]]
    tn <- try(wafc_tune(dgp[["x"]], dgp[["u"]], dgp[["y"]], rule = "gcv",
                        wavelet.table = wafc_pilot_table), silent = TRUE)
    el <- proc.time()[["elapsed"]] - t0
    if (inherits(tn, "try-error")) {
      for (lab in c("wafc.gcv", "wafc.gcv+max")) {
        rows[[length(rows) + 1L]] <- fail_row(cell, n, r, lab, dgp, active,
                                              tn)
      }
    } else {
      f <- tn[["fit"]]
      lam <- tn[["lambda"]]
      d_test <- wafc_design(test[["x"]], test[["u"]], spec = f[["design"]])
      d_grid <- wafc_design(matrix(1, n_grid, p), grid, spec = f[["design"]])
      cf <- wafc_raw_coef(f, s = lam)
      fh <- as.numeric(d_test[["Z"]] %*% cf[-1L, 1L]) + cf[1L, 1L]
      bh <- predict(f, newu = test[["u"]], s = lam, type = "beta")
      gh <- wafc_grid_components(f, grid, s = lam, design = d_grid)
      blk <- wafc_blocks(f, s = lam)[["nonzero"]] > 0L
      rows[[length(rows) + 1L]] <- one_row(
        cell, n, r, "wafc.gcv", dgp, test, grid,
        list(f_test = fh, beta_test = bh, time = el), active, blk, gh,
        list(J = tn[["J"]], lambda = lam,
             nzero = sum(cf[-1L, 1L][-f[["design"]][["unpenalized"]]] != 0)))
      gd <- tn[["guard"]]
      side <- side_add(side, list("gcv-guard" = data.frame(
        cell = cell[["name"]], n = n, rep = r, J = tn[["J"]], lambda = lam,
        excluded = gd[["excluded"]], points = gd[["points"]],
        decides = gd[["decides"]], J.unguarded = gd[["J.unguarded"]],
        lambda.unguarded = gd[["lambda.unguarded"]],
        df.unguarded = gd[["df.unguarded"]],
        df = tn[["nzero"]] + p, stringsAsFactors = FALSE)))
      th <- try(wafc_threshold(f, s = lam, rule = "max", c = 0.15,
                               refit = "none", y = dgp[["y"]]), silent = TRUE)
      if (inherits(th, "try-error")) {
        rows[[length(rows) + 1L]] <- fail_row(cell, n, r, "wafc.gcv+max", dgp,
                                              active, th)
      } else {
        ex <- th[["extra"]]
        rows[[length(rows) + 1L]] <- one_row(
          cell, n, r, "wafc.gcv+max", dgp, test, grid,
          list(f_test = predict(th, test[["x"]], test[["u"]]),
               beta_test = th[["beta"]](test[["u"]]),
               time = el + th[["time"]]),
          active, th[["blocks"]], wafc_grid_components(th, grid),
          list(J = ex[["J"]], lambda = ex[["lambda"]], nzero = ex[["nzero"]]))
        side <- side_add(side, list("thr-t" = data.frame(
          cell = cell[["name"]], n = n, rep = r, method = "wafc.gcv+max",
          t = ex[["t"]], c = ex[["c"]], stringsAsFactors = FALSE)))
      }
    }
  }

  ## 'gam.matched' is mgcv with the basis dimension of each smooth matched to
  ## the 2^J wavelet columns of a block at the J the WAFC with the LASSO
  ## selected on this replicate (wafc_k_matched(), step E2.4b), fitted by
  ## bam, which is what makes that dimension affordable (the note on
  ## wafc_fit_gam()). Its row carries that J, and its time is the time of
  ## the spline fit alone: the search that chose J is paid, and reported, in
  ## the column of the WAFC. Without a J to match, it fails as a row.
  ## 'klopp.free' is 'klopp' with penalize.levels = FALSE, and the two
  ## variants of step E2.5c, the one of step E2.5e and the one of step
  ## E2.5f keep the levels free as well.
  for (lab in setdiff(run_methods, c("wafc.lasso", "wafc.sglasso",
                                     "wafc.gcv"))) {
    from_start()
    mth <- sub(
      "\\.(matched|reml|gcv|cv|free|unit|merged|balanced|freecoarse)$", "",
      lab)
    own <- switch(lab,
                  klopp.free = list(penalize.levels = FALSE),
                  klopp.unit = list(penalize.levels = FALSE,
                                    chunk.weights = "unit"),
                  klopp.merged = list(penalize.levels = FALSE,
                                      merge.coarse = TRUE),
                  klopp.balanced = list(penalize.levels = FALSE,
                                        balanced = TRUE),
                  klopp.freecoarse = list(penalize.levels = FALSE,
                                          balanced = TRUE,
                                          free.coarse = TRUE),
                  gam.reml = list(k.select = "reml", engine = "bam"),
                  gam.gcv = list(k.select = "gcv", engine = gam_gcv_engine),
                  gam.cv = list(k.select = "cv", engine = "bam"),
                  list())
    if (lab == "gam.matched") {
      if (is.null(J_lasso)) {
        err <- try(stop("wafc.lasso failed, so there is no J to match"),
                   silent = TRUE)
        rows[[length(rows) + 1L]] <- fail_row(cell, n, r, lab, dgp, active,
                                              err)
        next
      }
      own <- list(k = wafc_k_matched(dgp[["u"]], J_lasso), engine = "bam")
    }
    f <- try(do.call(wafc_competitor,
                     c(list(mth, dgp[["x"]], dgp[["u"]], dgp[["y"]],
                            active = active, foldid = foldid,
                            wavelet.table = wafc_pilot_table), own)),
             silent = TRUE)
    if (inherits(f, "try-error")) {
      rows[[length(rows) + 1L]] <- fail_row(cell, n, r, lab, dgp, active, f)
      next
    }
    gh <- wafc_grid_components(f, grid)
    fc <- lab == "klopp.freecoarse"
    rows[[length(rows) + 1L]] <- one_row(
      cell, n, r, lab, dgp, test, grid,
      list(f_test = predict(f, test[["x"]], test[["u"]]),
           beta_test = f[["beta"]](test[["u"]]), time = f[["time"]]),
      active, if (fc) f[["extra"]][["blocks.fine"]] else f[["blocks"]], gh,
      list(J = if (lab == "gam.matched") J_lasso else f[["extra"]][["J"]],
           lambda = f[["extra"]][["lambda"]],
           nzero = if (fc) f[["extra"]][["nzero"]] else NULL))
    if (lab %in% thr_bases) {
      tr <- thr_rows(cell, n, r, lab, f, f[["time"]], dgp, test, grid,
                     active, foldid)
      rows <- c(rows, tr[["rows"]])
      side <- side_add(side, tr[["side"]])
    }
    ## the search over k of step E2.5h, one row per candidate fitted, with
    ## the dimensions chosen per modulator in the row of the candidate kept
    kt <- f[["extra"]][["k.table"]]
    if (!is.null(kt)) {
      kt[["chosen"]] <- kt[["k.used"]] ==
        paste(f[["extra"]][["k"]], collapse = ",")
      kt[["top"]] <- kt[["chosen"]] & f[["extra"]][["k.top"]]
      side <- side_add(side, list("gam-k" = cbind(
        data.frame(cell = cell[["name"]], n = n, rep = r, method = lab,
                   engine = f[["extra"]][["engine"]],
                   stringsAsFactors = FALSE), kt)))
    }
  }
  out <- do.call(rbind, rows)
  attr(out, "side") <- side
  out
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
    10000L * cell_index(cell) + 1000L * n_index(n) + r
  dgp <- draw_cell(cell, n, seed)
  test <- test_for(cell, dgp, seed)
  p <- cell[["p"]]
  grid <- grid_of(dgp, n_grid)
  active <- nzchar(dgp[["structure"]])
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:10, n))
  sp <- as.numeric(wafc_sprime(cell[["scenario"]]))
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

  ## One cross-validation serves the three rules that use it: cv.min and
  ## cv.1se read its two penalty levels, and the QUT takes its resolution.
  ## This part used to run the same cv.wafc() three times, once inside
  ## wafc_tune() for each of the first two rules and once for the QUT, on the
  ## same data and folds, which gives the same object three times: the
  ## numbers are unchanged and the part costs about a third. The time of
  ## each of the three rules is still the time of that cross-validation,
  ## plus, for the QUT, its own simulation and refit.
  t0 <- proc.time()[["elapsed"]]
  cv <- try(cv.wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], foldid = foldid,
                    wavelet.table = wafc_pilot_table), silent = TRUE)
  el_cv <- proc.time()[["elapsed"]] - t0
  if (inherits(cv, "try-error")) {
    for (rule in c("cv.min", "cv.1se", "qut")) {
      rows[[length(rows) + 1L]] <- fail_row(cell, n, r, rule, dgp, active, cv)
    }
  } else {
    fcv <- cv[["wafc.fit"]]
    rows[[length(rows) + 1L]] <- record("cv.min", fcv, cv[["lambda.min"]],
                                        el_cv, cv[["J.min"]])
    rows[[length(rows) + 1L]] <- record("cv.1se", fcv, cv[["lambda.1se"]],
                                        el_cv, cv[["J.min"]])
  }

  for (rule in c("bic", "ebic", "theory")) {
    t0 <- proc.time()[["elapsed"]]
    tn <- try(wafc_tune(dgp[["x"]], dgp[["u"]], dgp[["y"]], rule = rule,
                        foldid = foldid, s = sp, sigma = dgp[["sigma"]],
                        wavelet.table = wafc_pilot_table), silent = TRUE)
    if (inherits(tn, "try-error")) {
      rows[[length(rows) + 1L]] <- fail_row(cell, n, r, rule, dgp, active, tn)
      next
    }
    el <- proc.time()[["elapsed"]] - t0
    rows[[length(rows) + 1L]] <- record(rule, tn[["fit"]], tn[["lambda"]], el,
                                        tn[["J"]])
  }

  ## QUT: the resolution is not part of the rule, so it is taken from the
  ## same cross-validation the other rules use, and only lambda changes.
  if (!inherits(cv, "try-error")) {
    t0 <- proc.time()[["elapsed"]]
    lq <- wafc_lambda_qut(fcv[["design"]], dgp[["y"]], nsim = 200L,
                          seed = seed + 5L)
    fq2 <- wafc(design = fcv[["design"]], y = dgp[["y"]],
                lambda = wafc_path_to(lq, fcv[["design"]], dgp[["y"]]))
    el <- el_cv + proc.time()[["elapsed"]] - t0
    rows[[length(rows) + 1L]] <- record("qut", fq2, lq, el, cv[["J.min"]])
  }

  ## The oracle of the grid.
  t0 <- proc.time()[["elapsed"]]
  best <- NULL
  for (Ji in wafc_J_grid(NULL, n)) {
    fj <- wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], J = Ji,
               wavelet.table = wafc_pilot_table)
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
  ## exact evaluation here, and not the table of D31: this is the
  ## exception of that decision, a quantity read down to 1e-12
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
    10000L * cell_index(cell) + 1000L * n_index(n) + r
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
      fit <- try(wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], J = J, eps = e,
                      wavelet.table = wafc_pilot_table), silent = TRUE)
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
    10000L * cell_index(cell) + 1000L * n_index(n) + r
  dgp <- draw_cell(cell, n, seed)
  test <- test_for(cell, dgp, seed)
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:10, n))
  Jmax <- max(wafc_J_grid(NULL, n))
  rows <- list()
  for (lo in 1:2) {
    t0 <- proc.time()[["elapsed"]]
    cv <- try(cv.wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], J = lo:Jmax,
                      foldid = foldid, wavelet.table = wafc_pilot_table),
              silent = TRUE)
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
## Part "jgrid"
## ---------------------------------------------------------------------------

## The grid of cv.wall, 2:ceiling(log2(n)/2), against a grid that goes on
## to 'deep' levels. The question is not academic: the oracle of the grid
## of the part 'lambda' is at the top of the short grid in every replicate
## of every cell with components, so the short grid is truncated exactly
## where the choice is being made, and the resolution E1.3 says the
## inhomogeneous components need is above it.
run_jgrid <- function(cell, n, r, deep = 8L) {
  seed <- seed0 + 600000L +
    10000L * cell_index(cell) + 1000L * n_index(n) + r
  dgp <- draw_cell(cell, n, seed)
  test <- test_for(cell, dgp, seed)
  p <- cell[["p"]]
  grid <- grid_of(dgp, n_grid)
  active <- nzchar(dgp[["structure"]])
  set.seed(seed + 77L)
  foldid <- sample(rep_len(1:10, n))
  ## the short grid is the rule of cv.wall written out, since the default
  ## of cv.wafc() is now the deep one (decision D34)
  Jshort <- max(2L, as.integer(ceiling(log2(n) / 2)))
  rows <- list()
  for (lab in c("short", "deep")) {
    Jtop <- if (lab == "short") Jshort else deep
    t0 <- proc.time()[["elapsed"]]
    cv <- try(cv.wafc(dgp[["x"]], dgp[["u"]], dgp[["y"]], J = 2:Jtop,
                      foldid = foldid, wavelet.table = wafc_pilot_table),
              silent = TRUE)
    if (inherits(cv, "try-error")) next
    el <- proc.time()[["elapsed"]] - t0
    f <- cv[["wafc.fit"]]
    lam <- cv[["lambda.min"]]
    d_test <- wafc_design(test[["x"]], test[["u"]], spec = f[["design"]])
    d_grid <- wafc_design(matrix(1, n_grid, p), grid, spec = f[["design"]])
    cf <- wafc_raw_coef(f, s = lam)
    fh <- as.numeric(d_test[["Z"]] %*% cf[-1L, 1L]) + cf[1L, 1L]
    gh <- wafc_grid_components(f, grid, s = lam, design = d_grid)
    ise <- ise_components(gh, dgp, grid)
    rows[[length(rows) + 1L]] <- data.frame(
      cell = cell[["name"]], n = n, rep = r, grid = lab, Jtop = Jtop,
      J = cv[["J.min"]], lambda = lam, cvm = cv[["cvm.min"]],
      rmse_f = sqrt(mean((fh - test[["f"]])^2)),
      ise = sum(ise), ise_active = sum(ise[active]),
      nzero = sum(cf[-1L, 1L][-f[["design"]][["unpenalized"]]] != 0),
      time = el, stringsAsFactors = FALSE)
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
    ## a cell may carry its own budget, and then the smaller of the two
    ## counts: the cap of a part still applies to it
    nrep <- if (is.null(cell[["reps"]])) reps else min(reps, cell[["reps"]])
    for (n in ns_used) {
      for (r in seq_len(nrep)) {
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
  ## Failures are counted where they are printed, and not only in the
  ## messages: a job that fails as a whole leaves no row (a worker that dies
  ## returns a "try-error", which is not a data frame either), and a method
  ## that fails inside a job leaves a row with its message in 'error'.
  ok <- vapply(res, is.data.frame, TRUE)
  if (any(!ok)) {
    cat(sprintf("   %d of %d job(s) failed as a whole and left no row\n",
                sum(!ok), length(jobs)))
  }
  ## the side tables a job may carry (step E2.5g) are bound apart from the
  ## rows, and the rows are bound without them
  side <- list()
  for (z in res[ok]) {
    side <- side_add(side, attr(z, "side"))
  }
  out <- do.call(rbind, lapply(res[ok], function(z) {
    attr(z, "side") <- NULL
    z
  }))
  if (length(side)) attr(out, "side") <- side
  if (!is.null(out[["error"]]) && any(!is.na(out[["error"]]))) {
    bad <- out[!is.na(out[["error"]]), c("method", "error")]
    cat("   failed fits, kept as rows with no numbers:\n")
    for (m in unique(bad[["method"]])) {
      cat(sprintf("     %-14s %3d: %s\n", m, sum(bad[["method"]] == m),
                  substr(bad[["error"]][bad[["method"]] == m][1L], 1L, 100L)))
    }
  }
  out
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
  side <- attr(res, "side")
  attr(res, "side") <- NULL
  saveRDS(res, out_file("competitors"))
  for (nm in names(side)) saveRDS(side[[nm]], out_file(nm))
  res[["method"]] <- factor(res[["method"]], levels = method_levels)
  cat("\nmedians by cell, n and method (rmse_f out of sample, ISE of the",
      "components, blocks kept of the active and of the null ones, the J",
      "selected or matched, s):\n")
  print(med(res, c("rmse_f", "mse_beta", "ise", "ise_active", "ise_null",
                   "n_true", "n_false", "J", "time"),
            c("cell", "n", "method")),
        row.names = FALSE, digits = 3)
}
if ("competitors" %in% parts && "wafc.lasso" %in% run_methods) {
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
if ("competitors" %in% parts) {
  cat("\nstructure recovered: fraction of replicates with S-hat = S, and",
      "mean false positives and false negatives (blocks):\n")
  ok <- !is.na(res[["n_true"]])
  st <- res[ok, ]
  st[["exact"]] <- st[["n_true"]] == st[["n_active"]] & st[["n_false"]] == 0L
  st[["fp"]] <- st[["n_false"]]
  st[["fn"]] <- st[["n_active"]] - st[["n_true"]]
  a <- stats::aggregate(st[c("exact", "fp", "fn")],
                        st[c("cell", "n", "method")], FUN = mean)
  print(a[do.call(order, unname(as.list(a[c("cell", "n", "method")]))), ],
        row.names = FALSE, digits = 3)
}

if ("lambda" %in% parts) {
  res <- sweep_part(run_lambda, "lambda rules")
  saveRDS(res, out_file("lambda"))
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
  saveRDS(gm, out_file("gram"))
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
  saveRDS(res, out_file("margin"))
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
  saveRDS(at_cv, out_file("margin-atcv"))
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
  saveRDS(res, out_file("j1"))
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

if ("jgrid" %in% parts) {
  res <- sweep_part(run_jgrid, "depth of the grid of J",
                    cells_used = cells[vapply(cells, `[[`, "", "name") %in%
                                         c("smooth", "inhomogeneous")],
                    reps = R_small)
  saveRDS(res, out_file("jgrid"))
  cat("\nmedians by cell, n and grid:\n")
  print(med(res, c("Jtop", "J", "lambda", "nzero", "cvm", "rmse_f", "ise",
                   "ise_active", "time"), c("cell", "n", "grid")),
        row.names = FALSE, digits = 3)
  w <- stats::reshape(res[c("cell", "n", "rep", "grid", "J", "rmse_f", "ise",
                            "time")],
                      idvar = c("cell", "n", "rep"), timevar = "grid",
                      direction = "wide")
  w[["r_rmse"]] <- w[["rmse_f.deep"]] / w[["rmse_f.short"]]
  w[["r_ise"]] <- w[["ise.deep"]] / w[["ise.short"]]
  w[["r_time"]] <- w[["time.deep"]] / w[["time.short"]]
  cat("\nthe deep grid against the short one, paired, and how often the deep",
      "grid goes past the top of the short one:\n")
  tab <- stats::aggregate(
    cbind(r_rmse, r_ise, r_time, past = J.deep > J.short,
          same = J.deep == J.short) ~ cell + n, data = w,
    FUN = function(v) stats::median(v))
  print(tab[order(tab[["cell"]], tab[["n"]]), ], row.names = FALSE, digits = 3)
  cat("\nfraction of replicates in which the deep grid wins:\n")
  print(stats::aggregate(cbind(wins = r_rmse < 1, past = J.deep > J.short) ~
                           cell + n, data = w, FUN = mean),
        row.names = FALSE, digits = 3)
}

cat("\nOK\n")
