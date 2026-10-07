## E4.3 -- what the pilot of the simulation study decides.
##
##     Rscript wafc/scripts/11-pilot-read.R [--dir=wafc/cache/e43] [--top-n=2000]
##         [--arms-n=1000] [--workers=8] [--trigger=0.2] [--j-factor=2] [--hours=10]
##
## Run from the root of wafc-draft, after the round of
## wafc/scripts/11-pilot.yaml (its header lists the four steps). Reads the
## cached units of the round, <dir>/pilot/units/ (the core cells and the
## arms) and <dir>/scale/units/ (config/scale.yaml), and what each step
## left in <dir>/steps/ (<step>.<stamp>.log, the log of 01_simulate.R, and
## <step>.<stamp>.time, the report of /usr/bin/time -v; a step resumed
## after a stop has one pair per run, summed here). Prints, and writes to
## <dir>/read/:
##
##   (0) checks: the cells of the pilot are those of config/core.yaml and
##       config/arms.yaml, the entries it inherits are the study's, no seed
##       of the pilot is a seed of the production, the commit of the WAFC
##       code the units were fitted with, and the failed fits;
##   (i) the median time of every method in every core cell at n = top-n
##       and in every arm at n = arms-n (times.csv has every cell and n),
##       and the projection of the cost of the production (E4.4), core,
##       arms and scale, in processor hours and in hours on `workers`
##       workers (projection.csv, projection-designs.csv);
##  (ii) the cost of the scale arm against the ~10 h of D60(c), in both
##       readings of the 10 h;
## (iii) the fraction of the replicates at the top of the effective grid
##       by method and cell at n = top-n, against the trigger of D60(e),
##       and the entry "<top-n>" of config/study.yaml the rule asks for
##       (top.csv, rule-<top-n>.yaml);
##  (iv) the peak memory of a process and the processor time against the
##       elapsed time, by step (steps.csv);
##       and the health of the fits at n = top-n: cuts of the engines at
##       the chosen penalty, and the basis check of mgcv.
##
## Every time is the elapsed time of a unit (one method on one replicate,
## the drawing of the data and the reading of the fit included), which is
## what the production pays; the methods run single threaded (reference
## BLAS), and the ratio of processor to elapsed time by step checks it.
## The projection prices every unit of the production at the mean time of
## its (cell, n, method) in the pilot (the mean, not the median, because a
## total is a sum) and schedules the units of each design as run_study()
## of the compendium does: in the order of unit_order(), each to the first
## worker free (mclapply with mc.preschedule = FALSE); the designs run one
## after the other, as in run_all.R. A (cell, n, method) the pilot did not
## run is priced at the nearest n it did run and counted as a proxy; the
## round has none, the smoke has many.
##
## The rule (D60(e)). A method whose choice is at the top of its grid in
## more than `trigger` of the replicates of some cell with components at
## n = top-n gets a larger grid there, in every cell (the rule is by n, not
## by cell, D41): one more level of J (2:9) for the methods on the grid of
## J (wafc, lasso, klopp, oracle, and bsgl, whose dimensions are 2^J), and
## k = 120 added for the gam methods. The cell without components is shown
## and left out of the trigger, as D34 counted "replicates with a
## component"; the fraction pooled over the cells is shown too. The cost
## of the larger grids is an estimate, not a measurement: for the gam, the
## time of k = 120 extrapolated from those of k = 40 and k = 80 in the
## same search (side table gam_k), as a power law in k; for the grids of
## J, the unit time times `j-factor` (2 by default, because the level
## J = 9 has as many columns as the levels below it together; a
## conjecture).

## --- options ------------------------------------------------------------

opt <- list(dir = "wafc/cache/e43", compendium = "wafc-studies",
            pilot = "wafc/scripts/11-pilot.yaml", "top-n" = "2000",
            "arms-n" = "1000", workers = "8", trigger = "0.2",
            "j-factor" = "2", hours = "10")
for (a in commandArgs(trailingOnly = TRUE)) {
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1L]]
  if (!startsWith(a, "--") || length(kv) != 2L || !(kv[1L] %in% names(opt))) {
    stop("unknown argument '", a, "'; options: ",
         paste0("--", names(opt), collapse = ", "), call. = FALSE)
  }
  opt[[kv[1L]]] <- kv[2L]
}
top_n <- as.integer(opt[["top-n"]])
arms_n <- as.integer(opt[["arms-n"]])
workers <- as.integer(opt[["workers"]])
trigger <- as.numeric(opt[["trigger"]])
j_factor <- as.numeric(opt[["j-factor"]])
hours_max <- as.numeric(opt[["hours"]])
dir <- opt[["dir"]]
comp <- opt[["compendium"]]
outdir <- file.path(dir, "read")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## The compendium's own reading of a configuration, of its units and of
## their order of dispatch: the projection prices the units the
## production will run, in the order it will run them.
for (f in c("config.R", "seeds.R", "run.R", "methods.R", "metrics.R",
            "aggregate.R")) {
  source(file.path(comp, "R", f))
}
study_file <- file.path(comp, "config", "study.yaml")
prod <- lapply(c(core = "core", arms = "arms", scale = "scale"), function(d)
  study_config(file.path(comp, "config", paste0(d, ".yaml")), study_file))
pcfg <- study_config(opt[["pilot"]], study_file)

j_methods <- c("wafc", "lasso", "klopp", "oracle", "bsgl")
k_methods <- c("gam.reml", "gam.gcv")
method_order <- prod[["core"]][["methods"]]
core_cells <- names(prod[["core"]][["cells"]])
arm_cells <- names(prod[["arms"]][["cells"]])

section <- function(title) cat("\n", strrep("=", 76), "\n", title, "\n",
                               strrep("=", 76), "\n", sep = "")
write_csv <- function(x, name) {
  if (is.null(x) || nrow(x) == 0L) return(invisible(NULL))
  utils::write.csv(x, file.path(outdir, name), row.names = FALSE)
}
fmt <- function(x, d = 1L) formatC(x, format = "f", digits = d)

## --- reading the units --------------------------------------------------

read_units <- function(path) {
  f <- list.files(file.path(path, "units"), pattern = "^rep[0-9]+\\.rds$",
                  recursive = TRUE, full.names = TRUE)
  if (length(f) == 0L) return(NULL)
  us <- lapply(f, readRDS)
  first <- function(v) { v <- v[!is.na(v)]; if (length(v)) v[1L] else NA }
  units <- do.call(rbind, lapply(us, function(u) {
    r <- u[["rows"]]
    st <- u[["settings"]]
    data.frame(design = r[["design"]][1L], cell = st[["cell"]][["name"]],
               n = as.integer(st[["n"]]), method = st[["method"]],
               rep = as.integer(st[["rep"]]), elapsed = u[["elapsed"]],
               failed = sum(!is.na(r[["error"]])), rows = nrow(r),
               active = r[["n_active"]][1L],
               top = as.logical(first(r[["top"]])),
               J = as.integer(first(r[["J"]])), k = as.integer(first(r[["k"]])),
               kcheck_low = as.integer(first(r[["kcheck_n_low"]])),
               code = if (is.null(u[["code"]])) NA_character_ else u[["code"]],
               seeds = paste(unlist(st[["seeds"]]), collapse = ","),
               stringsAsFactors = FALSE)
  }))
  side <- function(nm) bind_rows(lapply(us, function(u) u[["side"]][[nm]]))
  list(units = units, rows = bind_rows(lapply(us, `[[`, "rows")),
       gam_k = side("gam_k"), convergence = side("convergence"))
}

pil <- read_units(file.path(dir, "pilot"))
sca <- read_units(file.path(dir, "scale"))
if (is.null(pil) && is.null(sca)) {
  stop("no unit under ", file.path(dir, "pilot"), " or ",
       file.path(dir, "scale"), ".", call. = FALSE)
}
all_units <- rbind(pil[["units"]], sca[["units"]])
all_rows <- bind_rows(list(pil[["rows"]], sca[["rows"]]))

## --- the top of the grids (iii), computed first: the projection with
## --- the raised grids depends on it ------------------------------------

wilson <- function(x, n, z = 1.96) {
  if (n == 0L) return(c(NA_real_, NA_real_))
  p <- x / n
  c0 <- p + z^2 / (2 * n)
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(c0 - h, c0 + h) / (1 + z^2 / n)
}
chosen_table <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_character_)
  tb <- table(v)
  paste(names(tb), tb, sep = ":", collapse = " ")
}
top_rows <- function(u, n, cells) {
  u <- u[u[["n"]] == n & u[["cell"]] %in% cells &
           u[["method"]] %in% c(j_methods, k_methods), , drop = FALSE]
  if (nrow(u) == 0L) return(NULL)
  key <- unique(u[, c("method", "cell")])
  out <- do.call(rbind, lapply(seq_len(nrow(key)), function(i) {
    s <- u[u[["method"]] == key[["method"]][i] &
             u[["cell"]] == key[["cell"]][i], , drop = FALSE]
    t <- s[["top"]]
    nv <- sum(!is.na(t))
    nt <- sum(t, na.rm = TRUE)
    ci <- wilson(nt, nv)
    dim <- if (key[["method"]][i] %in% k_methods) s[["k"]] else s[["J"]]
    data.frame(n = n, method = key[["method"]][i], cell = key[["cell"]][i],
               components = any(s[["active"]] > 0L, na.rm = TRUE),
               replicates = nrow(s), read = nv, at_top = nt,
               fraction = if (nv) nt / nv else NA_real_,
               ci_low = ci[1L], ci_high = ci[2L],
               chosen = chosen_table(dim), stringsAsFactors = FALSE)
  }))
  out[order(match(out[["method"]], method_order),
            match(out[["cell"]], cells)), , drop = FALSE]
}

top <- top_rows(all_units, top_n, core_cells)
rule <- NULL
if (!is.null(top)) {
  tc <- top[top[["components"]] & top[["read"]] > 0L, , drop = FALSE]
  rule <- do.call(rbind, lapply(intersect(c(j_methods, k_methods),
                                          unique(top[["method"]])), function(m) {
    s <- tc[tc[["method"]] == m, , drop = FALSE]
    if (nrow(s) == 0L) {
      return(data.frame(method = m, worst_cell = NA_character_,
                        worst = NA_real_, pooled = NA_real_, raise = NA,
                        borderline = NA, stringsAsFactors = FALSE))
    }
    w <- which.max(s[["fraction"]])
    worst <- s[["fraction"]][w]
    data.frame(method = m, worst_cell = s[["cell"]][w], worst = worst,
               pooled = sum(s[["at_top"]]) / sum(s[["read"]]),
               raise = worst > trigger,
               borderline = abs(worst - trigger) <= 0.1 + 1e-12,
               stringsAsFactors = FALSE)
  }))
}
raised <- if (is.null(rule)) character(0) else
  rule[["method"]][rule[["raise"]] %in% TRUE]

## The entry of config/study.yaml the rule asks for, in the form of the
## comment in that file.
entry <- if (length(raised) == 0L) {
  sprintf("  \"%d\": *tuning", top_n)
} else {
  w <- max(nchar(raised)) + 1L
  c(sprintf("  \"%d\":", top_n), "    <<: *tuning", "    methods:",
    vapply(raised, function(m) {
      tn <- method_tuning(prod[["core"]], top_n, m)
      g <- if (m %in% k_methods) {
        sprintf("{k: [%s]}", paste(c(tn[["k"]], 120L), collapse = ", "))
      } else {
        sprintf("{J: [%s]}", paste(c(tn[["J"]], max(tn[["J"]]) + 1L),
                                   collapse = ", "))
      }
      sprintf("      %-*s %s", w, paste0(m, ":"), g)
    }, "", USE.NAMES = FALSE))
}

## --- the projection -----------------------------------------------------

look <- stats::aggregate(elapsed ~ cell + n + method, all_units, mean)

## The time of every unit: the mean of its (cell, n, method) in the pilot,
## or of the nearest n the pilot ran for that (cell, method).
price <- function(units) {
  key <- unique(units[c("cell", "n", "method")])
  key[["time"]] <- NA_real_
  key[["proxy"]] <- FALSE
  for (i in seq_len(nrow(key))) {
    s <- look[look[["cell"]] == key[["cell"]][i] &
                look[["method"]] == key[["method"]][i], , drop = FALSE]
    if (nrow(s) == 0L) next
    j <- which(s[["n"]] == key[["n"]][i])
    if (length(j)) {
      key[["time"]][i] <- s[["elapsed"]][j]
    } else {
      key[["time"]][i] <- s[["elapsed"]][which.min(abs(log(s[["n"]] / key[["n"]][i])))]
      key[["proxy"]][i] <- TRUE
    }
  }
  m <- match(paste(units[["cell"]], units[["n"]], units[["method"]]),
             paste(key[["cell"]], key[["n"]], key[["method"]]))
  units[["time"]] <- key[["time"]][m]
  units[["proxy"]] <- key[["proxy"]][m]
  units
}

makespan <- function(d, w) {
  free <- numeric(w)
  for (x in d) { i <- which.min(free); free[i] <- free[i] + x }
  max(free)
}

project <- function(cfg, adjust = NULL) {
  u <- price(unit_order(study_units(cfg)))
  if (!is.null(adjust)) u[["time"]] <- adjust(u)
  ok <- !is.na(u[["time"]])
  by_n <- do.call(rbind, lapply(sort(unique(u[["n"]])), function(n) {
    s <- u[u[["n"]] == n, , drop = FALSE]
    data.frame(design = cfg[["design"]], n = n, units = nrow(s),
               unmeasured = sum(is.na(s[["time"]])),
               proxy = sum(s[["proxy"]]),
               cpu_h = sum(s[["time"]], na.rm = TRUE) / 3600)
  }))
  list(by_n = by_n, cpu_h = sum(u[["time"]][ok]) / 3600,
       wall_h = makespan(u[["time"]][ok], workers) / 3600, units = u)
}

pj <- lapply(prod, project)

## With the grids the rule asks for at n = top-n (an estimate; see the
## header).
extra_k <- NULL
proj_raised <- NULL
if (length(raised)) {
  gk <- pil[["gam_k"]]
  if (!is.null(gk) && any(raised %in% k_methods)) {
    gk <- gk[gk[["n"]] == top_n & gk[["k"]] %in% c(40L, 80L), , drop = FALSE]
    if (nrow(gk)) {
      gk <- stats::aggregate(time ~ method + cell + rep + k, gk, sum)
      w <- stats::reshape(gk, idvar = c("method", "cell", "rep"),
                          timevar = "k", direction = "wide")
      w <- w[is.finite(w[["time.40"]]) & is.finite(w[["time.80"]]) &
               w[["time.40"]] > 0, , drop = FALSE]
      w[["t120"]] <- w[["time.80"]] *
        1.5^(log(w[["time.80"]] / w[["time.40"]]) / log(2))
      extra_k <- stats::aggregate(cbind(time.80, t120) ~ method + cell, w, mean)
    }
  }
  adjust <- function(u) {
    d <- u[["time"]]
    at <- u[["n"]] == top_n
    for (m in intersect(raised, j_methods)) {
      i <- at & u[["method"]] == m
      d[i] <- d[i] * j_factor
    }
    for (m in intersect(raised, k_methods)) {
      i <- which(at & u[["method"]] == m)
      add <- if (is.null(extra_k)) rep(NA_real_, length(i)) else
        extra_k[["t120"]][match(paste(m, u[["cell"]][i]),
                                paste(extra_k[["method"]], extra_k[["cell"]]))]
      d[i] <- d[i] + ifelse(is.na(add), 0, add)
    }
    d
  }
  proj_raised <- project(prod[["core"]], adjust)
}

## --- (0) checks ---------------------------------------------------------

section("(0) Checks")

prod_cells <- c(prod[["core"]][["cells"]], prod[["arms"]][["cells"]])
drift <- character(0)
for (nm in names(pcfg[["cells"]])) {
  if (is.null(prod_cells[[nm]])) {
    drift <- c(drift, paste0(nm, " (not a cell of the production)"))
  } else if (!identical(pcfg[["cells"]][[nm]], prod_cells[[nm]])) {
    drift <- c(drift, nm)
  }
}
cat("cells of the pilot equal to the production's: ",
    if (length(drift)) paste("NO:", paste(drift, collapse = ", ")) else
      sprintf("yes, %d of %d", length(pcfg[["cells"]]),
              length(pcfg[["cells"]])), "\n", sep = "")
inherit <- c("code", "methods", "method_options", "tuning", "test_size",
             "grid_size", "curves", "labels", "reference")
diff_keys <- inherit[!vapply(inherit, function(k)
  identical(pcfg[[k]], prod[["core"]][[k]]), TRUE)]
cat("entries inherited from the study file and equal to the production's: ",
    if (length(diff_keys)) paste("NO:", paste(diff_keys, collapse = ", "))
    else paste(inherit, collapse = ", "), "\n", sep = "")
sd_keep <- setdiff(names(prod[["core"]][["seeds"]]), "master")
cat("seeds: master ", pcfg[["seeds"]][["master"]], " (the production's ",
    prod[["core"]][["seeds"]][["master"]], "); keys, sample sizes and radix ",
    "equal to the production's: ",
    identical(pcfg[["seeds"]][sd_keep], prod[["core"]][["seeds"]][sd_keep]),
    "\n", sep = "")

## Every seed the production draws, against the seeds the units of the
## pilot recorded; and those against the seeds the configurations of the
## pilot give, which says the units came from them.
all_seeds <- function(cfg, units) {
  key <- unique(units[, c("cell", "n", "rep")])
  unlist(lapply(seq_len(nrow(key)), function(i)
    unlist(unit_seeds(cfg, cfg[["cells"]][[key[["cell"]][i]]],
                      key[["n"]][i], key[["rep"]][i]))))
}
recorded <- function(u) {
  if (is.null(u)) return(integer(0))
  as.integer(unlist(strsplit(u[["seeds"]], ",", fixed = TRUE)))
}
prod_seeds <- unique(unlist(lapply(prod, function(cfg)
  all_seeds(cfg, study_units(cfg)))))
pilot_seeds <- unique(c(recorded(pil[["units"]]), recorded(sca[["units"]])))
expected <- unique(c(
  if (!is.null(pil)) all_seeds(pcfg, pil[["units"]]),
  if (!is.null(sca)) all_seeds(prod[["scale"]], sca[["units"]])))
cat(sprintf("seeds the units recorded that the production also draws: %d of %d (the production draws %d)\n",
            sum(pilot_seeds %in% prod_seeds), length(pilot_seeds),
            length(prod_seeds)))
cat("seeds the units recorded equal to those of the pilot's configurations: ",
    setequal(pilot_seeds, expected), "\n", sep = "")
codes <- table(all_units[["code"]], useNA = "ifany")
cat("WAFC code of the units:\n")
for (i in seq_along(codes)) {
  cat(sprintf("  %5d  %s\n", codes[[i]], names(codes)[i]))
}
nfail <- sum(all_units[["failed"]] > 0L)
cat(sprintf("units read: %d (%d in pilot/, %d in scale/); with a failed row: %d\n",
            nrow(all_units), NROW(pil[["units"]]), NROW(sca[["units"]]), nfail))
if (nfail) {
  bad <- all_rows[!is.na(all_rows[["error"]]),
                  c("cell", "n", "rep", "method", "error")]
  bad[["error"]] <- substr(gsub("\\s+", " ", bad[["error"]]), 1L, 70L)
  print(bad, row.names = FALSE)
}

## --- (i) times and the projection ---------------------------------------

section("(i) Median time of a unit, in seconds, by method and cell")

wide_time <- function(u, n, cells) {
  u <- u[u[["n"]] == n & u[["cell"]] %in% cells, , drop = FALSE]
  if (nrow(u) == 0L) return(NULL)
  cells <- intersect(cells, unique(u[["cell"]]))
  ms <- intersect(method_order, unique(u[["method"]]))
  m <- matrix(NA_real_, length(ms), length(cells), dimnames = list(ms, cells))
  r <- m
  for (i in ms) for (j in cells) {
    v <- u[["elapsed"]][u[["method"]] == i & u[["cell"]] == j]
    if (length(v)) { m[i, j] <- stats::median(v); r[i, j] <- length(v) }
  }
  list(median = m, reps = r)
}
show_time <- function(w, title) {
  if (is.null(w)) { cat("\n", title, ": no unit\n", sep = ""); return(invisible()) }
  m <- w[["median"]]
  out <- rbind(m, "sum of medians" = colSums(m, na.rm = TRUE))
  cat("\n", title, " (replicates per entry: ",
      paste(sort(unique(as.vector(w[["reps"]][!is.na(w[["reps"]])]))),
            collapse = ", "), ")\n", sep = "")
  print(round(out, 1), na.print = "-")
}
show_time(wide_time(all_units, top_n, core_cells),
          sprintf("core cells, n = %d", top_n))
show_time(wide_time(all_units, arms_n, arm_cells),
          sprintf("arms, n = %d", arms_n))

times <- stats::aggregate(elapsed ~ design + cell + n + method, all_units,
                          function(v) c(units = length(v),
                                        median = stats::median(v),
                                        mean = mean(v), max = max(v)))
times <- cbind(times[c("design", "cell", "n", "method")],
               as.data.frame(times[["elapsed"]]))
times <- times[order(times[["design"]], times[["cell"]], times[["n"]],
                     match(times[["method"]], method_order)), ]
write_csv(times, "times.csv")

section(sprintf("(i) Projection of the production (E4.4) on %d workers", workers))

by_n <- do.call(rbind, lapply(pj, `[[`, "by_n"))
pr <- by_n
pr[["cpu_h"]] <- fmt(pr[["cpu_h"]])
print(pr, row.names = FALSE)
tot <- data.frame(design = names(pj),
                  units = vapply(pj, function(p) nrow(p[["units"]]), 0L),
                  cpu_h = vapply(pj, `[[`, 0, "cpu_h"),
                  wall_h = vapply(pj, `[[`, 0, "wall_h"))
cat("\nBy design (processor hours; hours on", workers,
    "workers, the units of a design scheduled as run_study() does):\n")
print(transform(tot, cpu_h = fmt(cpu_h), wall_h = fmt(wall_h)),
      row.names = FALSE)
cat(sprintf("\ncore + arms:          %6.1f h of processor, %5.1f h on %d workers\n",
            sum(tot[["cpu_h"]][1:2]), sum(tot[["wall_h"]][1:2]), workers))
cat(sprintf("core + arms + scale:  %6.1f h of processor, %5.1f h on %d workers\n",
            sum(tot[["cpu_h"]]), sum(tot[["wall_h"]]), workers))
nprox <- sum(by_n[["proxy"]])
nunm <- sum(by_n[["unmeasured"]])
if (nprox + nunm > 0L) {
  cat(sprintf("INCOMPLETE: %d unit(s) priced at another n (proxy); %d unit(s) of a (cell, method) the pilot did not run, left out\n",
              nprox, nunm))
}
if (!is.null(proj_raised)) {
  if (!is.null(extra_k)) {
    cat("\nEstimated time of k = 120 at n =", top_n,
        "(power law through k = 40 and k = 80), seconds:\n")
    print(transform(extra_k, time.80 = fmt(time.80), t120 = fmt(t120)),
          row.names = FALSE)
  }
  rn <- proj_raised[["by_n"]]
  cat(sprintf("\nWith the grids the rule asks for at n = %d (%s; an estimate: J times %.1f, k = 120 extrapolated):\n",
              top_n, paste(raised, collapse = ", "), j_factor))
  cat(sprintf("  core:                %6.1f h of processor (n = %d: %.1f h, against %.1f h), %5.1f h on %d workers\n",
              proj_raised[["cpu_h"]], top_n, rn[["cpu_h"]][rn[["n"]] == top_n],
              by_n[["cpu_h"]][by_n[["design"]] == "core" & by_n[["n"]] == top_n],
              proj_raised[["wall_h"]], workers))
  cat(sprintf("  core + arms:         %6.1f h of processor, %5.1f h on %d workers\n",
              proj_raised[["cpu_h"]] + tot[["cpu_h"]][2L],
              proj_raised[["wall_h"]] + tot[["wall_h"]][2L], workers))
  cat(sprintf("  core + arms + scale: %6.1f h of processor, %5.1f h on %d workers\n",
              proj_raised[["cpu_h"]] + sum(tot[["cpu_h"]][2:3]),
              proj_raised[["wall_h"]] + sum(tot[["wall_h"]][2:3]), workers))
}
write_csv(rbind(cbind(by_n, grids = "current"),
                if (!is.null(proj_raised))
                  cbind(proj_raised[["by_n"]], grids = "raised")),
          "projection.csv")
write_csv(rbind(cbind(tot, grids = "current"),
                if (!is.null(proj_raised))
                  data.frame(design = "core", units = nrow(proj_raised[["units"]]),
                             cpu_h = proj_raised[["cpu_h"]],
                             wall_h = proj_raised[["wall_h"]],
                             grids = "raised")),
          "projection-designs.csv")

## --- (ii) the scale -------------------------------------------------------

section("(ii) The scale arm at n = 1000 against ~10 h (D60(c))")

ps <- pj[["scale"]]
su <- all_units[all_units[["cell"]] == "scale", , drop = FALSE]
if (nrow(su)) {
  st <- stats::aggregate(elapsed ~ method, su, function(v)
    c(units = length(v), median = stats::median(v), mean = mean(v)))
  st <- cbind(st["method"], round(as.data.frame(st[["elapsed"]]), 1))
  print(st[order(match(st[["method"]], method_order)), ], row.names = FALSE)
}
cat(sprintf("\nscale, %d replicates: %.1f h of processor, %.1f h on %d workers\n",
            prod[["scale"]][["replicates"]], ps[["cpu_h"]], ps[["wall_h"]],
            workers))
cat(sprintf("against ~%.0f h: %s in processor hours, %s in hours on %d workers\n",
            hours_max, if (ps[["cpu_h"]] < hours_max) "below" else "ABOVE",
            if (ps[["wall_h"]] < hours_max) "below" else "ABOVE", workers))
miss <- setdiff(cell_methods(prod[["scale"]], prod[["scale"]][["cells"]][["scale"]]),
                unique(su[["method"]]))
if (length(miss)) {
  cat("INCOMPLETE: methods of the scale arm not measured:",
      paste(miss, collapse = ", "), "\n")
}

## --- (iii) the top of the grids -------------------------------------------

section(sprintf("(iii) Replicates at the top of the effective grid, n = %d", top_n))

show_top <- function(tp) {
  pr <- tp
  pr[["fraction"]] <- fmt(pr[["fraction"]], 2L)
  pr[["ci95"]] <- sprintf("[%s, %s]", fmt(tp[["ci_low"]], 2L),
                          fmt(tp[["ci_high"]], 2L))
  pr[["components"]] <- ifelse(pr[["components"]], "yes", "no")
  print(pr[c("method", "cell", "components", "at_top", "read", "fraction",
             "ci95", "chosen")], row.names = FALSE)
}
if (is.null(top)) {
  cat("no unit at n =", top_n, "\n")
} else {
  show_top(top)
  cat(sprintf("\nThe rule: the largest fraction over the cells with components, against %.0f%%:\n",
              100 * trigger))
  pr <- rule
  pr[["worst"]] <- fmt(pr[["worst"]], 2L)
  pr[["pooled"]] <- fmt(pr[["pooled"]], 2L)
  pr[["raise"]] <- ifelse(is.na(rule[["raise"]]), "-",
                          ifelse(rule[["raise"]], "RAISE", "keep"))
  pr[["borderline"]] <- ifelse(rule[["borderline"]] %in% TRUE,
                               "within 0.1 of the trigger", "")
  print(pr, row.names = FALSE)
  write_csv(merge(top, rule[c("method", "raise", "borderline")],
                  by = "method", all.x = TRUE), "top.csv")
}
missing_rule <- setdiff(c(j_methods, k_methods), rule[["method"]])
if (length(missing_rule)) {
  cat("not measured at n = ", top_n, ": ", paste(missing_rule, collapse = ", "),
      "\n", sep = "")
}
cat("\nThe entry of config/study.yaml the rule asks for:\n\n")
cat(entry, sep = "\n")
writeLines(entry, file.path(outdir, sprintf("rule-%d.yaml", top_n)))

## The top elsewhere, for information: the grids of the arms and of the
## scale are not under the rule.
for (b in list(list(n = arms_n, cells = arm_cells, what = "the arms"),
               list(n = 1000L, cells = "scale",
                    what = "the scale arm (grid of J 2:7)"))) {
  tp <- top_rows(all_units, b[["n"]], b[["cells"]])
  if (is.null(tp)) next
  tp <- tp[tp[["components"]] & tp[["read"]] > 0L, , drop = FALSE]
  if (nrow(tp) == 0L) next
  agg <- stats::aggregate(cbind(at_top, read) ~ method, tp, sum)
  agg[["pooled"]] <- fmt(agg[["at_top"]] / agg[["read"]], 2L)
  agg[["worst"]] <- vapply(agg[["method"]], function(m)
    fmt(max(tp[["fraction"]][tp[["method"]] == m]), 2L), "")
  cat(sprintf("\nFor information, %s at n = %d:\n", b[["what"]], b[["n"]]))
  print(agg[order(match(agg[["method"]], method_order)), ], row.names = FALSE)
}

## --- (iv) memory and processor time, by step ------------------------------

section("(iv) Peak memory of a process, and processor against elapsed time, by step")

parse_time <- function(f) {
  x <- readLines(f, warn = FALSE)
  val <- function(label) {
    l <- grep(label, x, value = TRUE, fixed = TRUE)
    if (!length(l)) return(NA_character_)
    trimws(sub(label, "", l[1L], fixed = TRUE))
  }
  wall <- val("Elapsed (wall clock) time (h:mm:ss or m:ss):")
  ws <- if (is.na(wall)) NA_real_ else {
    p <- as.numeric(strsplit(wall, ":", fixed = TRUE)[[1L]])
    sum(p * 60^(rev(seq_along(p)) - 1L))
  }
  list(rss_gb = as.numeric(val("Maximum resident set size (kbytes):")) / 1024^2,
       user = as.numeric(val("User time (seconds):")),
       sys = as.numeric(val("System time (seconds):")),
       wall = ws, exit = as.integer(val("Exit status:")))
}
parse_log <- function(f) {
  x <- readLines(f, warn = FALSE)
  m <- regmatches(x, regexec("^\\[run\\] \\S+ n=\\d+ \\S+ rep \\d+: ([0-9.]+) s", x))
  secs <- as.numeric(vapply(m, function(v) if (length(v)) v[2L] else NA_character_, ""))
  hd <- grep("^\\[run\\] design", x, value = TRUE)
  list(fitted = sum(!is.na(secs)), unit_h = sum(secs, na.rm = TRUE) / 3600,
       failed = sum(grepl("UNIT FAILED", x, fixed = TRUE)),
       row_errors = sum(grepl("(row error:", x, fixed = TRUE)),
       cached = if (length(hd)) as.integer(sub(".* (\\d+) cached.*", "\\1", hd[1L]))
                else NA_integer_)
}
tf <- sort(list.files(file.path(dir, "steps"), pattern = "\\.time$",
                      full.names = TRUE))
runs <- do.call(rbind, lapply(tf, function(f) {
  tt <- parse_time(f)
  lf <- sub("\\.time$", ".log", f)
  lg <- if (file.exists(lf)) parse_log(lf) else
    list(fitted = NA, unit_h = NA, failed = NA, row_errors = NA, cached = NA)
  data.frame(step = sub("\\..*$", "", basename(f)), exit = tt[["exit"]],
             fitted = lg[["fitted"]], cached = lg[["cached"]],
             failed = lg[["failed"]], row_errors = lg[["row_errors"]],
             unit_h = lg[["unit_h"]], cpu_h = (tt[["user"]] + tt[["sys"]]) / 3600,
             wall_h = tt[["wall"]] / 3600, max_rss_gb = tt[["rss_gb"]],
             stringsAsFactors = FALSE)
}))
steps <- if (is.null(runs)) NULL else
  do.call(rbind, lapply(split(runs, factor(runs[["step"]], unique(runs[["step"]]))),
                        function(d) data.frame(
    step = d[["step"]][1L], runs = nrow(d), exit = max(d[["exit"]]),
    fitted = sum(d[["fitted"]]), cached_first = d[["cached"]][1L],
    failed = sum(d[["failed"]]), row_errors = sum(d[["row_errors"]]),
    unit_h = sum(d[["unit_h"]]), cpu_h = sum(d[["cpu_h"]]),
    wall_h = sum(d[["wall_h"]]), max_rss_gb = max(d[["max_rss_gb"]]),
    stringsAsFactors = FALSE)))
if (is.null(steps)) {
  cat("no <step>.time under", file.path(dir, "steps"), "\n")
} else {
  steps[["cpu_over_unit"]] <- ifelse(steps[["unit_h"]] > 0,
                                     steps[["cpu_h"]] / steps[["unit_h"]], NA)
  pr <- steps
  for (v in c("unit_h", "cpu_h", "wall_h", "max_rss_gb", "cpu_over_unit")) {
    pr[[v]] <- fmt(pr[[v]], 2L)
  }
  print(pr, row.names = FALSE)
  cat("\nmax_rss_gb: the largest resident set of one process of the step (GNU",
      "time reports the largest\nof the children it waited for, and",
      "mclapply waits for every worker); cpu_over_unit: the\nprocessor",
      "time of the step over the sum of its unit times, near 1 when the",
      "units run single\nthreaded (the setup of the R session adds a",
      "little).\n")
  write_csv(steps, "steps.csv")
}

## --- health of the fits at n = top-n ---------------------------------------

section(sprintf("Health of the fits at n = %d", top_n))

cv <- pil[["convergence"]]
if (!is.null(cv) && all(c("chosen", "cut.at.min") %in% names(cv))) {
  cv <- cv[cv[["n"]] == top_n & cv[["chosen"]] %in% TRUE, , drop = FALSE]
  if (nrow(cv)) {
    cv[["cut"]] <- cv[["cut.at.min"]] %in% TRUE
    h <- stats::aggregate(cut ~ method + cell, cv,
                          function(v) sprintf("%d/%d", sum(v), length(v)))
    cat("an engine cut at the chosen penalty of the chosen J (cut.at.min), fits:\n")
    print(stats::reshape(h, idvar = "method", timevar = "cell",
                         direction = "wide"), row.names = FALSE)
  }
}
kc <- all_units[all_units[["n"]] == top_n & all_units[["method"]] %in% k_methods &
                  !is.na(all_units[["kcheck_low"]]), , drop = FALSE]
if (nrow(kc)) {
  kc[["flag"]] <- kc[["kcheck_low"]] > 0L
  h <- stats::aggregate(flag ~ method + cell, kc,
                        function(v) sprintf("%d/%d", sum(v), length(v)))
  cat("\nmgcv::k.check flags some smooth (k-index < 1 and p < 0.05), fits:\n")
  print(stats::reshape(h, idvar = "method", timevar = "cell",
                       direction = "wide"), row.names = FALSE)
}
cat("\nwritten to ", outdir, "\n", sep = "")
