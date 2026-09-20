## wafc/scripts/05-sondagem-aplicacao.R -- step E6.1a: the empirical scouting
## of the candidate datasets for the application of E6. From the root of the
## repository:
##
##     Rscript wafc/scripts/05-sondagem-aplicacao.R [candidates] [n_max] [seed]
##
## Defaults: candidates = "all" (bike, beijing, housing), n_max = 0 (no
## subsampling), seed = 20260920. A candidate can also be given as a comma
## separated list, e.g. "bike,housing".
##
## The question this script answers is NOT which dataset the WAFC predicts
## best. It is the one of docs/aplicacao-candidatas.md: does the dataset
## carry an effect whose modulation is inhomogeneous? If beta_l(u) is smooth
## the spline is optimal (Xue and Yang, 2006, prove the univariate optimal
## rate under alpha_ls in C^{p+1}[0,1]) and the argument of the article
## collapses. So the script reports three things per candidate, and the
## verdict of the .md weighs the three together:
##
##   (1) out-of-sample prediction error of cv.wafc() against mgcv::gam with
##       s(u_m, by = x_l), on one train/test split, at two spline basis
##       sizes, because a WAFC win at k = 10 alone would be a win over the
##       competitor's basis size and not over its smoothness class;
##   (2) the levelwise energy of theta_hat by block, which is the direct
##       reading of the Besov exponent: with the B^{s'}_{2,2} convention
##       ||theta_{j.}||_2 is of order 2^{-j s'}, so the slope of
##       log2 ||theta_{j.}||_2 on j estimates -s'. A candidate whose active
##       blocks read s' close to 1/2 is the inhomogeneous case of D27; one
##       that reads s' >= 3/2 is the smooth case, where the spline wins;
##   (3) a localization index of the estimated component, the share of the
##       total variation carried by the largest 5 percent of the increments
##       on the plotting grid. A smooth function spreads its variation and
##       scores near 0.05 to 0.15; a component with a corner or a jump
##       concentrates it. The same index computed on the mgcv component is
##       the within-dataset floor, since that one is smooth by construction.
##
## The data are not versioned (docs/instrucoes.md, section 6, and .gitignore):
## they are downloaded into wafc/cache/data/ on first use, and every archive
## is checked against the SHA-256 recorded in wafc_sources below and in
## docs/aplicacao-candidatas.md. Nothing else in the repository reads them.
##
## Outputs, all under wafc/cache/ and none versioned: one PNG of the
## estimated components per candidate and 05-sondagem.rds with every number
## the .md quotes.

local({
  cand <- c("wafc/R/load.R", "../R/load.R", "R/load.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/R/load.R from ", getwd())
  source(hit[1L])
})

args <- commandArgs(trailingOnly = TRUE)
which_cand <- if (length(args) >= 1L) args[1L] else "all"
n_max <- if (length(args) >= 2L) as.integer(args[2L]) else 0L
seed <- if (length(args) >= 3L) as.integer(args[3L]) else 20260920L
stage <- if (length(args) >= 4L) args[4L] else "full"

cache_dir <- local({
  cand <- c("wafc/cache", "../cache", "cache")
  hit <- cand[dir.exists(dirname(cand)) | dir.exists(cand)]
  if (length(hit) == 0L) stop("Could not find wafc/cache/ from ", getwd())
  hit[1L]
})
data_dir <- file.path(cache_dir, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------------------
## Sources
## ---------------------------------------------------------------------------

## One entry per candidate: where the archive comes from, the SHA-256 of the
## file as downloaded on 2026-09-20, and the member of the archive the
## preparation reads. The checksum is the reproducibility contract: the UCI
## archive rewrites its zip files now and then, and a silent change of the
## data would be a silent change of every number below.
wafc_sources <- list(
  bike = list(
    url = "https://archive.ics.uci.edu/static/public/275/bike+sharing+dataset.zip",
    archive = "bike.zip",
    sha256 = "b70182d0d0508e9abbb79306ce5c0cec34869000f8220175ac83d11dbe845401",
    unpack = "bike",
    member = "bike/hour.csv"),
  beijing = list(
    url = "https://archive.ics.uci.edu/static/public/501/beijing+multi+site+air+quality+data.zip",
    archive = "beijing.zip",
    sha256 = "b04da438b2f331ac0ffd45aebdfec0d20d2367feb5f6948c4b1f7ce1191e33c4",
    unpack = "beijing",
    member = paste0("beijing/prsa/PRSA_Data_20130301-20170228/",
                    "PRSA_Data_Dongsi_20130301-20170228.csv"),
    inner = "beijing/PRSA2017_Data_20130301-20170228.zip",
    inner_dir = "beijing/prsa"),
  housing = list(
    url = "http://lib.stat.cmu.edu/datasets/houses.zip",
    archive = "houses.zip",
    sha256 = "8b18f0a01cf9c99a65174d18fa582aa31971dfe55a26ad794f3299937c3708d7",
    unpack = "houses",
    member = "houses/cadata.txt"))

## Downloads the archive if it is not there, checks the digest, unpacks it
## once. Returns the path of the file the preparation reads.
wafc_fetch <- function(key) {
  src <- wafc_sources[[key]]
  arc <- file.path(data_dir, src[["archive"]])
  if (!file.exists(arc)) {
    cat(sprintf("  downloading %s\n", src[["url"]]))
    utils::download.file(src[["url"]], arc, mode = "wb", quiet = TRUE)
  }
  got <- wafc_sha256(arc)
  if (!is.na(src[["sha256"]]) && !identical(got, src[["sha256"]])) {
    warning("SHA-256 of ", src[["archive"]], " is ", got,
            ", not the ", src[["sha256"]], " recorded in this script and in ",
            "docs/aplicacao-candidatas.md. The numbers below are not the ",
            "ones of the document.", call. = FALSE)
  }
  member <- file.path(data_dir, src[["member"]])
  if (!file.exists(member)) {
    utils::unzip(arc, exdir = file.path(data_dir, src[["unpack"]]))
    if (!is.null(src[["inner"]])) {
      utils::unzip(file.path(data_dir, src[["inner"]]),
                   exdir = file.path(data_dir, src[["inner_dir"]]))
    }
  }
  if (!file.exists(member)) {
    stop("Could not find ", member, " after unpacking ", arc, ".",
         call. = FALSE)
  }
  list(path = member, sha256 = got, url = src[["url"]])
}

## SHA-256 of a file, by the tools available: digest if it is installed,
## otherwise the system utility. The scouting does not add a dependency to
## wafc/R/load.R for a checksum.
wafc_sha256 <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(path, algo = "sha256", file = TRUE))
  }
  out <- tryCatch(system2("sha256sum", shQuote(path), stdout = TRUE),
                  warning = function(w) NA_character_,
                  error = function(e) NA_character_)
  if (length(out) == 0L || is.na(out[1L])) return(NA_character_)
  sub("\\s.*$", "", out[1L])
}

## ---------------------------------------------------------------------------
## Preparation of each candidate
## ---------------------------------------------------------------------------
##
## Every preparation returns the mapping of docs/aplicacao-candidatas.md in
## the notation of E1.1: the response y, the linear covariates x (column 1 is
## the constant, so that beta_1(u) is the baseline surface) and the
## modulating covariates u. It also returns the units, used only by the
## plots and the tables.

prep_bike <- function() {
  f <- wafc_fetch("bike")
  b <- utils::read.csv(f[["path"]])
  ## Days since the first day of the series: the seasonal modulator. Kept on
  ## the day scale and not folded into the year, because the two years of
  ## the series differ in level and the fold would confound them.
  day <- as.integer(as.Date(b[["dteday"]]) - as.Date("2011-01-01"))
  ## log count: the response is a count over three orders of magnitude and
  ## the additive-in-u model is a model for the log rate. Reported as such.
  y <- log(b[["cnt"]])
  x <- cbind(one = 1, temp = b[["temp"]], hum = b[["hum"]],
             wind = b[["windspeed"]])
  u <- cbind(hour = b[["hr"]], day = day)
  list(key = "bike", label = "Bike sharing (Capital bikeshare, hourly)",
       y = y, x = x, u = u, source = f,
       yname = "log(count)",
       note = paste("hour is discrete with 24 distinct values, so a block",
                    "on it has rank at most 23 whatever J is"))
}

prep_beijing <- function() {
  f <- wafc_fetch("beijing")
  d <- utils::read.csv(f[["path"]])
  ## Relative humidity from temperature and dew point (Magnus form, the
  ## coefficients of Alduchov and Eskridge). The dataset gives the dew
  ## point, and humidity is the modulator the scouting wants.
  mg <- function(t) exp(17.625 * t / (243.04 + t))
  rh <- 100 * mg(d[["DEWP"]]) / mg(d[["TEMP"]])
  keep <- stats::complete.cases(d[["PM2.5"]], d[["TEMP"]], d[["PRES"]],
                                d[["DEWP"]], d[["WSPM"]]) &
    is.finite(rh) & rh > 0 & rh <= 100 & d[["PM2.5"]] > 0
  d <- d[keep, ]
  rh <- rh[keep]
  y <- log(d[["PM2.5"]])
  x <- cbind(one = 1, wind = d[["WSPM"]], temp = d[["TEMP"]],
             pres = d[["PRES"]] - 1000)
  u <- cbind(hour = d[["hour"]], rh = rh)
  list(key = "beijing", label = "Beijing air quality, Dongsi site (hourly)",
       y = y, x = x, u = u, source = f,
       yname = "log(PM2.5)",
       note = paste("one of the twelve sites; hourly series, so the",
                    "residuals are autocorrelated and the split below is",
                    "random, not by time"))
}

prep_housing <- function() {
  f <- wafc_fetch("housing")
  ## cadata.txt is the StatLib flat file: 27 lines of prose and then nine
  ## columns, in the order the file documents.
  v <- scan(f[["path"]], what = numeric(), skip = 27L, quiet = TRUE)
  m <- matrix(v, ncol = 9L, byrow = TRUE)
  colnames(m) <- c("value", "income", "age", "rooms", "bedrooms", "pop",
                   "hh", "lat", "lon")
  y <- log(m[, "value"])
  x <- cbind(one = 1,
             income = m[, "income"],
             rooms = log(m[, "rooms"] / m[, "hh"]),
             age = m[, "age"] / 10)
  u <- cbind(lat = m[, "lat"], lon = m[, "lon"])
  list(key = "housing", label = "California block groups, 1990 census",
       y = y, x = x, u = u, source = f,
       yname = "log(median house value)",
       note = paste("the response is top-coded at 500001 (",
                    sum(m[, "value"] >= 500001), " of ", nrow(m),
                    " block groups), and additivity in latitude and",
                    "longitude is a strong assumption on a spatial field",
                    sep = ""))
}

preps <- list(bike = prep_bike, beijing = prep_beijing,
              housing = prep_housing)

## ---------------------------------------------------------------------------
## Diagnostics
## ---------------------------------------------------------------------------

## Levelwise energy of theta_hat, block by block, and the Besov exponent it
## reads. With the B^{s'}_{2,2} convention ||theta_{j.}||_2 is of order
## 2^{-j s'}, so the slope of log2 ||theta_{j.}||_2 on j estimates -s'.
##
## Two readings are returned because the LASSO shrinks, and shrinkage falls
## hardest on the fine levels, which biases s' upwards, that is, towards the
## smooth case that would sink the argument: 'lasso' is the fit as it comes,
## 'relaxed' is the least squares refit on the selected support, which is
## the conservative reading of the same support.
wafc_level_energy <- function(fit, s, y, relaxed = TRUE) {
  design <- fit[["design"]]
  p <- design[["p"]]
  q <- design[["q"]]
  J <- design[["J"]]
  j0 <- design[["j0"]]
  b <- wafc_raw_coef(fit, s = s)[-1L, 1L]

  if (relaxed) {
    nz <- which(b != 0)
    unp <- design[["unpenalized"]]
    keep <- sort(unique(c(unp, nz)))
    if (length(keep) > 0L && length(keep) < length(b)) {
      Z <- as.matrix(design[["Z"]][, keep, drop = FALSE])
      ## No extra intercept column when a constant linear covariate is
      ## already among the unpenalized ones: the two would be collinear and
      ## the refit would be dropped for rank deficiency.
      add_int <- length(design[["constant"]]) == 0L
      if (add_int) Z <- cbind(`(Intercept)` = 1, Z)
      qrz <- qr(Z, LAPACK = FALSE)
      if (qrz[["rank"]] == ncol(Z)) {
        cf <- qr.coef(qrz, y)
        b[] <- 0
        b[keep] <- if (add_int) cf[-1L] else cf
      }
    }
  }

  ## Column j of the basis of U_m in the order wafc_colnames() writes them:
  ## levels j0 .. J-1, with 2^j translations each.
  lev_of <- function(Jm) rep(seq.int(j0, Jm - 1L), times = 2^seq.int(j0, Jm - 1L))

  out <- list()
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      idx <- design[["blocks"]][[wafc_block_name(design, l, m)]]
      bl <- b[idx]
      lv <- lev_of(J[m])
      nrm <- tapply(bl, lv, function(z) sqrt(sum(z^2)))
      nz <- tapply(bl, lv, function(z) sum(z != 0))
      j <- as.integer(names(nrm))
      tot <- sum(nrm^2)
      fine <- if (tot > 0) sum(nrm[j >= J[m] - 2L]^2) / tot else NA_real_
      ## slope over the levels that carry mass; needs at least three of them
      pos <- which(nrm > 0)
      shat <- NA_real_
      if (length(pos) >= 3L) {
        shat <- -unname(stats::coef(stats::lm(log2(nrm[pos]) ~ j[pos]))[2L])
      }
      out[[length(out) + 1L]] <- data.frame(
        x = design[["xnames"]][l], u = design[["unames"]][m],
        norm = sqrt(tot), nonzero = sum(bl != 0),
        fine.share = fine, s.hat = shat,
        levels = paste(sprintf("%.3f", nrm), collapse = "/"),
        nz.levels = paste(nz, collapse = "/"),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

## Share of the total variation of g on the grid carried by the largest 5
## percent of the increments. Near 0.05 to 0.15 for a function whose
## variation is spread, large for one with a corner or a jump.
wafc_localization <- function(g, frac = 0.05) {
  d <- abs(diff(g))
  tot <- sum(d)
  if (!is.finite(tot) || tot <= 0) return(NA_real_)
  k <- max(1L, ceiling(frac * length(d)))
  sum(sort(d, decreasing = TRUE)[seq_len(k)]) / tot
}

rmse <- function(a, b) sqrt(mean((a - b)^2))

## Resumable cache by unit, the scheme E4.1 asks of the compendium: the
## cross-validation over (J, lambda) costs hours on these sample sizes, and
## a rerun that adds one competitor must not pay for it again. Delete
## wafc/cache/05-cache-*.rds to force a recomputation.
cache_get <- function(name, expr) {
  path <- file.path(cache_dir, sprintf("05-cache-%s.rds", name))
  if (file.exists(path)) {
    cat(sprintf("  [cache] %s\n", basename(path)))
    return(readRDS(path))
  }
  v <- expr
  saveRDS(v, path)
  v
}

## The index has no absolute meaning, so it is anchored on the shapes of
## wafc_component(), which are the vocabulary D27 uses: 'sine' and 'cubic'
## are the smooth case (s' = 3/2 in the scenario of the pilot) and 'bumps',
## 'blocks' and 'heavisine' the inhomogeneous one (s' = 1/2).
wafc_localization_reference <- function(G = 512L) {
  gr <- seq(0, 1, length.out = G)
  nm <- c("sine", "cubic", "heavisine", "bumps", "blocks")
  v <- vapply(nm, function(z) wafc_localization(wafc_component(z)(gr)), 0)
  v
}

## The gam of wafc_competitor() takes one basis size for every smooth, which
## caps every modulator at the coarsest one: the hour of the day has 24
## distinct values, so no smooth of that fit can exceed k = 23 however many
## distinct values latitude or the day index have. That cap is what makes
## the comparison above one of dimension and not one of basis, and the
## comparison of dimension is the one that decides whether the wavelet basis
## is doing anything a spline could not.
##
## This fits the same model with the basis size chosen per modulator,
## matched to the 2^J of the resolution level the cross-validation selected,
## and truncated at what the data admit. It is fitted by mgcv::bam with
## fREML and discretized covariates, because at these dimensions gam() does
## not finish. That is a change of fitting algorithm plus a binning of the
## covariates, not a change of estimator, and it is labelled as such
## wherever the number appears; on the candidate where both are affordable
## the two agree to the fourth decimal.
fit_gam_matched <- function(x, u, y, k, xn, un) {
  p <- ncol(x); q <- ncol(u)
  dat <- as.data.frame(cbind(x, u))
  names(dat) <- c(xn, un)
  dat[["y"]] <- y
  terms <- character(0)
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      terms <- c(terms, sprintf("s(%s, k = %d, by = %s)", un[m], k[m], xn[l]))
    }
  }
  fo <- stats::as.formula(paste("y ~ 0 +", paste(c(xn, terms), collapse = " + ")))
  fit <- mgcv::bam(fo, data = dat, method = "fREML", discrete = TRUE,
                   select = TRUE)
  cc <- stats::coef(fit)[xn]
  term_at <- function(newu) {
    nd <- as.data.frame(cbind(matrix(1, nrow(newu), p), newu))
    names(nd) <- c(xn, un)
    stats::predict(fit, newdata = nd, type = "terms")
  }
  list(fit = fit, k = k, cc = cc,
       predict = function(newx, newu) {
         tm <- term_at(newu)
         b <- matrix(rep(cc, each = nrow(newu)), nrow(newu), p)
         i <- 0L
         for (l in seq_len(p)) for (m in seq_len(q)) {
           i <- i + 1L
           b[, l] <- b[, l] + tm[, wafc_gam_term(colnames(tm), i)]
         }
         as.numeric(rowSums(newx * b))
       },
       g = function(grid) {
         tm <- term_at(grid)
         g <- vector("list", p * q); dim(g) <- c(p, q)
         i <- 0L
         for (l in seq_len(p)) for (m in seq_len(q)) {
           i <- i + 1L
           g[[l, m]] <- as.numeric(tm[, wafc_gam_term(colnames(tm), i)])
           g[[l, m]] <- g[[l, m]] - mean(g[[l, m]])
         }
         g
       })
}

## ---------------------------------------------------------------------------
## One candidate
## ---------------------------------------------------------------------------

run_one <- function(key) {
  cat(sprintf("\n=== %s ===\n", key))
  dat <- preps[[key]]()
  y <- dat[["y"]]
  x <- dat[["x"]]
  u <- dat[["u"]]
  n_full <- length(y)

  set.seed(seed)
  if (n_max > 0L && n_max < n_full) {
    take <- sort(sample.int(n_full, n_max))
    y <- y[take]; x <- x[take, , drop = FALSE]; u <- u[take, , drop = FALSE]
  }
  n <- length(y)
  p <- ncol(x); q <- ncol(u)
  tr <- sample.int(n, floor(0.7 * n))
  te <- setdiff(seq_len(n), tr)
  cat(sprintf("  %s\n  n = %d (of %d), p = %d, q = %d; train %d, test %d\n",
              dat[["label"]], n, n_full, p, q, length(tr), length(te)))
  cat(sprintf("  note: %s\n", dat[["note"]]))
  cat(sprintf("  sd(y) = %.4f\n", stats::sd(y)))

  ## ---- WAFC, defaults of cv.wafc() ---------------------------------------
  ## One evaluation table for the whole candidate, passed explicitly, which
  ## is D31: the basis is fixed and evaluated by table because the same
  ## basis is rebuilt once per candidate resolution level and once per fold.
  ## Nothing else in the call departs from the defaults of cv.wafc().
  wtab <- WaveBased::wtable(family = "Daublets", filter.size = 8L,
                            prec.wavelet = 30L, check = FALSE)
  t0 <- proc.time()[["elapsed"]]
  cv <- cache_get(paste0("cv-", key),
                  cv.wafc(x[tr, ], u[tr, ], y[tr], trace = TRUE,
                          wavelet.table = wtab))
  t_wafc <- proc.time()[["elapsed"]] - t0
  cat(sprintf("  cv.wafc: J = %d, lambda.min = %.6f, %d nonzero, %.1f s\n",
              cv[["J.min"]], cv[["lambda.min"]], cv[["nzero.min"]], t_wafc))
  print(cv[["cvtab"]])

  pred <- list()
  for (s in c("lambda.min", "lambda.1se")) {
    pred[[paste0("wafc.", s)]] <-
      as.numeric(predict(cv, x[te, ], u[te, ], s = s))
  }

  ## ---- mgcv, the competitor of L2d, at two basis sizes -------------------
  ## mgcv refuses a smooth whose basis is larger than the number of distinct
  ## values of its covariate, which binds whenever a modulator is discrete
  ## (the hour of the day has 24 values). The large size is therefore the
  ## largest one the data admit, capped at 30.
  k_max <- min(vapply(seq_len(q), function(m) length(unique(u[tr, m])), 0L)) - 1L
  k_big <- max(10L, min(30L, k_max))
  cat(sprintf("  gam basis sizes: 10 and %d (cap from the data: %d)\n",
              k_big, k_max))
  gams <- list()
  for (k in unique(c(10L, k_big))) {
    t0 <- proc.time()[["elapsed"]]
    g <- tryCatch(cache_get(sprintf("gam-%s-%d", key, k),
                            wafc_competitor("gam", x[tr, ], u[tr, ], y[tr],
                                            k = k)),
                  error = function(e) {
                    cat(sprintf("  gam k = %d failed: %s\n", k,
                                conditionMessage(e)))
                    NULL
                  })
    if (is.null(g)) next
    tg <- proc.time()[["elapsed"]] - t0
    gams[[as.character(k)]] <- g
    pred[[paste0("gam.k", k)]] <- predict(g, x[te, ], u[te, ])
    cat(sprintf("  gam k = %2d: %d of %d blocks kept, %.1f s\n",
                k, sum(g[["blocks"]]), length(g[["blocks"]]), tg))
  }

  ## ---- mgcv at the dimension of the selected WAFC sieve -------------------
  ## 2^J per modulator, truncated at the number of distinct values: the
  ## comparison that separates "the wavelet basis helps" from "the WAFC was
  ## simply allowed more columns".
  k_match <- pmin(2^cv[["J.min"]],
                  vapply(seq_len(q), function(m) length(unique(u[tr, m])), 0L) - 1L)
  gm <- tryCatch(cache_get(sprintf("gamm-%s", key),
                           fit_gam_matched(x[tr, ], u[tr, ], y[tr], k_match,
                                           colnames(x), colnames(u))),
                 error = function(e) {
                   cat(sprintf("  matched gam failed: %s\n",
                               conditionMessage(e)))
                   NULL
                 })
  if (!is.null(gm)) {
    pred[["gam.matched"]] <- gm[["predict"]](x[te, ], u[te, ])
    cat(sprintf("  gam matched: k = (%s), by bam/fREML/discrete\n",
                paste(k_match, collapse = ", ")))
  }

  ## ---- out-of-sample error -----------------------------------------------
  err <- data.frame(method = names(pred),
                    rmse = vapply(pred, function(z) rmse(z, y[te]), 0),
                    stringsAsFactors = FALSE)
  err[["rel"]] <- err[["rmse"]] / min(err[["rmse"]])
  rownames(err) <- NULL
  cat("\n  out-of-sample RMSE:\n")
  print(err, digits = 5)

  ## ---- structure ----------------------------------------------------------
  fit <- cv[["wafc.fit"]]
  blk <- wafc_blocks(fit, s = cv[["lambda.min"]])
  cat("\n  nonzero coefficients by block (lambda.min):\n")
  print(blk[["nonzero"]])
  cat("  block norms:\n")
  print(round(blk[["norm"]], 4))

  lev <- wafc_level_energy(fit, s = cv[["lambda.min"]], y = y[tr],
                           relaxed = TRUE)
  lev_raw <- wafc_level_energy(fit, s = cv[["lambda.min"]], y = y[tr],
                               relaxed = FALSE)
  lev[["s.hat.lasso"]] <- lev_raw[["s.hat"]]
  lev[["fine.share.lasso"]] <- lev_raw[["fine.share"]]
  cat("\n  levelwise energy of theta_hat, refit on the selected support\n")
  cat("  (norm/nonzero/fine.share/s.hat relaxed; s.hat.lasso is the shrunk reading)\n")
  print(lev[, c("x", "u", "norm", "nonzero", "fine.share", "s.hat",
                "s.hat.lasso", "levels")], digits = 3, row.names = FALSE)

  ## ---- components on a grid, and how localized they are -------------------
  G <- 512L
  grid <- matrix(0, G, q)
  for (m in seq_len(q)) {
    rg <- range(u[, m])
    grid[, m] <- seq(rg[1L], rg[2L], length.out = G)
  }
  g_wafc <- wafc_grid_components(fit, grid, s = cv[["lambda.min"]])
  g_big <- gams[[as.character(k_big)]]
  g_gam <- if (!is.null(g_big)) wafc_grid_components(g_big, grid) else NULL
  g_mat <- if (!is.null(gm)) gm[["g"]](grid) else NULL

  loc <- NULL
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      loc <- rbind(loc, data.frame(
        x = colnames(x)[l], u = colnames(u)[m],
        wafc = wafc_localization(g_wafc[[l, m]]),
        gam = if (is.null(g_gam)) NA_real_ else wafc_localization(g_gam[[l, m]]),
        gam.matched = if (is.null(g_mat)) NA_real_
                      else wafc_localization(g_mat[[l, m]]),
        sd.wafc = stats::sd(g_wafc[[l, m]]),
        sd.gam = if (is.null(g_gam)) NA_real_ else stats::sd(g_gam[[l, m]]),
        sd.gam.matched = if (is.null(g_mat)) NA_real_
                         else stats::sd(g_mat[[l, m]]),
        stringsAsFactors = FALSE))
    }
  }
  cat("\n  localization index (share of total variation in the top 5% of increments):\n")
  print(loc, digits = 3, row.names = FALSE)

  ## ---- figure -------------------------------------------------------------
  png_path <- file.path(cache_dir, sprintf("05-%s.png", key))
  grDevices::png(png_path, width = 1100, height = 260 * p, res = 110)
  op <- graphics::par(mfrow = c(p, q), mar = c(4, 4, 2.5, 1))
  for (l in seq_len(p)) {
    for (m in seq_len(q)) {
      yy <- g_wafc[[l, m]]
      if (!is.null(g_gam)) yy <- c(yy, g_gam[[l, m]])
      if (!is.null(g_mat)) yy <- c(yy, g_mat[[l, m]])
      graphics::plot(grid[, m], g_wafc[[l, m]], type = "l", lwd = 2,
                     ylim = range(yy, finite = TRUE),
                     xlab = colnames(u)[m], ylab = "",
                     main = sprintf("g[%s, %s]", colnames(x)[l],
                                    colnames(u)[m]))
      if (!is.null(g_gam)) {
        graphics::lines(grid[, m], g_gam[[l, m]], col = "firebrick", lty = 2,
                        lwd = 2)
      }
      if (!is.null(g_mat)) {
        graphics::lines(grid[, m], g_mat[[l, m]], col = "steelblue", lty = 3,
                        lwd = 2)
      }
      graphics::abline(h = 0, col = "grey70")
      if (l == 1L && m == 1L) {
        graphics::legend("topleft",
                         c("WAFC", sprintf("mgcv, k = %d", k_big),
                           sprintf("mgcv, k = (%s)",
                                   paste(k_match, collapse = ", "))),
                         lwd = 2, lty = c(1, 2, 3),
                         col = c("black", "firebrick", "steelblue"),
                         bty = "n", cex = 0.8)
      }
    }
  }
  graphics::par(op)
  grDevices::dev.off()
  cat(sprintf("\n  figure: %s\n", png_path))

  list(key = key, label = dat[["label"]], note = dat[["note"]],
       source = dat[["source"]], n = n, n_full = n_full, p = p, q = q,
       xnames = colnames(x), unames = colnames(u), yname = dat[["yname"]],
       sd_y = stats::sd(y), seed = seed, ntrain = length(tr),
       J.min = cv[["J.min"]], lambda.min = cv[["lambda.min"]],
       cvtab = cv[["cvtab"]], time.wafc = t_wafc,
       err = err, blocks = blk, levels = lev, localization = loc,
       k.big = k_big, k.matched = k_match,
       gam.blocks = lapply(gams, `[[`, "blocks"),
       gam.edf = lapply(gams, function(z) z[["extra"]][["edf"]]),
       png = png_path)
}

## ---------------------------------------------------------------------------
## The matched-dimension stage on its own
## ---------------------------------------------------------------------------
##
## The cross-validation over (J, lambda) costs hours at these sample sizes,
## and the run of 2026-09-20 was made before the matched-dimension gam
## existed, so its cv objects were not cached. This stage adds that one
## competitor to a finished run: it rebuilds the same split from the same
## seed, reads the selected J from 05-sondagem.rds, fits the matched gam and
## merges its error and its localization into the stored tables. A run from
## a clean cache never needs it, because stage "full" already fits it.
run_matched <- function(key, prev) {
  cat(sprintf("\n=== %s (matched-dimension stage) ===\n", key))
  dat <- preps[[key]]()
  y <- dat[["y"]]; x <- dat[["x"]]; u <- dat[["u"]]
  n_full <- length(y)
  set.seed(seed)
  if (n_max > 0L && n_max < n_full) {
    take <- sort(sample.int(n_full, n_max))
    y <- y[take]; x <- x[take, , drop = FALSE]; u <- u[take, , drop = FALSE]
  }
  n <- length(y); p <- ncol(x); q <- ncol(u)
  tr <- sample.int(n, floor(0.7 * n))
  te <- setdiff(seq_len(n), tr)
  stopifnot(identical(n, prev[["n"]]), identical(length(tr), prev[["ntrain"]]))

  k_match <- pmin(2^prev[["J.min"]],
                  vapply(seq_len(q), function(m) length(unique(u[tr, m])), 0L) - 1L)
  t0 <- proc.time()[["elapsed"]]
  gm <- cache_get(sprintf("gamm-%s", key),
                  fit_gam_matched(x[tr, ], u[tr, ], y[tr], k_match,
                                  colnames(x), colnames(u)))
  cat(sprintf("  gam matched: k = (%s) by bam/fREML/discrete, %.1f s\n",
              paste(k_match, collapse = ", "),
              proc.time()[["elapsed"]] - t0))

  err <- prev[["err"]]
  err <- err[err[["method"]] != "gam.matched", c("method", "rmse")]
  err <- rbind(err, data.frame(method = "gam.matched",
                               rmse = rmse(gm[["predict"]](x[te, ], u[te, ]),
                                           y[te])))
  err[["rel"]] <- err[["rmse"]] / min(err[["rmse"]])
  rownames(err) <- NULL
  cat("\n  out-of-sample RMSE:\n")
  print(err, digits = 5)

  G <- 512L
  grid <- matrix(0, G, q)
  for (m in seq_len(q)) {
    rg <- range(u[, m])
    grid[, m] <- seq(rg[1L], rg[2L], length.out = G)
  }
  g_mat <- gm[["g"]](grid)
  loc <- prev[["localization"]]
  loc[["gam.matched"]] <- vapply(seq_len(nrow(loc)), function(i) {
    l <- match(loc[["x"]][i], colnames(x)); m <- match(loc[["u"]][i], colnames(u))
    wafc_localization(g_mat[[l, m]])
  }, 0)
  loc[["sd.gam.matched"]] <- vapply(seq_len(nrow(loc)), function(i) {
    l <- match(loc[["x"]][i], colnames(x)); m <- match(loc[["u"]][i], colnames(u))
    stats::sd(g_mat[[l, m]])
  }, 0)
  cat("\n  localization index (share of total variation in the top 5% of increments):\n")
  print(loc, digits = 3, row.names = FALSE)

  prev[["err"]] <- err
  prev[["localization"]] <- loc
  prev[["k.matched"]] <- k_match
  prev
}

## ---------------------------------------------------------------------------
## One level above the default grid
## ---------------------------------------------------------------------------
##
## For two of the three candidates the cross-validation selected the largest
## J the default grid offers, with the curve still falling: the grid
## 2:ceiling(log2 n / 2) of cv.wafc() was binding, and the comparison above
## is therefore of a WAFC that was not allowed to go as deep as it wanted.
## This stage fits the level above the ceiling and reports what it changes.
## Five folds instead of ten, which is a departure from the defaults made
## for cost and labelled wherever the number appears: the point is the
## direction and the size of the change, not a tuned fit.
run_deep <- function(key, prev) {
  cat(sprintf("\n=== %s (one level above the default grid) ===\n", key))
  dat <- preps[[key]]()
  y <- dat[["y"]]; x <- dat[["x"]]; u <- dat[["u"]]
  n_full <- length(y)
  set.seed(seed)
  if (n_max > 0L && n_max < n_full) {
    take <- sort(sample.int(n_full, n_max))
    y <- y[take]; x <- x[take, , drop = FALSE]; u <- u[take, , drop = FALSE]
  }
  n <- length(y); q <- ncol(u)
  tr <- sample.int(n, floor(0.7 * n))
  te <- setdiff(seq_len(n), tr)
  stopifnot(identical(n, prev[["n"]]), identical(length(tr), prev[["ntrain"]]))

  Jtop <- max(prev[["cvtab"]][["J"]])
  Jnew <- Jtop + 1L
  wtab <- WaveBased::wtable(family = "Daublets", filter.size = 8L,
                            prec.wavelet = 30L, check = FALSE)
  t0 <- proc.time()[["elapsed"]]
  cvd <- cache_get(sprintf("deep-%s-%d", key, Jnew),
                   cv.wafc(x[tr, ], u[tr, ], y[tr], J = Jnew, nfolds = 5L,
                           trace = TRUE, wavelet.table = wtab))
  cat(sprintf("  J = %d, 5 folds: lambda.min = %.6f, %d nonzero, %.1f s\n",
              Jnew, cvd[["lambda.min"]], cvd[["nzero.min"]],
              proc.time()[["elapsed"]] - t0))

  r <- data.frame(
    s = c("lambda.min", "lambda.1se"),
    rmse = vapply(c("lambda.min", "lambda.1se"), function(ss) {
      rmse(as.numeric(predict(cvd, x[te, ], u[te, ], s = ss)), y[te])
    }, 0), stringsAsFactors = FALSE)
  base <- prev[["err"]]
  cat(sprintf("\n  test RMSE at J = %d (5 folds): %.5f / %.5f\n",
              Jnew, r[["rmse"]][1L], r[["rmse"]][2L]))
  cat(sprintf("  against J = %d (10 folds): %.5f / %.5f\n", prev[["J.min"]],
              base[["rmse"]][base[["method"]] == "wafc.lambda.min"],
              base[["rmse"]][base[["method"]] == "wafc.lambda.1se"]))
  cat(sprintf("  against the matched-dimension gam at J = %d: %.5f\n",
              prev[["J.min"]],
              base[["rmse"]][base[["method"]] == "gam.matched"]))

  ## The spline has to follow the sieve up, or the comparison stops being
  ## one of basis and becomes one of dimension again, in the other
  ## direction this time.
  k_deep <- pmin(2^Jnew,
                 vapply(seq_len(q), function(m) length(unique(u[tr, m])), 0L) - 1L)
  t0 <- proc.time()[["elapsed"]]
  gmd <- cache_get(sprintf("gamm-deep-%s-%d", key, Jnew),
                   fit_gam_matched(x[tr, ], u[tr, ], y[tr], k_deep,
                                   colnames(x), colnames(u)))
  r_gmd <- rmse(gmd[["predict"]](x[te, ], u[te, ]), y[te])
  cat(sprintf("  against the gam matched to J = %d, k = (%s): %.5f  (%.1f s)\n",
              Jnew, paste(k_deep, collapse = ", "), r_gmd,
              proc.time()[["elapsed"]] - t0))

  lv <- wafc_level_energy(cvd[["wafc.fit"]], s = cvd[["lambda.min"]],
                          y = y[tr], relaxed = TRUE)
  cat("\n  levelwise energy at the deeper level:\n")
  print(lv[, c("x", "u", "norm", "nonzero", "fine.share", "s.hat")],
        digits = 3, row.names = FALSE)

  prev[["deep"]] <- list(J = Jnew, nfolds = 5L, rmse = r,
                         cvm = cvd[["cvm.min"]], nzero = cvd[["nzero.min"]],
                         levels = lv, k.matched = k_deep,
                         rmse.gam.matched = r_gmd)
  prev
}

## ---------------------------------------------------------------------------
## The blocked split
## ---------------------------------------------------------------------------
##
## Only for bike, and only because the deeper level put the WAFC ahead there
## by a route that a random split over an hourly series rewards: the gain sat
## in the blocks modulated by the day index, whose components are spiky, and
## a random split leaves test hours in the same days as training hours, so a
## fine component in the day index can carry a day specific level from one
## to the other. That is memorising, not adaptation, and it has to be ruled
## out before the edge is believed.
##
## The split here holds out whole weeks, drawn at random: test hours never
## share a day with a training hour, and the range of the day index is still
## covered, so nothing is extrapolated and the comparison stays about
## interpolation. Five folds, and the folds themselves are blocked by week
## for the same reason.
run_blocked <- function(key, prev) {
  cat(sprintf("\n=== %s (blocked split, whole weeks held out) ===\n", key))
  dat <- preps[[key]]()
  y <- dat[["y"]]; x <- dat[["x"]]; u <- dat[["u"]]
  q <- ncol(u)
  day_col <- match("day", colnames(u))
  if (is.na(day_col)) {
    stop("run_blocked() needs a modulator named \"day\".", call. = FALSE)
  }
  n <- length(y)
  week <- as.integer(u[, day_col]) %/% 7L
  uw <- sort(unique(week))
  set.seed(seed)
  te_w <- sample(uw, ceiling(0.3 * length(uw)))
  te <- which(week %in% te_w)
  tr <- setdiff(seq_len(n), te)
  cat(sprintf("  %d weeks, %d held out; train %d, test %d rows\n",
              length(uw), length(te_w), length(tr), length(te)))
  ## folds blocked by week as well, so that the penalty is not chosen on a
  ## partition looser than the one it is judged by
  wtr <- week[tr]
  fold_of_week <- stats::setNames(
    sample(rep_len(seq_len(5L), length(unique(wtr)))), sort(unique(wtr)))
  foldid <- unname(fold_of_week[as.character(wtr)])

  Jgrid <- seq(2L, max(prev[["cvtab"]][["J"]]) + 1L)
  wtab <- WaveBased::wtable(family = "Daublets", filter.size = 8L,
                            prec.wavelet = 30L, check = FALSE)
  t0 <- proc.time()[["elapsed"]]
  cvb <- cache_get(sprintf("blocked-%s", key),
                   cv.wafc(x[tr, ], u[tr, ], y[tr], J = Jgrid,
                           foldid = foldid, trace = TRUE,
                           wavelet.table = wtab))
  cat(sprintf("  cv.wafc on the blocked folds: J = %d, %d nonzero, %.1f s\n",
              cvb[["J.min"]], cvb[["nzero.min"]],
              proc.time()[["elapsed"]] - t0))
  print(cvb[["cvtab"]])

  out <- data.frame(
    method = c("wafc.lambda.min", "wafc.lambda.1se"),
    rmse = vapply(c("lambda.min", "lambda.1se"), function(ss) {
      rmse(as.numeric(predict(cvb, x[te, ], u[te, ], s = ss)), y[te])
    }, 0), stringsAsFactors = FALSE)
  for (J in unique(c(cvb[["J.min"]], max(Jgrid)))) {
    k_b <- pmin(2^J,
                vapply(seq_len(q), function(m) length(unique(u[tr, m])), 0L) - 1L)
    gb <- cache_get(sprintf("blocked-gam-%s-%d", key, J),
                    fit_gam_matched(x[tr, ], u[tr, ], y[tr], k_b,
                                    colnames(x), colnames(u)))
    out <- rbind(out, data.frame(
      method = sprintf("gam.matched.J%d(k=%s)", J, paste(k_b, collapse = ",")),
      rmse = rmse(gb[["predict"]](x[te, ], u[te, ]), y[te])))
  }
  out[["rel"]] <- out[["rmse"]] / min(out[["rmse"]])
  rownames(out) <- NULL
  cat("\n  test RMSE on the held-out weeks:\n")
  print(out, digits = 5)

  lv <- wafc_level_energy(cvb[["wafc.fit"]], s = cvb[["lambda.min"]],
                          y = y[tr], relaxed = TRUE)
  cat("\n  levelwise energy under the blocked split:\n")
  print(lv[, c("x", "u", "norm", "nonzero", "fine.share", "s.hat")],
        digits = 3, row.names = FALSE)

  prev[["blocked"]] <- list(J = cvb[["J.min"]], Jgrid = Jgrid, err = out,
                            cvtab = cvb[["cvtab"]], levels = lv,
                            nweeks = length(uw), ntest = length(te))
  prev
}

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------

keys <- if (identical(which_cand, "all")) names(preps) else {
  k <- trimws(strsplit(which_cand, ",", fixed = TRUE)[[1L]])
  bad <- setdiff(k, names(preps))
  if (length(bad) > 0L) {
    stop("Unknown candidate(s): ", paste(bad, collapse = ", "),
         ". Known: ", paste(names(preps), collapse = ", "), ".", call. = FALSE)
  }
  k
}

cat(sprintf("WAFC application scouting (E6.1a): %s; n_max = %s, seed = %d\n",
            paste(keys, collapse = ", "),
            if (n_max > 0L) format(n_max) else "no cap", seed))

loc_ref <- wafc_localization_reference()
cat("\nlocalization index on the shapes of wafc_component(), for scale:\n")
print(round(loc_ref, 3))

out_path <- file.path(cache_dir, "05-sondagem.rds")

previous <- if (stage %in% c("matched", "deep", "blocked")) {
  if (!file.exists(out_path)) {
    stop("stage \"", stage, "\" needs a finished run in ", out_path, ".",
         call. = FALSE)
  }
  readRDS(out_path)
} else NULL

res <- if (is.null(previous)) list() else previous[["results"]]
for (k in keys) {
  res[[k]] <- tryCatch(
    if (identical(stage, "matched")) run_matched(k, previous[["results"]][[k]])
    else if (identical(stage, "deep")) run_deep(k, previous[["results"]][[k]])
    else if (identical(stage, "blocked")) run_blocked(k, previous[["results"]][[k]])
    else run_one(k),
    error = function(e) {
      cat(sprintf("  FAILED: %s\n", conditionMessage(e)))
      structure(list(key = k, error = conditionMessage(e)),
                class = "wafc_scout_error")
    })
}

saveRDS(list(results = res, seed = seed, n_max = n_max,
             localization.reference = loc_ref,
             sources = wafc_sources,
             sessionInfo = utils::sessionInfo()),
        out_path)

cat("\n=== summary ===\n")
for (k in names(res)) {
  r <- res[[k]]
  if (inherits(r, "wafc_scout_error")) {
    cat(sprintf("%-8s FAILED: %s\n", k, r[["error"]]))
    next
  }
  best_gam <- min(r[["err"]][["rmse"]][grepl("^gam", r[["err"]][["method"]])])
  w <- r[["err"]][["rmse"]][r[["err"]][["method"]] == "wafc.lambda.min"]
  act <- r[["levels"]][r[["levels"]][["nonzero"]] > 0, , drop = FALSE]
  cat(sprintf(paste0("%-8s n = %5d, J = %d | RMSE wafc %.4f vs gam %.4f ",
                     "(ratio %.3f) | active blocks %d/%d | median s.hat %.2f\n"),
              k, r[["n"]], r[["J.min"]], w, best_gam, w / best_gam,
              nrow(act), nrow(r[["levels"]]),
              stats::median(act[["s.hat"]], na.rm = TRUE)))
}
cat(sprintf("\nwritten: %s\n", out_path))
