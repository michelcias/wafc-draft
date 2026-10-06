## wafc/scripts/10-sondagem-aplicacao-b.R -- step E6.1b: the scouting of the
## application after decision D44. From the root of the repository:
##
##     Rscript wafc/scripts/10-sondagem-aplicacao-b.R [bases] [parts] [ncores] [nsplit] [n_max]
##
## Defaults: bases = "all", parts = "fit,report", ncores = 6, nsplit = 20,
## (parts: "fit", "extra" and "report"; see the part "extra" below),
## n_max = 0 (no subsampling). 'bases' is a comma separated list of the
## keys of wafc_bases below; "all" is every base whose data are already in
## wafc/cache/data/ or can be fetched from a source recorded here, except
## the base "marylebone.ukair" of step E6.1c, which is named explicitly.
##
## Environment: E61B_OUT (the folder of the outputs, default
## wafc/cache/e61b), E61B_REF (the commit the code is taken from when the
## snapshot does not exist yet, default b0ea096, the HEAD at which step
## E6.1b was catalogued; step E3.1 changed defaults after it, and every
## call here spells its arguments out), E61B_JGRID (the grid of J of
## the WAFC, default "2:8", the grid of decision D34), E61B_NFOLDS (folds
## of the cross-validation inside each training sample, default 10).
##
## What changed since step E6.1a (docs/aplicacao-candidatas.md, sections 1
## to 8) and why this is a second script and not a stage of the first.
## Decision D44 moved the thesis of the article from "the WAFC beats the
## spline by adaptation" to "the WAFC recovers the structure, with a
## prediction error competitive with the tuned spline", and changed the
## estimator: the WAFC is now the block LASSO in the balanced form (free
## levels, the weights sqrt(|G|) of grpreg, fine chunks of b_n to 2b_n - 1
## columns, the coarse levels of a block in one chunk; step E2.5e) followed
## by a threshold on the block norms with the rule cv1se (decision D45).
## The spline of the comparison is mgcv with the basis dimension chosen on
## the grid of decision D41 by REML and by GCV (decision D46). Three
## lessons of E6.1a are built in: the comparison is at a tuned dimension,
## never at k = 10; the partition is by block wherever the data depend
## (weeks of the series, tiles of the map; section 4.4 of the .md); and one
## partition is not evidence, so the partition is repeated and every number
## comes with a standard error between partitions.
##
## The structure is read from the same repetitions: the training sample of
## each partition is a subsample of 70 percent of the blocks, and the
## frequency with which a block (l, m) is kept by the thresholded fit over
## the partitions is its stability. Beside it, the same frequency for the
## smooths of mgcv with effective degrees of freedom above 0.1 (select =
## TRUE shrinks a null smooth towards zero but never to zero exactly) and
## the median effective degrees of freedom.
##
## The code comes from a snapshot. Step E3.1 edits wafc/R/ while this runs,
## so the functions are loaded from a copy of wafc/R extracted from a
## commit by git archive into <E61B_OUT>/snap/, once; the hash is recorded
## in snap/HEAD.txt and in every output, and the copy is checked against
## the commit, file by file, at every start. Nothing of wafc/R/ in the
## working tree is read.
##
## Every call to the methods spells out its arguments, the ones equal to
## the defaults included, so that a change of a default in E3.1 cannot move
## a number here.
##
## Outputs, none versioned: <E61B_OUT>/units/<base>-<split>.rds, one per
## partition (the run resumes from them), <E61B_OUT>/e61b-summary.rds with
## every table of the report, <E61B_OUT>/e61b-<base>.png with the
## components, and the log printed to the standard output.

## ---------------------------------------------------------------------------
## Paths, arguments, snapshot of the code
## ---------------------------------------------------------------------------

root <- local({
  cand <- c(".", "..", "../..")
  hit <- cand[file.exists(file.path(cand, "wafc", "R", "load.R"))]
  if (length(hit) == 0L) stop("Run from the root of wafc-draft.")
  normalizePath(hit[1L])
})
setwd(root)

args <- commandArgs(trailingOnly = TRUE)
which_bases <- if (length(args) >= 1L && nzchar(args[1L])) args[1L] else "all"
parts <- strsplit(if (length(args) >= 2L) args[2L] else "fit,report", ",")[[1L]]
ncores <- if (length(args) >= 3L) as.integer(args[3L]) else 6L
nsplit <- if (length(args) >= 4L) as.integer(args[4L]) else 20L
n_max <- if (length(args) >= 5L) as.integer(args[5L]) else 0L

out_dir <- Sys.getenv("E61B_OUT", file.path("wafc", "cache", "e61b"))
J_grid <- eval(parse(text = Sys.getenv("E61B_JGRID", "2:8")))
nfolds <- as.integer(Sys.getenv("E61B_NFOLDS", "10"))
unit_dir <- file.path(out_dir, "units")
dir.create(unit_dir, showWarnings = FALSE, recursive = TRUE)
extra_dir <- file.path(out_dir, "units-extra")
extra_methods <- strsplit(Sys.getenv("E61B_EXTRA", "gam.cv,gam.reml.wide"),
                          ",")[[1L]]
data_dir <- file.path("wafc", "cache", "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

## The snapshot: extracted once, from E61B_REF (b0ea096 by default), and
## checked against that commit at every start. A snapshot that does not
## match its commit stops the run: the numbers would not be the ones of the
## hash they carry.
snap_dir <- file.path(out_dir, "snap")
snap_head <- file.path(snap_dir, "HEAD.txt")
if (!file.exists(snap_head)) {
  ref <- Sys.getenv("E61B_REF", "b0ea096c5eacdd5c8c3ac4f1d8213da408901e39")
  hash <- system2("git", c("rev-parse", ref), stdout = TRUE)
  dir.create(snap_dir, showWarnings = FALSE, recursive = TRUE)
  st <- system(sprintf("git archive %s wafc/R | tar -x -C %s", hash,
                       shQuote(snap_dir)))
  if (st != 0L) stop("git archive failed.")
  writeLines(hash, snap_head)
}
code_hash <- readLines(snap_head)[1L]
local({
  fs <- sort(list.files(file.path(snap_dir, "wafc", "R"), pattern = "\\.R$"))
  ref <- system2("git", c("ls-tree", "--name-only", code_hash, "wafc/R/"),
                 stdout = TRUE)
  ref <- basename(ref)
  if (!identical(fs, sort(ref[grepl("\\.R$", ref)]))) {
    stop("The snapshot in ", snap_dir, " does not have the files of ",
         code_hash, ".")
  }
  for (f in fs) {
    a <- readLines(file.path(snap_dir, "wafc", "R", f), warn = FALSE)
    b <- system2("git", c("show", paste0(code_hash, ":wafc/R/", f)),
                 stdout = TRUE)
    if (!identical(a, b)) {
      stop("wafc/R/", f, " in the snapshot differs from ", code_hash, ".")
    }
  }
})
source(file.path(snap_dir, "wafc", "R", "load.R"))
if (!identical(normalizePath(wafc_r_dir),
               normalizePath(file.path(snap_dir, "wafc", "R")))) {
  stop("load.R resolved ", wafc_r_dir, ", not the snapshot.")
}
for (pk in c("grpreg", "mgcv")) {
  if (!requireNamespace(pk, quietly = TRUE)) stop("Needs the package ", pk, ".")
}
suppressPackageStartupMessages(library(parallel))

cat(sprintf("E6.1b: code of %s (snapshot in %s)\n", code_hash, snap_dir))
cat(sprintf("bases %s | parts %s | %d processes | %d partitions | J %s | %d folds%s\n",
            which_bases, paste(parts, collapse = ","), ncores, nsplit,
            paste(range(J_grid), collapse = ":"), nfolds,
            if (n_max > 0L) sprintf(" | n_max %d", n_max) else ""))

## ---------------------------------------------------------------------------
## Sources
## ---------------------------------------------------------------------------

## The three bases of step E6.1a keep the sources, digests and members of
## wafc/scripts/05-sondagem-aplicacao.R (copied here, so that this script
## does not source the other). The new candidates of step E6.1b (ii) are
## added with the digest of the file as downloaded; see
## docs/aplicacao-candidatas.md, sections 9 and following.
wafc_sources <- list(
  bike = list(
    url = "https://archive.ics.uci.edu/static/public/275/bike+sharing+dataset.zip",
    archive = "bike.zip",
    sha256 = "b70182d0d0508e9abbb79306ce5c0cec34869000f8220175ac83d11dbe845401",
    unpack = "bike", member = "bike/hour.csv"),
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
    unpack = "houses", member = "houses/cadata.txt"),
  ## The example data of openair 3.1.0 (Carslaw and Ropkins, 2012), read
  ## from the CRAN mirror on GitHub at the tag of that version; the file is
  ## the data/mydata.rda of the CRAN tarball. Downloaded on 2026-10-03.
  marylebone = list(
    url = "https://raw.githubusercontent.com/cran/openair/3.1.0/data/mydata.rda",
    archive = "marylebone/mydata.rda",
    sha256 = "f40da24c7855076874d04dfc13e612d53313d029fe3c24692bff73d19d92fc38",
    unpack = "marylebone", member = "marylebone/mydata.rda"),
  ## The 10-minute SCADA of 2017 of Kelmarsh wind farm (Plumley, 2022,
  ## Zenodo, doi:10.5281/zenodo.5841834, CC BY 4.0), the first full
  ## calendar year of the record; only the file of turbine 1 is extracted.
  ## The MD5 of the archive matches the one Zenodo publishes
  ## (c78263ee52ee0e48e2cb4bbaa1ba211a). Downloaded on 2026-10-03.
  kelmarsh = list(
    url = "https://zenodo.org/records/5841834/files/Kelmarsh_SCADA_2017_3083.zip?download=1",
    archive = "kelmarsh/Kelmarsh_SCADA_2017_3083.zip",
    sha256 = "528c942caec09bf02a25651cb13b6ecb13d4a09d287b2bb546f10adc614a14e3",
    unpack = "kelmarsh",
    member = "kelmarsh/Turbine_Data_Kelmarsh_1_2017-01-01_-_2018-01-01_228.csv",
    extract = "Turbine_Data_Kelmarsh_1_2017-01-01_-_2018-01-01_228.csv"))

wafc_sha256 <- function(path) {
  out <- tryCatch(system2("sha256sum", shQuote(path), stdout = TRUE),
                  warning = function(w) NA_character_,
                  error = function(e) NA_character_)
  if (length(out) == 0L || is.na(out[1L])) return(NA_character_)
  sub("\\s.*$", "", out[1L])
}

## The archive is downloaded only when it is not in wafc/cache/data/; the
## digest is checked and a mismatch is a warning, as in script 05.
wafc_fetch <- function(key) {
  src <- wafc_sources[[key]]
  arc <- file.path(data_dir, src[["archive"]])
  if (!file.exists(arc)) {
    dir.create(dirname(arc), showWarnings = FALSE, recursive = TRUE)
    cat(sprintf("  downloading %s\n", src[["url"]]))
    utils::download.file(src[["url"]], arc, mode = "wb", quiet = TRUE)
  }
  got <- wafc_sha256(arc)
  if (!is.na(src[["sha256"]]) && !identical(got, src[["sha256"]])) {
    warning("SHA-256 of ", src[["archive"]], " is ", got, ", not the ",
            src[["sha256"]], " recorded here.", call. = FALSE)
  }
  member <- file.path(data_dir, src[["member"]])
  if (!file.exists(member)) {
    utils::unzip(arc, files = src[["extract"]],
                 exdir = file.path(data_dir, src[["unpack"]]))
    if (!is.null(src[["inner"]])) {
      utils::unzip(file.path(data_dir, src[["inner"]]),
                   exdir = file.path(data_dir, src[["inner_dir"]]))
    }
  }
  if (!file.exists(member)) stop("Could not find ", member, ".")
  list(path = member, sha256 = got, url = src[["url"]])
}

## Step E6.1c (decision D57): the base "marylebone" again, from sources
## with a declared licence, as the base "marylebone.ukair"; several files,
## so it has its own list and fetch. The pollutants are the hourly series of
## the Automatic Urban and Rural Network at Marylebone Road (site MY1), in
## the yearly files that openair::importAURN() reads, published by Defra on
## uk-air.defra.gov.uk under the Open Government Licence (the data page of
## UK-AIR links version 2, the footer of the site states version 3); they
## carry no wind before about 2010. The wind is the 10 m wind of the ERA5
## reanalysis (Hersbach et al., 2020, QJRMS 146: 1999-2049; Copernicus
## Climate Change Service, CC BY 4.0, doi:10.24381/cds.adbb2d47), hourly
## and instantaneous, at the grid point nearest MY1 (51.5225 N, 0.1546 W;
## the API answers 51.5 N, 0.25 W, 35 m), served by the historical API of
## Open-Meteo (CC BY 4.0). Downloaded on 2026-10-05 with the permission of
## the author. The file of Open-Meteo is generated at each request, so a
## new download may differ in its digest without differing in its values.
ukair_sources <- list(
  dir = "marylebone.ukair",
  files = c(
    stats::setNames(lapply(1998:2005, function(y) {
      list(url = sprintf("https://uk-air.defra.gov.uk/openair/R_data/MY1_%d.RData", y))
    }), sprintf("MY1_%d.RData", 1998:2005)),
    list("era5-wind-my1-1998-2005.csv" = list(
      url = paste0("https://archive-api.open-meteo.com/v1/archive?",
                   "latitude=51.5225&longitude=-0.1546&start_date=1998-01-01&",
                   "end_date=2005-12-31&hourly=wind_speed_10m,wind_direction_10m&",
                   "models=era5&wind_speed_unit=ms&timezone=GMT&format=csv")))),
  sha256 = c(
    MY1_1998.RData = "503ee40cef24944022f7f78424d2092c81934753b8a429d83cc91cc74fa85412",
    MY1_1999.RData = "4007f98e51d3efdcdc39ef1b97e0119b36d25e04a6bcf65eec4a5ae0c45e085b",
    MY1_2000.RData = "f157a23b0a2e6221c2da10603c5a325fd35fea8e29dcd3428b27f9c4eb915715",
    MY1_2001.RData = "63b8283f629327bc36fb1fc34693b84f0b9989331d516867940217926e6a48fd",
    MY1_2002.RData = "6db6ebae8d683d69879c27f6102e267df38d800bb926a47756102d37b3589c7f",
    MY1_2003.RData = "2b32e149445b83d49125c3b8f3046ba8a1efb95a22fbfe09ef5a0b993a37ec6b",
    MY1_2004.RData = "fe88b12e189772938251e5a1031ffc0f6a7ecf33236f90c3913bb663591c3651",
    MY1_2005.RData = "a3ed2c09c8feee9085bf54863e4757ed13a8220693634a30da9f281787e91324",
    "era5-wind-my1-1998-2005.csv" =
      "e026ea3fa1ca311e599e87b52b9415499f419d7f87d05983974d5f9d64f25641"))

wafc_fetch_ukair <- function() {
  dir <- file.path(data_dir, ukair_sources[["dir"]])
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  fs <- names(ukair_sources[["files"]])
  got <- vapply(fs, function(f) {
    p <- file.path(dir, f)
    if (!file.exists(p)) {
      cat(sprintf("  downloading %s\n", ukair_sources[["files"]][[f]][["url"]]))
      utils::download.file(ukair_sources[["files"]][[f]][["url"]], p,
                           mode = "wb", quiet = TRUE)
    }
    h <- wafc_sha256(p)
    if (!identical(h, ukair_sources[["sha256"]][[f]])) {
      warning("SHA-256 of ", f, " is ", h, ", not the ",
              ukair_sources[["sha256"]][[f]], " recorded here.", call. = FALSE)
    }
    h
  }, "")
  list(path = dir, sha256 = got,
       url = vapply(ukair_sources[["files"]], `[[`, "", "url"))
}

## ---------------------------------------------------------------------------
## The bases
## ---------------------------------------------------------------------------
##
## Every preparation returns y, x (column 1 the constant), u, the block of
## each observation for the partition ('block', integer), a grid of each
## modulator for the components, and labels. The mapping of the three bases
## of E6.1a is the one of docs/aplicacao-candidatas.md, section 2.

## A grid of G points between the 0.5 and 99.5 percent quantiles of each
## modulator, the range in which every partition has data.
u_grid <- function(u, G = 256L) {
  apply(u, 2L, function(v) {
    r <- stats::quantile(v, c(0.005, 0.995), names = FALSE)
    seq(r[1L], r[2L], length.out = G)
  })
}

prep_bike <- function() {
  f <- wafc_fetch("bike")
  b <- utils::read.csv(f[["path"]])
  day <- as.integer(as.Date(b[["dteday"]]) - as.Date("2011-01-01"))
  y <- log(b[["cnt"]])
  x <- cbind(one = 1, temp = b[["temp"]], hum = b[["hum"]],
             wind = b[["windspeed"]])
  u <- cbind(hour = b[["hr"]], day = day)
  ## the block is the week, the unit of section 4.4 of the .md
  list(key = "bike", label = "Bike sharing (Capital bikeshare, hourly)",
       y = y, x = x, u = u, block = day %/% 7L, block.unit = "week",
       source = f, yname = "log(count)")
}

## The Dongsi site of the Beijing multi-site data, read once for the two
## mappings.
read_beijing <- function() {
  f <- wafc_fetch("beijing")
  d <- utils::read.csv(f[["path"]])
  mg <- function(t) exp(17.625 * t / (243.04 + t))
  rh <- 100 * mg(d[["DEWP"]]) / mg(d[["TEMP"]])
  keep <- stats::complete.cases(d[["PM2.5"]], d[["TEMP"]], d[["PRES"]],
                                d[["DEWP"]], d[["WSPM"]]) &
    is.finite(rh) & rh > 0 & rh <= 100 & d[["PM2.5"]] > 0
  d <- d[keep, ]
  d[["rh"]] <- rh[keep]
  d[["date"]] <- as.Date(sprintf("%04d-%02d-%02d", d[["year"]], d[["month"]],
                                 d[["day"]]))
  list(d = d, f = f)
}

prep_beijing <- function() {
  r <- read_beijing()
  d <- r[["d"]]
  y <- log(d[["PM2.5"]])
  x <- cbind(one = 1, wind = d[["WSPM"]], temp = d[["TEMP"]],
             pres = d[["PRES"]] - 1000)
  u <- cbind(hour = d[["hour"]], rh = d[["rh"]])
  wk <- as.integer(d[["date"]] - as.Date("2013-03-01")) %/% 7L
  list(key = "beijing", label = "Beijing air quality, Dongsi (hourly)",
       y = y, x = x, u = u, block = wk, block.unit = "week",
       source = r[["f"]], yname = "log(PM2.5)")
}

## The same site with the day of the year as the first modulator (step
## E6.1b (ii)). The winter heating season of Beijing and the North China
## Plain runs usually from 15 November to 15 March (Liang et al., 2015,
## Proc. R. Soc. A 471: 20150257, section 7), a nominal date (the effective
## dates move with the temperature; ibid., p. 16) at which the emissions,
## and so the level of PM2.5 and its response to temperature, change; the
## day of the year is periodic, which the periodized basis
## matches. Day of the year 1 is 1 January; 15 November is day 319 (320 in
## a leap year) and 15 March day 74 (75).
prep_beijing_heat <- function() {
  r <- read_beijing()
  d <- r[["d"]]
  y <- log(d[["PM2.5"]])
  x <- cbind(one = 1, wind = d[["WSPM"]], temp = d[["TEMP"]],
             pres = d[["PRES"]] - 1000)
  doy <- as.integer(format(d[["date"]], "%j")) - 1 + d[["hour"]] / 24
  u <- cbind(doy = doy, rh = d[["rh"]])
  wk <- as.integer(d[["date"]] - as.Date("2013-03-01")) %/% 7L
  out <- list(key = "beijing.heat",
              label = "Beijing air quality, Dongsi, by day of the year",
              y = y, x = x, u = u, block = wk, block.unit = "week",
              source = r[["f"]], yname = "log(PM2.5)",
              marks = list(doy = c(heat.on = 318, heat.off = 73)))
  out[["grid"]] <- cbind(doy = seq(0, 365, length.out = 256L),
                         rh = u_grid(u)[, "rh"])
  out
}

prep_housing <- function() {
  f <- wafc_fetch("housing")
  v <- scan(f[["path"]], what = numeric(), skip = 27L, quiet = TRUE)
  m <- matrix(v, ncol = 9L, byrow = TRUE)
  colnames(m) <- c("value", "income", "age", "rooms", "bedrooms", "pop",
                   "hh", "lat", "lon")
  y <- log(m[, "value"])
  x <- cbind(one = 1, income = m[, "income"],
             rooms = log(m[, "rooms"] / m[, "hh"]), age = m[, "age"] / 10)
  u <- cbind(lat = m[, "lat"], lon = m[, "lon"])
  ## The block is a tile of a quarter of a degree in latitude and longitude:
  ## the region of section 6 of the .md. A held-out tile leaves its
  ## latitude covered by the tiles to its east and west and its longitude
  ## by the ones to its north and south, so the additive model is not asked
  ## to extrapolate.
  tile <- as.integer(factor(paste(floor(u[, "lat"] / 0.25),
                                  floor(u[, "lon"] / 0.25))))
  list(key = "housing", label = "California block groups, 1990 census",
       y = y, x = x, u = u, block = tile, block.unit = "0.25-degree tile",
       source = f, yname = "log(median house value)")
}

## Marylebone Road, London, hourly, 1998-01-01 to 2005-06-23 (the mydata of
## openair). The model is the oxidant relation of Clapp and Jenkin (2001):
## the total oxidant OX = NO2 + O3 at a roadside site is close to linear in
## NOx, with intercept the regional background of oxidant and slope the
## local contribution to it, an estimate of the fraction of NOx emitted
## directly as NO2 (Carslaw and Beevers, 2005, section 2.4) that also
## carries the NO2 formed by 2NO + O2 and from HONO and is lower at night
## (Clapp and Jenkin, 2001, section 2.3). At Marylebone Road, Carslaw (2005,
## Atmospheric Environment 39: 4793-4802, sections 3.1 and 3.2, Fig. 3(a))
## estimated that fraction at about 10 percent by volume from 1997 to 2002,
## rising through 2002 and 2003 to about 23 percent at the end of 2003, a
## change he links to the particle filters fitted to the London buses (his
## 5 to 6 percent in 1997 and 17 percent in 2003 are the mean of the London
## sites, not this one); so beta_nox(u) estimates the primary NO2 fraction,
## modulated by the date (the documented change) and by the wind speed
## (dispersion and the mixing of ozone into the street).
## The wind direction is not a modulator: it is recorded in steps of 10
## degrees (38 values), the discrete case section 5 of the .md rules out.
prep_marylebone <- function() {
  f <- wafc_fetch("marylebone")
  e <- new.env()
  load(f[["path"]], envir = e)
  d <- as.data.frame(e[["mydata"]])
  keep <- stats::complete.cases(d[["ws"]], d[["nox"]], d[["no2"]], d[["o3"]])
  d <- d[keep, ]
  t0 <- as.POSIXct("1998-01-01 00:00:00", tz = "GMT")
  day <- as.numeric(difftime(d[["date"]], t0, units = "days"))
  y <- d[["no2"]] + d[["o3"]]
  x <- cbind(one = 1, nox = d[["nox"]] / 100)
  u <- cbind(day = day, ws = d[["ws"]])
  list(key = "marylebone", label = "Marylebone Road, London (hourly oxidant)",
       y = y, x = x, u = u, block = as.integer(floor(day)) %/% 7L,
       block.unit = "week", source = f, yname = "NO2 + O3 (ppb)",
       ## 2003-01-01, the year of the step Carslaw (2005) reports
       marks = list(day = 1826))
}

## Kelmarsh wind farm, turbine 1 (a Senvion MM92 of 2050 kW), 10-minute
## SCADA of 2017. The power curve has documented thresholds in the wind
## speed: the cut-in speed, below which the turbine does not produce, and
## the rated speed, above which the pitch control holds the power at the
## rated value. The density of the air raises the power below the rated
## speed and does not move it above, which is why the power curve is
## normalized to the density only below rated (IEC 61400-12-1); Lee, Ding,
## Genton and Xie (2015, JASA 110: 56-67) estimate the curve with the
## density as a covariate. With the ambient temperature as the inverse
## proxy of the density (no pressure is recorded), beta_temp(u) should be
## negative below the rated speed and drop to zero above it: the effect of
## a covariate that switches off at a known value of the modulator. The
## wind direction, recorded continuously, carries the wakes of the other
## five turbines, in sectors fixed by the layout. Intervals with downtime
## or curtailment (lost production to downtime and curtailment above zero,
## available capacity below the rated power) are removed, which is the
## usual filter of the power curve; the turbulence intensity is not used,
## because the standard deviation of the wind speed is missing in 73
## percent of the intervals.
prep_kelmarsh <- function() {
  f <- wafc_fetch("kelmarsh")
  d <- utils::read.csv(f[["path"]], skip = 9L, check.names = FALSE)
  ws <- d[["Wind speed (m/s)"]]
  wd <- d[["Wind direction (\u00b0)"]]
  pw <- d[["Power (kW)"]]
  tp <- d[["Nacelle ambient temperature (\u00b0C)"]]
  lost <- d[["Lost Production to Downtime and Curtailment Total (kWh)"]]
  cap <- d[["Available Capacity for Production (kW)"]]
  keep <- stats::complete.cases(ws, wd, pw, tp, lost, cap) & lost == 0 &
    cap >= 2050
  tm <- as.POSIXct(d[["# Date and time"]], tz = "UTC")[keep]
  day <- as.numeric(difftime(tm, as.POSIXct("2017-01-01", tz = "UTC"),
                             units = "days"))
  y <- pw[keep] / 1000
  x <- cbind(one = 1, temp = tp[keep] / 10)
  u <- cbind(ws = ws[keep], wd = wd[keep])
  list(key = "kelmarsh", label = "Kelmarsh wind farm, turbine 1 (10 min, 2017)",
       y = y, x = x, u = u, block = as.integer(floor(day)) %/% 7L,
       block.unit = "week", source = f, yname = "power (MW)",
       ## cut-in 3 m/s and rated 12.5 m/s of the MM92 data sheet [VERIFICAR]
       marks = list(ws = c(cut.in = 3, rated = 12.5)))
}

## Marylebone Road from UK-AIR and ERA5 (step E6.1c; sources above), with
## the mapping of prep_marylebone(): Y = NO2 + O3, X = (1, NOx / 100),
## U = (days since 1998-01-01, wind speed), the block the week. UK-AIR
## reports ug/m3 at 20 C and 1013 mb, NOx as NO2; the oxidant relation is
## molar, so the three are converted to ppb, the unit of the mydata of
## openair, with the factors of Defra for the UK Air Quality Archive (1 ppb
## of NO2, and of NOx as NO2, is 1.9125 ug/m3; of O3, 1.9957; "Conversion
## Factors Between ppb and ug m-3 and ppm and mgm-3", uk-air.defra.gov.uk,
## report cat06 0502160851). The date of UK-AIR is the beginning of the
## hour of the mean (the first hour of a year is 00:00); the wind of ERA5
## is instantaneous, so the hour gets the speed at its middle, by linear
## interpolation between t and t + 1 h. The window is the one of the mydata
## of openair, 1998-01-01 00:00 to 2005-06-23 12:00 GMT, so that the
## numbers compare with those of the base "marylebone" of step E6.1b; the
## files cover 1998 to 2005. The wind direction of ERA5 is read and not
## used, as in "marylebone".
prep_marylebone_ukair <- function() {
  f <- wafc_fetch_ukair()
  a <- do.call(rbind, lapply(1998:2005, function(yr) {
    e <- new.env()
    load(file.path(f[["path"]], sprintf("MY1_%d.RData", yr)), envir = e)
    as.data.frame(e[[sprintf("MY1_%d", yr)]])[, c("date", "NOXasNO2", "NO2",
                                                  "O3")]
  }))
  w <- utils::read.csv(file.path(f[["path"]], "era5-wind-my1-1998-2005.csv"),
                       skip = 3L, check.names = FALSE)
  wt <- as.POSIXct(w[[1L]], format = "%Y-%m-%dT%H:%M", tz = "GMT")
  ws <- stats::approx(as.numeric(wt), w[[2L]],
                      xout = as.numeric(a[["date"]]) + 1800)[["y"]]
  t0 <- as.POSIXct("1998-01-01 00:00:00", tz = "GMT")
  t1 <- as.POSIXct("2005-06-23 12:00:00", tz = "GMT")
  keep <- a[["date"]] >= t0 & a[["date"]] <= t1 &
    stats::complete.cases(ws, a[["NOXasNO2"]], a[["NO2"]], a[["O3"]])
  a <- a[keep, ]
  ws <- ws[keep]
  day <- as.numeric(difftime(a[["date"]], t0, units = "days"))
  y <- a[["NO2"]] / 1.9125 + a[["O3"]] / 1.9957
  x <- cbind(one = 1, nox = a[["NOXasNO2"]] / 1.9125 / 100)
  u <- cbind(day = day, ws = ws)
  list(key = "marylebone.ukair",
       label = "Marylebone Road, London, UK-AIR and ERA5 (hourly oxidant)",
       y = y, x = x, u = u, block = as.integer(floor(day)) %/% 7L,
       block.unit = "week", source = f, yname = "NO2 + O3 (ppb)",
       ## 2003-01-01, the year of the step Carslaw (2005) reports
       marks = list(day = 1826))
}

preps <- list(bike = prep_bike, beijing = prep_beijing,
              housing = prep_housing, beijing.heat = prep_beijing_heat,
              marylebone = prep_marylebone, kelmarsh = prep_kelmarsh,
              marylebone.ukair = prep_marylebone_ukair)

## The order fixes the seeds: the index of a base in this vector, and not
## its position in the command line, enters the seed of its partitions.
## A new base goes at the end, so that the seeds of the others do not move.
base_order <- c("bike", "beijing", "housing", "beijing.heat", "marylebone",
                "kelmarsh", "marylebone.ukair")

## One preparation per process, kept for the units of the same base.
data_cache <- new.env()
get_data <- function(key) {
  if (!exists(key, envir = data_cache, inherits = FALSE)) {
    d <- preps[[key]]()
    if (n_max > 0L && n_max < length(d[["y"]])) {
      set.seed(20261003L)
      take <- sort(sample.int(length(d[["y"]]), n_max))
      d[["y"]] <- d[["y"]][take]
      d[["x"]] <- d[["x"]][take, , drop = FALSE]
      d[["u"]] <- d[["u"]][take, , drop = FALSE]
      d[["block"]] <- d[["block"]][take]
    }
    if (is.null(d[["grid"]])) d[["grid"]] <- u_grid(d[["u"]])
    assign(key, d, envir = data_cache)
  }
  get(key, envir = data_cache, inherits = FALSE)
}

## ---------------------------------------------------------------------------
## One partition
## ---------------------------------------------------------------------------

## 30 percent of the blocks to the test sample, drawn at random, and the
## blocks of the training sample spread over the folds, so that the
## cross-validation inside the training sample is blocked as well (section
## 4.4 of the .md: a penalty chosen on folds that split a week would learn
## the level of that week).
make_split <- function(block, s, key) {
  set.seed(20261003L + 100000L * match(key, base_order) + s)
  ub <- sort(unique(block))
  te_b <- sample(ub, round(0.3 * length(ub)))
  te <- which(block %in% te_b)
  tr <- which(!(block %in% te_b))
  tb <- sort(unique(block[tr]))
  fold_of <- stats::setNames(sample(rep_len(seq_len(nfolds), length(tb))), tb)
  list(tr = tr, te = te, foldid = unname(fold_of[as.character(block[tr])]),
       ntest.blocks = length(te_b), nblocks = length(ub))
}

## The evaluation table of the basis, built once per process (decision D31).
wtab <- NULL
get_wtab <- function() {
  if (is.null(wtab)) {
    wtab <<- WaveBased::wtable(family = "Daublets", filter.size = 8L,
                               prec.wavelet = 30L, check = FALSE)
  }
  wtab
}

## Norm of each resolution level of each block, from the coefficients of a
## fit in the coordinates of its design: the levelwise energy of step E6.1a
## (section 3 of the .md), read on the block LASSO.
level_norms <- function(b, des) {
  out <- list()
  for (l in seq_len(des[["p"]])) {
    for (m in seq_len(des[["q"]])) {
      idx <- des[["blocks"]][[wafc_block_name(des, l, m)]]
      Jm <- des[["J"]][m]
      lv <- rep(0:(Jm - 1L), times = 2^(0:(Jm - 1L)))
      out[[wafc_block_name(des, l, m)]] <- as.numeric(tapply(b[idx], lv,
                                                             function(z) sqrt(sum(z^2))))
    }
  }
  out
}

## Share of the total variation of g on the grid carried by the largest 5
## percent of the increments (step E6.1a, section 3 of the .md).
localization <- function(g, frac = 0.05) {
  d <- abs(diff(g))
  tot <- sum(d)
  if (!is.finite(tot) || tot <= 0) return(NA_real_)
  k <- max(1L, ceiling(frac * length(d)))
  sum(sort(d, decreasing = TRUE)[seq_len(k)]) / tot
}

rmse <- function(a, b) sqrt(mean((a - b)^2))

run_unit <- function(key, s) {
  path <- file.path(unit_dir, sprintf("%s-%02d.rds", key, s))
  if (file.exists(path)) return(path)
  t_unit <- proc.time()[["elapsed"]]
  d <- get_data(key)
  sp <- make_split(d[["block"]], s, key)
  tr <- sp[["tr"]]
  te <- sp[["te"]]
  xtr <- d[["x"]][tr, , drop = FALSE]
  utr <- d[["u"]][tr, , drop = FALSE]
  ytr <- d[["y"]][tr]
  xte <- d[["x"]][te, , drop = FALSE]
  ute <- d[["u"]][te, , drop = FALSE]
  yte <- d[["y"]][te]
  grid <- d[["grid"]]
  res <- list()
  comp <- list()
  err <- list()
  keep_row <- function(name, obj, extra = list()) {
    res[[name]] <<- c(list(method = name,
                           rmse = rmse(yte, predict(obj, xte, ute)),
                           rmse.train = rmse(ytr, predict(obj)),
                           time = obj[["time"]],
                           blocks = obj[["blocks"]]), extra)
    g <- wafc_grid_components(obj, grid)
    if (!is.null(g)) comp[[name]] <<- g
  }
  try_fit <- function(name, expr) {
    tryCatch(expr, error = function(e) {
      err[[name]] <<- conditionMessage(e)
      cat(sprintf("  [%s %02d] %s failed: %s\n", key, s, name,
                  conditionMessage(e)))
      NULL
    })
  }

  ## The WAFC of decision D44, every argument written out: the block LASSO
  ## in the balanced form with free levels and the weights of grpreg, on
  ## the design of the WAFC (periodized Daublets 8, j0 = 0, rescaled to
  ## [0, 1] with eps = 0, decisions D15, D35), J on the grid of D34 by
  ## cross-validation on the blocked folds.
  kp <- try_fit("wafc.block", {
    t0 <- proc.time()[["elapsed"]]
    o <- wafc_competitor("klopp", x = xtr, u = utr, y = ytr, active = NULL,
                         J = J_grid, block.size = NULL,
                         penalize.levels = FALSE, chunk.weights = "sqrt",
                         merge.coarse = FALSE, balanced = TRUE,
                         free.coarse = FALSE, nfolds = nfolds,
                         foldid = sp[["foldid"]], j0 = 0L,
                         family = "Daublets", filter.size = 8L,
                         prec.wavelet = 30L, wavelet.filter = NULL,
                         boundary = "periodic", rescale = TRUE, eps = 0,
                         use.table = "auto", wavelet.table = get_wtab(),
                         sparse = "auto")
    o[["time"]] <- proc.time()[["elapsed"]] - t0
    o
  })
  if (!is.null(kp)) {
    ex <- kp[["extra"]]
    keep_row("wafc.block", kp,
             list(J = ex[["J"]], lambda = ex[["lambda"]], cve = ex[["cve"]],
                  block.size = ex[["block.size"]], conv = ex[["conv"]]))
    ## The folds of the threshold, fitted once for the two rules (D45).
    t0 <- proc.time()[["elapsed"]]
    ff <- try_fit("folds", wafc_threshold_folds(kp, s = NULL, y = NULL,
                                                foldid = NULL))
    t_ff <- proc.time()[["elapsed"]] - t0
    for (rule in c("cv1se", "cv")) {
      nm <- paste0("wafc.block+", rule)
      th <- if (is.null(ff)) NULL else
        try_fit(nm, wafc_threshold(kp, t = NULL, rule = rule, c = 0.15,
                                   refit = "none", s = NULL, y = NULL,
                                   foldid = NULL, truth = NULL,
                                   fold.fits = ff, gate = "none",
                                   alpha = 0.05, nsim = 200L,
                                   gate.seed = NULL, gate.test = NULL))
      if (is.null(th)) next
      th[["time"]] <- th[["time"]] + t_ff + kp[["time"]]
      keep_row(nm, th,
               list(J = th[["extra"]][["J"]], t = th[["extra"]][["t"]],
                    c = th[["extra"]][["c"]], norm = th[["extra"]][["norm"]],
                    levels = if (rule == "cv1se")
                      level_norms(th[["coef"]][-1L], th[["design"]]) else NULL,
                    candidates = th[["extra"]][["candidates"]]))
    }
  }

  ## The spline of decision D46, the basis dimension on the grid of D41
  ## (5, 10, 20, 40, 80, truncated at the distinct values of each
  ## modulator minus one), by REML and by GCV, on the engine of the pilot
  ## of E2.5h (bam; under GCV, bam without discretization).
  for (crit in c("reml", "gcv")) {
    nm <- paste0("gam.", crit)
    g <- try_fit(nm, wafc_competitor("gam", x = xtr, u = utr, y = ytr,
                                     active = NULL, k = NULL,
                                     k.select = crit, select = TRUE,
                                     edf.tol = 0.1, engine = "bam",
                                     nfolds = 10L, foldid = NULL))
    if (is.null(g)) next
    keep_row(nm, g, list(k = g[["extra"]][["k"]],
                         k.top = g[["extra"]][["k.top"]],
                         edf = g[["extra"]][["edf"]],
                         k.table = g[["extra"]][["k.table"]]))
  }
  li <- try_fit("linear", wafc_competitor("linear", x = xtr, u = utr, y = ytr,
                                          active = NULL))
  if (!is.null(li)) keep_row("linear", li)

  out <- list(key = key, split = s, code = code_hash, rows = res,
              components = comp, errors = err, grid = grid,
              ntrain = length(tr), ntest = length(te),
              ntest.blocks = sp[["ntest.blocks"]], nblocks = sp[["nblocks"]],
              sd.test = stats::sd(yte), J.grid = J_grid, nfolds = nfolds,
              n_max = n_max, time = proc.time()[["elapsed"]] - t_unit,
              mem = sum(gc()[, 2L]))
  saveRDS(out, path)
  cat(sprintf("  [%s %02d] done in %.0f s (J = %s, rmse wafc+cv1se %s, gam.reml %s)\n",
              key, s, out[["time"]],
              if (is.null(res[["wafc.block"]])) "-" else res[["wafc.block"]][["J"]],
              if (is.null(res[["wafc.block+cv1se"]])) "-" else
                sprintf("%.4f", res[["wafc.block+cv1se"]][["rmse"]]),
              if (is.null(res[["gam.reml"]])) "-" else
                sprintf("%.4f", res[["gam.reml"]][["rmse"]])))
  path
}

## ---------------------------------------------------------------------------
## Bases to run
## ---------------------------------------------------------------------------

## "all" stays the six bases of step E6.1b, so that its reproduction does
## not move; the base of step E6.1c is named on the command line.
keys <- if (identical(which_bases, "all")) setdiff(names(preps), "marylebone.ukair") else
  strsplit(which_bases, ",")[[1L]]
bad <- setdiff(keys, names(preps))
if (length(bad) > 0L) stop("Unknown base(s): ", paste(bad, collapse = ", "))

## ---------------------------------------------------------------------------
## Part "fit"
## ---------------------------------------------------------------------------

if ("fit" %in% parts) {
  ## The heaviest bases first, so that the last units to finish are short.
  weight <- c(beijing = 4, beijing.heat = 4, marylebone = 4, kelmarsh = 4,
              housing = 2, bike = 1, marylebone.ukair = 4)
  todo <- expand.grid(s = seq_len(nsplit), key = keys,
                      stringsAsFactors = FALSE)
  todo <- todo[order(-weight[todo[["key"]]], todo[["s"]]), ]
  todo <- todo[!file.exists(file.path(unit_dir, sprintf("%s-%02d.rds",
                                                        todo[["key"]],
                                                        todo[["s"]]))), ]
  cat(sprintf("\n[fit] %d units to run\n", nrow(todo)))
  t0 <- proc.time()[["elapsed"]]
  if (nrow(todo) > 0L) {
    ## the data and the table of the basis are prepared in the parent, so
    ## that every forked unit inherits them instead of rebuilding them
    for (k in unique(todo[["key"]])) invisible(get_data(k))
    invisible(get_wtab())
    run <- function(i) {
      tryCatch(run_unit(todo[["key"]][i], todo[["s"]][i]),
               error = function(e) {
                 cat(sprintf("  [%s %02d] FAILED: %s\n", todo[["key"]][i],
                             todo[["s"]][i], conditionMessage(e)))
                 NA_character_
               })
    }
    got <- if (ncores > 1L) {
      mclapply(seq_len(nrow(todo)), run, mc.cores = ncores,
               mc.preschedule = FALSE)
    } else lapply(seq_len(nrow(todo)), run)
    cat(sprintf("[fit] %d of %d units written, %.1f min\n",
                sum(!is.na(unlist(got))), nrow(todo),
                (proc.time()[["elapsed"]] - t0) / 60))
  }
  cat("fit exit\n")
}

## ---------------------------------------------------------------------------
## Part "extra": the check of the splines
## ---------------------------------------------------------------------------
##
## Added after the first reading of the run (2026-10-05). In beijing.heat
## the WAFC predicted 8 percent better than gam.reml in 20 of 20
## partitions, with k at the top of the grid of D41 (80) in every one and
## about 70 effective degrees of freedom in the smooths of the day of the
## year. Two readings fit that, and they say opposite things about the
## WAFC: the grid caps the spline (the lesson of E6.1a on dimension), or
## REML, which assumes independent errors, undersmooths an hourly series
## whose residuals are correlated (Opsomer, Wang and Yang, 2001), while the
## WAFC chooses (J, lambda) by cross-validation on folds of whole weeks.
## Two splines separate them, on the same partitions (make_split() is
## deterministic): "gam.cv" chooses k on the grid of D41 by the squared
## error on the blocked folds of the WAFC (k.select = "cv" of
## wafc_fit_gam(), REML inside each fold, step E2.5h), and "gam.reml.wide"
## is REML on the grid 80, 160, 320. E61B_EXTRA restricts the methods. The
## units go to <E61B_OUT>/units-extra/, and the report joins them to the
## units of the part "fit" when they exist.

run_extra <- function(key, s) {
  path <- file.path(extra_dir, sprintf("%s-%02d.rds", key, s))
  if (file.exists(path)) return(path)
  t_unit <- proc.time()[["elapsed"]]
  d <- get_data(key)
  sp <- make_split(d[["block"]], s, key)
  tr <- sp[["tr"]]
  te <- sp[["te"]]
  xtr <- d[["x"]][tr, , drop = FALSE]
  utr <- d[["u"]][tr, , drop = FALSE]
  ytr <- d[["y"]][tr]
  res <- list()
  comp <- list()
  err <- list()
  for (nm in extra_methods) {
    a <- switch(nm,
                gam.cv = list(k = NULL, k.select = "cv"),
                gam.reml.wide = list(k = c(80L, 160L, 320L),
                                     k.select = "reml"),
                stop("Unknown extra method ", nm, "."))
    g <- tryCatch(
      wafc_competitor("gam", x = xtr, u = utr, y = ytr, active = NULL,
                      k = a[["k"]], k.select = a[["k.select"]],
                      select = TRUE, edf.tol = 0.1, engine = "bam",
                      nfolds = nfolds, foldid = sp[["foldid"]]),
      error = function(e) {
        err[[nm]] <<- conditionMessage(e)
        cat(sprintf("  [%s %02d] %s failed: %s\n", key, s, nm,
                    conditionMessage(e)))
        NULL
      })
    if (is.null(g)) next
    res[[nm]] <- list(method = nm,
                      rmse = rmse(d[["y"]][te],
                                  predict(g, d[["x"]][te, , drop = FALSE],
                                          d[["u"]][te, , drop = FALSE])),
                      rmse.train = rmse(ytr, predict(g)), time = g[["time"]],
                      blocks = g[["blocks"]], k = g[["extra"]][["k"]],
                      k.top = g[["extra"]][["k.top"]],
                      edf = g[["extra"]][["edf"]],
                      k.table = g[["extra"]][["k.table"]])
    comp[[nm]] <- wafc_grid_components(g, d[["grid"]])
  }
  out <- list(key = key, split = s, code = code_hash, rows = res,
              components = comp, errors = err, n_max = n_max,
              time = proc.time()[["elapsed"]] - t_unit)
  saveRDS(out, path)
  cat(sprintf("  [%s %02d] extra done in %.0f s (%s)\n", key, s, out[["time"]],
              paste(sprintf("%s %.4f k %s", names(res),
                            vapply(res, `[[`, 0, "rmse"),
                            vapply(res, function(r) paste(r[["k"]], collapse = ","), "")),
                    collapse = "; ")))
  path
}

if ("extra" %in% parts) {
  dir.create(extra_dir, showWarnings = FALSE, recursive = TRUE)
  todo <- expand.grid(s = seq_len(nsplit), key = keys,
                      stringsAsFactors = FALSE)
  todo <- todo[!file.exists(file.path(extra_dir, sprintf("%s-%02d.rds",
                                                         todo[["key"]],
                                                         todo[["s"]]))), ]
  cat(sprintf("\n[extra] %d units to run (%s)\n", nrow(todo),
              paste(extra_methods, collapse = ", ")))
  t0 <- proc.time()[["elapsed"]]
  if (nrow(todo) > 0L) {
    for (k in unique(todo[["key"]])) invisible(get_data(k))
    run <- function(i) {
      tryCatch(run_extra(todo[["key"]][i], todo[["s"]][i]),
               error = function(e) {
                 cat(sprintf("  [%s %02d] extra FAILED: %s\n",
                             todo[["key"]][i], todo[["s"]][i],
                             conditionMessage(e)))
                 NA_character_
               })
    }
    got <- if (ncores > 1L) {
      mclapply(seq_len(nrow(todo)), run, mc.cores = ncores,
               mc.preschedule = FALSE)
    } else lapply(seq_len(nrow(todo)), run)
    cat(sprintf("[extra] %d of %d units written, %.1f min\n",
                sum(!is.na(unlist(got))), nrow(todo),
                (proc.time()[["elapsed"]] - t0) / 60))
  }
  cat("extra exit\n")
}

## ---------------------------------------------------------------------------
## Part "report"
## ---------------------------------------------------------------------------

if ("report" %in% parts) {
  methods <- c("wafc.block+cv1se", "wafc.block+cv", "wafc.block", "gam.reml",
               "gam.gcv", "linear")
  summary_all <- list()
  for (key in keys) {
    fs <- sort(list.files(unit_dir, pattern = sprintf("^%s-[0-9]+\\.rds$",
                                                      gsub(".", "\\.", key,
                                                           fixed = TRUE)),
                          full.names = TRUE))
    if (length(fs) == 0L) next
    U <- lapply(fs, readRDS)
    U <- U[vapply(U, function(z) identical(z[["n_max"]], n_max), TRUE)]
    if (length(U) == 0L) next
    ## the splines of the part "extra", when they were run
    U <- lapply(U, function(z) {
      ep <- file.path(extra_dir, sprintf("%s-%02d.rds", z[["key"]],
                                         z[["split"]]))
      if (file.exists(ep)) {
        ex <- readRDS(ep)
        z[["rows"]] <- c(z[["rows"]], ex[["rows"]])
        z[["components"]] <- c(z[["components"]], ex[["components"]])
      }
      z
    })
    K <- length(U)
    mth <- c(methods, intersect(c("gam.cv", "gam.reml.wide"),
                                unique(unlist(lapply(U, function(z)
                                  names(z[["rows"]]))))))
    gams <- grep("^gam\\.", mth, value = TRUE)
    d <- get_data(key)
    cat(sprintf("\n=== %s: %s ===\n", key, d[["label"]]))
    cat(sprintf("  n = %d, %d partitions, test %d to %d obs (%d of %d %ss), code %s\n",
                length(d[["y"]]), K, min(vapply(U, `[[`, 0L, "ntest")),
                max(vapply(U, `[[`, 0L, "ntest")), U[[1L]][["ntest.blocks"]],
                U[[1L]][["nblocks"]], d[["block.unit"]],
                paste(unique(vapply(U, `[[`, "", "code")), collapse = ",")))

    ## Prediction: the RMSE on the test sample of each partition, its mean
    ## and the standard error between partitions; the ratio to the best of
    ## the two splines of the same partition, and the paired difference to
    ## that spline with its standard error, naive and corrected for the
    ## overlap of the training samples by the factor of Nadeau and Bengio
    ## (2003), (1/K + n_test/n_train) / (1/K).
    ## The best spline is the best of every spline measured: the two of
    ## D46, and the two of the part "extra" when they were run.
    R <- sapply(mth, function(m) vapply(U, function(z) {
      r <- z[["rows"]][[m]]
      if (is.null(r)) NA_real_ else r[["rmse"]]
    }, 0))
    R <- matrix(R, nrow = K, dimnames = list(NULL, mth))
    best_gam <- apply(R[, gams, drop = FALSE], 1L, min, na.rm = TRUE)
    which_best <- names(which.min(colMeans(R[, gams, drop = FALSE],
                                           na.rm = TRUE)))
    rho <- mean(vapply(U, function(z) z[["ntest"]] / z[["ntrain"]], 0))
    nb <- sqrt((1 / K + rho) / (1 / K))
    pred <- data.frame(
      method = mth,
      mean = colMeans(R, na.rm = TRUE),
      se = apply(R, 2L, stats::sd, na.rm = TRUE) / sqrt(colSums(!is.na(R))),
      ratio.best.gam = colMeans(R / best_gam, na.rm = TRUE),
      diff.vs = which_best,
      diff = colMeans(R - R[, which_best], na.rm = TRUE),
      diff.se = apply(R - R[, which_best], 2L, stats::sd, na.rm = TRUE) /
        sqrt(K),
      wins = colMeans(R < R[, which_best], na.rm = TRUE),
      row.names = NULL)
    pred[["diff.se.nb"]] <- pred[["diff.se"]] * nb
    pred[["R2"]] <- 1 - (pred[["mean"]] / mean(vapply(U, `[[`, 0,
                                                     "sd.test")))^2
    cat("\n  prediction, RMSE on the held-out blocks (mean over partitions):\n")
    print(format(pred, digits = 4), row.names = FALSE)
    cat(sprintf("  best spline by mean: %s; Nadeau-Bengio factor on the SE %.2f\n",
                which_best, nb))

    ## Choices: J of the block LASSO, k of the splines, t of the threshold.
    pick <- function(m, what) vapply(U, function(z) {
      r <- z[["rows"]][[m]]
      if (is.null(r) || is.null(r[[what]])) NA_character_ else
        paste(r[[what]], collapse = ",")
    }, "")
    choices <- list(J = table(pick("wafc.block", "J")),
                    k.reml = table(pick("gam.reml", "k")),
                    k.gcv = table(pick("gam.gcv", "k")),
                    k.top.reml = table(pick("gam.reml", "k.top")),
                    k.top.gcv = table(pick("gam.gcv", "k.top")))
    for (m in setdiff(gams, c("gam.reml", "gam.gcv"))) {
      choices[[paste0("k.", sub("^gam\\.", "", m))]] <- table(pick(m, "k"))
      choices[[paste0("k.top.", sub("^gam\\.", "", m))]] <-
        table(pick(m, "k.top"))
    }
    cat("\n  choices over the partitions:\n")
    for (nm in names(choices)) {
      cat(sprintf("    %-10s %s\n", nm,
                  paste(sprintf("%s (%d)", names(choices[[nm]]),
                                as.integer(choices[[nm]])), collapse = "  ")))
    }

    ## Stability of the structure: the frequency with which each block is
    ## kept over the partitions, by method; for the splines, kept means an
    ## effective degrees of freedom above 0.1, and the median edf is beside.
    bl <- function(m) {
      A <- lapply(U, function(z) z[["rows"]][[m]][["blocks"]])
      A <- A[!vapply(A, is.null, TRUE)]
      if (length(A) == 0L) return(NULL)
      Reduce(`+`, lapply(A, function(a) a * 1)) / length(A)
    }
    edf_med <- function(m) {
      A <- lapply(U, function(z) z[["rows"]][[m]][["edf"]])
      A <- A[!vapply(A, is.null, TRUE)]
      if (length(A) == 0L) return(NULL)
      apply(simplify2array(A), c(1L, 2L), stats::median)
    }
    nrm_med <- local({
      A <- lapply(U, function(z) z[["rows"]][["wafc.block+cv1se"]][["norm"]])
      A <- A[!vapply(A, is.null, TRUE)]
      apply(simplify2array(A), c(1L, 2L), stats::median)
    })
    xn <- colnames(d[["x"]])
    un <- colnames(d[["u"]])
    stab <- expand.grid(x = xn, u = un, stringsAsFactors = FALSE)
    stab <- stab[order(match(stab[["x"]], xn)), ]
    at <- cbind(match(stab[["x"]], xn), match(stab[["u"]], un))
    for (m in c("wafc.block+cv1se", "wafc.block+cv", "wafc.block", gams)) {
      f <- bl(m)
      stab[[m]] <- if (is.null(f)) NA_real_ else f[at]
    }
    stab[["norm.cv1se"]] <- nrm_med[at]
    for (m in gams) {
      e <- edf_med(m)
      stab[[paste0("edf.", m)]] <- if (is.null(e)) NA_real_ else e[at]
    }
    ## Shape: the localization index of the median component of the WAFC
    ## and of the REML spline, and where the largest increment of the
    ## median component is.
    comp_med <- function(m) {
      A <- lapply(U, function(z) z[["components"]][[m]])
      A <- A[!vapply(A, is.null, TRUE)]
      if (length(A) == 0L) return(NULL)
      out <- A[[1L]]
      for (i in seq_along(out)) {
        out[[i]] <- apply(sapply(A, function(a) a[[i]]), 1L, stats::median)
      }
      out
    }
    cw <- comp_med("wafc.block+cv1se")
    cg <- comp_med("gam.reml")
    G <- d[["grid"]]
    stab[["loc.wafc"]] <- vapply(seq_len(nrow(at)), function(i)
      localization(cw[[at[i, 1L], at[i, 2L]]]), 0)
    stab[["loc.gam"]] <- vapply(seq_len(nrow(at)), function(i)
      localization(cg[[at[i, 1L], at[i, 2L]]]), 0)
    stab[["jump.at.wafc"]] <- vapply(seq_len(nrow(at)), function(i) {
      g <- cw[[at[i, 1L], at[i, 2L]]]
      if (all(g == 0)) return(NA_real_)
      k <- which.max(abs(diff(g)))
      (G[k, at[i, 2L]] + G[k + 1L, at[i, 2L]]) / 2
    }, 0)
    stab[["range.wafc"]] <- vapply(seq_len(nrow(at)), function(i)
      diff(range(cw[[at[i, 1L], at[i, 2L]]])), 0)
    stab[["range.gam"]] <- vapply(seq_len(nrow(at)), function(i)
      diff(range(cg[[at[i, 1L], at[i, 2L]]])), 0)
    ## A modulator with few distinct values (the hour, 24) leaves a block
    ## of the WAFC undetermined between them when 2^J - 1 exceeds that
    ## number: the fit at the observed values is identified, the curve on
    ## a continuous grid and the norm of the block are not (section 2 of
    ## the .md, on bike). The shape statistics of those blocks are not
    ## reported.
    disc <- un[vapply(seq_along(un), function(m)
      length(unique(d[["u"]][, m])) <= 30L, TRUE)]
    for (v in c("loc.wafc", "jump.at.wafc", "range.wafc")) {
      stab[[v]][stab[["u"]] %in% disc] <- NA_real_
    }
    cat("\n  stability of the structure (frequency kept over the partitions):\n")
    print(format(stab, digits = 3), row.names = FALSE)

    ## The time and the memory of a unit.
    tm <- sapply(mth, function(m) stats::median(vapply(U, function(z) {
      r <- z[["rows"]][[m]]
      if (is.null(r)) NA_real_ else r[["time"]]
    }, 0), na.rm = TRUE))
    cat(sprintf("\n  median seconds: %s\n  unit: median %.0f s, max %.0f s\n",
                paste(sprintf("%s %.0f", names(tm), tm), collapse = ", "),
                stats::median(vapply(U, `[[`, 0, "time")),
                max(vapply(U, `[[`, 0, "time"))))
    errs <- unlist(lapply(U, `[[`, "errors"))
    if (length(errs) > 0L) {
      cat("  errors:\n")
      print(table(names(errs)))
    }

    ## Figure: the components of the WAFC (median and the range over the
    ## partitions) with the median component of the REML spline.
    png(file.path(out_dir, sprintf("e61b-%s.png", key)),
        width = 300 * length(un), height = 220 * length(xn))
    op <- par(mfrow = c(length(xn), length(un)), mar = c(3, 3, 2, 1),
              mgp = c(1.8, 0.6, 0))
    for (l in seq_along(xn)) {
      for (m in seq_along(un)) {
        A <- sapply(U, function(z) {
          a <- z[["components"]][["wafc.block+cv1se"]]
          if (is.null(a)) rep(NA_real_, nrow(G)) else a[[l, m]]
        })
        lo <- apply(A, 1L, min, na.rm = TRUE)
        hi <- apply(A, 1L, max, na.rm = TRUE)
        yl <- if (un[m] %in% disc) range(cg[[l, m]], na.rm = TRUE) else
          range(c(lo, hi, cg[[l, m]]), na.rm = TRUE)
        plot(G[, m], cw[[l, m]], type = "n", ylim = yl, xlab = un[m],
             ylab = "", main = sprintf("g[%s, %s]  kept %.2f / %.2f", xn[l],
                                       un[m],
                                       stab[["wafc.block+cv1se"]][stab[["x"]] == xn[l] & stab[["u"]] == un[m]],
                                       stab[["gam.reml"]][stab[["x"]] == xn[l] & stab[["u"]] == un[m]]))
        polygon(c(G[, m], rev(G[, m])), c(lo, rev(hi)), col = "grey85",
                border = NA)
        lines(G[, m], cw[[l, m]], lwd = 2)
        lines(G[, m], cg[[l, m]], col = "firebrick", lty = 2, lwd = 2)
        mk <- d[["marks"]][[un[m]]]
        if (!is.null(mk)) abline(v = mk, col = "steelblue", lty = 3)
      }
    }
    par(op)
    dev.off()

    summary_all[[key]] <- list(prediction = pred, choices = choices,
                               stability = stab, components.wafc = cw,
                               components.gam = cg, grid = G, K = K,
                               nb.factor = nb, time = tm, code = code_hash,
                               rmse = R, label = d[["label"]],
                               source = d[["source"]],
                               block.unit = d[["block.unit"]])
  }
  saveRDS(list(summary = summary_all, code = code_hash, J.grid = J_grid,
               nfolds = nfolds, nsplit = nsplit, n_max = n_max,
               sessionInfo = utils::sessionInfo()),
          file.path(out_dir, "e61b-summary.rds"))
  cat(sprintf("\nwritten: %s\nreport exit\n",
              file.path(out_dir, "e61b-summary.rds")))
}
