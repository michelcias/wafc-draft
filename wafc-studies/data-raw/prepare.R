## data-raw/prepare.R -- the data of the two applications, from the sources.
##
##   Rscript data-raw/prepare.R [--dest=outputs/data-raw] [--check]
##
## Run from the root of the compendium, after data-raw/fetch.R. Writes
## data/marylebone.rds and data/beijing-dongsi.rds from the source files in
## <dest>; with --check, rebuilds them and compares them with the versioned
## files value by value, writing nothing. The response, the covariates and
## the modulators of each application are formed from these data in
## R/application_data.R; this script only reads, converts, aligns and
## removes incomplete records.

if (!file.exists(file.path("data-raw", "sources.yaml"))) {
  stop("run from the root of the compendium.", call. = FALSE)
}

opt <- local({
  out <- list(dest = file.path("outputs", "data-raw"), check = FALSE)
  for (a in commandArgs(trailingOnly = TRUE)) {
    if (identical(a, "--check")) {
      out[["check"]] <- TRUE
      next
    }
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1L]]
    if (!startsWith(a, "--") || length(kv) != 2L || kv[1L] != "dest") {
      stop("unrecognized argument '", a, "'; use --dest= or --check.",
           call. = FALSE)
    }
    out[["dest"]] <- kv[2L]
  }
  out
})
src <- yaml::read_yaml(file.path("data-raw", "sources.yaml"))

## Marylebone Road (site MY1 of the Automatic Urban and Rural Network),
## hourly, from 1998-01-01 00:00 to 2005-06-23 12:00 GMT.
##
## UK-AIR reports ug/m3 at 20 C and 1013 mb, NOx as NO2. The oxidant
## relation the application uses is molar, so the three pollutants are
## converted to ppb with the factors of Defra for the UK Air Quality Archive:
## 1 ppb of NO2, and of NOx as NO2, is 1.9125 ug/m3, and 1 ppb of O3 is
## 1.9957 ug/m3 (Defra, 2005, "Conversion Factors Between ppb and ug m-3
## and ppm and mgm-3"). The date of UK-AIR is the beginning of the hour of
## the mean. The wind of ERA5 is instantaneous, so each hour gets the speed
## at its middle, by linear interpolation between t and t + 1 h. The wind
## direction is read and not used. Hours with any of the four variables
## missing are removed.
prepare_marylebone <- function(dir) {
  a <- do.call(rbind, lapply(1998:2005, function(yr) {
    e <- new.env()
    load(file.path(dir, sprintf("MY1_%d.RData", yr)), envir = e)
    as.data.frame(e[[sprintf("MY1_%d", yr)]])[, c("date", "NOXasNO2", "NO2",
                                                  "O3")]
  }))
  w <- utils::read.csv(file.path(dir, "era5-wind-my1-1998-2005.csv"),
                       skip = 3L, check.names = FALSE)
  wt <- as.POSIXct(w[[1L]], format = "%Y-%m-%dT%H:%M", tz = "GMT")
  ws <- stats::approx(as.numeric(wt), w[[2L]],
                      xout = as.numeric(a[["date"]]) + 1800)[["y"]]
  t0 <- as.POSIXct("1998-01-01 00:00:00", tz = "GMT")
  t1 <- as.POSIXct("2005-06-23 12:00:00", tz = "GMT")
  keep <- a[["date"]] >= t0 & a[["date"]] <= t1 &
    stats::complete.cases(ws, a[["NOXasNO2"]], a[["NO2"]], a[["O3"]])
  a <- a[keep, ]
  out <- data.frame(date = a[["date"]],
                    nox = a[["NOXasNO2"]] / 1.9125,
                    no2 = a[["NO2"]] / 1.9125,
                    o3 = a[["O3"]] / 1.9957,
                    ws = ws[keep])
  rownames(out) <- NULL
  attr(out, "units") <- c(nox = "ppb (NOx as NO2)", no2 = "ppb", o3 = "ppb",
                          ws = "m/s, 10 m, ERA5")
  attr(out, "source") <- paste("data-raw/sources.yaml (marylebone); licence and",
                               "attribution in data/README.md")
  out
}

## Beijing, site Dongsi, hourly. The relative humidity is computed from the
## temperature and the dew point by the Magnus form (coefficients 17.625 and
## 243.04 C). Hours with PM2.5, temperature, pressure, dew point or wind
## speed missing, a relative humidity outside (0, 100] or a PM2.5 of zero
## are removed.
prepare_beijing <- function(dir) {
  d <- utils::read.csv(file.path(dir, src[["beijing"]][["files"]][[1L]][[
    "extract"]][[1L]][["member"]]))
  mg <- function(t) exp(17.625 * t / (243.04 + t))
  rh <- 100 * mg(d[["DEWP"]]) / mg(d[["TEMP"]])
  keep <- stats::complete.cases(d[["PM2.5"]], d[["TEMP"]], d[["PRES"]],
                                d[["DEWP"]], d[["WSPM"]]) &
    is.finite(rh) & rh > 0 & rh <= 100 & d[["PM2.5"]] > 0
  d <- d[keep, ]
  out <- data.frame(date = as.Date(sprintf("%04d-%02d-%02d", d[["year"]],
                                           d[["month"]], d[["day"]])),
                    hour = d[["hour"]], pm25 = d[["PM2.5"]],
                    temp = d[["TEMP"]], pres = d[["PRES"]],
                    dewp = d[["DEWP"]], wspm = d[["WSPM"]], rh = rh[keep])
  rownames(out) <- NULL
  attr(out, "units") <- c(pm25 = "ug/m3", temp = "C", pres = "hPa",
                          dewp = "C", wspm = "m/s", rh = "percent")
  attr(out, "source") <- paste("data-raw/sources.yaml (beijing); licence and",
                               "attribution in data/README.md")
  out
}

built <- list(
  marylebone = prepare_marylebone(file.path(opt[["dest"]],
                                            src[["marylebone"]][["dir"]])),
  "beijing-dongsi" = prepare_beijing(file.path(opt[["dest"]],
                                               src[["beijing"]][["dir"]])))

ok <- TRUE
for (nm in names(built)) {
  path <- file.path("data", paste0(nm, ".rds"))
  if (opt[["check"]]) {
    same <- file.exists(path) && identical(readRDS(path), built[[nm]])
    message(sprintf("%-15s %d rows: %s", nm, nrow(built[[nm]]),
                    if (same) "identical to the versioned file" else
                      "DIFFERS from the versioned file"))
    ok <- ok && same
  } else {
    dir.create("data", showWarnings = FALSE)
    saveRDS(built[[nm]], path, compress = "xz")
    message(sprintf("%-15s %d rows written to %s", nm, nrow(built[[nm]]),
                    path))
  }
}
if (!ok) quit(status = 1L)
