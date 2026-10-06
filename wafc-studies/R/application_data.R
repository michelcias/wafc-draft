## R/application_data.R -- the two data applications.
##
## Each application reads one file of data/ (built by data-raw/prepare.R)
## and forms the response y, the linear covariates x (column 1 the
## constant), the modulating covariates u, and the block of each
## observation, the unit of the partitions: a week, because the series are
## hourly and their errors are correlated over hours and days, so a test
## sample of hours scattered among the training hours would reward a method
## for learning the level of a week. The folds of the cross-validation
## inside a training sample are blocked the same way (R/application.R).
##
##   marylebone.ukair  Marylebone Road, London; the structure of the
##                     coefficients (the article)
##   beijing.heat      Beijing, site Dongsi; the prediction (the
##                     supplementary material)

#' The applications
#'
#' For each one: a label, the file of data/, the function that forms the
#' model from it, the labels of the plots, the marks drawn on the plots
#' of a modulator, and the modulators whose grid has a fixed range (a
#' periodic modulator is shown over its whole period).
application_bases <- list(

  ## Marylebone Road, London, hourly, 1998-01-01 to 2005-06-23 (UK-AIR and
  ## ERA5; data-raw/sources.yaml). The model is the oxidant relation of
  ## Clapp and Jenkin (2001): at a roadside site the total oxidant
  ## OX = NO2 + O3 is close to linear in NOx, with intercept the regional
  ## background of oxidant and slope the local contribution to it. The
  ## slope is an estimate of the fraction of NOx emitted directly as NO2,
  ## the primary fraction (Carslaw and Beevers, 2005, section 2.4): it also
  ## carries the NO2 formed by the thermal reaction 2NO + O2 and from HONO,
  ## and it is lower at night (Clapp and Jenkin, 2001, section 2.3). With
  ## NOx in units of 100 ppb, beta_nox(u) reads directly as a percentage.
  ## The modulators are the date (days since 1998-01-01) and the wind
  ## speed. Carslaw (2005, sections 3.1 and 3.2, Fig. 3(a)) estimated the
  ## primary fraction at Marylebone Road at about 10 percent by volume from
  ## 1997 to 2002, rising through 2002 and 2003 to about 23 percent at the
  ## end of 2003, a change he relates to the particle filters fitted to the
  ## London buses; the step of beta_nox along the date is compared with
  ## that rise, not with the mean of the London sites. The wind enters
  ## both coefficients: dispersion moves the background, and the mixing of
  ## ozone into the street can move the slope. The wind direction is not a
  ## modulator.
  marylebone.ukair = list(
    label = "Marylebone Road, London (hourly oxidant)",
    file = "marylebone.rds",
    yname = "NO2 + O3 (ppb)",
    xlab = c(one = "Background", nox = "NOx"),
    betalab = c(one = "Oxidant background (ppb)",
                nox = "Slope (percent of NOx)"),
    ulab = c(day = "Date", ws = "Wind speed (m/s)"),
    build = function(d) {
      t0 <- as.POSIXct("1998-01-01 00:00:00", tz = "GMT")
      day <- as.numeric(difftime(d[["date"]], t0, units = "days"))
      list(y = d[["no2"]] + d[["o3"]],
           x = cbind(one = 1, nox = d[["nox"]] / 100),
           u = cbind(day = day, ws = d[["ws"]]),
           block = as.integer(floor(day)) %/% 7L)
    },
    ## the date axis is printed as calendar dates from this origin
    origin = list(day = as.Date("1998-01-01")),
    ## 2003-01-01, the year of the rise reported by Carslaw (2005)
    marks = list(day = 1826)),

  ## Beijing, site Dongsi, hourly, March 2013 to February 2017 (UCI;
  ## data-raw/sources.yaml). Response log PM2.5; linear covariates the wind
  ## speed, the temperature and the pressure (centred at 1000 hPa); the
  ## modulators are the day of the year (0 on 1 January, with the hour as a
  ## fraction of the day) and the relative humidity. The day of the year is
  ## periodic, as the periodized wavelet basis is, and its grid covers the
  ## whole year. The winter heating season of Beijing usually runs from 15
  ## November to 15 March, and the weather-adjusted level of PM2.5 rises at
  ## its start and falls at its end (Liang et al., 2015, section 7): the
  ## emissions, and with them the response of PM2.5 to the weather, change
  ## there. The effective dates move with the temperature from year to
  ## year, so the marks on the plots, day 318 and day 73 of the year, are
  ## the nominal dates.
  beijing.heat = list(
    label = "Beijing, Dongsi, by day of the year (hourly PM2.5)",
    file = "beijing-dongsi.rds",
    yname = "log PM2.5",
    xlab = c(one = "Level", wind = "Wind speed", temp = "Temperature",
             pres = "Pressure"),
    betalab = c(one = "Level (log PM2.5)", wind = "Per m/s",
                temp = "Per degree C", pres = "Per hPa"),
    ulab = c(doy = "Day of the year", rh = "Relative humidity (%)"),
    build = function(d) {
      doy <- as.integer(format(d[["date"]], "%j")) - 1 + d[["hour"]] / 24
      list(y = log(d[["pm25"]]),
           x = cbind(one = 1, wind = d[["wspm"]], temp = d[["temp"]],
                     pres = d[["pres"]] - 1000),
           u = cbind(doy = doy, rh = d[["rh"]]),
           block = as.integer(d[["date"]] - as.Date("2013-03-01")) %/% 7L)
    },
    fixed_range = list(doy = c(0, 365)),
    marks = list(doy = c(318, 73)))
)

#' The data of one application
#'
#' Reads the file of data/, forms the model, and adds the grid on which the
#' components are read: `grid_size` points over each modulator, between
#' the quantiles `grid_quantiles` of its values (the range in which every
#' training sample of the partitions has data), or over the fixed range of
#' a periodic modulator.
#'
#' @return A list with `y`, `x`, `u`, `block`, `grid`, the entry of
#'   application_bases, and the MD5 digest of the data file.
application_data <- function(cfg, base) {
  def <- application_bases[[base]]
  if (is.null(def)) stop("unknown application '", base, "'.", call. = FALSE)
  path <- file.path(cfg[["data_dir"]], def[["file"]])
  if (!file.exists(path)) {
    stop("data file not found: ", path, " (data-raw/prepare.R builds it).",
         call. = FALSE)
  }
  m <- def[["build"]](readRDS(path))
  G <- as.integer(cfg[["grid_size"]])
  qs <- unlist(cfg[["grid_quantiles"]])
  grid <- vapply(colnames(m[["u"]]), function(v) {
    r <- def[["fixed_range"]][[v]]
    if (is.null(r)) r <- stats::quantile(m[["u"]][, v], qs, names = FALSE)
    seq(r[1L], r[2L], length.out = G)
  }, numeric(G))
  c(m, list(grid = grid, base = base, def = def,
            digest = unname(tools::md5sum(path))))
}
