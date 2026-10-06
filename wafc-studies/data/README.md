# Data of the applications

Two data sets, read by `R/application_data.R`. Both are built from public
sources by `data-raw/` and versioned here, so that the applications run
without downloading anything. The response, the covariates and the
modulators of each application are formed from these files in
`R/application_data.R`.

| File | Rows | Size | Source | Licence |
|---|---|---|---|---|
| `marylebone.rds` | 61 280 hours | 0.35 MB | UK-AIR (Defra), site MY1, 1998-01-01 00:00 to 2005-06-23 12:00 GMT; ERA5 wind through Open-Meteo | OGL (pollutants); CC BY 4.0 (wind) |
| `beijing-dongsi.rds` | 34 287 hours | 0.33 MB | UCI Machine Learning Repository, data set 501, site Dongsi, 2013-03-01 to 2017-02-28 | CC BY 4.0 |

## marylebone.rds

Hourly means at Marylebone Road, London (site MY1 of the Automatic Urban
and Rural Network), with the hourly 10 m wind speed of the ERA5 reanalysis
at the nearest grid point. Columns:

- `date`: beginning of the hour of the mean (POSIXct, GMT);
- `nox`, `no2`, `o3`: NOx (as NO2), NO2 and O3 in ppb, converted from the
  ug/m3 at 20 C and 1013 mb of UK-AIR with the factors of Defra (1 ppb of
  NO2, and of NOx as NO2, is 1.9125 ug/m3; of O3, 1.9957 ug/m3);
- `ws`: wind speed at 10 m (m/s), the instantaneous values of ERA5
  interpolated linearly to the middle of the hour.

Hours with any of the four variables missing are left out.

Attribution:

> (c) Crown 2026 copyright Defra via uk-air.defra.gov.uk, licenced under
> the Open Government Licence (OGL).

> Generated using or contains modified Copernicus Climate Change Service
> information 2026. Neither the European Commission nor ECMWF is
> responsible for any use that may be made of the Copernicus information or
> data it contains.

> Hourly ERA5 wind data were retrieved through the Open-Meteo historical
> weather API (Zippenfenig, 2024; https://open-meteo.com/), licensed under
> CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/); the
> instantaneous values were linearly interpolated to the middle of each
> hour.

References: Hersbach et al. (2023), "ERA5 hourly data on single levels
from 1940 to present", Copernicus Climate Change Service (C3S) Climate
Data Store, doi:10.24381/cds.adbb2d47; Zippenfenig (2024), "Open-Meteo.com
Weather API", Zenodo, doi:10.5281/zenodo.7970649; Defra (2005),
"Conversion Factors Between ppb and ug m-3 and ppm and mgm-3", UK-AIR.

## beijing-dongsi.rds

Hourly observations at the site Dongsi of the Beijing Multi-Site Air
Quality data. Columns: `date` (Date), `hour` (0 to 23), `pm25` (ug/m3),
`temp` (C), `pres` (hPa), `dewp` (dew point, C), `wspm` (wind speed, m/s),
and `rh`, the relative humidity (percent) computed from the temperature
and the dew point by the Magnus form. Hours with PM2.5, temperature,
pressure, dew point or wind speed missing, a relative humidity outside
(0, 100] or a PM2.5 of zero are left out.

Attribution:

> Beijing Multi-Site Air Quality data, UCI Machine Learning Repository
> (https://doi.org/10.24432/C5RK5G), licensed under CC BY 4.0
> (https://creativecommons.org/licenses/by/4.0/); only the site Dongsi is
> used.

Reference: Zhang, Guo, Dong, He, Xu and Chen (2017), "Cautionary tales on
air-quality improvement in Beijing", *Proceedings of the Royal Society A*
473, 20170457.

## Rebuilding the data

From the root of the compendium:

```bash
Rscript data-raw/fetch.R
Rscript data-raw/prepare.R --check
```

`fetch.R` downloads the source files listed in `data-raw/sources.yaml` to
`outputs/data-raw/` and checks each against the SHA-256 recorded there;
`--from=DIR` copies them from a local folder instead. `prepare.R --check`
rebuilds the two files and compares them with these, value by value;
without `--check` it overwrites them. The wind file of Open-Meteo is
generated at each request, so a new download may differ in its digest
without differing in its values: `fetch.R` warns, and `prepare.R --check`
is the test that matters.
