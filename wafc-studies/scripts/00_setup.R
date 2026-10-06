## scripts/00_setup.R -- the session every other script starts from.
##
## Checks that R runs from the root of the compendium, that the locked
## environment is the one in use, and sources the functions of R/. It
## computes nothing and is safe to source again. The WAFC code itself is
## sourced by the scripts that fit, from the path in the study
## configuration (`code`).

if (!file.exists(file.path("config", "study.yaml")) ||
    !file.exists(file.path("R", "run.R"))) {
  stop("run from the root of the compendium (the folder with config/ and ",
       "R/).", call. = FALSE)
}

if (file.exists("renv.lock")) {
  if (!requireNamespace("renv", quietly = TRUE)) {
    stop("renv.lock is present but the package 'renv' is not installed: ",
         "install.packages(\"renv\"), then renv::restore().", call. = FALSE)
  }
  if (is.null(renv::project())) {
    message("[setup] the renv project is not active; start R from the root ",
            "of the compendium so that .Rprofile activates it.")
  }
}

for (pkg in c("yaml", "WaveBased", "mgcv", "grpreg", "glmnet", "Matrix",
              "VCBART")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("package '", pkg, "' is not installed; run renv::restore().",
         call. = FALSE)
  }
}

for (f in sort(list.files("R", pattern = "\\.[Rr]$", full.names = TRUE))) {
  source(f)
}
