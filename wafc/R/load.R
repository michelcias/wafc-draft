## wafc/R/load.R -- dependencies and loader of the WAFC code.
##
## Loads every wafc/R/*.R file into the calling session and declares the
## packages the code depends on. From the root of the repository:
##
##     source("wafc/R/load.R")
##
## Decision D4 (docs/ESTADO.md): the method lives in this folder, not in
## WaveBased, which is used only as a dependency for the wavelet bases
## (wbasis(), wtable()). Every new dependency is declared here and in
## docs/CONTINUAR.md in the same round (docs/instrucoes.md, section 6).

## Packages the code in wafc/R/ requires. 'suggests' are needed by scripts
## and by the variants of the estimator, not by the functions loaded here.
## VCBART is the engine of wafc_fit_vcbart(); step E2.4 used it and did not
## declare it, because wafc/R/load.R was outside its catalogue.
wafc_depends <- c("WaveBased", "glmnet", "Matrix")
wafc_suggests <- c("testthat", "sparsegl", "grpreg", "mgcv", "VCBART")

## Directory holding this file, resolved from the source() frame when
## available and from the working directory otherwise.
wafc_r_dir <- local({
  ofile <- NULL
  for (i in rev(seq_len(sys.nframe()))) {
    cand <- get0("ofile", envir = sys.frame(i), inherits = FALSE)
    if (is.character(cand) && length(cand) == 1L && nzchar(cand)) {
      ofile <- cand
      break
    }
  }
  if (!is.null(ofile)) {
    dirname(normalizePath(ofile, mustWork = TRUE))
  } else if (file.exists(file.path("wafc", "R", "design.R"))) {
    normalizePath(file.path("wafc", "R"))
  } else if (file.exists(file.path("R", "design.R"))) {
    normalizePath("R")
  } else if (file.exists(file.path("..", "R", "design.R"))) {
    normalizePath(file.path("..", "R"))
  } else {
    stop("Could not locate wafc/R/. Call source(\"wafc/R/load.R\") from the root of the repository.")
  }
})

## Attaches the dependencies, with an informative error naming every
## missing package at once instead of failing on the first one.
wafc_attach <- function(pkgs = wafc_depends) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, TRUE, quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("The following package(s) are required by wafc/ and are not installed: ",
         paste(missing, collapse = ", "),
         ". See docs/CONTINUAR.md, section 2.", call. = FALSE)
  }
  for (pk in pkgs) {
    suppressPackageStartupMessages(library(pk, character.only = TRUE))
  }
  invisible(pkgs)
}

## Reports the packages of 'wafc_suggests' that are not installed, without
## failing: the functions in wafc/R/ do not need them.
wafc_check_suggests <- function(pkgs = wafc_suggests) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, TRUE, quietly = TRUE)]
  if (length(missing) > 0L) {
    message("wafc: optional package(s) not installed: ",
            paste(missing, collapse = ", "), ".")
  }
  invisible(missing)
}

wafc_attach()

local({
  files <- sort(list.files(wafc_r_dir, pattern = "\\.[rR]$", full.names = TRUE))
  files <- files[basename(files) != "load.R"]
  for (f in files) source(f, local = FALSE, echo = FALSE)
  invisible(files)
})
