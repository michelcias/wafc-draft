## E8.3 -- fixed-width digest of a results CSV (context saver; base R only).
##
##     Rscript wafc/scripts/15-digest.R <csv> [--cells a,b] [--n 500,1000]
##         [--methods x,y] [--cols c1,c2,...] [--ratio <ref method>]
##         [--by <column>] [--digits 3]
##
## Prints one header line and one row per (design, cell, n, method) left after
## the filters, at most 40 rows. Above 40 rows it prints one line with the row
## count and the distinct cell, n and method values, and exits with status 2.
## --ratio M adds "ratio": the first --cols value over the value of method M
## in the same (design, cell, n); the reference row is found before the
## filters, so M need not be in --methods. --by names the column that holds
## the method labels (default "method"). Without --cols, a summary.csv shows
## rmse_f_mean, exact_support, time_median; any other CSV lists its columns.

MAXROWS <- 40L

usage <- function() {
  cat("Usage: Rscript wafc/scripts/15-digest.R <csv> [options]\n",
      "  --cells a,b      keep these values of column 'cell'\n",
      "  --n 500,1000     keep these values of column 'n'\n",
      "  --methods x,y    keep these values of column 'method' (or --by)\n",
      "  --cols c1,c2     numeric columns shown (default for summary.csv:\n",
      "                   rmse_f_mean,exact_support,time_median)\n",
      "  --ratio M        add ratio = first col / same col of method M\n",
      "  --by <column>    column holding the method labels (default method)\n",
      "  --digits 3       decimals shown (default 3)\n",
      "Prints at most ", MAXROWS, " rows; more rows -> one summary line, status 2.\n",
      "Without --cols on a non-summary CSV, lists the columns and exits.\n", sep = "")
}

die <- function(msg, status = 1L) {
  cat("digest: ", msg, "\n", sep = "")
  quit(save = "no", status = status)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || any(args %in% c("--help", "-h"))) {
  usage()
  quit(save = "no", status = 0L)
}

## Parse: first non-option is the CSV; options take one value each.
path <- NULL
opt <- list()
known <- c("cells", "n", "methods", "cols", "ratio", "by", "digits")
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (grepl("^--", a)) {
    key <- sub("^--", "", a)
    val <- NULL
    if (grepl("=", key, fixed = TRUE)) {
      val <- sub("^[^=]*=", "", key)
      key <- sub("=.*$", "", key)
    } else {
      if (i == length(args)) die(paste0("option --", key, " needs a value"))
      i <- i + 1L
      val <- args[[i]]
    }
    if (!(key %in% known)) die(paste0("unknown option --", key))
    opt[[key]] <- val
  } else if (is.null(path)) {
    path <- a
  } else {
    die(paste0("unexpected argument ", a))
  }
  i <- i + 1L
}
if (is.null(path)) die("no CSV given")
if (!file.exists(path)) die(paste0("file not found: ", path))

splitv <- function(x) strsplit(x, ",", fixed = TRUE)[[1L]]

dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
nlines <- length(readLines(path, warn = FALSE))
if (nrow(dat) != nlines - 1L)
  die(paste0("read ", nrow(dat), " rows but the file has ", nlines - 1L,
             " data lines (embedded newlines?)"))

digits <- if (is.null(opt[["digits"]])) 3L else suppressWarnings(as.integer(opt[["digits"]]))
if (is.na(digits) || digits < 0L) die("--digits must be a non-negative integer")

mcol <- if (is.null(opt[["by"]])) "method" else opt[["by"]]
if (!(mcol %in% names(dat))) die(paste0("column '", mcol, "' not found; columns: ",
                                        paste(names(dat), collapse = ", ")))

## Columns shown.
if (is.null(opt[["cols"]])) {
  if (!("rmse_f_mean" %in% names(dat))) {
    cat("digest: no --cols given; available columns:\n",
        paste(names(dat), collapse = ", "), "\n", sep = "")
    quit(save = "no", status = 1L)
  }
  cols <- intersect(c("rmse_f_mean", "exact_support", "time_median"), names(dat))
} else {
  cols <- splitv(opt[["cols"]])
}
miss <- setdiff(cols, names(dat))
if (length(miss) > 0L) die(paste0("column(s) not found: ", paste(miss, collapse = ", ")))
for (cc in cols) dat[[cc]] <- suppressWarnings(as.numeric(dat[[cc]]))

## Key columns present in the CSV.
keys <- intersect(c("design", "cell", "n"), names(dat))
rowkey <- if (length(keys) > 0L) do.call(paste, c(unname(as.list(dat[keys])), sep = "\r")) else rep("", nrow(dat))

## Ratio before filtering, so the reference method may be filtered out.
if (!is.null(opt[["ratio"]])) {
  ref <- dat[dat[[mcol]] == opt[["ratio"]], , drop = FALSE]
  if (nrow(ref) == 0L) die(paste0("reference method '", opt[["ratio"]], "' not found in ", mcol))
  refval <- ref[[cols[[1L]]]][match(rowkey, rowkey[dat[[mcol]] == opt[["ratio"]]])]
  dat[["ratio"]] <- dat[[cols[[1L]]]] / refval
}

## Filters.
filt <- function(d, col, vals, label) {
  if (is.null(vals)) return(d)
  if (!(col %in% names(d))) die(paste0("filter --", label, ": column '", col, "' does not exist"))
  d[as.character(d[[col]]) %in% splitv(vals), , drop = FALSE]
}
dat <- filt(dat, "cell", opt[["cells"]], "cells")
dat <- filt(dat, "n", opt[["n"]], "n")
dat <- filt(dat, mcol, opt[["methods"]], "methods")

## Order by design, cell, n (numeric), method.
ordcols <- lapply(intersect(c("design", "cell", "n", mcol), names(dat)),
                  function(k) if (k == "n") as.numeric(dat[[k]]) else dat[[k]])
if (nrow(dat) > 1L && length(ordcols) > 0L) dat <- dat[do.call(order, ordcols), , drop = FALSE]

## Size guard.
if (nrow(dat) > MAXROWS) {
  dv <- function(k) if (k %in% names(dat)) paste(unique(dat[[k]]), collapse = ",") else "-"
  cat("digest: ", nrow(dat), " rows > ", MAXROWS, "; narrow with filters. cell: ",
      dv("cell"), " | n: ", dv("n"), " | method: ", dv(mcol), "\n", sep = "")
  quit(save = "no", status = 2L)
}
if (nrow(dat) == 0L) die("no rows left after the filters")

## Fixed-width print.
shown <- c(keys, mcol, cols, if (!is.null(opt[["ratio"]])) "ratio")
shown <- unique(shown)
out <- dat[shown]
for (k in shown) {
  v <- out[[k]]
  out[[k]] <- if (k %in% c(cols, "ratio")) {
    ifelse(is.na(v), "NA", formatC(v, format = "f", digits = digits))
  } else as.character(v)
}
w <- vapply(shown, function(k) max(nchar(c(k, out[[k]]))), integer(1L))
pad <- function(s, k) if (k %in% c(cols, "ratio", "n")) formatC(s, width = w[[k]]) else formatC(s, width = -w[[k]])
cat(paste(vapply(shown, function(k) pad(k, k), ""), collapse = "  "), "\n", sep = "")
for (r in seq_len(nrow(out)))
  cat(paste(vapply(shown, function(k) pad(out[[k]][[r]], k), ""), collapse = "  "), "\n", sep = "")
