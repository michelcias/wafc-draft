## data-raw/fetch.R -- the source files of the two data applications.
##
##   Rscript data-raw/fetch.R [--dest=outputs/data-raw] [--from=DIR]
##                            [--sets=marylebone,beijing]
##
## Run from the root of the compendium. Every file of data-raw/sources.yaml
## is placed at <dest>/<dir>/<file> and checked against its SHA-256; the
## member of an archive the preparation reads is extracted beside it and
## checked as well. A file already in <dest> with the right digest is kept.
## With --from, the files are copied from a local folder, searched by name
## in it and its subfolders, instead of being downloaded; a copy whose
## digest differs from the recorded one is refused.
##
## The data in data/ are built from these files by data-raw/prepare.R, and
## are versioned; this script is needed only to rebuild them. The sources,
## their licences and the attribution each one asks for are in
## data-raw/sources.yaml and data/README.md.

if (!file.exists(file.path("data-raw", "sources.yaml"))) {
  stop("run from the root of the compendium.", call. = FALSE)
}

opt <- local({
  out <- list(dest = file.path("outputs", "data-raw"), from = NULL,
              sets = NULL)
  for (a in commandArgs(trailingOnly = TRUE)) {
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1L]]
    if (!startsWith(a, "--") || length(kv) != 2L ||
        !kv[1L] %in% names(out)) {
      stop("unrecognized argument '", a, "'; use --dest=, --from= or ",
           "--sets=.", call. = FALSE)
    }
    out[[kv[1L]]] <- kv[2L]
  }
  out
})

src <- yaml::read_yaml(file.path("data-raw", "sources.yaml"))
sets <- if (is.null(opt[["sets"]])) names(src) else
  strsplit(opt[["sets"]], ",", fixed = TRUE)[[1L]]
bad <- setdiff(sets, names(src))
if (length(bad) > 0L) stop("unknown data set(s): ", paste(bad, collapse = ", "),
                           call. = FALSE)

sha <- function(path) unname(tools::sha256sum(path))

## Where a file comes from: a copy from the local folder, or the download.
obtain <- function(f, path) {
  if (!is.null(opt[["from"]])) {
    hit <- list.files(opt[["from"]], recursive = TRUE, full.names = TRUE)
    hit <- hit[basename(hit) == f[["file"]]]
    hit <- hit[vapply(hit, sha, "") == f[["sha256"]]]
    if (length(hit) == 0L) {
      stop("no copy of ", f[["file"]], " with the recorded digest under ",
           opt[["from"]], ".", call. = FALSE)
    }
    file.copy(hit[1L], path, overwrite = TRUE)
    message("  copied   ", f[["file"]], " from ", hit[1L])
  } else {
    message("  download ", f[["url"]])
    utils::download.file(f[["url"]], path, mode = "wb", quiet = TRUE)
  }
}

check <- function(path, want, volatile = FALSE) {
  got <- sha(path)
  if (identical(got, want)) return(TRUE)
  msg <- sprintf("SHA-256 of %s is %s, not the recorded %s.", path, got, want)
  if (volatile) {
    warning(msg, " The file is generated at each request; ",
            "data-raw/prepare.R --check compares its values.", call. = FALSE)
    return(FALSE)
  }
  stop(msg, call. = FALSE)
}

for (s in sets) {
  d <- src[[s]]
  dir <- file.path(opt[["dest"]], d[["dir"]])
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  message(s, ": ", dir)
  for (f in d[["files"]]) {
    path <- file.path(dir, f[["file"]])
    if (!file.exists(path) || !identical(sha(path), f[["sha256"]])) {
      obtain(f, path)
    }
    check(path, f[["sha256"]], isTRUE(f[["volatile"]]))
    for (e in f[["extract"]]) {
      inner <- file.path(dir, e[["archive"]])
      if (!file.exists(inner) || !identical(sha(inner), e[["sha256"]])) {
        utils::unzip(path, files = e[["archive"]], exdir = dir)
      }
      check(inner, e[["sha256"]])
      member <- file.path(dir, e[["member"]])
      if (!file.exists(member) ||
          !identical(sha(member), e[["member_sha256"]])) {
        utils::unzip(inner, files = e[["member"]], exdir = dir)
      }
      check(member, e[["member_sha256"]])
    }
  }
}
message("Source files in ", opt[["dest"]], ".")
