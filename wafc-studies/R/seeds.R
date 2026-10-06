## R/seeds.R -- one seed per (cell, sample size, replicate, stream).
##
## The single source of randomness is `seeds$master` of config/study.yaml.
## Every other seed is master + index, where the index is an injective
## mixed-radix encoding of the position of the cell's seed key in
## `seeds$keys`, the position of the sample size in `seeds$sample_sizes`,
## the replicate and the stream, so two units of the study never share a
## seed unless they are meant to. Nothing else calls set.seed().
##
## Three streams per replicate: the training sample, the test sample and
## the folds. The methods of a replicate start from the random stream the
## folds leave behind, so what one method draws does not depend on which
## other methods were run (only VCBART and the permutation test of
## mgcv::k.check() draw anything). A cell's seed key is its name unless it
## names another one: an arm of the study takes the key of the core cell it
## varies, and with it that cell's data.

seed_streams <- c(data = 0L, test = 1L, folds = 2L)

#' Seeds of one replicate
#'
#' @param cfg The configuration of study_config().
#' @param cell A cell of the configuration (a list with `seed_key`).
#' @param n The sample size.
#' @param rep The replicate, between 1 and `seeds$rep_radix - 1`.
#' @return A named list of integer seeds, one per stream of `seed_streams`.
unit_seeds <- function(cfg, cell, n, rep) {
  sd <- cfg[["seeds"]]
  k0 <- match(cell[["seed_key"]], unlist(sd[["keys"]])) - 1L
  n0 <- match(n, unlist(sd[["sample_sizes"]])) - 1L
  if (is.na(k0)) {
    stop("seed key '", cell[["seed_key"]], "' is not in seeds$keys of the ",
         "study configuration; append it there.", call. = FALSE)
  }
  if (is.na(n0)) {
    stop("sample size ", n, " is not in seeds$sample_sizes of the study ",
         "configuration; append it there.", call. = FALSE)
  }
  radix <- as.double(sd[["rep_radix"]])
  if (rep < 1 || rep >= radix) {
    stop("replicate ", rep, " is outside 1, ..., ", radix - 1, ".",
         call. = FALSE)
  }
  nn <- length(sd[["sample_sizes"]])
  ns <- length(seed_streams)
  idx <- ((as.double(k0) * nn + n0) * radix + rep) * ns + seed_streams
  seed <- as.double(sd[["master"]]) + idx
  if (any(seed >= .Machine$integer.max)) {
    stop("seed overflow: lower seeds$master or seeds$rep_radix.",
         call. = FALSE)
  }
  as.list(stats::setNames(as.integer(seed), names(seed_streams)))
}
