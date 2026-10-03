## wafc/tests/helper-slow.R -- the switch of the slow tests (step E3.1).
##
## testthat sources every helper-*.R file of the directory before the tests.
## The default suite has to run in under a minute (plano-projeto.md, E3.1),
## so the tests that hold an earlier form against a frozen copy of it, redo
## by hand the search of a competitor over its grid, or need a large n run
## only when asked for:
##
##     WAFC_SLOW_TESTS=1 Rscript -e 'testthat::test_dir("wafc/tests")'
##
## Without the variable they are skipped, each with the reason printed, so
## no test is lost: the whole suite is the default one plus these. Each such
## test calls skip_slow() as its first line, which is how they are found.

wafc_slow_tests <- function() identical(Sys.getenv("WAFC_SLOW_TESTS"), "1")

skip_slow <- function() {
  testthat::skip_if_not(wafc_slow_tests(),
                        "slow test; set WAFC_SLOW_TESTS=1 to run it")
}
