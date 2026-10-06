## scripts/01_simulate.R -- fit the units of one design.
##
##   Rscript scripts/01_simulate.R --config=config/core.yaml --workers=8
##
## Options in R/cli.R. Every unit (cell, n, method, replicate) is cached
## under <out>/units/ as soon as it finishes, and a second run resumes
## where the first stopped. The configuration and the output directory are
## arguments, so a run with another configuration writes elsewhere and
## leaves the outputs of the study alone.

source(file.path("scripts", "00_setup.R"))
cli <- parse_cli()
cfg <- study_config(cli[["config"]], cli[["study"]])
out <- run_dir(cli, cfg)

if (cli[["list"]]) {
  u <- study_units(cfg, cells = cli[["cells"]], sample_sizes = cli[["sizes"]],
                   methods = cli[["methods"]], reps = cli[["reps"]])
  print(stats::aggregate(rep ~ cell + n + method, u, length))
  message(nrow(u), " unit(s); output directory ", out)
} else {
  load_wafc_code(cfg)
  run_study(cfg, out, workers = run_workers(cli, cfg),
            refresh = cli[["refresh"]], cells = cli[["cells"]],
            sample_sizes = cli[["sizes"]], methods = cli[["methods"]],
            reps = cli[["reps"]])
}
