## scripts/02_aggregate.R -- the tables of one design.
##
##   Rscript scripts/02_aggregate.R --config=config/core.yaml
##
## Reads the cached units under <out>/units/ and writes <out>/tables/ (the
## list of tables is in R/aggregate.R). Takes the options of
## scripts/01_simulate.R that select units; it fits nothing.

source(file.path("scripts", "00_setup.R"))
cli <- parse_cli()
cfg <- study_config(cli[["config"]], cli[["study"]])
load_wafc_code(cfg)
aggregate_study(cfg, run_dir(cli, cfg), cells = cli[["cells"]],
                sample_sizes = cli[["sizes"]], methods = cli[["methods"]],
                reps = cli[["reps"]])
