## run_all.R -- the whole simulation study, in one call.
##
##   Rscript run_all.R            (from the root of the compendium)
##
## For each design listed in `designs` of config/study.yaml: fits every
## unit (scripts/01_simulate.R) and writes its tables
## (scripts/02_aggregate.R), under outputs/<design>/. Every unit is cached
## as it finishes, so an interrupted run is resumed by running this file
## again. The full study takes days of processor time; README.md gives the
## cost and how to spread it over workers.

source(file.path("scripts", "00_setup.R"))
study <- file.path("config", "study.yaml")
designs <- unlist(yaml::read_yaml(study)[["designs"]])
for (d in designs) {
  cfg <- study_config(file.path("config", paste0(d, ".yaml")), study)
  out <- file.path("outputs", d)
  message("\n>>> design ", d)
  load_wafc_code(cfg)
  run_study(cfg, out, workers = run_workers(list(), cfg))
  aggregate_study(cfg, out)
}
message("\nDone. Tables under outputs/<design>/tables/.")
