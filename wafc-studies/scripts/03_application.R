## scripts/03_application.R -- the two data applications.
##
##   Rscript scripts/03_application.R --workers=3
##
## Fits every unit of config/application.yaml (each method on the whole
## sample and on each partition by block of each application; R/application.R)
## and writes the tables and figures under <out>/tables/ and <out>/figures/
## (R/application_report.R). Every unit is cached under <out>/units/ as soon
## as it finishes, and a second run resumes where the first stopped.
##
## Options (R/cli.R and application_cli()): --config (default
## config/application.yaml), --study, --out (default outputs/application),
## --workers, --methods, --refresh, --list, and --bases=a,b,
## --splits=FROM:TO (0 is the whole sample) and --parts=fit,report.

source(file.path("scripts", "00_setup.R"))
cli <- application_cli()
cfg <- application_config(cli[["config"]], cli[["study"]])
out <- run_dir(cli, cfg)

if (cli[["list"]]) {
  u <- application_units(cfg, bases = cli[["bases"]], splits = cli[["splits"]],
                         methods = cli[["methods"]])
  print(stats::aggregate(split ~ base + method, u, length))
  message(nrow(u), " unit(s); output directory ", out)
} else {
  load_wafc_code(cfg)
  if ("fit" %in% cli[["parts"]]) {
    run_application(cfg, out, workers = run_workers(cli, cfg),
                    refresh = cli[["refresh"]], bases = cli[["bases"]],
                    splits = cli[["splits"]], methods = cli[["methods"]])
  }
  if ("report" %in% cli[["parts"]]) {
    application_report(cfg, out, bases = cli[["bases"]])
  }
}
