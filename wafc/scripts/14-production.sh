#!/bin/bash
## ============================================================
## wafc/scripts/14-production.sh -- the production of the simulation study (E4.4)
## ============================================================
## Runs the three designs of the compendium (core, arms, scale) with the
## study's own configuration and master seed, each fitted by
## scripts/01_simulate.R and aggregated by scripts/02_aggregate.R, under
## outputs/<design>/ of wafc-studies/. The data applications are not run
## here (they are done; run_all.R would refit them).
##
##   bash wafc/scripts/14-production.sh [workers]     (from the root of wafc-draft)
##
## Workers default to 10, one per physical core of the machine of E4.4;
## the time of a unit is a measure of the study, so the BLAS is held to one
## thread and the hyperthreads are left idle. Every unit is cached as it
## finishes, so a stopped run is resumed by running this script again.
##
## Marks in wafc-studies/outputs/production.log ("production end" at the
## end); the PID of this script in wafc-studies/outputs/production.pid; the
## machine in wafc-studies/outputs/machine.txt; GNU time records of each
## step in wafc-studies/outputs/steps/.
## ------------------------------------------------------------

set -u
W=${1:-10}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
cd "$REPO/wafc-studies" || exit 1
O=outputs
mkdir -p "$O/steps"
L=$O/production.log
echo $$ > "$O/production.pid"

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

if [ -n "$(git -C "$REPO" status --porcelain -- wafc/R wafc-studies/R wafc-studies/config wafc-studies/scripts)" ]; then
  echo "uncommitted changes in the code or the configuration; commit them first" | tee -a "$L" >&2
  exit 3
fi

{
  echo "# machine of the production, $(date -Iseconds)"
  echo "# wafc-draft at $(git -C "$REPO" rev-parse HEAD)"
  lscpu
  echo
  free -g
  echo
  Rscript -e 'cat(R.version.string, "\n"); print(extSoftVersion()["BLAS"]); print(La_library())'
} > "$O/machine.txt" 2>&1

step() {
  local name=$1; shift
  echo "$name start $(date -Iseconds)" >> "$L"
  /usr/bin/time -v -o "$O/steps/$name.$(date +%Y%m%dT%H%M%S).time" "$@" \
    > "$O/steps/$name.$(date +%Y%m%dT%H%M%S).log" 2>&1
  local code=$?
  echo "$name end $(date -Iseconds) exit=$code" >> "$L"
  return $code
}

echo "production start $(date -Iseconds) workers=$W" >> "$L"
for d in core arms scale; do
  step "fit-$d" Rscript scripts/01_simulate.R --config=config/$d.yaml --workers=$W
  step "aggregate-$d" Rscript scripts/02_aggregate.R --config=config/$d.yaml
done
echo "production end $(date -Iseconds)" >> "$L"
