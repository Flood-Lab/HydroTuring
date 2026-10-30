#!/bin/bash
# Evaluate one model on every probe, or the ones named, on the gate seeds.
#
#   sbatch --chdir=<repo> --job-name=ht-<model> scripts/hpc/run_model.sbatch <model> [probe ...]
#
# Probes run side by side, as many as fit the CPUs and memory the manifest asks
# for per case, since Apptainer cannot enforce those itself and a model that runs
# out of either comes back as an ERROR that says nothing about it. The image is
# read from node-local disk: on the shared filesystem a model takes several times
# longer to start, enough to cross a probe's time budget under load.
#
# Each probe writes its own report; the rows are gathered, in probe order, into
# $HT_RUNS/<model>/result.csv, ready to append to models/result.csv once they are
# worth archiving.
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=3700M
#SBATCH --time=23:30:00
#SBATCH --output=%x-%j.out
set -u
model="$1"; shift
HT_REPO="$PWD"; . "$HT_REPO/scripts/hpc/env.sh"
out="$HT_RUNS/$model"; mkdir -p "$out/probes"
if [ "$#" -gt 0 ]; then
  probes=("$@")
else
  mapfile -t probes < <(find probes -mindepth 3 -maxdepth 3 -name probe.yaml | sed -e 's#^probes/##' -e 's#/probe.yaml$##' | sort)
fi
# Build first and once, then copy: the probes would otherwise queue on the lock.
read -r image parallel < <(python - "$model" <<'P'
import os, sys
from dataclasses import replace
from hydroturing import registry
from hydroturing.runner.apptainer_runner import build
model = registry.find_model(sys.argv[1])
image = "-"
if model.runner != "subprocess":
    image = build(replace(model, runner="apptainer"))
resources = model.resources or {}
tasks = int(os.environ.get("SLURM_NTASKS", 8))
gigabytes = tasks * 3.6
print(image, max(1, min(tasks // max(1, int(resources.get("cpu", 2))),
                        int(gigabytes // max(1, resources.get("memory_gb", 4))))))
P
)
if [ "$image" != "-" ] && [ -n "${SLURM_SCRATCH:-}" ]; then
  cp "$image" "$SLURM_SCRATCH/" && export HT_APPTAINER_IMAGES="$SLURM_SCRATCH"
fi
echo "image: $image; $parallel probes side by side"
one() {
  local probe="$1" tag="${1//\//__}" runner=()
  [ "$(grep -E '^runner:' "models/$model/model.yaml" | awk '{print $2}')" = "subprocess" ] || runner=(--runner apptainer)
  rm -f "$out/probes/$tag.csv"
  python -m hydroturing run --model "$model" --probe "$probe" --gate-seeds "${runner[@]}" \
    --json "$out/probes/$tag.json" --csv "$out/probes/$tag.csv" --workdir "$TMPDIR/$model-$SLURM_JOB_ID/$tag" \
    > "$out/probes/$tag.txt" 2>&1
  echo "$probe exit=$?"
  rm -rf "$TMPDIR/$model-$SLURM_JOB_ID/$tag"
}
export -f one; export model out
printf '%s\n' "${probes[@]}" | xargs -P "$parallel" -I{} bash -c 'one {}'
python - "$out" <<'P'
import csv, glob, os, sys
out = sys.argv[1]
rows, header = [], None
for path in sorted(glob.glob(os.path.join(out, "probes", "*.csv"))):
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        header = reader.fieldnames
        rows += list(reader)
with open(os.path.join(out, "result.csv"), "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=header, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
tally = {}
for row in rows:
    key = (row["verdict"], row["reason"])
    tally[key] = tally.get(key, 0) + 1
print(sorted(tally.items()))
for row in rows:
    if row["verdict"] != "PASS" and row["verdict"] != "N/A":
        print(row["probe"], row["verdict"], row["reason"], row["detail"][:220])
P
