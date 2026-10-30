# Running on a cluster

A shared cluster usually has no Docker daemon, and no root to start one. Apptainer
(formerly Singularity) puts a model behind the same isolation without either:

```bash
ht verify-adapter --model cwatm --runner apptainer
ht run --model cwatm --runner apptainer --gate-seeds --csv models/result.csv
```

A model's manifest still says `runner: docker`, and its `Dockerfile` is still the
only recipe. `src/hydroturing/runner/apptainer_runner.py` translates it, stage by
stage, into Apptainer definition files, and a Dockerfile instruction it cannot
translate is an error rather than an image that quietly differs. The build needs
`apptainer build --fakeroot` to work, which needs unprivileged user namespaces, and
it needs network access to the registries. The image goes under
`HT_APPTAINER_IMAGES`, named by the model's version and a digest of its build
context, so a changed adapter rebuilds and an unchanged one does not.

## What holds, and what does not

| Docker | Apptainer | |
| --- | --- | --- |
| `--network none` | `--net --network none` | the same: loopback only |
| `--read-only` | a SIF is read-only | the same |
| `--tmpfs /tmp` (64 MB) | `--containall` | a tmpfs limited by the site's `sessiondir max size`, 64 MB by default |
| only `/io` mounted | `--containall --no-home --no-mount bind-paths,home,cwd,hostfs` | the same: without `--no-mount` a site's bind paths would expose its project filesystems, and the probes with them |
| `--user uid:gid` | always the invoking user | the same |
| `--cap-drop ALL`, `no-new-privileges` | `--no-privs` | the same for a non-root user |
| `--memory`, `--cpus`, `--pids-limit` | **not applied** | Apptainer needs cgroup delegation a batch node does not grant, so the scheduler's limits on the job are the cap |

`model.yaml`'s `resources` therefore do not bind on this runner. Ask the scheduler
for what the model needs. A CPU-bound model that spawns one thread per core sees
the node's cores unless its own environment says otherwise, so a Dockerfile should
keep setting `OMP_NUM_THREADS` and its kin, as the ones here do.

## On a Slurm cluster

`scripts/hpc/` has the two jobs this takes, set up for CU Boulder's Alpine and easy
to change for elsewhere. `env.sh` keeps everything that grows on scratch; it expects
a virtual environment with the harness's dependencies at `$HT_VENV`.

```bash
# Build the image and check the adapter. A build needs network access; Alpine's CPU
# nodes have it, but check that yours do.
sbatch --chdir="$PWD" scripts/hpc/build_model.sbatch cwatm

# Evaluate on every probe, side by side, one task per probe.
sbatch --chdir="$PWD" --job-name=ht-cwatm scripts/hpc/run_model.sbatch cwatm
```

The evaluation writes one report per probe under `$HT_RUNS/<model>/probes/`, and
gathers the archive rows in probe order into `$HT_RUNS/<model>/result.csv`, for
appending to `models/result.csv`. Do not run either on a login node: Apptainer is
often not installed there, and a model evaluation is what those nodes ask you not to
run.

A row records the runner it came from, in the `runner` column of the archive, so an
Apptainer evaluation can be told from a Docker one.
