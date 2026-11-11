#!/bin/bash
# Build one model's image, then check its adapter against the contract.
#
#   sbatch --chdir=<repo> scripts/hpc/build_model.sbatch <model> [probe]
#
# Alpine's CPU nodes reach the network, which a build needs. Its compile partition
# would too, but caps a user at 8 CPUs and 4 jobs, too few for several models at once.
#SBATCH --job-name=ht-build
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
set -u
model="$1"; probe="${2:-}"
HT_REPO="$PWD"; . "$HT_REPO/scripts/hpc/env.sh"
python -m hydroturing verify-adapter --model "$model" --runner apptainer ${probe:+--probe "$probe"}
