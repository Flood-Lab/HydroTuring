#!/bin/bash
# Build one model's image, then check its adapter against the contract.
#
#   sbatch --chdir=<repo> scripts/hpc/build_model.sbatch <model> [probe]
#
# On the compile partition: its nodes reach the network, which a build needs and
# an evaluation does not.
#SBATCH --job-name=ht-build
#SBATCH --partition=acompile
#SBATCH --qos=compile
#SBATCH --ntasks=2
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
set -u
model="$1"; probe="${2:-}"
HT_REPO="$PWD"; . "$HT_REPO/scripts/hpc/env.sh"
python -m hydroturing verify-adapter --model "$model" --runner apptainer ${probe:+--probe "$probe"}
