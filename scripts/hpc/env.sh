# Sourced by the batch scripts here. Everything that grows lives on scratch:
# a cluster's home is small, and an image cache is large and can be rebuilt.
# Override any of these before submitting.
export HT_SCRATCH="${HT_SCRATCH:-/scratch/alpine/$USER/hydroturing}"
export HT_VENV="${HT_VENV:-$HT_SCRATCH/venv}"
export HT_RUNS="${HT_RUNS:-$HT_SCRATCH/runs}"
export HT_APPTAINER_IMAGES="${HT_APPTAINER_IMAGES:-$HT_SCRATCH/sif}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$HT_SCRATCH/cache/apptainer}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$HT_SCRATCH/tmp}"
export TMPDIR="${TMPDIR_HT:-$HT_SCRATCH/tmp}"
export UV_CACHE_DIR="${UV_CACHE_DIR:-$HT_SCRATCH/cache/uv}"
export UV_PYTHON_INSTALL_DIR="${UV_PYTHON_INSTALL_DIR:-$HT_SCRATCH/cache/uv-python}"
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$HT_RUNS" "$HT_APPTAINER_IMAGES" "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
. "$HT_VENV/bin/activate"
export PYTHONPATH="${HT_REPO:-$PWD}/src${PYTHONPATH:+:$PYTHONPATH}"
