# LISFLOOD (EC JRC) for HydroTuring.
#
# The model code is EUPL-1.2. The pinned source is fetched at build time for
# its LICENSE; the model itself is the PyPI release of the same tag. The
# container runs with no network and the adapter never downloads anything.
#
# Why amd64: LISFLOOD 5.0.0 imports PCRaster at every level (its dynamic
# framework, the LDD and routing operators), and conda-forge publishes
# PCRaster for linux-64, osx-64, osx-arm64 and win-64 but not for
# linux-aarch64. The image is therefore pinned to linux/amd64: native on the
# CI runner, emulated on an Apple-silicon host. It builds on either.

FROM --platform=linux/amd64 debian:13-slim AS source
# ec-jrc/lisflood-code at tag v5.0.0.
ARG REVISION=025cff03f8304cbe55cdbef3879312cc5950bc91
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src \
 && cd /src \
 && git remote add origin https://github.com/ec-jrc/lisflood-code.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD \
 && cat VERSION

FROM --platform=linux/amd64 mambaorg/micromamba:2.9.0-debian13 AS env
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
# The binary stack at the versions of the upstream environment.yml for this
# tag (python 3.12, pcraster 4.4.2, gdal 3.12.4, numpy 2.2.6, numba 0.65.1).
# One extraction and download thread: the parallel extraction deadlocked
# under amd64 emulation on an arm64 build host.
ENV MAMBA_EXTRACT_THREADS=1 MAMBA_DOWNLOAD_THREADS=1
RUN micromamba install -y -n base -c conda-forge --override-channels \
        python=3.12 pcraster=4.4.2 libgdal-core=3.12.4 \
        numpy=2.2.6 numba=0.65.1 numexpr=2.14.1 scipy=1.17.1 \
        pandas=2.3.3 xarray netcdf4=1.7.4 cftime=1.6.5 pyproj=3.7.2 \
        nine future pip \
 && micromamba clean --all --yes
# LISFLOOD 5.0.0 from PyPI. --no-deps: its requirement list pins the test
# tooling and GDAL's Python bindings, neither of which the model imports.
# Its setup.py calls gdal-config (shipped by libgdal-core) for the GDAL
# version string, so the environment's bin directory must be on PATH.
ENV PATH=/opt/conda/bin:$PATH
RUN pip install --no-deps lisflood-model==5.0.0 \
 && python -c "import lisflood; print('lisflood', lisflood.__version__)"

FROM env
USER root
COPY --from=source /src/LICENSE /model/LICENSE
WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py
# NUMBA_DISABLE_JIT: LISFLOOD's jitted loops run over pixels, and the adapter
# runs one. Compiling them costs tens of seconds per run on a read-only
# container that cannot keep a cache; as Python they close the same budget.
ENV PATH=/opt/conda/bin:$PATH \
    GDAL_DATA=/opt/conda/share/gdal PROJ_DATA=/opt/conda/share/proj \
    OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMBA_DISABLE_JIT=1 \
    NUMEXPR_MAX_THREADS=1 NUMEXPR_NUM_THREADS=1
# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
