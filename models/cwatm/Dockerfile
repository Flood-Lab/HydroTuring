# CWatM, the Community Water Model (IIASA), for HydroTuring.
#
# The network is available at build time only. The pinned source is baked in
# here; the container runs with no network and the adapter downloads nothing.
# CWatM is GPL-3.0; the licence travels with the image at /opt/cwatm/LICENSE.

# python:3.11-slim (Python 3.11.16, Debian trixie) pinned by the index digest
# of the image the archived evaluation was built on.
FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534 AS source
# iiasa/CWatM at release 1.11 (tag 1.11 points at this commit).
ARG REVISION=5baaadd8c4379f4e82956b31c1d4b90d3ebf7e00
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates g++ \
 && rm -rf /var/lib/apt/lists/*
# Only the model package and the licence: the repository's tutorials, toolkit
# and MODFLOW binaries (half a gigabyte) are not needed to run it.
RUN git init -q /src \
 && cd /src \
 && git remote add origin https://github.com/iiasa/CWatM.git \
 && git sparse-checkout set --cone cwatm \
 && git fetch -q --depth 1 --filter=blob:none origin "${REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test -f LICENSE
# CWatM loads its routing and runoff-concentration kernels from a C++ library
# it ships as source (t5.cpp) and as prebuilt x86-64 binaries. Compile the
# shipped source for the build architecture so the image runs on arm64 too.
RUN cd /src/cwatm/hydrological_modules/routing_reservoirs \
 && g++ -O3 -fPIC -shared -static-libstdc++ -static-libgcc -o t5_linux.so t5.cpp \
 && rm -f t5.dll t5cyg.so t5_mac.so

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
# rasterio's wheel links the system expat.
RUN apt-get update \
 && apt-get install -y --no-install-recommends libexpat1 \
 && rm -rf /var/lib/apt/lists/*
# The imports CWatM's daily water-balance path needs. flopy/xmipy (MODFLOW),
# xarray/metpy (pySnowClim) and openpyxl (Excel settings) are not used here.
RUN pip install "numpy==1.26.4" "scipy==1.13.1" "netCDF4==1.7.2" "pandas==2.2.3" "rasterio==1.4.4"
COPY --from=source /src/cwatm /opt/cwatm/cwatm
COPY --from=source /src/LICENSE /opt/cwatm/LICENSE
RUN python -m compileall -q /opt/cwatm/cwatm
ENV PYTHONPATH=/opt/cwatm PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
