# RVIC 1.1.1 from UW-Hydro, for HydroTuring.
#
# The tag rvic.1.1.1, installed unmodified from its source. RVIC predates recent NumPy
# and netCDF4: `np.float` and a real datetime from netCDF4.num2date both need the
# versions pinned below, which are the last that keep them on a Python with wheels.
# RVIC is GPL-3.0; LICENSE.txt is copied into the image, which is never pushed anywhere.
# matplotlib is left out: RVIC imports it optionally and uses it for diagnostic plots
# only. CDO is left out: it is needed only for REMAP=True, and the grid here is not remapped.
# Pinned to linux/amd64, as models/lisflood/Dockerfile is:
# the pinned NumPy-era stack (netCDF4 1.4.2 and friends) has no linux/arm64 wheels,
# and building them from source needs HDF5 headers this image does not carry.
# Native on an x86-64 host, emulated on an arm64 one.
FROM --platform=linux/amd64 python:3.7-slim-bookworm@sha256:b53f496ca43e5af6994f8e316cf03af31050bf7944e0e4a308ad86c001cf028b

ENV PYTHONUNBUFFERED=1
ENV OMP_NUM_THREADS=1
ENV OPENBLAS_NUM_THREADS=1

# The commit the tag rvic.1.1.1 points at.
ARG REVISION=afde3841d9b8de243c88a9e539cb3bf13b0fe313

RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates git gcc libc6-dev \
 && rm -rf /var/lib/apt/lists/*

RUN git init -q /src && cd /src \
 && git remote add origin https://github.com/UW-Hydro/RVIC.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git -c filter.lfs.smudge= -c filter.lfs.required=false checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${REVISION}"

RUN pip install --no-cache-dir numpy==1.16.6 scipy==1.2.3 pandas==0.24.2 cftime==1.0.4.2 netCDF4==1.4.2 \
 && pip install --no-cache-dir --no-deps /src \
 && python -c "import rvic.parameters, rvic.convolution, rvic_convolution" 2>&1 | tail -3 \
 && mkdir -p /model && cp /src/LICENSE.txt /model/LICENSE

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
