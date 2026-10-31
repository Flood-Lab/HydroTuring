# SHUD for HydroTuring.
#
# SHUD (MIT) is compiled unmodified from the pinned tag with upstream's own
# Makefile against CVODE 6.0.0, the release its ./configure asks for, whose
# tarball is checked against the sha256 recorded here before ./configure finds
# it and builds it. The network is available at build time only.
#
# Upstream's Makefile compiles its optional BoomerAMG solver wrapper into every
# binary, so `make shud` needs Hypre, MPI and OpenBLAS headers and libraries
# even though the default solver path, which is the one used here, never calls
# them (SHUD_LINSOL unset selects SPGMR). They come from Debian 13's archive
# at build time. The image is one stage because those libraries are needed at
# load time and the build tools are small beside them.

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534 AS build
# Tag cpu-accel-v1.1.1, the revision named in issue #106.
ARG REVISION=c7404df0ebaffd67ea267bce32fff36ba36f6cfe
ARG CVODE_SHA256=c968eb5ee891618a423fe49dc82349bfd1fee5c904cd968fff3aaf33f433c48f
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1 PYTHONDONTWRITEBYTECODE=1
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates git curl cmake make gcc g++ libc6-dev libhypre-dev libopenmpi-dev libopenblas-dev \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src && cd /src \
 && git remote add origin https://github.com/SHUD-System/SHUD.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${REVISION}"
# ./configure uses a cvode.tar.gz it finds beside itself, and otherwise downloads the same URL.
RUN cd /src \
 && curl -sSL -o cvode.tar.gz https://github.com/LLNL/sundials/releases/download/v6.0.0/cvode-6.0.0.tar.gz \
 && echo "${CVODE_SHA256}  cvode.tar.gz" | sha256sum -c - \
 && ./configure \
 && make HYPRE_INCDIR=/usr/include/hypre HYPRE_LIBDIR=/usr/lib/x86_64-linux-gnu shud \
 && test -x shud
RUN pip install numpy==2.1.3

FROM build
ENV OMP_NUM_THREADS=1 OMPI_MCA_btl=self
RUN mkdir -p /model /model/upstream \
 && cp /src/shud /model/shud \
 && cp /src/LICENSE /model/LICENSE \
 && cp -r /src/input/ccw /model/upstream/ccw \
 && rm -rf /src/input/ccw/forcing.csv /src/input/qhh /src/input/heihe /src/cvode.tar.gz /src/cvode-6.0.0 /src/InstallSundials/builddir /src/InstallSundials/example
WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
