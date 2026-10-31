# mizuRoute v3.1.1 from ESCOMP, for HydroTuring.
#
# Built from the pinned tag commit, unmodified, with upstream's own Makefile, which
# in turn builds ParallelIO at the tag named in upstream's .gitmodules. Only the
# stand-alone executable is used; ht_adapter.py writes its control file and inputs
# and reads its history file. The network is available at build time only.
# mizuRoute is Apache-2.0; ParallelIO is Apache-2.0. Both LICENSE files are in the image.
#
# One stage: the executable links dynamically against Debian's MPI and netCDF, and
# keeping one stage keeps those libraries certain to be the ones it was linked with.
FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534

ENV PYTHONUNBUFFERED=1
ENV OMP_NUM_THREADS=1
ENV OMPI_MCA_btl_vader_single_copy_mechanism=none
ENV OMPI_MCA_rmaps_base_oversubscribe=1

# The tag commit of v3.1.1, as named in issue #132.
ARG REVISION=28514dea3fd28a2f820c7ecaa9014f820c1bb220
# pio2_6_6, the fxtag in upstream's .gitmodules at that commit.
ARG PIO_REVISION=1d516b82768ef9d2cade524724c38ab63fbe08e4

RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates git cmake make gcc g++ gfortran perl \
      libnetcdf-dev libnetcdff-dev mpich libmpich-dev \
 && rm -rf /var/lib/apt/lists/*

RUN git init -q /src && cd /src \
 && git remote add origin https://github.com/ESCOMP/mizuRoute.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${REVISION}" \
 && mkdir -p /src/libraries/parallelio && cd /src/libraries/parallelio \
 && git init -q . \
 && git remote add origin https://github.com/NCAR/ParallelIO.git \
 && git fetch -q --depth 1 origin "${PIO_REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${PIO_REVISION}"

# MODE=fast is upstream's optimised build. No pnetcdf and no OpenMP: the network here
# is two reaches, and upstream documents netCDF-only ParallelIO as sufficient for a
# small domain. GPTL stays on because the driver's perf_mod is part of it. None of
# those options changes the routing equations.
# Upstream's Makefile and ParallelIO's CMake look under one prefix for lib/ and include/,
# which Debian's multiarch layout does not have, so /opt/netcdf is that prefix in links.
# The Makefile lists netCDF before the static ParallelIO libraries that need it, which
# this linker's --as-needed drops; LIBNETCDF says the same thing with that turned off.
RUN mkdir -p /opt/netcdf \
 && ln -s /usr/include /opt/netcdf/include \
 && ln -s /usr/lib/x86_64-linux-gnu /opt/netcdf/lib \
 && cd /src/route/build \
 && export CC=mpicc \
 && make FC=gnu FC_EXE=mpif90 EXE=route_runoff F_MASTER=/src/route/ NCDF_PATH=/opt/netcdf \
      MODE=fast isOpenMP=no isPIO=yes isGPTL=yes PIO_FILESYSTEM_HINTS= \
      LIBNETCDF="-L/opt/netcdf/lib -Wl,--no-as-needed -lnetcdff -lnetcdf" \
 && test -x /src/route/bin/route_runoff

RUN pip install --no-cache-dir numpy==2.2.6 netCDF4==1.7.2

RUN mkdir -p /model \
 && cp /src/route/bin/route_runoff /model/route_runoff \
 && cp /src/route/ancillary_data/param.nml.default /model/param.nml.default \
 && cp /src/LICENSE.txt /model/LICENSE \
 && cp /src/libraries/parallelio/LICENSE /model/LICENSE.parallelio

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
