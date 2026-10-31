# LISFLOOD-FP v8.1, CPU FV1 solver, for HydroTuring.
#
# The Zenodo archive named in issue #77, checked against its published MD5, compiled
# with upstream's CMake. The one departure from `config.default.cmake` is upstream's own
# documented switch: _NUMERIC_MODE=1, double precision, because a record of 3.8e6 s in
# single precision resolves time to about 0.25 s and a bed of 40 m to about 4e-6 m.
# The source is not modified. README.md states what the archive says about its licence.
# Never pushed to a registry.
FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534 AS build
ARG ARCHIVE_URL=https://zenodo.org/api/records/6912932/files/LISFLOOD-FP%20v8.1.zip/content
ARG ARCHIVE_MD5=77fa78a9c074d7baf9fbb2157fac4090
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl unzip cmake make g++ libnuma-dev libnetcdf-dev \
 && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /src && cd /src \
 && curl -fsSL -o lisflood.zip "${ARCHIVE_URL}" \
 && echo "${ARCHIVE_MD5}  lisflood.zip" | md5sum -c - \
 && unzip -q lisflood.zip -x 'LISFLOOD-FP-trunk/testing/*' 'LISFLOOD-FP-trunk/windep/*' \
 && rm lisflood.zip
# Not in the archive's own build/, which holds a developer's stale CMakeCache.txt.
RUN cd /src/LISFLOOD-FP-trunk \
 && sed 's/_NUMERIC_MODE=0/_NUMERIC_MODE=1/' config.default.cmake > config.double.cmake \
 && cmake -S . -B /src/out -D_CONFIG=config.double.cmake -DCMAKE_BUILD_TYPE=Release \
 && cmake --build /src/out --target lisflood -j4 \
 && test -x /src/out/lisflood

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534
ENV PYTHONUNBUFFERED=1
ENV OMP_NUM_THREADS=1
RUN apt-get update \
 && apt-get install -y --no-install-recommends libnuma1 libgomp1 libnetcdf-dev \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/out/lisflood /model/lisflood
COPY --from=build /src/LISFLOOD-FP-trunk/LICENCE-3RD-PARTY.txt /model/LICENCE-3RD-PARTY.txt
WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
