# CFE (Conceptual Functional Equivalent) from NOAA-OWP, for HydroTuring.
#
# Upstream's shared library is built from the pinned commit, unmodified, with
# upstream's own CMake target. ht_cfe_driver.c links against it and only feeds
# it and reads it. The network is available at build time only; the container
# runs without one and nothing is fetched at run time. CFE is Apache-2.0; the
# licence text is copied into the image beside the adapter.
#
# Both stages use the same python:3.11-slim index, pinned by digest, so the
# runtime has the glibc the library was linked against. That index is Debian 13, whose
# archive the build tools resolve from at build time (the evaluated image used cmake
# 3.31.6 and gcc 14.2); none of them is in the final image.

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534 AS build
# NOAA-OWP/cfe at the revision named in issue #127.
ARG REVISION=f9182dfb3d81c407c66c5b782114173893c30cf2
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates git cmake make gcc g++ libc6-dev \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src && cd /src \
 && git remote add origin https://github.com/NOAA-OWP/cfe.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${REVISION}"
COPY ht_cfe_driver.c /src/ht_cfe_driver.c
# -DNGEN=ON is upstream's option for the BMI shared library, and the only one
# that needs neither submodule. Without any option upstream's CMake refuses to
# configure. The test executables it declares are not built.
RUN cmake -S /src -B /src/build -DNGEN=ON -DCMAKE_BUILD_TYPE=Release \
 && cmake --build /src/build --target cfebmi -j4 \
 && gcc -O2 -I/src/include /src/ht_cfe_driver.c -L/src/build -lcfebmi -lm \
      -Wl,-rpath,/model -o /src/ht_cfe_driver

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534
ENV PYTHONUNBUFFERED=1
COPY --from=build /src/build/libcfebmi.so.1.0.0 /model/libcfebmi.so.1.0.0
COPY --from=build /src/ht_cfe_driver /model/ht_cfe_driver
COPY --from=build /src/LICENSE /model/LICENSE
WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
