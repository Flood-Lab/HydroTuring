# SUMMA (Structure for Unifying Multiple Modeling Alternatives) for HydroTuring.
#
# The model is compiled here from the pinned release and nothing is fetched
# at run time: the container runs with no network. SUMMA is GPL-3.0; the
# licence text is copied into the image beside the adapter.

FROM debian:bookworm-slim AS build
# CH-Earth/summa, release v4.0.0 (2026-09-04), the revision named in issue #30.
ARG SUMMA_REPO=https://github.com/CH-Earth/summa.git
ARG SUMMA_TAG=v4.0.0
ARG SUMMA_COMMIT=f787fa5e63d3c67030721a85f2f244d22a8d691e
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates git cmake make gfortran gcc g++ \
      libnetcdff-dev libnetcdf-dev libopenblas-dev liblapack-dev \
 && rm -rf /var/lib/apt/lists/*
# The tag is fetched so SUMMA's build stamps its own version string, and the
# checkout is refused unless the tag still points at the pinned commit.
RUN git init -q /src/summa && cd /src/summa \
 && git remote add origin "${SUMMA_REPO}" \
 && git fetch -q --depth 1 origin "refs/tags/${SUMMA_TAG}:refs/tags/${SUMMA_TAG}" \
 && git checkout -q "${SUMMA_TAG}" \
 && test "$(git rev-parse HEAD)" = "${SUMMA_COMMIT}"
# The release's default build: CMake, Release flags, no SUNDIALS (the
# homegrown backward-Euler solver needs none), no NextGen, no OpenWQ.
RUN cmake -B /src/cmake_build -S /src/summa/build -DCMAKE_BUILD_TYPE=Release \
 && cmake --build /src/cmake_build --target all -j4 \
 && /src/summa/bin/summa.exe -v

FROM debian:bookworm-slim
ENV PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
RUN apt-get update && apt-get install -y --no-install-recommends \
      libgfortran5 libnetcdff7 libopenblas0-serial \
      python3 python3-numpy python3-netcdf4 \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/summa/bin/summa.exe /opt/summa/bin/summa.exe
COPY --from=build /src/cmake_build/libsumma.so /usr/local/lib/libsumma.so
# The only complete SUMMA setup the release ships: the decisions, the default
# parameter tables and the Noah-MP lookup tables of the CAMELS basin 01073000
# test case. The adapter reads them from here; nothing in them is edited.
COPY --from=build /src/summa/test_ngen/gauge_01073000/settings/SUMMA /model/summa_settings
COPY --from=build /src/summa/COPYING /model/LICENSE
COPY --from=build /src/summa/LICENSE.txt /model/LICENSE.txt
RUN ldconfig && chmod -R a+rX /model /opt/summa

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
