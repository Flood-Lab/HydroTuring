# mHM v5.13.4 for HydroTuring.
#
# mHM and FORCES (the library it builds against) are compiled unmodified from
# the pinned commits with upstream's own CMake. The network is available at
# build time only; the container runs without one and nothing is fetched at run
# time. GPL-3.0-or-later (mHM), LGPL-3.0 (FORCES); the licences are copied into
# the image, which is never pushed to a registry.
#
# The morphological inputs of upstream's test domain are baked in from the same
# commit, and derive_domain.py cuts the sub-basin the adapter runs on from them
# at build time. No model input is committed to this repository.

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534 AS build
# Tag v5.13.4.
ARG REVISION=5f2de41bf244e9bd0d4b7a27806e0c96ff8186db
# Tag v0.10.0, the version named by version_forces.txt at that tag.
ARG FORCES_REVISION=1e0194777d35ddc4b0158aab16522dcaabc4c378
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates git cmake make gcc gfortran libc6-dev libnetcdf-dev libnetcdff-dev \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src && cd /src \
 && git remote add origin https://github.com/mhm-ufz/mHM.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${REVISION}" \
 && test "$(cat version_forces.txt)" = "v0.10.0"
# Upstream's CMake takes a checkout at ./forces in place of downloading the tag.
RUN git init -q /src/forces && cd /src/forces \
 && git remote add origin https://github.com/mhm-ufz/FORCES.git \
 && git fetch -q --depth 1 origin "${FORCES_REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${FORCES_REVISION}"
# Upstream's own build, as a subdirectory, with one more executable linked to the
# library it makes. ht_mhm_driver.f90 says what that executable adds.
COPY CMakeLists.txt /src/ht/CMakeLists.txt
COPY ht_mhm_driver.f90 /src/ht/ht_mhm_driver.f90
RUN cmake -S /src/ht -B /src/build -DCMAKE_BUILD_TYPE=Release -DBUILD_MHM_DRIVER=ON -DCMAKE_WITH_OpenMP=OFF -DCMAKE_WITH_MPI=OFF \
 && cmake --build /src/build --target mhm ht_mhm_driver -j4 \
 && find /src/build -type f -executable \( -name mhm -o -name ht_mhm_driver \) -exec cp {} /src/ \; \
 && strip /src/mhm /src/ht_mhm_driver

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
RUN apt-get update \
 && apt-get install -y --no-install-recommends libnetcdff7 libgfortran5 \
 && rm -rf /var/lib/apt/lists/*
RUN pip install numpy==2.1.3 netCDF4==1.7.2
COPY --from=build /src/mhm /model/mhm
COPY --from=build /src/ht_mhm_driver /model/ht_mhm_driver
COPY --from=build /src/LICENSE.md /model/LICENSE.md
COPY --from=build /src/COPYING /model/COPYING
COPY --from=build /src/COPYING.LESSER /model/COPYING.LESSER
COPY --from=build /src/test_domain /model/upstream/test_domain
COPY derive_domain.py /model/derive_domain.py
COPY midrange_parameters.py /model/midrange_parameters.py
COPY --from=build /src/mhm_parameter.nml /model/upstream/mhm_parameter_root.nml
RUN python3 /model/derive_domain.py /model/domain
# Only README.md's sensitivity runs read these; a case never does.
RUN mkdir -p /model/sensitivity \
 && python3 /model/derive_domain.py /model/sensitivity/8km 92 65 8000 \
 && python3 /model/derive_domain.py /model/sensitivity/monthly_lai 92 65 4000 monthly \
 && python3 /model/derive_domain.py /model/sensitivity/mosel_24km 32 169 24000 \
 && python3 /model/midrange_parameters.py /model/upstream/test_domain/mhm_parameter.nml /model/sensitivity/midrange.nml
WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
