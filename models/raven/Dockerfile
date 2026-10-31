# Raven v4.12, in its HBV-EC emulation, for HydroTuring.
#
# The executable is upstream's own CMake target, built unmodified at the commit
# named in issue #116. Upstream's Salmon River HBV-EC emulation files are copied
# into the image, byte for byte, and ht_adapter.py applies a named list of edits
# to them for each case. The network is available at build time only; nothing is
# fetched at run time. Raven is Artistic-2.0; the licence travels with the image.
#
# Both stages use the same python:3.11-slim index (Debian 13), pinned by digest. The
# C++ runtime is linked statically, because the slim runtime image does not carry it.
# NetCDF is optional to upstream's CMake and is installed for its double-precision output:
# Raven's CSV files carry six significant figures (ht_adapter.py). netcdf-bin brings the
# library and ncdump into the runtime, so the adapter stays standard library.

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534 AS build
# CSHS-CWRA/RavenHydroFramework, v4.12, released 2026-03-05.
ARG REVISION=81011f061e53ef1289ebcab36dc3622cea05c76d
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates git cmake make g++ libnetcdf-dev \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src && cd /src \
 && git remote add origin https://github.com/CSHS-CWRA/RavenHydroFramework.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${REVISION}"
RUN cmake -S /src -B /src/build -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXE_LINKER_FLAGS="-static-libstdc++ -static-libgcc" \
 && cmake --build /src/build --target Raven -j4 \
 && strip /src/build/Raven

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
RUN apt-get update \
 && apt-get install -y --no-install-recommends netcdf-bin \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/build/Raven /model/Raven
COPY --from=build /src/LICENSE /model/LICENSE
COPY --from=build /src/benchmarking/_InputFiles/Salmon_HBV/raven-hbv-salmon.rvi /model/template/raven-hbv-salmon.rvi
COPY --from=build /src/benchmarking/_InputFiles/Salmon_HBV/raven-hbv-salmon.rvp /model/template/raven-hbv-salmon.rvp
COPY --from=build /src/benchmarking/_InputFiles/Salmon_HBV/raven-hbv-salmon.rvh /model/template/raven-hbv-salmon.rvh
COPY --from=build /src/benchmarking/_InputFiles/Salmon_HBV/raven-hbv-salmon.rvc /model/template/raven-hbv-salmon.rvc
COPY --from=build /src/benchmarking/_InputFiles/Salmon_HBV/raven-hbv-salmon.rvt /model/template/raven-hbv-salmon.rvt
WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
