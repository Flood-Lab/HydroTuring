# HBV-EC as emulated by Raven v4.12, for HydroTuring.
#
# The executable is upstream's own CMake target, built unmodified at the commit
# named in issue #116. The configuration is Raven's HBV-EC tutorial model,
# Alouette2 from RavenTutorialFiles.zip, which issue #207 asked for in place of the
# benchmarking files used before. Its four model files are copied into the image,
# byte for byte, and ht_adapter.py applies a named list of edits to them for each
# case. The network is available at build time only; nothing is fetched at run time.
# Raven is Artistic-2.0; the licence travels with the image.
#
# Every stage uses the same python:3.11-slim index (Debian 13), pinned by digest. The
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

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534 AS tutorial
# The Raven tutorial archive (last modified 2025-01-30). Its URL carries no version, so
# the four Alouette2 model files are pinned by content: if the URL serves other files,
# the build fails rather than evaluating a different configuration. The download and
# the check are one step, so changing a pin refetches the archive.
ARG TUTORIAL_URL=https://raven.uwaterloo.ca/files/RavenTutorialFiles.zip
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl unzip \
 && rm -rf /var/lib/apt/lists/*
RUN curl -sSL -o /tmp/tutorial.zip "${TUTORIAL_URL}" \
 && unzip -q /tmp/tutorial.zip 'RavenTutorialFiles/Alouette2/Alouette2.rv[iphc]' -d /tmp/tutorial \
 && rm /tmp/tutorial.zip \
 && mkdir /template && mv /tmp/tutorial/RavenTutorialFiles/Alouette2/Alouette2.rv? /template/ \
 && cd /template \
 && echo "3b1ac60ff2414b8926b2fcb7252d90be2ee17c66c22f71a925eb51df6a3f7559  Alouette2.rvi" | sha256sum -c - \
 && echo "a305fcf4a8eb6daa9b09234a3b082d66a045ebfd9b8e001af8a42d16f6d3a0e0  Alouette2.rvp" | sha256sum -c - \
 && echo "c434ff7010541fc21418998cc5f6c4bde42f67dd4c6065ff7b52b124a6c4985c  Alouette2.rvh" | sha256sum -c - \
 && echo "e3c23271e122179516dfdca23807fafce831cfa4a82ed2c9c05449b8e87771e9  Alouette2.rvc" | sha256sum -c -

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
RUN apt-get update \
 && apt-get install -y --no-install-recommends netcdf-bin \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/build/Raven /model/Raven
COPY --from=build /src/LICENSE /model/LICENSE
COPY --from=tutorial /template/ /model/template/
WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
