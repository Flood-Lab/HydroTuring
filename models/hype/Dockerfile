# HYPE 5.36.0 from SMHI, for HydroTuring.
#
# The executable is upstream's makefile target, built unmodified from the versioned
# release archive named in issue #113, with gfortran, as the release's README says.
# The parameter and class files come from upstream's demo2 for the same release and are
# copied unchanged; ht_adapter.py applies a named list of edits for each case. Both
# archives are checked against the checksum recorded when the evaluation was made. The
# network is available at build time only. HYPE is LGPL-3.0; both licence texts travel
# with the image.
#
# Both stages use the same python:3.11-slim index (Debian 13), pinned by digest, so the
# runtime has the libgfortran the executable was linked against.

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534 AS build
ARG SRC_URL=https://sourceforge.net/projects/hype/files/release_hype_5_36_0/hype_5_36_0_src.tgz/download
ARG SRC_SHA256=d90bc3670f341bb80105de94a23531775497ab472dbbc3560cd25b2afb03f43a
ARG DEMO_URL=https://sourceforge.net/projects/hype/files/release_hype_5_36_0/demo2.zip/download
ARG DEMO_SHA256=117fdd721e3ea8a29ed57b4c6fcde5db4c2b6cfdbea6864096da14281af59024
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl unzip make gfortran \
 && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /src && cd /src \
 && curl -fsSL -o hype.tgz "${SRC_URL}" \
 && echo "${SRC_SHA256}  hype.tgz" | sha256sum -c - \
 && curl -fsSL -o demo2.zip "${DEMO_URL}" \
 && echo "${DEMO_SHA256}  demo2.zip" | sha256sum -c - \
 && tar xzf hype.tgz \
 && unzip -q demo2.zip
RUN cd /src/hype_5_36_0_src \
 && make hype \
 && strip hype

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
RUN apt-get update \
 && apt-get install -y --no-install-recommends libgfortran5 \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/hype_5_36_0_src/hype /model/hype
COPY --from=build /src/hype_5_36_0_src/lgpl-3_0.txt /model/lgpl-3_0.txt
COPY --from=build /src/hype_5_36_0_src/gpl-3_0.txt /model/gpl-3_0.txt
COPY --from=build /src/demo2/par.txt /model/template/par.txt
COPY --from=build /src/demo2/GeoClass.txt /model/template/GeoClass.txt
COPY --from=build /src/demo2/GeoData.txt /model/template/GeoData.txt
COPY --from=build /src/demo2/info.txt /model/template/info.txt
WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
