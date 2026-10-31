# VIC 5.1.0, classic driver, for HydroTuring.
#
# Compiled from the pinned tag with upstream's own Makefile, passing LOG_LVL, which the
# Makefile documents, to stop it logging warnings, and the compiler flags below. The sample
# parameters come from the commit the 5.1.0 tag pins as its samples/data submodule.
# The network is available at build time only. VIC is MIT, the sample data GPL-2.0;
# both licences are copied into the image.

FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171 AS build
ARG VIC_COMMIT=14a371a8dbee7ea834152f8ac65e2d6d54a9c30a
ARG SAMPLE_COMMIT=3f2dcb445e3e985456589a29711a2fee59f8faaf
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates git make gcc libc6-dev \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src/VIC && cd /src/VIC \
 && git remote add origin https://github.com/UW-Hydro/VIC.git \
 && git fetch -q --depth 1 origin "refs/tags/5.1.0:refs/tags/5.1.0" \
 && git checkout -q 5.1.0 \
 && test "$(git rev-parse HEAD)" = "${VIC_COMMIT}" \
 && test "$(git ls-tree HEAD samples/data | cut -d' ' -f3 | cut -f1)" = "${SAMPLE_COMMIT}"
RUN git init -q /src/sample && cd /src/sample \
 && git remote add origin https://github.com/UW-Hydro/VIC_sample_data.git \
 && git fetch -q --depth 1 origin "${SAMPLE_COMMIT}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${SAMPLE_COMMIT}"
# -fcommon because VIC's headers define globals (LOG_DEST, funcd) that GCC 10 and later
# refuse to link by default. -O3 is the Makefile's own "normal optimized" line, which it
# ships commented out in favour of -g. Upstream's source is as tagged.
RUN cd /src/VIC/vic/drivers/classic && make LOG_LVL=30 CC="gcc -fcommon -O3" all

FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171
ENV PYTHONUNBUFFERED=1
RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/VIC/vic/drivers/classic/vic_classic.exe /opt/vic/vic_classic.exe
COPY --from=build /src/VIC/LICENSE.txt /model/LICENSE
COPY --from=build /src/sample/LICENSE /model/LICENSE.sample_data
COPY --from=build /src/sample/classic/Stehekin/parameters /model/stehekin
RUN chmod -R a+rX /model /opt/vic
WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
