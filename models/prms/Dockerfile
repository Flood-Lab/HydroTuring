# PRMS 6.0.0 for HydroTuring: the official USGS Linux executable, unmodified.
#
# The zip is pinned by SHA-256 because USGS publishes it at a bare URL with no commit
# to name. The build checks it, runs the bundled Sagehen example and refuses to finish
# unless the water budget totals match upstream's own recorded output. Nothing is
# fetched at run time. PRMS is CC0 per the release README, which the image carries.
#
# Both stages use one bookworm-slim index by digest, so the runtime has the glibc and
# libgfortran the executable was checked against.

FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171 AS build
ARG PRMS_URL=https://water.usgs.gov/water-resources/software/PRMS/prms_6.0.0_linux.zip
ARG PRMS_SHA256=450c5b65077d343e74cbe31a78a51ef8685a20f846394120f6d882af44d00423
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl unzip libgfortran5 \
 && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /src && cd /src \
 && curl -sSfL -o prms.zip "${PRMS_URL}" \
 && test "$(sha256sum prms.zip | cut -d' ' -f1)" = "${PRMS_SHA256}" \
 && unzip -q prms.zip && chmod +x prms_6.0.0_linux/bin/prms
RUN cd /src/prms_6.0.0_linux/projects/sagehen && mkdir -p output \
 && ../../bin/prms ./control/climate_hru.control > /dev/null \
 && test "$(grep 'Total for run' output/sagehen_cbh.out)" = "$(grep 'Total for run' output-test/climate_hru/sagehen_cbh.out)"

FROM debian:bookworm-slim@sha256:88200866dfff7ea7f5cbcb6ec7c8a701889efe6fe859fe64d6990e4b07ea4171
ENV PYTHONUNBUFFERED=1
RUN apt-get update \
 && apt-get install -y --no-install-recommends libgfortran5 python3 \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/prms_6.0.0_linux/bin/prms /opt/prms/bin/prms
COPY --from=build /src/prms_6.0.0_linux/Readme_PRMS_6.0.0.pdf /model/Readme_PRMS_6.0.0.pdf
RUN chmod -R a+rX /model /opt/prms
WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
