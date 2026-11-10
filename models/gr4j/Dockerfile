# GR4J (airGR, INRAE) for HydroTuring.
#
# The network is available at build time only. airGR is built from the CRAN mirror at a
# pinned commit and jsonlite from CRAN's archive at a pinned checksum; the container runs
# with no network and installs nothing.
#
# rocker/r-ver:4.4.2, by digest and without its tag: Apptainer's image references do not
# accept a tag and a digest together, and the digest is what pins it.

FROM rocker/r-ver@sha256:df26749182af64d5263bf64149d51a427b476ed28c4e046997143be3f97fdd7c AS source
# cran/airGR at the commit of airGR 1.7.9 (CRAN, 2026-07-08), as issue #96 names it.
ARG REVISION=4bb401c6529e20116db357369c9f2bc34f959992
ARG JSONLITE_SHA256=89f130e0e1163328c01decd54e7712b5ebf3d0a667da0052833722cb9a6e90b0
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates curl \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src \
 && cd /src \
 && git remote add origin https://github.com/cran/airGR.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${REVISION}" \
 && rm -rf /src/.git
RUN mkdir -p /dl \
 && curl -fsSL -o /dl/jsonlite.tar.gz https://cran.r-project.org/src/contrib/Archive/jsonlite/jsonlite_1.8.9.tar.gz \
 && echo "${JSONLITE_SHA256}  /dl/jsonlite.tar.gz" | sha256sum -c -

FROM rocker/r-ver@sha256:df26749182af64d5263bf64149d51a427b476ed28c4e046997143be3f97fdd7c
COPY --from=source /src /opt/airGR
COPY --from=source /dl/jsonlite.tar.gz /opt/jsonlite.tar.gz
# GPL-2: the source the image was built from stays in it, under /opt/airGR.
RUN R CMD INSTALL --no-docs /opt/jsonlite.tar.gz \
 && R CMD INSTALL --no-docs /opt/airGR \
 && Rscript --vanilla -e 'stopifnot(packageVersion("airGR") == "1.7.9")'

WORKDIR /model
COPY ht_adapter.R /model/ht_adapter.R
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
