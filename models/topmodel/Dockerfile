# TOPMODEL (ICHydro/topmodel, R package, Buytaert) for HydroTuring.
#
# The network is available at build time only. The package is built from a pinned commit and
# jsonlite from CRAN's archive at a pinned checksum; the container runs with no network.
#
# rocker/r-ver:4.4.2, by digest and without its tag: Apptainer's image references do not
# accept a tag and a digest together, and the digest is what pins it.

FROM rocker/r-ver@sha256:df26749182af64d5263bf64149d51a427b476ed28c4e046997143be3f97fdd7c AS source
# ICHydro/topmodel at "update to 0.7.5" (2022-11-15). The repository has no tags.
ARG REVISION=f60752ced41ffac47778ffb6616ae01cae56456e
ARG JSONLITE_SHA256=89f130e0e1163328c01decd54e7712b5ebf3d0a667da0052833722cb9a6e90b0
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates curl \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src \
 && cd /src \
 && git remote add origin https://github.com/ICHydro/topmodel.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${REVISION}" \
 && rm -rf /src/.git
RUN mkdir -p /dl \
 && curl -fsSL -o /dl/jsonlite.tar.gz https://cran.r-project.org/src/contrib/Archive/jsonlite/jsonlite_1.8.9.tar.gz \
 && echo "${JSONLITE_SHA256}  /dl/jsonlite.tar.gz" | sha256sum -c -

FROM rocker/r-ver@sha256:df26749182af64d5263bf64149d51a427b476ed28c4e046997143be3f97fdd7c
COPY --from=source /src /opt/topmodel-src
COPY --from=source /dl/jsonlite.tar.gz /opt/jsonlite.tar.gz
# GPL (>= 2): the source the image was built from stays in it, under /opt/topmodel-src.
RUN R CMD INSTALL --no-docs /opt/jsonlite.tar.gz \
 && R CMD INSTALL --no-docs /opt/topmodel-src/topmodel \
 && Rscript --vanilla -e 'stopifnot(packageVersion("topmodel") == "0.7.5")'

# The driver is the package's own core files, unchanged, plus ht_driver.c, which only copies
# out the stores c_topmodel() frees. topmodel_init.c stays out: it registers the package's
# symbols and turns dynamic lookup off.
WORKDIR /opt/topmodel
COPY ht_driver.c /opt/topmodel/ht_driver.c
COPY check_driver.R /opt/topmodel/check_driver.R
RUN cp /opt/topmodel-src/topmodel/src/topmodel.h /opt/topmodel/ \
 && for f in c_topmodel core_topmodel get_f misc output param_init startroutines NS; do cp /opt/topmodel-src/topmodel/src/$f.c /opt/topmodel/; done \
 && R CMD SHLIB -o ht_topmodel.so ht_driver.c c_topmodel.c core_topmodel.c get_f.c misc.c output.c param_init.c startroutines.c NS.c \
 && Rscript --vanilla /opt/topmodel/check_driver.R

WORKDIR /model
COPY ht_adapter.R /model/ht_adapter.R
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
