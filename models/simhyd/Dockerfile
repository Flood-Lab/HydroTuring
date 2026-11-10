# SIMHYD as MARRMoT v2.1.2 implements it (m_18_simhyd_7p_3s), for HydroTuring.
#
# The model code is MARRMoT's own MATLAB, unmodified, run under GNU Octave 6.4.0
# with the optim package, the combination its README says it was tested with. The
# network is available at build time only; the container runs without one.
# MARRMoT is GPL-3.0; the licence text is copied into the image beside the adapter.
#
# Ubuntu 22.04 is the release whose archive carries Octave 6.4.0. Its index is pinned
# by digest and the two packages that decide the numerics by version, so a build
# fails rather than quietly changes interpreter or solver.

FROM ubuntu:22.04@sha256:076e7bf82319619f28dd0aa9c02ab34b055ad37100395b4f49a2ffba04fb5e11 AS source
# wknoben/MARRMoT, tag v2.1.2, the release named in issue #164.
ARG MARRMOT_TAG=v2.1.2
ARG MARRMOT_COMMIT=a95b925669385e6c747637dec45c848601177ea4
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates git \
 && rm -rf /var/lib/apt/lists/*
# The build is refused unless the tag still points at the pinned commit.
RUN git init -q /src && cd /src \
 && git remote add origin https://github.com/wknoben/MARRMoT.git \
 && git fetch -q --depth 1 origin "refs/tags/${MARRMOT_TAG}:refs/tags/${MARRMOT_TAG}" \
 && git checkout -q "${MARRMOT_TAG}" \
 && test "$(git rev-parse HEAD)" = "${MARRMOT_COMMIT}"

FROM ubuntu:22.04@sha256:076e7bf82319619f28dd0aa9c02ab34b055ad37100395b4f49a2ffba04fb5e11
ENV PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MARRMOT_ROOT=/opt/MARRMoT
RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 octave=6.4.0-2 octave-optim=1.6.1-2build1 \
 && rm -rf /var/lib/apt/lists/*
COPY --from=source /src/MARRMoT/Functions /opt/MARRMoT/Functions
COPY --from=source /src/MARRMoT/Models /opt/MARRMoT/Models
COPY --from=source /src/LICENSE /model/LICENSE
WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py
COPY simhyd_run.m /model/simhyd_run.m
RUN octave-cli --norc --no-history --eval "assert(strcmp(version(),'6.4.0')); pkg load optim; disp(version())"

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
