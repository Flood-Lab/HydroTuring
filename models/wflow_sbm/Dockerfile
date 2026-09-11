# wflow_sbm (Wflow.jl, Deltares) for HydroTuring.
#
# Wflow.jl is fetched at a pinned commit at build time and installed with exactly the
# dependency versions of the Manifest.toml committed next to this file (Wflow's own
# Manifest at that commit plus JSON, for the adapter). The adapter and the Wflow code it
# calls are compiled here, into a depot that is read-only at run time; the container runs
# with no network and never installs or compiles anything.

FROM julia:1.12.6 AS source
# Deltares/Wflow.jl tag v1.0.4.
ARG WFLOW_REVISION=82df72031511339d50fd9142fa159d0ec13e73c5
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src \
 && cd /src \
 && git remote add origin https://github.com/Deltares/Wflow.jl.git \
 && git fetch -q --depth 1 origin "${WFLOW_REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${WFLOW_REVISION}"

FROM julia:1.12.6
COPY --from=source /src/Wflow /opt/wflow/Wflow
COPY --from=source /src/LICENSE /opt/wflow/LICENSE

# A generic CPU target, so the compiled caches load on any CPU of this architecture rather
# than only on the one that built the image. The same value is kept at run time.
ENV JULIA_DEPOT_PATH=/opt/julia-depot \
    JULIA_CPU_TARGET=generic \
    JULIA_NUM_PRECOMPILE_TASKS=2 \
    JULIA_PKG_PRECOMPILE_AUTO=0

WORKDIR /opt/wflow_sbm
COPY Project.toml Manifest.toml ./
COPY src ./src
COPY ht_adapter.jl ./
RUN julia --project=/opt/wflow_sbm -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); using WflowSbmAdapter' \
 && rm -rf /opt/julia-depot/registries /opt/julia-depot/logs \
 && chmod -R a+rX /opt/julia-depot /opt/wflow /opt/wflow_sbm

# At run time the image's depot is read-only and the container runs as an arbitrary user:
# anything Julia wants to write goes to the first depot, on the /tmp tmpfs.
ENV JULIA_DEPOT_PATH=/tmp/julia-depot:/opt/julia-depot

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
