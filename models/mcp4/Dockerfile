# MCPModel4 (YuanHWang/5MCPmodelnet_HydroTuring) for HydroTuring.
#
# The network is available at build time only. The pinned source and the
# checkpoint it ships are baked in here; the container itself runs with no
# network, and the adapter never downloads anything. Code and weights are
# Apache-2.0.

# python:3.10-slim, the version the repository was tested on, at a fixed digest.
FROM python:3.10-slim@sha256:6ff506466f8b1e981719b468cc0023b8ef000b983cdec0eafe119d01f379e31d AS source
ARG REVISION=5944b5f23adb7265a84fd10d62a8bd5501e25208
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src \
 && cd /src \
 && git remote add origin https://github.com/YuanHWang/5MCPmodelnet_HydroTuring.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${REVISION}"

FROM python:3.10-slim@sha256:6ff506466f8b1e981719b468cc0023b8ef000b983cdec0eafe119d01f379e31d
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
# The versions the repository's requirements.txt names for what the model code
# imports. torch comes from the CPU index: the evaluation script is CPU by
# default and the store loops gain nothing from a GPU. The base image's pip is too
# old for that index's package names, and torch's own dependencies are pinned too.
RUN pip install pip==25.3 \
 && pip install --no-deps torch==2.9.0+cpu --index-url https://download.pytorch.org/whl/cpu \
 && pip install numpy==2.2.2 typing-extensions==4.16.0 filelock==4.1.0 fsspec==2026.9.0 \
      jinja2==3.1.6 markupsafe==3.0.4 mpmath==1.3.0 networkx==3.4.2 sympy==1.14.0

# Only what inference reads. Leaf River data and training scripts stay out.
COPY --from=source /src/MCPBRNN_lib_tools /model/upstream/MCPBRNN_lib_tools
COPY --from=source /src/MCPmodel4/best_model_epoch2707.pt /model/upstream/MCPmodel4/best_model_epoch2707.pt
COPY --from=source /src/LICENSE /model/upstream/LICENSE

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py
ENV MCP_UPSTREAM=/model/upstream OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
