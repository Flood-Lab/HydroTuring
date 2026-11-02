# PyTOPKAPI 0.4.0 for HydroTuring.
#
# The network is available at build time only. The pinned source is baked in
# here; the container runs with no network and the adapter downloads nothing.

# Pinned to linux/amd64, as models/lisflood/Dockerfile is:
# h5py 3.11.0 has no linux/arm64 wheel for Python 3.8, and building it from source
# needs HDF5 and pkg-config, which this image does not carry.
# Native on an x86-64 host, emulated on an arm64 one.
FROM --platform=linux/amd64 python:3.8.20-slim@sha256:1d52838af602b4b5a831beb13a0e4d073280665ea7be7f69ce2382f29c5a613f AS source
# Tag v0.4.0.
ARG REVISION=f7b953eebcb201403f80970846b429bb1081739c
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src \
 && cd /src \
 && git remote add origin https://github.com/sahg/PyTOPKAPI.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD

# Python 3.8 and networkx 1.11 are not a preference. pretreatment.sort_cell
# takes np.array() of nx.topological_sort(), which is a list only before
# networkx 2.0, and networkx 1.11 itself does not import on Python 3.9+.
# PyTOPKAPI's raster pre-processing (parameter_utils.create_file) is the only
# code that needs GDAL and it is not on this path, so GDAL is not installed.
# matplotlib is here because `import pytopkapi` imports its plotting tools.
FROM --platform=linux/amd64 python:3.8.20-slim@sha256:1d52838af602b4b5a831beb13a0e4d073280665ea7be7f69ce2382f29c5a613f
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
RUN pip install numpy==1.24.4 scipy==1.10.1 h5py==3.11.0 networkx==1.11 tqdm==4.66.5 matplotlib==3.7.5 pandas==2.0.3

COPY --from=source /src/pytopkapi /model/pytopkapi
COPY --from=source /src/LICENSE /model/LICENSE

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py
ENV MPLBACKEND=Agg MPLCONFIGDIR=/tmp OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONWARNINGS=ignore

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
