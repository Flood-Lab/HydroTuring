# δHBV 2.0 (mhpi/dhbv2) for HydroTuring.
#
# The network is available at build time only. The pinned source and the
# published weights are baked in here; the container itself runs with no
# network, and the adapter never downloads anything.
#
# The code and the weights are under the PSU Non-Commercial Software
# License. This image is for the non-commercial evaluation the license
# permits; do not push it to a public registry.

FROM python:3.11-slim AS source
# mhpi/dhbv2 at the revision named in the submission (tag v0.5.4).
ARG REVISION=6d4bb686ab00ac81e3ea2150b0a18eb29ed99938
ARG WEIGHTS_URL=https://mhpi-spatial.s3.us-east-2.amazonaws.com/mhpi-release/models/owp/dhbv_2.zip
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates curl unzip \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src \
 && cd /src \
 && git remote add origin https://github.com/mhpi/dhbv2.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD
# The published daily weights: config.yaml, normalization_statistics.json
# and hbv_2_ep100.pt, from the MHPI release bucket (public, unsigned).
RUN curl -sSL -o /tmp/dhbv_2.zip "${WEIGHTS_URL}" \
 && unzip -q /tmp/dhbv_2.zip -d /weights \
 && rm /tmp/dhbv_2.zip \
 && ls /weights/dhbv_2

FROM python:3.11-slim
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
# CPU-only torch from the PyTorch index; the framework (dmg) and the physics
# (hydrodl2) pinned to the versions dhbv2 v0.5.4 resolves against today.
RUN pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cpu \
 && pip install "numpy<2" pandas pyyaml pydantic bmipy dmg==1.4.3 hydrodl2==1.3.5

# The module itself from the release wheel of the same tag (v0.5.4 is the
# commit fetched above; the package takes its version from that tag, which
# a one-commit fetch does not carry, so the wheel is the reproducible way in).
RUN pip install --no-deps dhbv2==0.5.4

COPY --from=source /weights/dhbv_2 /model/dhbv_2
COPY --from=source /src/LICENSE /model/LICENSE

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py
ENV DHBV_MODEL_DIR=/model/dhbv_2 DHBV_THREADS=2 OMP_NUM_THREADS=2

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
