# Tank Hydrologic Model (nzahasan/tank-model) for HydroTuring.
#
# The network is available at build time only. The pinned source is baked in
# here; the container runs with no network and the adapter downloads nothing.

FROM python:3.12.15-slim@sha256:a0f57064cbdcd6285fa45f51237ab0b35ea8ebaed06a1f5973e89e19dfb86704 AS source
# Tag 3.0.0b3. tank_core is vendored from this commit rather than built, so the
# image needs neither hatchling nor upstream's plotting and CLI dependencies.
ARG REVISION=aea105ca32d0957fee76fbcfe3a5b8993938c9a1
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src \
 && cd /src \
 && git remote add origin https://github.com/nzahasan/tank-model.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD

FROM python:3.12.15-slim@sha256:a0f57064cbdcd6285fa45f51237ab0b35ea8ebaed06a1f5973e89e19dfb86704
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
# The versions upstream's uv.lock resolves at this commit, for what tank_core imports.
RUN pip install numpy==2.5.1 pandas==3.0.3 scipy==1.18.0

COPY --from=source /src/tank_core /model/tank_core
COPY --from=source /src/LICENSE /model/LICENSE

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
