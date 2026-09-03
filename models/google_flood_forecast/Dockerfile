# Google Flood Hub model (OpenHydroNet) for HydroTuring.
#
# The network is available at build time only. The pinned source and the
# published weights are baked in here; the container itself runs with no
# network, and the adapter never downloads anything.

FROM python:3.12-slim AS source
# google-research/flood-forecasting at the revision named in the submission.
ARG REVISION=828dfc5a66f4e6f2c86b09d500e4547c4d92ed25
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src \
 && cd /src \
 && git remote add origin https://github.com/google-research/flood-forecasting.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD

FROM python:3.12-slim
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
# CPU-only torch from the PyTorch index, the rest from PyPI. These are the
# imports the model's inference path needs; training-only dependencies
# (gcsfs, tensorboard, bokeh, ...) are left out.
RUN pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cpu \
 && pip install numpy pandas xarray h5netcdf h5py dask pydantic "ruamel.yaml" more-itertools

COPY --from=source /src /opt/flood-forecasting
RUN pip install --no-deps /opt/flood-forecasting

# Only the run directory the adapter needs: configuration, feature scaler
# and the epoch-55 weights. The optimizer state and the second (NSE-filtered)
# model are not copied.
COPY --from=source \
     /src/pretrained-models/google-floodhub-settings-55-epochs/config.yml \
     /src/pretrained-models/google-floodhub-settings-55-epochs/scaler.nc \
     /src/pretrained-models/google-floodhub-settings-55-epochs/model_epoch055.pt \
     /model/run/

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py
ENV GFF_RUN_DIR=/model/run GFF_THREADS=8 OMP_NUM_THREADS=8

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
