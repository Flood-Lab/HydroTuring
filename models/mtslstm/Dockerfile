# MTS-LSTM (neuralhydrology at 9626578, Gauch et al. 2021) for HydroTuring.
#
# The network is available at build time only. The pinned source and the published
# weights are baked in here; the container itself runs with no network, and the
# adapter never downloads anything. Nothing is trained.

# python:3.12-slim at a fixed digest.
FROM python:3.12.15-slim@sha256:a0f57064cbdcd6285fa45f51237ab0b35ea8ebaed06a1f5973e89e19dfb86704 AS source
ARG REVISION=9626578d6013aa08864caf29d9520321cdcafc24
ARG RECORD=https://zenodo.org/api/records/18953458/files
ARG RUN=mtslstm/mtslstm_seed333_1208_060947
ARG SCALER=smtslstm/smtslstm_noregularization_seed555_1408_051556/train_data/train_data_scaler.p
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates curl \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src \
 && cd /src \
 && git remote add origin https://github.com/neuralhydrology/neuralhydrology.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${REVISION}"
# The published archives are 22 and 6 GB. GNU tar stops at the last member it was
# named, and these sit near the front of their streams, so a few minutes of download
# is enough. The whole archives are not checked: the members are, by checksum.
#
# The MTS-LSTM archive's train_data directories are empty, so this is not the run's
# own scaler: it comes from the sMTS-LSTM archive, whose runs read the same data over
# the same period but not with the same hourly sequence length (72 here, 336 there).
# What justifies it is a measurement, in README.md: with it, the checkpoint below
# reproduces its own published simulation to single precision.
RUN mkdir -p /dl /run_dir && cd /dl \
 && curl -sSL --retry 5 "${RECORD}/models-mtslstm.tar.gz/content" \
    | tar -xz --occurrence=1 "${RUN}/config.yml" "${RUN}/model_epoch030.pt" \
 ; cd /dl \
 && curl -sSL --retry 5 "${RECORD}/models-smtslstm.tar.gz/content" \
    | tar -xz --occurrence=1 "${SCALER}" \
 ; cd /dl \
 && echo "c16f5454769101186c19ce3876a090c9cd587189cad282a3f526a7416aa72d8f  ${RUN}/config.yml" | sha256sum -c - \
 && echo "387f8b67b381571cd02a27f9206317e1e352ddaca7fe2c212735ff7d3e5cddfa  ${RUN}/model_epoch030.pt" | sha256sum -c - \
 && echo "256c5285f44e0dd91f67759aac9de44c2acf9fd494223d054ac6f5ca0c390da6  ${SCALER}" | sha256sum -c - \
 && mv "${RUN}/config.yml" "${RUN}/model_epoch030.pt" "${SCALER}" /run_dir/ \
 && chmod 644 /run_dir/* \
 && rm -rf /dl

# The scaler's pickle holds 2020 xarray objects that no longer open on a current stack,
# so it is read here under the versions in the paper's environment.yml. numpy, pandas
# and xarray are its pins; python 3.8 is the oldest slim image with those wheels.
FROM python:3.8.20-slim@sha256:1d52838af602b4b5a831beb13a0e4d073280665ea7be7f69ce2382f29c5a613f AS scaler
ENV PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
RUN pip install numpy==1.18.5 pandas==1.0.5 xarray==0.16.0 setuptools==65.5.0
COPY --from=source /run_dir/train_data_scaler.p /work/train_data_scaler.p
COPY export_scaler.py /work/export_scaler.py
RUN cd /work && python3 export_scaler.py train_data_scaler.p scaler.json

FROM python:3.12.15-slim@sha256:a0f57064cbdcd6285fa45f51237ab0b35ea8ebaed06a1f5973e89e19dfb86704
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
# CPU torch from the PyTorch index. The rest is what the package imports on the way
# to the network; training-only dependencies such as tensorboard are left out.
RUN pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cpu \
 && pip install numpy==1.26.4 pandas==2.2.3 xarray==2026.9.0 numba==0.68.1 scipy==1.17.1 \
      matplotlib==3.11.2 tqdm==4.70.1 ruamel.yaml==0.19.1

# Only what inference reads.
COPY --from=source /src/neuralhydrology /model/upstream/neuralhydrology
COPY --from=source /src/LICENSE /model/upstream/LICENSE
COPY --from=source /run_dir/config.yml /model/run/config.yml
COPY --from=source /run_dir/model_epoch030.pt /model/run/model_epoch030.pt
COPY --from=scaler /work/scaler.json /model/run/scaler.json

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py
ENV PYTHONPATH=/model/upstream MTSLSTM_RUN_DIR=/model/run MTSLSTM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=1 MPLCONFIGDIR=/tmp

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
