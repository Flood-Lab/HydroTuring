# EA-LSTM (kratzert/ealstm_regional_modeling) for HydroTuring.
#
# The network is available at build time only. The pinned source and the published
# weights are baked in here; the container itself runs with no network, and the
# adapter never downloads anything. Nothing is trained.

# python:3.12-slim at a fixed digest.
FROM python:3.12.15-slim@sha256:a0f57064cbdcd6285fa45f51237ab0b35ea8ebaed06a1f5973e89e19dfb86704 AS source
ARG REVISION=d11815824793737f28f125d55477e42f8900a78e
ARG ARCHIVE_MD5=224985a9c4d0932e607f4c9fe8d61ee8
ARG RUN=run_1906_1004_seed111
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates curl unzip \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src \
 && cd /src \
 && git remote add origin https://github.com/kratzert/ealstm_regional_modeling.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test "$(git rev-parse HEAD)" = "${REVISION}"
# The published archive, checked against the checksum HydroShare states for it. It
# holds 48 runs and their simulations; only one run's weights, configuration and
# attribute table are kept.
RUN mkdir -p /dl /run_dir \
 && curl -sSfL --retry 5 -o /dl/runs.zip "https://www.hydroshare.org/resource/83ea5312635e44dc824eeb99eda12f06/data/contents/Kratzert_et_al_2019_HESSD_LSTM_4_REGIONAL_MODELING.zip" \
 && echo "${ARCHIVE_MD5}  /dl/runs.zip" | md5sum -c - \
 && unzip -q -j /dl/runs.zip \
      "Kratzert_et_al_2019_HESSD_LSTM_4_REGIONAL_MODELING/runs/${RUN}/cfg.json" \
      "Kratzert_et_al_2019_HESSD_LSTM_4_REGIONAL_MODELING/runs/${RUN}/model_epoch30.pt" \
      "Kratzert_et_al_2019_HESSD_LSTM_4_REGIONAL_MODELING/runs/${RUN}/attributes.db" \
      -d /run_dir \
 && rm -rf /dl

FROM python:3.12.15-slim@sha256:a0f57064cbdcd6285fa45f51237ab0b35ea8ebaed06a1f5973e89e19dfb86704
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
# CPU torch from the PyTorch index. The rest is what upstream's main.py imports.
RUN pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cpu \
 && pip install numpy==1.26.4 pandas==2.2.3 h5py==3.16.0 numba==0.68.1 tqdm==4.70.1

# Only what inference reads.
COPY --from=source /src/main.py /model/upstream/main.py
COPY --from=source /src/papercode /model/upstream/papercode
COPY --from=source /src/data/basin_list.txt /model/upstream/data/basin_list.txt
COPY --from=source /src/LICENSE /model/upstream/LICENSE
COPY --from=source /run_dir /model/run

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py
ENV EALSTM_UPSTREAM=/model/upstream EALSTM_RUN_DIR=/model/run EALSTM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=1

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
