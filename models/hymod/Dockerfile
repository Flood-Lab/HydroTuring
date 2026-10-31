# HYMOD, from OuyangWenyu/hydromodel, for HydroTuring.
#
# The network is available at build time only. The pinned source is baked in
# here; the container runs with no network and the adapter downloads nothing.
# hydromodel is GPL-3.0; the licence travels with the image.

# python:3.11-slim (Python 3.11.16, Debian trixie), the digest models/cwatm pins.
FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534 AS source
# hydromodel v0.4.0, the commit named in the submission (issue #133).
ARG REVISION=89d7a8ed1d72ce4fffbbd9897490b089382ecbac
RUN apt-get update \
 && apt-get install -y --no-install-recommends git ca-certificates \
 && rm -rf /var/lib/apt/lists/*
RUN git init -q /src \
 && cd /src \
 && git remote add origin https://github.com/OuyangWenyu/hydromodel.git \
 && git fetch -q --depth 1 origin "${REVISION}" \
 && git checkout -q FETCH_HEAD \
 && test -f LICENSE
# Only the model and the two modules it imports, byte for byte, under empty
# package files. The package's own __init__ files import hydroutils for a cache
# path, and through it async-retriever, cartopy and the rest of a data stack
# that forward simulation never touches. Upstream's uv.lock cannot be installed
# either: it resolves hydrodataset from a sibling checkout outside the repository.
RUN mkdir -p /vendor/hydromodel/models \
 && touch /vendor/hydromodel/__init__.py /vendor/hydromodel/models/__init__.py \
 && cp /src/hydromodel/models/hymod.py /src/hydromodel/models/model_config.py /src/hydromodel/models/param_utils.py /vendor/hydromodel/models/

FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
# The versions upstream's uv.lock records at this commit.
RUN pip install "numpy==2.2.6" "scipy==1.16.0" "numba==0.61.2" "llvmlite==0.44.0" "pyyaml==6.0.2"
COPY --from=source /vendor/hydromodel /opt/hydromodel/hydromodel
COPY --from=source /src/LICENSE /opt/hydromodel/LICENSE
ENV PYTHONPATH=/opt/hydromodel PYTHONDONTWRITEBYTECODE=1 NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
