FROM python:3.12-slim

ARG MF6_VERSION=6.7.0
ARG MF6_SHA256=f83e743bac1b3009836a44e6e4aa4b4c522ac2e2c7d291e280db47f68353ebdc

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates curl unzip \
    && rm -rf /var/lib/apt/lists/* \
    && curl -fsSL -o /tmp/mf6.zip \
       "https://github.com/MODFLOW-ORG/modflow6/releases/download/${MF6_VERSION}/mf6.7.0_linux.zip" \
    && echo "${MF6_SHA256}  /tmp/mf6.zip" | sha256sum -c - \
    && mkdir -p /opt/modflow6 \
    && unzip -q /tmp/mf6.zip -d /opt/modflow6 \
    && find /opt/modflow6 -type f -name mf6 -exec chmod +x {} \; \
    && rm /tmp/mf6.zip

RUN pip install --no-cache-dir flopy==3.11.0 numpy

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py
ENV MODFLOW6_ROOT=/opt/modflow6
# model.yaml supplies the entrypoint so the image follows the repository's
# standard model packaging convention.
ENTRYPOINT []