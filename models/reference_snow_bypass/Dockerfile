# Reference models run via the subprocess runner in CI, so this image is not
# on the critical path. It is kept as the worked example a submitter copies.
FROM python:3.11-slim
WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py
ENV PYTHONUNBUFFERED=1
ENTRYPOINT []
