FROM python:3.11-slim

# No network at run time. Bake in weights, parameter files and lookup tables.
# RUN pip install --no-cache-dir your-model==1.2.3

WORKDIR /model
COPY ht_adapter.py /model/ht_adapter.py
ENV PYTHONUNBUFFERED=1

# Leave the entrypoint empty. model.yaml supplies the argv.
ENTRYPOINT []
