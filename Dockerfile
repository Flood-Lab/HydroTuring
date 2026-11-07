# Builds nothing of the benchmark itself: this image holds the harness and a
# Docker client, and running it does the checking (scripts/check_all.py).
#
# The models cannot share one image. They need R, Julia, and Python from 3.7 to
# 3.12, and two are pinned to linux/amd64, so each keeps its own Dockerfile,
# and a Dockerfile cannot build another Dockerfile. The container therefore
# drives the host's daemon through its socket:
#
#   docker build -t hydroturing-check .
#   mkdir -p "$HOME/.hydroturing-check"
#   docker run --rm \
#     -v /var/run/docker.sock:/var/run/docker.sock \
#     -v "$HOME/.hydroturing-check:$HOME/.hydroturing-check" \
#     -e TMPDIR="$HOME/.hydroturing-check" \
#     hydroturing-check
#
# The harness mounts each case's /io from its temp directory, and the daemon
# resolves that path on the host, so TMPDIR must be the same path on both
# sides. On Docker Desktop it must also be a shared directory, such as one
# under /Users. Model images land in the host's image store under the tags the
# harness uses, so a later `ht run` reuses them.
#
# The tests run in the container, not during `docker build`, so a failing test
# does not stop the models from being checked. Arguments after the image name
# go to check_all.py: --skip-tests, --skip-models, --models NAME..., --jobs N.

# Same index the models pin, and the Python CI tests on.
FROM python:3.11-slim@sha256:9534e5a8e315485d4061ed659af0fd78a284c015f9b73661b41d6bab25604534

COPY --from=docker:28.5.2-cli@sha256:625d9431a9f54c5a2bc90f24f0e1c3d55b1349fd857dd85035f98c2c9acbdd4d \
     /usr/local/bin/docker /usr/local/bin/docker
COPY --from=docker:28.5.2-cli@sha256:625d9431a9f54c5a2bc90f24f0e1c3d55b1349fd857dd85035f98c2c9acbdd4d \
     /usr/local/libexec/docker/cli-plugins/docker-buildx /usr/local/libexec/docker/cli-plugins/docker-buildx

WORKDIR /repo

# Dependencies first, so a source change does not reinstall them.
COPY pyproject.toml README.md ./
RUN mkdir -p src && pip install --no-cache-dir -e '.[dev]'

COPY . .
RUN pip install --no-cache-dir --no-deps -e .

ENV PYTHONUNBUFFERED=1 PYTHONUTF8=1

ENTRYPOINT ["python3", "scripts/check_all.py"]
