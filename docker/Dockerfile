FROM python:3.10-slim-bookworm

ENV LANG=C.UTF-8 \
  LC_ALL=C.UTF-8 \
  UV_COMPILE_BYTECODE=1 \
  UV_LINK_MODE=copy

# Build toolchain plus htslib/zlib, which the _cppext extension links against.
# HTSLIB_ROOT tells CMakeLists where to find them outside a conda environment.
RUN apt-get update && \
  apt-get install -y --no-install-recommends \
  build-essential \
  cmake \
  ninja-build \
  libhts-dev \
  zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

ENV HTSLIB_ROOT=/usr

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

WORKDIR /src
COPY . /src

# Allow installing dev dependencies to run tests.
# The previous form lost the ARG interpolation ("if [ == 'true' ]") and used
# --no-dev, which Poetry 2.x removed, so this stage had been failing.
ARG INSTALL_DEV=false
RUN if [ "$INSTALL_DEV" = "true" ]; then \
      uv sync --frozen --group dev ; \
    else \
      uv sync --frozen --no-dev ; \
    fi

ENV PATH="/src/.venv/bin:${PATH}"

RUN mkdir -p /workspace
WORKDIR /workspace
