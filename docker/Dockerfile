# Integrating poetry with docker <https://stackoverflow.com/questions/72465421/how-to-use-poetry-with-docker>
FROM --platform=linux/amd64 python:3.10-slim AS builder

ENV POETRY_VERSOIN=1.8.3
ENV POETRY_HOME=/opt/poetry
ENV POETRY_VENV=/opt/poetry-venv

ENV POETRY_CACHE_DIR=/opt/.cache

FROM builder as poetry-base

RUN python3 -m venv $POETRY_VENV \
    && $POETRY_VENV/bin/pip install -U pip setuptools \
    && $POETRY_VENV/bin/pip install poetry==${POETRY_VERSOIN}

FROM builder as app

RUN apt-get update \
    && apt-get install -y \
    gcc g++ wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/list/*

COPY --from=poetry-base ${POETRY_VENV} ${POETRY_VENV}

ENV PATH="${PATH}:${POETRY_VENV}/bin"

# Create directory structure and copy files
WORKDIR /opt/mypackage/autoafids
ARG CACHEBUST=1

# For some strange reason, I need to copy the README file; otherwise, 
# the container complains that Poetry didn't install the package.
COPY poetry.lock pyproject.toml README.md ./
COPY autoafids /opt/mypackage/autoafids/autoafids

# Install the pipeline.
RUN poetry install --no-interaction --no-cache --without dev \
    && poetry cache clear pypi --all

# Run the pipeline
ENTRYPOINT ["poetry","run","autoafids"]