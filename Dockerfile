# --- Stage 1: Build & Environment Prep ---
FROM ghcr.io/prefix-dev/pixi:latest AS build

WORKDIR /src
COPY . /src

# Create the environment based on your lockfile and install dependencies
# This creates an isolated, self-contained environment in /src/.pixi/envs/default
RUN pixi install --locked

# Run your autoafids pre-caching tasks inside the pixi environment
RUN pixi run ./autoafids/run.py test_data/bids_T1w test_out participant --participant-label 002 --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    pixi run ./autoafids/run.py test_data/bids_T1w test_out participant --participant-label 002 --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs --detect_with_nnlm && \
    pixi run ./autoafids/run.py test_data/bids_T2w test_out participant --modality T2w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    pixi run ./autoafids/run.py test_data/bids_ct test_out participant --modality ct --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    pixi run ./autoafids/run.py test_data/bids_T1w test_out participant --participant-label 002 --stereotaxy STN --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    pixi run ./autoafids/run.py test_data/bids_T1w test_out participant --participant-label 002 --fidqc --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs

# --- Stage 2: Tiny Production Image ---
# Use a lightweight glibc base image instead of a heavy Conda/Miniforge image
FROM ubuntu:24.04 AS production

WORKDIR /src

# Copy your source code, pre-cached conda environments, and the built Pixi environment
COPY --from=build /src /src

# Add the pixi environment binaries straight to the PATH so you don't have to "activate" it
ENV PATH="/src/.pixi/envs/default/bin:$PATH"
ENV SNAKEMAKE_PROFILE=/src/autoafids/workflow/profiles/docker-conda
ENV PYTHONNOUSERSITE=1
ENV SNAKEMAKE_PROFILE=/src/autoafids/autoafids/workflow/profiles/docker-conda

ENTRYPOINT ["/src/entrypoint.sh"]