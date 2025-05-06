FROM condaforge/miniforge3:latest

WORKDIR /src/

# Copy your code
COPY . /src/

# Disable user site packages
ENV PYTHONNOUSERSITE=1

# Use bash for the following RUNs
SHELL ["/bin/bash", "-c"]

# ---- ONE SINGLE RUN ----
RUN set -e && \
    conda install -n base -c conda-forge mamba -y && \
    mamba create -y -n snakebids-env -c conda-forge -c bioconda snakebids unzip && \
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate snakebids-env && \
    ./autoafids/run.py test_data/bids_T1w test_out participant --participant-label 001 --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./autoafids/run.py test_data/bids_T2w test_out participant --modality T2w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./autoafids/run.py test_data/bids_ct test_out participant --modality ct --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./autoafids/run.py test_data/bids_T1w test_out participant --participant-label 001 --stereotaxy STN --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./autoafids/run.py test_data/bids_T1w test_out participant --participant-label 001 --fidqc --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    conda clean --all -y && \
    rm -rf /opt/conda/pkgs /root/.cache && \
    cp /src/entrypoint_quick.sh

# Set snakemake profile
ENV SNAKEMAKE_PROFILE=/src/profiles/docker-conda

# Set entrypoint
ENTRYPOINT ["/src/entrypoint.sh"]

