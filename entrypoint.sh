#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate snakebids-env
exec /src/autoafids/run.py "$@" 