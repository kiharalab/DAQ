#!/bin/bash

ENV_NAME=${1:-daq}

conda create -n "$ENV_NAME" python=3.10 gcc gxx -c conda-forge -y
conda activate "$ENV_NAME"

pip3 install -r requirements.txt -r torch-requirements.txt
