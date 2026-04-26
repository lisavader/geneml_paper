#!/bin/bash

# Build Docker images

# ANNEVO
docker build -t annevo docker_files/annevo

# AUGUSTUS
git clone --branch v3.5.0 https://github.com/Gaius-Augustus/Augustus.git
docker build -t augustus Augustus

# BRAKER3
docker pull teambraker/braker3:v3.0.7.6

# geneML
docker build -t geneml docker_files/geneml

# Helixer
docker pull gglyptodon/helixer-docker:helixer_v0.3.6_cuda_12.2.2-cudnn8_1
docker run --name helixer-setup gglyptodon/helixer-docker:helixer_v0.3.6_cuda_12.2.2-cudnn8_1 \
python3 /home/helixer_user/Helixer/scripts/fetch_helixer_models.py --lineage fungi
docker commit helixer-setup helixer-with-models
docker rm helixer-setup
