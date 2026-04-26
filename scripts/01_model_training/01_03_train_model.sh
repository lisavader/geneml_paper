#!/bin/bash

# Train regular model
python3 geneML/trainer/train_model.py \
--num-epochs 30 \
--context-length 800 \
--dataset-name final \
--train-path datafile/dataset_all_800_800_final.h5

# Train benchmarking model
python3 geneML/trainer/train_model.py \
--num-epochs 30 \
--context-length 800 \
--dataset-name benchmarking \
--train-path datafile_benchmarking/dataset_all_800_800_benchmarking.h5
