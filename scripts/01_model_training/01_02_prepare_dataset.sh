#!/bin/bash

ln -s geneML/trainer trainer

mkdir -p annotations

# Build annotation TSV files
for f in training_genomes/*.gff; do
    python3 trainer/build_annotation_tsv_from_gff.py \
    $(basename $f .gff) \
    $f \
    annotations
done

# Build data files with validated genes
for f in training_genomes/*.gff; do
    bash trainer/create_datafile.sh \
    $(basename $f .gff) \
    datafile \
    annotations \
    training_genomes
done

# Reduce data files for benchmarking
mkdir -p datafile_benchmarking
for f in datafile/*.h5; do
  name=$(basename $f .h5)
  echo $name | grep -qf data/exclude_benchmarking.txt \
  || ln -s $(realpath $f) datafile_benchmarking/$name.h5
done

# Create dataset for regular model
python3 trainer/create_dataset.py \
--data-dir datafile \
--multigenome --mode all \
--suffix final \
--CL_max 800 --SL 800 \
--max-genes-per-genome 10000

# Create dataset for benchmarking model
python3 trainer/create_dataset.py \
--data-dir datafile_benchmarking \
--multigenome --mode all \
--suffix benchmarking \
--CL_max 800 --SL 800 \
--max-genes-per-genome 10000
