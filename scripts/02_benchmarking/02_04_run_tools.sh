#!/bin/bash

# Build results directories
for tool in annevo augustus braker3 geneml helixer_with_UTR; do
    dir=benchmarking_results/$tool
    mkdir -p $dir
    chmod 777 $dir
done

DOCKER_OPTIONS="--rm --cpuset-cpus=0-7 --cpuset-mems=0 -v $(pwd):/data"

# ANNEVO scripts
rm -f annevo_scripts
for a in $(ls benchmarking_dataset); do
    echo "docker run $DOCKER_OPTIONS --shm-size=8g \
annevo python annevo/annotation.py \
--genome /data/benchmarking_dataset/$a/chromosomes.fna \
--model_path annevo/ANNEVO_model/ANNEVO_Fungi.pt \
--output /data/benchmarking_results/annevo/$a.gff3 \
--threads 8" \
    >> annevo_scripts
done

# AUGUSTUS scripts
rm -f augustus_scripts
for a in $(ls benchmarking_dataset); do
    species_model=$(awk -F',' -v val="$a" '$1==val {print $2}' data/augustus_models.csv)
    echo "docker run $DOCKER_OPTIONS \
augustus augustus \
/data/benchmarking_dataset/$a/chromosomes.fna \
--species=$species_model \
--gff3=on \
> $(pwd)/benchmarking_results/augustus/$a.gff3" \
    >> augustus_scripts
done

# BRAKER3 scripts
rm -f braker3_scripts
for a in $(ls benchmarking_dataset); do
    wd=$(echo benchmarking_results/braker3/$a)
    echo "mkdir -p $wd; chmod 777 $wd; \
docker run $DOCKER_OPTIONS \
teambraker/braker3:v3.0.7.6 braker.pl \
--genome=/data/benchmarking_dataset/$a/chromosomes.fna \
--prot_seq=/data/orthodb/odb12v2_og_aa_filtered.faa \
--workingdir=/data/$wd \
--threads=8 \
--fungus \
--gff3" \
    >> braker3_scripts
done

# geneML scripts
rm -f geneml_scripts
for a in $(ls benchmarking_dataset); do
    echo "docker run $DOCKER_OPTIONS \
geneml geneml \
/data/benchmarking_dataset/$a/chromosomes.fna \
-o /data/benchmarking_results/geneml/$a.gff3 \
--model /data/models/benchmarking_final.keras \
--context-length 800 \
--cores 8 \
--cpu-only" \
    >> geneml_scripts
done

# Helixer scripts
rm -f helixer_scripts
for a in $(ls benchmarking_dataset); do
    echo "docker run $DOCKER_OPTIONS \
helixer-with-models Helixer.py \
--fasta-path /data/benchmarking_dataset/$a/chromosomes.fna \
--gff-output-path /data/benchmarking_results/helixer_with_UTR/$a.gff3 \
--lineage fungi" \
    >> helixer_scripts
done

cat annevo_scripts augustus_scripts braker3_scripts geneml_scripts helixer_scripts > all_scripts

# Run with parallel for job logging
parallel --joblog benchmarking_log.txt -j1 < all_scripts
