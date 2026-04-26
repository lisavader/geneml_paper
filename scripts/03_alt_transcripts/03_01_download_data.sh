#!/bin/bash

DATA_PATH=http://123.56.75.179/sites/default/files/wwwroot/fgbase

wget -P isoseq_data \
$DATA_PATH/Genome-sequence.fasta.zip \
"$DATA_PATH/Transcript-annotation-5'collapsed.gff3.zip"

for f in isoseq_data/*.zip; do
    unzip "$f" -d isoseq_data
done
