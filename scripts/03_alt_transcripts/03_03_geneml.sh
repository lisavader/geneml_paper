#!/bin/bash

# Run geneML
for genome in alt_transcript_data/*/*.fna; do
    geneml $genome -o $(dirname $genome)/geneml.gff3 -c8;
done
