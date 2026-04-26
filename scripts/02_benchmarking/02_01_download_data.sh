#!/bin/bash

## Download dataset
#datasets download genome accession GCA_964035595.1 GCF_000143535.2 GCF_000146045.2 GCF_000240135.3 \
#GCF_000835755.1 GCF_009017415.1 GCF_021901695.1 GCF_026210795.1 GCF_033473495.1 \
#--include genome,gff3

unzip ncbi_dataset.zip

mkdir -p benchmarking_dataset
mv ncbi_dataset/data/GC* benchmarking_dataset/

for genome in benchmarking_dataset/*; do
    # Step 1: Trim headers and remove mitochondrial sequences from the FASTA file
    awk '
        NR == FNR {
            # Read mitochondrial IDs into an array
            ids[$1]; next
        }
        /^>/ {
            # Process header lines
            id = $1; sub(/^>/, "", id);
            drop = (id in ids)  # Check if the ID is in the mitochondrial list
            if (!drop) print ">" id  # Print the ID if not in the list
            next
        }
        !drop {
            # Print non-header lines if the sequence is not mitochondrial
            print
        }
    ' data/mitochondrial_ids.txt $genome/*genomic.fna > $genome/chromosomes.fna

    # Step 2: Remove mitochondrial features from the GFF file
    awk '
        NR == FNR {
            # Read mitochondrial IDs into an array
            ids[$1]; next
        }
        /^#/ {
            # Print comment lines as-is
            print; next
        }
        !($1 in ids) {
            # Print lines where the first column is not in the mitochondrial list
            print
        }
    ' data/mitochondrial_ids.txt $genome/genomic.gff > $genome/chromosomes.gff

    # Step 3: Remove non-CDS features
    python3 scripts/python/simplify_gff.py $genome/chromosomes.gff \
    $genome/chromosomes_simplified.gff --accession $(basename $genome)

done
