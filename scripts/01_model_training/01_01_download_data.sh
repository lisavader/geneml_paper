#!/bin/bash

python3 scripts/python/select_training_set.py \
data/fungi_summary.json \
--stats training_set_stats.tsv \
--paths training_set_paths.txt \
--file-types genome,gff \

# Upon manual inspection, genomes with bad genus annotations were excluded
grep -v -E 'GCA_000412225.2|GCA_002572855.1|GCA_043643485.1|GCA_000825705.1' training_set_paths.txt \
> training_set_paths_clean.txt

while read url; do
    wget $url -P training_genomes;
done < training_set_paths_clean.txt

gunzip training_genomes/*.gz
