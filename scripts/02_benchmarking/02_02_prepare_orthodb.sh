#!/bin/bash

# Filter out sequences from genera that are part of the benchmark
python3 scripts/python/filter_orthodb.py \
databases/orthodb/odb12v2_og_aa_fasta \
databases/orthodb/odb12v2_species.tab \
databases/orthodb/odb12v2_levels.tab \
databases/orthodb/odb12v2_level2species.tab \
--output databases/orthodb/odb12v2_og_aa_filtered.faa \
--include fungi \
--exclude aspergillus,botrytis,cercospora,fusarium,saccharomyces,cryptococcus,puccinia,somion,rhizophagus
