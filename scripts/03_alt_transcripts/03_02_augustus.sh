#!/bin/bash

# Run AUGUSTUS with alternative transcript support
# Note: AUGUSTUS does not come with gene models for G. lingzhi and Z. tritici

## Fusarium graminearum
DIR=alt_transcript_data/F_graminearum

docker run -v $(pwd):/data --rm augustus \
augustus /data/$DIR/PH-1_YL1.fna \
--species=fusarium_graminearum \
--gff3=on \
--alternatives-from-sampling=true \
--maxtracks=5 \
> $DIR/augustus.gff3

## Coprinopsis cinerea
DIR=alt_transcript_data/C_cinerea

docker run -v $(pwd):/data --rm augustus \
augustus /data/$DIR/CopciAB_new_jgi_20220113.fna \
--species=coprinus_cinereus \
--gff3=on \
--alternatives-from-sampling=true \
--maxtracks=5 \
> $DIR/augustus.gff3
