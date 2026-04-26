#!/bin/bash

# Run AUGUSTUS with alternative transcript support
docker run -v $(pwd):/data --rm augustus \
augustus /data/isoseq_data/Genome-sequence.fasta/PH-1_YL1.fasta \
--species=fusarium_graminearum \
--gff3=on \
--alternatives-from-sampling=true \
--maxtracks=5 \
> isoseq_data/PH-1_YL1_augustus.gff3
