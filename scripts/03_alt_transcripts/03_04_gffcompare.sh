#!/bin/bash

python3 scripts/python/simplify_gff.py isoseq_data/"Transcript-annotation-5'collapsed.gff3" \
isoseq_data/transcript_annotation_simplified.gff3

for tool in augustus geneml; do
    gffcompare -r isoseq_data/transcript_annotation_simplified.gff3 \
    isoseq_data/PH-1_YL1_$tool.gff3 -T --no-exon-merge --strict-match -e 0
    mkdir -p gffcompare_isoseq
    mv gffcmp.stats gffcompare_isoseq/$tool.stats
done

#clean up
rm gffcmp.*
