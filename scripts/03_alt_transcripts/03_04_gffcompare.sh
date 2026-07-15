#!/bin/bash

# Remove UTRs from GFFs (we only compare the CDS)
python3 scripts/python/simplify_gff.py \
"alt_transcript_data/F_graminearum/Transcript-annotation-5'collapsed.gff3" \
alt_transcript_data/F_graminearum/simplified_annotation.gff3 \
--write-transcript-counts alt_transcript_data/F_graminearum/transcript_counts.tsv

python3 scripts/python/simplify_gff.py \
alt_transcript_data/Z_tritici/z.tritici.IP0323.reannot.isoforms.gff3 \
alt_transcript_data/Z_tritici/simplified_annotation.gff3 \
--contig-map data/contig_mapping_Z_tritici.json \
--write-transcript-counts alt_transcript_data/Z_tritici/transcript_counts.tsv

python3 scripts/python/simplify_gff.py \
alt_transcript_data/G_lingzhi/GL102-53.Annotations.gff3 \
alt_transcript_data/G_lingzhi/simplified_annotation.gff3 \
--write-transcript-counts alt_transcript_data/G_lingzhi/transcript_counts.tsv

python3 scripts/python/simplify_gff.py \
alt_transcript_data/C_cinerea/CopciAB_new_jgi_20220113.gtf \
alt_transcript_data/C_cinerea/simplified_annotation.gff3 \
--gtf \
--split-gene-ids \
--write-transcript-counts alt_transcript_data/C_cinerea/transcript_counts.tsv

for tool in augustus geneml; do
    for result in alt_transcript_data/*/$tool.gff3; do
        dir=$(dirname $result)
        gffcompare -r $dir/simplified_annotation.gff3 $result \
        -T --no-exon-merge --strict-match -e 0
        mkdir -p $dir/gffcmp_$tool
        mv gffcmp* $dir/gffcmp_$tool
    done
done
