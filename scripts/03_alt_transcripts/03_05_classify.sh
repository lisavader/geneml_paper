
# Add transcript classes to the reference GFF3 files
for ref in alt_transcript_data/*/simplified_annotation.gff3; do
    geneml-classify-gff3 $ref $(dirname $ref)/simplified_annotation_w_classes.gff3
done

# Count exact matches per transcript class
for tool in augustus geneml; do
    for result in alt_transcript_data/*/gffcmp_$tool/gffcmp.tracking; do
        species=$(echo $result | cut -d'/' -f2)
        python scripts/python/count_matches_per_transcript_class.py \
        alt_transcript_data/$species/simplified_annotation_w_classes.gff3 \
        $result \
        alt_transcript_data/$species/alt_transcript_counts_$tool.tsv
    done
done
