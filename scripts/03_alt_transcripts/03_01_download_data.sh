#!/bin/bash

## Fusarium graminearum
DIR=alt_transcript_data/F_graminearum
mkdir -p $DIR

DATA_PATH=http://123.56.75.179/sites/default/files/wwwroot/fgbase
wget -P $DIR \
$DATA_PATH/Genome-sequence.fasta.zip \
"$DATA_PATH/Transcript-annotation-5'collapsed.gff3.zip"

for f in $DIR/*.zip; do
    unzip "$f" -d $DIR
    rm $f
done

mv $DIR/Genome-sequence.fasta/PH-1_YL1.fasta $DIR/PH-1_YL1.fna
rm -r $DIR/Genome-sequence.fasta

## Zymoseptoria tritici
DIR=alt_transcript_data/Z_tritici
mkdir -p $DIR

wget -O $DIR/z.tritici.IP0323.reannot.isoforms.gff3 \
https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/GHHLA0
wget -P $DIR \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/219/625/GCF_000219625.1_MYCGR_v2.0/GCF_000219625.1_MYCGR_v2.0_genomic.fna.gz
gunzip $DIR/GCF_000219625.1_MYCGR_v2.0_genomic.fna.gz

## Ganoderma lingzhi
DIR=alt_transcript_data/G_lingzhi
mkdir -p $DIR

wget -O $DIR/Ganoderma_lucidum_GL102-53.genome.zip \
http://www.gpgenome.com/download/346
unzip alt_transcript_data/G_lingzhi/Ganoderma_lucidum_GL102-53.genome.zip -d $DIR

mv alt_transcript_data/G_lingzhi/GL102-53.genome.fasta alt_transcript_data/G_lingzhi/GL102-53.genome.fna
rm alt_transcript_data/G_lingzhi/Ganoderma_lucidum_GL102-53.genome.zip

## Coprinopsis cinerea
DIR=alt_transcript_data/C_cinerea
mkdir -p $DIR

DATA_PATH=https://mushroomdb.brc.hu/files/genome_and_genes
wget -P $DIR \
$DATA_PATH/CopciAB_new_jgi_20220113.fasta.gz \
$DATA_PATH/CopciAB_new_jgi_20220113.gtf.gz
gunzip alt_transcript_data/C_cinerea/*.gz

mv alt_transcript_data/C_cinerea/CopciAB_new_jgi_20220113.fasta alt_transcript_data/C_cinerea/CopciAB_new_jgi_20220113.fna
