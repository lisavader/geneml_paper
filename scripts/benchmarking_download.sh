## datasets version: 18.4.1
## Download dataset
datasets download genome accession GCA_964035595.1 GCF_000143535.2 GCF_000146045.2 GCF_000240135.3 GCF_000835755.1 GCF_009017415.1 GCF_021901695.1 GCF_026210795.1 GCF_033473495.1 --include genome,gff3
mv ncbi_dataset/data/GC* benchmarking_dataset/

## Remove mitochondrial sequences and annotations
for genome in $(ls benchmarking_dataset); do
    awk 'NR==FNR {ids[$1]; next} /^>/ {id=$1; sub(/^>/,"",id); drop=(id in ids)} !drop' mitochondrial_ids.txt benchmarking_dataset/$genome/*genomic.fna > benchmarking_dataset/$genome/chromosomes.fna
    awk 'NR==FNR {ids[$1]; next} /^#/ {print; next} !($1 in ids)' mitochondrial_ids.txt benchmarking_dataset/$genome/genomic.gff > benchmarking_dataset/$genome/chromosomes.gff
done
