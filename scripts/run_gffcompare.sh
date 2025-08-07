for genome in $(ls)
do
cd $genome
gffcompare -T -r ../../benchmarking_dataset/$genome/genomic.gff *.gff
cd ..
done
