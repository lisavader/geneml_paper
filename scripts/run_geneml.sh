for genome in $(ls benchmarking_dataset)
do for path in $(cat model_paths.txt)
do model=$(basename $path .keras)
if [ ! -f geneml/$genome/${genome}_${model}.gff ]
then mkdir -p geneml/$genome
geneml benchmarking_dataset/$genome/*.fna geneml/$genome/${genome}_${model}.gff --model $path 2> geneml/$genome/${genome}_${model}_stderr.txt
fi
done
done
