import argparse
import json
import statistics

from collections import defaultdict
from dataclasses import dataclass
from ete3 import NCBITaxa


@dataclass
class GenomeInfo:
    """ Basic genome statistics
    """
    accession: str
    assembly_name: str
    scaffolds: int
    scaffold_n50: int
    genome_size: int
    genes: int
    gc_content: float
    warnings: list

    def in_refseq(self):
        return self.accession.startswith("GCF")

    @property
    def full_name(self):
        return "_".join([self.accession, self.assembly_name])

    @property
    def download_path(self):
        prefix, numerical_id = self.accession.split("_")
        numbers = numerical_id.split(".")[0]
        split_id = "/".join([numbers[:3],numbers[3:6],numbers[6:9]])
        return "/".join(["https://ftp.ncbi.nlm.nih.gov/genomes/all", prefix, split_id, self.full_name])

    def get_genome_path(self):
        file_name = "_".join([self.full_name, "genomic.fna.gz"])
        return "/".join([self.download_path, file_name])

    def get_gff_path(self):
        file_name = "_".join([self.full_name, "genomic.gff.gz"])
        return "/".join([self.download_path, file_name])

    def write_stats(self):
        return "\t".join([str(value) for value in self.__dict__.values()])


def remove_problematic(genomes):
    clean = []
    preselected = [genome for genome in genomes if not genome.warnings]
    if not preselected:
        return None
    median_size = statistics.median([genome.genome_size for genome in preselected])
    median_genes = statistics.median([genome.genes for genome in preselected])
    for genome in preselected:
        if genome.scaffold_n50 >= 100000 and \
           2/3*median_size < genome.genome_size < (1+1/3)*median_size and \
           2/3*median_genes < genome.genes < (1+1/3)*median_genes:
            clean.append(genome)
    return clean

def select_best(genomes, boost_gc=False):
    if not genomes:
        return None
    selected = []
    lower_gc_boundary = 44
    upper_gc_boundary = 58
    gc_boosted_per_genus = 5

    sorted_genomes = sorted(genomes, key=lambda g: (g.in_refseq(), g.scaffold_n50), reverse=True)
    selected = [sorted_genomes[0]]
    if boost_gc:
        for genome in sorted_genomes[1:]:
            if not lower_gc_boundary < genome.gc_content < upper_gc_boundary and \
                len(selected) < gc_boosted_per_genus:
                selected.append(genome)
    return selected

def main(summary, boost_gc, stats, paths):
    genomes_by_genus = defaultdict(list)
    selected_genomes = []
    with open(summary, 'r') as stream_in:
        data = json.load(stream_in)
        for report in data["reports"]:
            accession = report["accession"]
            assembly_name = report["assembly_info"]["assembly_name"]
            assembly_name = assembly_name.replace(",","").replace(" ","_").replace("__","_")
            scaffold_n50 = report["assembly_stats"]["scaffold_n50"]
            scaffolds = report["assembly_stats"]["number_of_scaffolds"]
            genome_size = report["assembly_stats"]["total_sequence_length"]
            genes = report["annotation_info"]["stats"]["gene_counts"]["total"]
            gc_content = report["assembly_stats"]["gc_percent"]
            try:
                warnings = report["assembly_info"]["atypical"]["warnings"]
            except KeyError:
                warnings = []
            genome_stats = GenomeInfo(accession, assembly_name, int(scaffolds), int(scaffold_n50), int(genome_size), int(genes), float(gc_content), warnings)
            genus = report["organism"]["organism_name"].split(" ")[0].strip("[]")
            genomes_by_genus[genus].append(genome_stats)

    genera = []
    for genus, genomes in genomes_by_genus.items():
        genomes = remove_problematic(genomes)
        best_genomes = select_best(genomes, boost_gc)
        if best_genomes:
            genera.append(genus)
            for genome in best_genomes:
                selected_genomes.append(genome)

    phylum_counter = defaultdict(int)
    for taxids in NCBITaxa().get_name_translator(genera).values():
        lineage = NCBITaxa().get_rank(NCBITaxa().get_lineage(taxids[0]))
        for taxid, rank in lineage.items():
            if rank == "phylum":
                phylum_counter[taxid] += 1
                if taxid == 6040:
                    print(taxids)
    print(phylum_counter)

    if stats:
        with open(stats, "w", encoding="utf-8") as s:
            s.write("\t".join(["accession", "assembly_name", "scaffolds", "scaffold_n50", "genome_size", "genes", "gc_content", "warnings"])+"\n")
            for genome in selected_genomes:
                s.write(genome.write_stats()+"\n")

    if paths:
        with open(paths, "w", encoding="utf-8") as s:
            for genome in selected_genomes:
                s.write(genome.get_genome_path()+"\n")
                s.write(genome.get_gff_path()+"\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("summary", type=str,
                        help="NCBI datasets summary .json file")
    parser.add_argument("--boost_gc", action="store_true",
                        help="Boost genomes with high or low GC content")
    parser.add_argument("--stats", type=str, nargs='?', default=None,
                        help="Output genome stats")
    parser.add_argument("--paths", type=str, nargs='?', default=None,
                        help="Output genome download paths")
    args = parser.parse_args()
    main(**vars(args))
