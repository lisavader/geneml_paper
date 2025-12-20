import argparse
import json
import statistics

from collections import defaultdict
from dataclasses import dataclass


@dataclass
class GenomeInfo:
    """ Basic genome statistics
    """
    accession: str
    species: str
    assembly_name: str
    contigs: int
    contig_n50: int
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

    def get_path(self, file_type):
        file_extensions = {
            "cds" : "cds_from_genomic.fna.gz",
            "genome" : "genomic.fna.gz",
            "gff" : "genomic.gff.gz",
            "protein" : "protein.faa.gz",
            "rna" : "rna_from_genomic.fna.gz",
        }
        try:
            ext = file_extensions[file_type]
        except KeyError:
            raise ValueError(f"Unknown file type: {file_type}")

        file_name = "_".join([self.full_name, ext])
        return "/".join([self.download_path, file_name])

    def write_stats(self):
        return "\t".join([str(value) for value in self.__dict__.values()][:-1]) # exclude warnings


def remove_problematic(genomes, contig_n50_threshold):
    clean = []
    preselected = [genome for genome in genomes if not genome.warnings]
    if not preselected:
        return None
    for genome in preselected:
        genes_per_mb = genome.genes / (genome.genome_size / 1_000_000)
        if genome.contig_n50 >= contig_n50_threshold and \
            200 <= genes_per_mb <= 800:
            clean.append(genome)
    return clean

def select_best(genomes, boost_gc=False):
    if not genomes:
        return None
    selected = []
    lower_gc_boundary = 44
    upper_gc_boundary = 58
    gc_boosted_per_genus = 5

    sorted_genomes = sorted(genomes, key=lambda g: (g.in_refseq(), g.contig_n50), reverse=True)
    selected = [sorted_genomes[0]]
    if boost_gc:
        for genome in sorted_genomes[1:]:
            if not lower_gc_boundary < genome.gc_content < upper_gc_boundary and \
                len(selected) < gc_boosted_per_genus:
                selected.append(genome)
    return selected

def main(summary, boost_gc, stats, paths, file_types, exclude_species, contig_n50_threshold):
    genomes_by_genus = defaultdict(list)
    selected_genomes = []
    species_list = []
    if exclude_species:
        species_list = exclude_species.split(',')
    with open(summary, 'r') as stream_in:
        data = json.load(stream_in)
        for report in data["reports"]:
            accession = report["accession"]
            species = report["organism"]["organism_name"]
            if any(s in species for s in species_list):
                print(f"Skipped {accession} of species {species}")
                continue
            assembly_name = report["assembly_info"]["assembly_name"]
            assembly_name = assembly_name.replace(",","").replace(" ","_").replace("__","_")
            contig_n50 = report["assembly_stats"]["contig_n50"]
            contigs = report["assembly_stats"]["number_of_contigs"]
            genome_size = report["assembly_stats"]["total_sequence_length"]
            genes = report["annotation_info"]["stats"]["gene_counts"]["total"]
            gc_content = report["assembly_stats"]["gc_percent"]
            try:
                warnings = report["assembly_info"]["atypical"]["warnings"]
            except KeyError:
                warnings = []
            genus = species.split(" ")[0].strip("[]")
            genome_stats = GenomeInfo(accession, species, assembly_name, int(contigs), int(contig_n50), int(genome_size), int(genes), float(gc_content), warnings)
            genomes_by_genus[genus].append(genome_stats)

    genera = []
    for genus, genomes in genomes_by_genus.items():
        genomes = remove_problematic(genomes, contig_n50_threshold)
        best_genomes = select_best(genomes, boost_gc)
        if best_genomes:
            genera.append(genus)
            for genome in best_genomes:
                selected_genomes.append(genome)

    if stats:
        with open(stats, "w", encoding="utf-8") as s:
            s.write("\t".join(["accession", "species", "assembly_name", "contigs", "contig_n50", "genome_size", "genes", "gc_content"]) + "\n")
            for genome in selected_genomes:
                s.write(genome.write_stats()+"\n")

    if paths:
        with open(paths, "w", encoding="utf-8") as s:
            for genome in selected_genomes:
                if file_types == "all":
                    files = ["cds", "genome", "gff", "protein", "rna"]
                else:
                    files = file_types.split(',')
                for file_type in files:
                    s.write(genome.get_path(file_type)+"\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("summary", type=str,
                        help="NCBI datasets summary .json file")
    parser.add_argument("--boost_gc", action="store_true",
                        help="Boost genomes with high or low GC content")
    parser.add_argument("--stats", type=str, nargs='?', default=None,
                        help="Output genome stats")
    parser.add_argument("--paths", type=str, nargs='?', default=None,
                        help="Output NCBI download paths")
    parser.add_argument("--file_types", type=str, nargs='?', default="genome",
                        help="Which download paths to generate, comma separated. \
                            Choose from 'cds', 'genome', 'gff', 'protein', 'rna', 'all'. (default: %(default)s)")
    parser.add_argument("--exclude-species", type=str, nargs='?', default="",
                        help="Species names to exclude from dataset, comma separated")
    parser.add_argument("--contig_n50_threshold", type=int, nargs='?', default=50000,
                        help="Minimum contig N50 threshold in bp (default: %(default)s)")
    args = parser.parse_args()
    main(**vars(args))
