import argparse
import re

from collections import defaultdict

def parse_gff3(filename):
    features_by_record = defaultdict(list)
    with open(filename) as stream:
        for line in stream:
            if line.startswith("#") or not line.strip():
                continue
            record, source, feature, start, end, score, strand, phase, attributes = line.strip().split("\t")
            start, end = int(start), int(end)

            try:
                parent = re.search("(?<=Parent=).+?(?=;)", attributes).group(0)
            except AttributeError:
                parent = None

            if feature in ["region", "gene", "mRNA", "exon", "CDS"]:
                features_by_record[record].append({
                    "feature" : feature,
                    "start" : start,
                    "end" : end,
                    "strand" : strand,
                    "parent" : parent,
                })
        return features_by_record

def add_intergenic(features):
    gene_coordinates = sorted([[f["start"],f["end"]] for f in features if f["feature"] == "gene"])
    if not gene_coordinates:
        return features

    merged = [gene_coordinates[0]]
    for gene in gene_coordinates[1:]:
        #If there is a distance between the start of this gene and the end of the previous
        if gene[0] - merged[-1][1] > 1:
            merged.append(gene)
        else:
            merged[-1] = [merged[-1][0], max(merged[-1][1], gene[1])]

    regions = [feature for feature in features if feature["feature"] == "region"]
    assert len(regions) == 1
    region = regions[0]
    assert region["start"] == 1

    starts = [1]+[p[1]+1 for p in merged]
    ends = [p[0]-1 for p in merged]+[region["end"]]
    for start, end in zip(starts, ends):
        if start < end:
            features.append({
                "feature" : "intergenic",
                "start" : start,
                "end" : end,
                "strand" : None,
                "parent" : None,
            })
    return features

def main(gff_file, output):
    features_by_record = parse_gff3(gff_file)
    for record, features in features_by_record.items():
        features = add_intergenic(features)
    with open(output, 'w') as stream:
        header = ["record","feature","start","end","strand","parent"]
        stream.write('\t'.join(header)+'\n')
        for record, feature_list in features_by_record.items():
            for feature in feature_list:
                items = [record] + list(feature.values())
                stream.write('\t'.join([str(item) for item in items])+'\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gff_file")
    parser.add_argument("output")
    args = parser.parse_args()
    main(**vars(args))
