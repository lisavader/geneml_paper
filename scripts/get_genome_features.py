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
            assert end >= start

            feature_id = re.search("(?<=ID=).+?(?=;)", attributes).group(0)
            try:
                parent = re.search("(?<=Parent=).+?(?=;)", attributes).group(0)
            except AttributeError:
                parent = None

            if feature in ["region", "gene", "mRNA", "exon", "CDS"]:
                features_by_record[record].append({
                    "feature" : feature,
                    "id" : feature_id,
                    "start" : start,
                    "end" : end,
                    "strand" : strand,
                    "parent" : parent,
                })
        return features_by_record

def add_intergenic(features):
    gene_features = [f for f in features if f["feature"] == "gene"]
    if not gene_features:
        return features

    gene_features = sorted(gene_features, key=lambda f: f["start"])
    [region] = [feature for feature in features if feature["feature"] == "region"]
    assert region["start"] == 1
    intergenic_start = 1

    for gene in gene_features:
        intergenic_end = gene["start"]-1
        if intergenic_end >= intergenic_start:
            features.append({
                "feature" : "intergenic",
                "id" : None,
                "start" : intergenic_start,
                "end" : intergenic_end,
                "strand" : None,
                "parent" : None,
            })
        intergenic_start = gene["end"]+1

    if region["end"] >= intergenic_start:
        features.append({
            "feature" : "intergenic",
            "id" : None,
            "start" : intergenic_start,
            "end" : region["end"],
            "strand" : None,
            "parent" : None,
        })
    return features

def add_introns(features):
    mRNA_features = [f for f in features if f["feature"] == "mRNA"]
    for feature in mRNA_features:
        exons = [f for f in features if f["feature"] == "exon" and f["parent"] == feature["id"]]
        exons = sorted(exons, key=lambda e: e["start"])
        if len(exons) > 1:
            (strand,) = set(exon["strand"] for exon in exons)
            intron_start = exons[0]["end"]+1
            for exon in exons[1:]:
                intron_end = exon["start"]-1
                assert intron_end >= intron_start
                features.append({
                    "feature" : "intron",
                    "id" : None,
                    "start" : intron_start,
                    "end" : intron_end,
                    "strand" : strand,
                    "parent" : feature["id"],
                })
                intron_start = exon["end"]+1
    return features

def main(gff_file, output):
    features_by_record = parse_gff3(gff_file)
    for record, features in features_by_record.items():
        features = add_intergenic(features)
        features = add_introns(features)
    with open(output, 'w') as stream:
        header = ["record","feature","id","start","end","strand","parent"]
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
