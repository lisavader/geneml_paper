import argparse
import re

def parse_gff3(filename):
    lines_by_id = {}
    with open(filename, 'r') as stream:
        for line in stream:
            if line.startswith("#") or not line.strip():
                continue
            record, source, feature, start, end, score, strand, phase, attributes = line.strip().split("\t")
            if feature != "CDS":
                continue

            feature_id = re.search("(?<=ID=cds-).+?(?=;|$)", attributes).group(0)
            try:
                parent = re.search("(?<=Parent=).+?(?=;|$)", attributes).group(0)
            except AttributeError:
                parent = None

            features.append({
                    "feature" : feature,
                    "id" : feature_id,
                    "parent" : parent,
                    "line" : line,
                    "record" : record,
                    "start" : start,
                })
        return features

def select_from_accessions(features, accessions):
    selected_lines = set()
    parents = set()
    sorted_features = sorted(features, key=lambda f: FEATURE_ORDER[f["feature"]],
                             reverse=True)
    for feature in sorted_features:
        if feature["feature"] == "CDS":
            protein_id = feature["id"].replace("cds-","")
            if protein_id in accessions:
                selected_lines.add(feature["line"])
                parents.add(feature["parent"])
        elif feature["id"] in parents:
            selected_lines.add(feature["line"])
            if feature["parent"]:
                parents.add(feature["parent"])
    return selected_lines

def main(accessions, gff_in, gff_out):
    with open(accessions, 'r') as s:
        accs = [line.strip() for line in s]
    features = parse_gff3(gff_in)
    lines = select_from_accessions(features, accs)
    if not lines:
        raise ValueError("No valid features found.")
    with open(gff_out, 'w') as s:
        s.writelines(lines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("accessions", type=str, help="NCBI protein accessions to select for")
    parser.add_argument("gff_in", type=str, help="Input .gff file")
    parser.add_argument("gff_out", type=str, help="Output .gff file")
    args = parser.parse_args()
    main(**vars(args))
