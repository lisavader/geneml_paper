import argparse

from Bio import SeqIO

def get_gff_line(seqid, feature_type, start, end, strand, feature_id, parent_id = None):
    attr = f"ID={feature_id}"
    if parent_id:
        attr = ";".join([attr,f"Parent={parent_id}"])

    gff_line = f"{seqid}\tGenBank\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{attr}\n"
    return gff_line

def get_feature_id(feature_type, start, end):
    return "_".join([feature_type, str(start), str(end)])

def main(gbk_path, gff3_path):
    with open(gff3_path, "w") as gff_out:
        gff_lines = {}
        for record in SeqIO.parse(gbk_path, "genbank"):
            seqid = record.id
            feature_id = None
            last_id = None
            for feature in record.features:
                strand = "+" if feature.location.strand == 1 else "-"
                start = int(feature.location.start) + 1
                end = int(feature.location.end)
                feature_id = get_feature_id(feature.type, start, end)

                if feature.type in ["source", "subregion", "region"] or feature_id in gff_lines:
                    continue

                if feature.type in ["gene", "mRNA"]:
                    last_id = None
                    for feature in ["gene", "mRNA"]:
                        feature_id = get_feature_id(feature, start, end)
                        gff_lines[feature_id] = get_gff_line(seqid, feature, start, end, strand, feature_id, last_id)
                        last_id = feature_id

                elif feature.type == "CDS":
                    for location in feature.location.parts:
                        exon_start = int(location.start) + 1
                        exon_end = int(location.end)
                        feature_id = get_feature_id("CDS", start, end)
                        gff_lines[feature_id] = get_gff_line(seqid, "CDS", exon_start, exon_end, strand, feature_id, last_id)

                else:
                    gff_lines[feature_id] = get_gff_line(seqid, feature.type, start, end, strand, feature_id, last_id)

        for line in gff_lines.values():
            gff_out.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gbk_path")
    parser.add_argument("gff3_path")
    args = parser.parse_args()
    main(**vars(args))
