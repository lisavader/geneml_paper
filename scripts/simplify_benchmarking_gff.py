import argparse
import gffutils

FEATURE_ORDER = {"gene": 0, "mRNA": 1, "CDS": 2}

def filter_features(gff_db):
    filtered_features = []

    transcripts = gff_db.features_of_type("mRNA")
    for transcript in transcripts:
        cds_list = list(gff_db.children(transcript, featuretype="CDS"))
        if not cds_list:
            gff_db.delete(transcript.id)
            continue
        filtered_features.extend(cds_list)

        # update transcript coordinates
        transcript.start = min(c.start for c in cds_list)
        transcript.end   = max(c.end for c in cds_list)

    genes = gff_db.features_of_type("gene")
    for gene in genes:
        gene_transcripts = list(gff_db.children(gene, featuretype="mRNA"))
        if not gene_transcripts:
            continue
        filtered_features.extend(gene_transcripts)

        # update gene coordinates
        gene.start = min(t.start for t in gene_transcripts)
        gene.end   = max(t.end for t in gene_transcripts)
        filtered_features.append(gene)

    return filtered_features


def main(input, output, accession):
    gff_db = gffutils.create_db(input, dbfn=":memory:", merge_strategy="create_unique")
    filtered_features = filter_features(gff_db)

    gff_header = [
        '##gff-version 3',
        f'#!genome-build-accession NCBI_Assembly: {accession}',
        '##Note: UTRs removed; features span CDS only',
    ]
    with open(output, "w", encoding="utf-8") as out_f:
        for line in gff_header:
            out_f.write(line + "\n")
        for feature in sorted(filtered_features, key=lambda x: (
            x.seqid, x.start, FEATURE_ORDER.get(x.featuretype, 99))):
            out_f.write(str(feature) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str,
                        help="Input GFF file to clean")
    parser.add_argument("output", type=str,
                        help="Output GFF file")
    parser.add_argument("accession", type=str,
                        help="Genome accession")
    args = parser.parse_args()
    main(**vars(args))
