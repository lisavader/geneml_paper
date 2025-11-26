import argparse
import gffutils


def format_gff_feature(feature):
    # Build a GFF3 line from a gffutils feature, preserving input format
    attr_str = ";".join([f"{k}={','.join(v)}" for k, v in feature.attributes.items()])
    fields = [
        feature.seqid,
        feature.source,
        feature.featuretype,
        str(feature.start),
        str(feature.end),
        str(feature.score) if feature.score is not None else '.',
        feature.strand if feature.strand is not None else '.',
        feature.frame if feature.frame is not None else '.',
        attr_str
    ]
    return "\t".join(fields)


def main(input, output, accession):
    gff_db = gffutils.create_db(input, dbfn=":memory:", merge_strategy="create_unique")

    gff_header = [
        '##gff-version 3',
        f'#!genome-build-accession NCBI_Assembly: {accession}',
        '##Note: UTRs removed; features span CDS only',
    ]
    with open(output, "w", encoding="utf-8") as out_f:
        for line in gff_header:
            out_f.write(line + "\n")

        # Group by sequence_id, then for each gene: gene, mRNA, CDSes
        genes = sorted(list(gff_db.features_of_type("gene")), key=lambda g: (g.seqid, g.start))
        for gene in genes:
            mrnas = sorted(list(gff_db.children(gene, featuretype="mRNA")), key=lambda m: m.start)
            # Update gene coordinates to span only CDS regions
            cds_starts = []
            cds_ends = []
            mrna_cdses_dict = {}
            for mrna in mrnas:
                cdses = sorted(list(gff_db.children(mrna, featuretype="CDS")), key=lambda c: c.start)
                mrna_cdses_dict[mrna.id] = cdses
                if cdses:
                    cds_starts.extend([cds.start for cds in cdses])
                    cds_ends.extend([cds.end for cds in cdses])
            if cds_starts and cds_ends:
                gene.start = min(cds_starts)
                gene.end = max(cds_ends)
            out_f.write(format_gff_feature(gene) + "\n")
            for mrna in mrnas:
                cdses = mrna_cdses_dict[mrna.id]
                # Update mRNA coordinates to span only CDS regions
                if cdses:
                    mrna.start = min(cds.start for cds in cdses)
                    mrna.end = max(cds.end for cds in cdses)
                out_f.write(format_gff_feature(mrna) + "\n")
                for cds in cdses:
                    out_f.write(format_gff_feature(cds) + "\n")


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
