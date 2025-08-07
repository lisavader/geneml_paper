import argparse
import glob
import os

from collections import defaultdict

from antismash.common.secmet.record import Record

COLUMN_NAMES = ["record",
                "core",
                "additional",
                "transport",
                "regulatory",
                "resistance",
                "other",
                "antismash_domains",
                "pfam_domains",
                ]

def get_counts(record):
    counts = defaultdict(int)
    for cds in record.get_cds_features():
        counts[cds.gene_function.name] += 1
    counts["antismash_domain"] = len(record.get_antismash_domains())
    counts["pfam_domain"] = len(record.get_pfam_domains())
    return counts

def main(input_dir, output_tsv):
    with open(output_tsv, 'w') as out:
        out.write("\t".join(COLUMN_NAMES)+"\n")
        for file in glob.glob(input_dir+"/*.gbk"):
            file_name = os.path.basename(file).rstrip("*.gbk")
            record = Record.from_genbank(file, taxon="fungi")[0]
            counts = get_counts(record)
            stats = [file_name,
                     counts.get("CORE", 0),
                     counts.get("ADDITIONAL", 0),
                     counts.get("TRANSPORT", 0),
                     counts.get("REGULATORY", 0),
                     counts.get("RESISTANCE", 0),
                     counts.get("OTHER", 0),
                     counts["antismash_domain"],
                     counts["pfam_domain"],
                     ]
            out.write("\t".join([str(item) for item in stats])+"\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir")
    parser.add_argument("output_tsv")
    args = parser.parse_args()
    main(**vars(args))
