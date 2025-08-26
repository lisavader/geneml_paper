import argparse
import glob
import pandas as pd

from pfam2go import pfam2go

from antismash.common.serialiser import AntismashResults

def get_pfam_counts(record):
    counts = {}
    pfams = record.get_pfam_domains()
    for pfam in pfams:
        try:
            counts[pfam.identifier] += 1
        except KeyError:
            counts[pfam.identifier] = 1
    return counts

def get_go_table(record, file_name):
    pfam_counts = get_pfam_counts(record)
    print(pfam_counts.keys())
    go_df = pfam2go(pfam_counts.keys())
    go_df["count"] = go_df["Pfam accession"].apply(lambda x: pfam_counts[x])
    go_df["record"] = record.id
    go_df["file"] = file_name
    return go_df[["file","record","Pfam accession","GO accession","name","count"]]

def main(results_dir, go_table):
        go_tables = []
        path = ''.join([results_dir,"/*/*.json"])
        for file in glob.glob(path):
            as_results = AntismashResults.from_file(file)
            records = as_results.records
            for record in records:
                go_tables.append(get_go_table(record, as_results.input_file))
        pd.concat(go_tables).to_csv(go_table, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir")
    parser.add_argument("go_table")
    args = parser.parse_args()
    main(**vars(args))
