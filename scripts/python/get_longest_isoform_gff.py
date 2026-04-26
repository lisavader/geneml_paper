#!/usr/bin/env python3
"""
Extract the longest isoform per gene from a GFF3 file.
Only gene, mRNA, and CDS features are considered and written to output.
Longest isoform is defined as the greatest total CDS length.

Usage:
    python extract_longest_isoform.py input.gff3 output.gff3
"""

import argparse
import shutil
from collections import defaultdict


def parse_attributes(attr_string):
    attrs = {}
    for part in attr_string.strip().rstrip(";").split(";"):
        part = part.strip()
        if "=" in part:
            key, _, value = part.partition("=")
            attrs[key.strip()] = value.strip()
    return attrs


def main():
    parser = argparse.ArgumentParser(
        description="Extract the longest isoform per gene from a GFF3 file."
    )
    parser.add_argument("input",  help="Input GFF3 file")
    parser.add_argument("output", help="Output GFF3 file")
    args = parser.parse_args()

    headers = []
    genes = {}                        # gene_id -> line
    mrnas = {}                        # mrna_id -> line
    cds_by_mrna = defaultdict(list)   # mrna_id -> [cds_line, ...]
    mrnas_by_gene = defaultdict(list) # gene_id -> [mrna_id, ...]

    with open(args.input) as fh:
        for line in fh:
            if line.startswith("#"):
                headers.append(line)
                continue
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                continue
            ftype = parts[2]
            attrs = parse_attributes(parts[8])

            if ftype == "gene":
                gid = attrs.get("ID", "")
                if gid:
                    genes[gid] = line

            elif ftype == "mRNA":
                mid    = attrs.get("ID", "")
                parent = attrs.get("Parent", "").split(",")[0]
                if mid and parent:
                    mrnas[mid] = line
                    mrnas_by_gene[parent].append(mid)

            elif ftype == "CDS":
                parent = attrs.get("Parent", "").split(",")[0]
                if parent:
                    cds_by_mrna[parent].append(line)

    #If mRNA features are missing there is only a single isoform per gene,
    #so output equals input
    if not mrnas:
        print(f"Warning: {args.input} has no mRNA features, returning input with no changes")
        shutil.copy(args.input, args.output)
        exit(0)

    def cds_length(mrna_id):
        return sum(
            int(c.split("\t")[4]) - int(c.split("\t")[3]) + 1
            for c in cds_by_mrna.get(mrna_id, [])
        )

    with open(args.output, "w") as out:
        for h in headers:
            out.write(h)

        for gene_id, mrna_ids in mrnas_by_gene.items():
            best = max(mrna_ids, key=cds_length)

            if gene_id in genes:
                out.write(genes[gene_id])
            out.write(mrnas[best])
            for cds_line in cds_by_mrna.get(best, []):
                out.write(cds_line)


if __name__ == "__main__":
    main()
