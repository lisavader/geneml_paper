#!/usr/bin/env python3
"""
Trim gene and mRNA spans to their CDS extent, dropping UTRs.
Only gene, mRNA, and CDS features are written to output.
Duplicate mRNAs (identical CDS coordinates) are skipped.

Usage:
    python trim_to_cds.py input.gff3 output.gff3
    python trim_to_cds.py input.gff3 output.gff3 --accession GCF_000001405.40
"""

import argparse
import re
from collections import defaultdict


def parse_attributes(attr_string, gtf=False):
    attrs = {}
    if gtf:
        for key, value in re.findall(r'(\S+)\s+"([^"]+)"', attr_string):
            attrs[key.strip()] = value.strip()
    else:
        for part in attr_string.strip().rstrip(";").split(";"):
            part = part.strip()
            if "=" in part:
                key, _, value = part.partition("=")
                attrs[key.strip()] = value.strip()
    return attrs


def format_gff3_attributes(attrs):
    return ";".join(f"{key}={value}" for key, value in attrs.items()) + ";"


def main():
    parser = argparse.ArgumentParser(
        description="Trim gene/mRNA spans to CDS extent, removing UTRs."
    )
    parser.add_argument("input",  help="Input GFF3 file")
    parser.add_argument("output", help="Output GFF3 file")
    parser.add_argument("--accession", default=None,
                        help="NCBI assembly accession for header (optional)")
    parser.add_argument("--filter-mRNA", help="String to filter mRNA features (optional)")
    parser.add_argument("--gtf", action="store_true",
                        help="Parse input as GTF instead of GFF3")
    args = parser.parse_args()

    genes = {}                        # gene_id -> list of 9 columns
    mrnas = {}                        # mrna_id -> list of 9 columns
    cds_by_mrna = defaultdict(list)   # mrna_id -> [cols, ...]
    mrnas_by_gene = defaultdict(list) # gene_id -> [mrna_id, ...]
    orphan_mrnas = []                 # mrna_ids with no gene parent

    # ---- Parse -------------------------------------------------------
    with open(args.input) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            ftype = cols[2]
            attrs = parse_attributes(cols[8], gtf=args.gtf)

            if ftype == "gene":
                gid = attrs.get("gene_id", "") if args.gtf else attrs.get("ID", "")
                if gid:
                    if args.gtf:
                        gene_cols = cols[:]
                        gene_attrs = {"ID": gid}
                        for key, value in attrs.items():
                            if key != "gene_id":
                                gene_attrs[key] = value
                        gene_cols[8] = format_gff3_attributes(gene_attrs)
                        genes[gid] = gene_cols
                    else:
                        genes[gid] = cols

            elif ftype == ("transcript" if args.gtf else "mRNA"):
                if args.filter_mRNA and args.filter_mRNA not in cols[8]:
                    continue
                mid = attrs.get("transcript_id", "") if args.gtf else attrs.get("ID", "")
                parent = (
                    attrs.get("gene_id", "")
                    if args.gtf
                    else attrs.get("Parent", "")
                ).split(",")[0]
                if mid:
                    if args.gtf:
                        mrna_cols = cols[:]
                        mrna_cols[2] = "mRNA"
                        mrna_attrs = {"ID": mid}
                        if parent:
                            mrna_attrs["Parent"] = parent
                        for key, value in attrs.items():
                            if key not in {"gene_id", "transcript_id"}:
                                mrna_attrs[key] = value
                        mrna_cols[8] = format_gff3_attributes(mrna_attrs)
                        mrnas[mid] = mrna_cols
                    else:
                        mrnas[mid] = cols
                    if parent:
                        mrnas_by_gene[parent].append(mid)
                    else:
                        orphan_mrnas.append(mid)

            elif ftype == "CDS":
                parent = (
                    attrs.get("transcript_id", "")
                    if args.gtf
                    else attrs.get("Parent", "")
                ).split(",")[0]
                if parent:
                    if args.gtf:
                        cds_cols = cols[:]
                        cds_attrs = {"Parent": parent}
                        for key, value in attrs.items():
                            if key != "Parent":
                                cds_attrs[key] = value
                        cds_cols[8] = format_gff3_attributes(cds_attrs)
                        cds_by_mrna[parent].append(cds_cols)
                    else:
                        cds_by_mrna[parent].append(cols)

    def cds_extent(mrna_id):
        """Return (min_start, max_end) across all CDS of an mRNA, or None."""
        cds = cds_by_mrna.get(mrna_id, [])
        if not cds:
            return None
        starts = [int(c[3]) for c in cds]
        ends   = [int(c[4]) for c in cds]
        return min(starts), max(ends)

    def cds_signature(mrna_id):
        """(start, end, strand) of the CDS extent for duplicate detection."""
        cds = cds_by_mrna.get(mrna_id, [])
        strand = cds[0][6] if cds else "."
        return (min(int(c[3]) for c in cds), max(int(c[4]) for c in cds), strand)

    def write_mrna_block(out, mrna_id, seen_signatures):
        """Write mRNA + CDS lines, trimmed to CDS extent. Returns True if written."""
        extent = cds_extent(mrna_id)
        if extent is None:
            return False
        sig = cds_signature(mrna_id)
        if sig in seen_signatures:
            return False
        seen_signatures.add(sig)

        mrna_cols = mrnas[mrna_id][:]
        mrna_cols[3] = str(extent[0])
        mrna_cols[4] = str(extent[1])
        out.write("\t".join(mrna_cols) + "\n")
        for cds_cols in cds_by_mrna[mrna_id]:
            out.write("\t".join(cds_cols) + "\n")
        return True

    # ---- Write -------------------------------------------------------
    with open(args.output, "w") as out:
        # Header
        out.write("##gff-version 3\n")
        if args.accession:
            out.write(f"#!genome-build-accession NCBI_Assembly: {args.accession}\n")
        out.write("##Note: UTRs removed; features span CDS only\n")

        seen_signatures = set()

        # Genes with mRNA children
        for gene_id, mrna_ids in mrnas_by_gene.items():
            # Collect trimmed extents to update gene span
            extents = [cds_extent(mid) for mid in mrna_ids]
            extents = [e for e in extents if e is not None]
            if not extents:
                continue

            if gene_id in genes:
                gene_cols = genes[gene_id][:]
                gene_cols[3] = str(min(e[0] for e in extents))
                gene_cols[4] = str(max(e[1] for e in extents))
                out.write("\t".join(gene_cols) + "\n")

            for mrna_id in mrna_ids:
                write_mrna_block(out, mrna_id, seen_signatures)

        # Orphan mRNAs (no gene parent)
        for mrna_id in orphan_mrnas:
            write_mrna_block(out, mrna_id, seen_signatures)


if __name__ == "__main__":
    main()
