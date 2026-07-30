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
import json
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
    parser.add_argument("--contig-map", default=None,
                        help="JSON file mapping old contig IDs to new contig IDs")
    parser.add_argument("--split-gene-ids", action="store_true",
                        help="Split gene_id values on '.' and use first part only (for GTF)")
    parser.add_argument("--write-transcript-counts",
                        help="File to write the number of transcripts per locus (optional)")
    args = parser.parse_args()

    contig_map = {}
    if args.contig_map:
        with open(args.contig_map) as fh:
            contig_map = json.load(fh)
        if not isinstance(contig_map, dict):
            raise SystemExit("--contig-map must be a JSON object mapping contig IDs")

    genes = {}                        # gene_id -> list of 9 columns
    mrnas = {}                        # mrna_id -> list of 9 columns
    cds_by_mrna = defaultdict(list)   # mrna_id -> [cols, ...]
    start_codons_by_mrna = defaultdict(list) # mrna_id -> [cols, ...]
    stop_codons_by_mrna = defaultdict(list)   # mrna_id -> [cols, ...]
    mrnas_by_gene = defaultdict(list) # gene_id -> [mrna_id, ...]
    orphan_mrnas = []                 # mrna_ids with no gene parent

    transcript_counts = defaultdict(int)  # transcript_number -> count of loci with that many transcripts

    # ---- Parse -------------------------------------------------------
    with open(args.input) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            if contig_map:
                old_contig = cols[0]
                if old_contig not in contig_map:
                    raise SystemExit(f"Contig '{old_contig}' not found in mapping")
                cols[0] = contig_map[old_contig]
            ftype = cols[2]
            attrs = parse_attributes(cols[8], gtf=args.gtf)

            if ftype == "gene":
                gid = attrs.get("gene_id", "") if args.gtf else attrs.get("ID", "")
                if gid:
                    if args.split_gene_ids:
                        gid = gid.split(".")[0]
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
                parent = attrs.get("Parent")
                if not parent:
                    parent = attrs.get("gene_id")
                parent = parent.split(",")[0] if parent else None

                if mid:
                    mrna_cols = cols[:]

                    if args.gtf:
                        mrna_cols[2] = "mRNA"
                        mrna_attrs = {"ID": mid}
                        for key, value in attrs.items():
                            if key not in {"gene_id", "transcript_id"}:
                                mrna_attrs[key] = value
                    else:
                        mrna_attrs = attrs

                    if parent:
                        if args.split_gene_ids:
                            parent = parent.split(".")[0]
                        mrna_attrs["Parent"] = parent
                        mrnas_by_gene[parent].append(mid)
                    else:
                        orphan_mrnas.append(mid)

                    mrna_cols[8] = format_gff3_attributes(mrna_attrs)
                    mrnas[mid] = mrna_cols

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

            elif ftype == "start_codon":
                parent = (
                    attrs.get("transcript_id", "")
                    if args.gtf
                    else attrs.get("Parent", "")
                ).split(",")[0]
                if parent:
                    start_codons_by_mrna[parent].append(cols)

            elif ftype == "stop_codon":
                parent = (
                    attrs.get("transcript_id", "")
                    if args.gtf
                    else attrs.get("Parent", "")
                ).split(",")[0]
                if parent:
                    stop_codons_by_mrna[parent].append(cols)

    def adjusted_cds_rows(mrna_id):
        """Return CDS rows in original GFF order with start/stop codons merged into the ends."""
        cds = cds_by_mrna.get(mrna_id, [])
        if not cds:
            return []

        strand = cds[0][6]

        # Preserve original GFF order
        merged = [row[:] for row in cds]

        # Identify transcript first/last CDS without reordering
        if strand == "-":
            first = max(merged, key=lambda c: int(c[4]))  # highest genomic end
            last  = min(merged, key=lambda c: int(c[3]))  # lowest genomic start
        else:
            first = min(merged, key=lambda c: int(c[3]))  # lowest genomic start
            last  = max(merged, key=lambda c: int(c[4]))  # highest genomic end

        # Merge start codon into first CDS
        if start_codons_by_mrna.get(mrna_id):
            start_start = min(int(c[3]) for c in start_codons_by_mrna[mrna_id])
            start_end = max(int(c[4]) for c in start_codons_by_mrna[mrna_id])

            if strand == "-":
                first[4] = str(max(int(first[4]), start_end))
            else:
                first[3] = str(min(int(first[3]), start_start))

        # Merge stop codon into last CDS
        if stop_codons_by_mrna.get(mrna_id):
            stop_start = min(int(c[3]) for c in stop_codons_by_mrna[mrna_id])
            stop_end = max(int(c[4]) for c in stop_codons_by_mrna[mrna_id])

            if strand == "-":
                last[3] = str(min(int(last[3]), stop_start))
            else:
                last[4] = str(max(int(last[4]), stop_end))

        return merged

    def cds_extent(mrna_id):
        """Return (min_start, max_end) across merged CDS rows of an mRNA, or None."""
        cds = adjusted_cds_rows(mrna_id)
        if not cds:
            return None
        starts = [int(c[3]) for c in cds]
        ends   = [int(c[4]) for c in cds]
        return min(starts), max(ends)

    def cds_signature(mrna_id):
        cds = adjusted_cds_rows(mrna_id)
        if not cds:
            return None

        return (
            cds[0][6],  # strand
            tuple((int(c[3]), int(c[4])) for c in cds)
        )

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
        for cds_cols in adjusted_cds_rows(mrna_id):
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

            unique_mrna_count = 0
            for mrna_id in mrna_ids:
                written = write_mrna_block(out, mrna_id, seen_signatures)
                if args.write_transcript_counts and written:
                    unique_mrna_count += 1

            if args.write_transcript_counts and unique_mrna_count > 0:
                transcript_counts[unique_mrna_count] += 1

        # Orphan mRNAs (no gene parent)
        for mrna_id in orphan_mrnas:
            written = write_mrna_block(out, mrna_id, seen_signatures)
            if args.write_transcript_counts and written:
                transcript_counts[1] += 1

    # Write transcript counts per locus if requested
    if args.write_transcript_counts:
        with open(args.write_transcript_counts, "w") as out:
            out.write("Transcripts_per_locus\tCount\n")
            for transcript_count, locus_count in sorted(transcript_counts.items()):
                out.write(f"{transcript_count}\t{locus_count}\n")

if __name__ == "__main__":
    main()
