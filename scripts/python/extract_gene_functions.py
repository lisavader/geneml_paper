import argparse
from collections import defaultdict

from goatools.obo_parser import GODag

MGE_BIN_ID = 1
UNKNOWN_BIN_ID = 100


def get_pfam_to_gos(pfam2go_file) -> dict[str, list[str]]:
    pfam_to_gos = defaultdict(list)
    with open(pfam2go_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("!") or not line.strip():
                continue
            line = line.strip()
            pfam_id = line.split(" ")[0].split(":")[1]
            go_id = line.split(" ; ")[1]
            pfam_to_gos[pfam_id].append(go_id)
    return pfam_to_gos


def get_anchor_to_bin(go_bins_file) -> tuple[dict[str, int], dict[str, int]]:
    BP_anchor_to_bin = {}
    MF_anchor_to_bin = {}
    with open(go_bins_file, "r", encoding="utf-8") as f:
        next(f)  # Skip header
        for line in f:
            cols = line.strip().split("\t")
            bin_id = int(cols[0])
            BP_anchors = cols[2].split(",")
            MF_anchors = cols[3].split(",")
            for anchor in BP_anchors:
                BP_anchor_to_bin[anchor] = bin_id
            for anchor in MF_anchors:
                MF_anchor_to_bin[anchor] = bin_id
    return BP_anchor_to_bin, MF_anchor_to_bin


def get_bin_from_go(go_id, BP_anchor_to_bin, MF_anchor_to_bin, go_dag) -> int:
    anchors = []
    gos = {go_id}
    if go_id in go_dag:
        parents = go_dag[go_id].get_all_parents()
        gos.update(parents)
    # First check Biological Process anchors
    for go in gos:
        if go in BP_anchor_to_bin:
            anchors.append(BP_anchor_to_bin[go])
    # If none found, check Molecular Function anchors
    if not anchors:
        for go in gos:
            if go in MF_anchor_to_bin:
                anchors.append(MF_anchor_to_bin[go])
    if anchors:
        return min(anchors)   # Return the most specific bin (lowest number)
    return UNKNOWN_BIN_ID


def get_gene_to_pfams(transcript_annotations) -> dict[str, list[str]]:
    gene_to_pfams: dict[str, list[str]] = {}
    with open(transcript_annotations, "r", encoding="utf-8") as s:
        next(s)  # Skip header
        for line in s:
            fields = line.strip().split("\t")
            gene_id = fields[0]
            pfam_ids = fields[2].split(",") if fields[2] else []
            if gene_id in gene_to_pfams:
                gene_to_pfams[gene_id].extend(pfam_ids)
            else:
                gene_to_pfams[gene_id] = pfam_ids
    return gene_to_pfams


def write_gene_functions(gene_to_pfams, output_file, pfam_to_gos, BP_anchor_to_bin,
                         MF_anchor_to_bin, go_dag, pfams_mge) -> None:
    with open(output_file, "w", encoding="utf-8") as outfile:
        outfile.write("gene_id\tpfam_ids\tgo_ids\tbin_id\n")
        for gene_id, pfam_ids in gene_to_pfams.items():
            go_ids = set()
            bin_ids = set()
            for pfam_id in pfam_ids:
                base_pfam_id = pfam_id.split(".")[0]  # Remove version number if present
                if base_pfam_id in pfams_mge:
                    bin_id = MGE_BIN_ID  # MGE-specific bin
                    bin_ids.add(bin_id)
                    break
                if base_pfam_id in pfam_to_gos:
                    for go_id in pfam_to_gos[base_pfam_id]:
                        go_ids.add(go_id)
                        bin_id = get_bin_from_go(go_id, BP_anchor_to_bin, MF_anchor_to_bin, go_dag)
                        bin_ids.add(bin_id)
            pfam_ids_str = ",".join(sorted(pfam_ids))
            go_ids_str = ",".join(sorted(go_ids)) if go_ids else "NA"
            bin_id = min(bin_ids) if bin_ids else UNKNOWN_BIN_ID  # Choose the most specific bin
            outfile.write(f"{gene_id}\t{pfam_ids_str}\t{go_ids_str}\t{bin_id}\n")


def main(transcript_annotations, output, go_obo, pfam2go, go_bins, pfams_mge):
    # Data loading
    with open(pfams_mge, "r", encoding="utf-8") as f:
        pfams_mge_list = [line.strip() for line in f if line.strip()]
    go_dag = GODag(go_obo)
    pfam_to_gos = get_pfam_to_gos(pfam2go)
    BP_anchor_to_bin, MF_anchor_to_bin = get_anchor_to_bin(go_bins)

    # Processing of pfam hits
    gene_to_pfams = get_gene_to_pfams(transcript_annotations)
    write_gene_functions(gene_to_pfams, output, pfam_to_gos, BP_anchor_to_bin, MF_anchor_to_bin,
                         go_dag, pfams_mge_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--transcript_annotations", type=str,
                        help="Custom transcript annotation file")
    parser.add_argument("--output", type=str,
                        help="Output file with gene functions")
    parser.add_argument("--go-obo", type=str,
                        help="GO OBO file")
    parser.add_argument("--pfam2go", type=str,
                        help="Mapping file of PFAMs to GOs")
    parser.add_argument("--go-bins", type=str,
                        help="Mapping file of GOs to bins")
    parser.add_argument("--pfams-mge", type=str,
                        help="PFAM annotations for MGE genes")
    args = parser.parse_args()
    main(**vars(args))
