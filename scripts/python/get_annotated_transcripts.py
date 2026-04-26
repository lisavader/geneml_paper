import argparse
from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum
import gffutils
import re


class DomainAnnotation(Enum):
    COMPLETE = "Complete"
    PARTIAL = "Partial"
    NONE = "None"


class DomainChange(Enum):
    UNCHANGED = "Unchanged"
    GAIN = "Gain"
    LOSS = "Loss"
    EXTENDED = "Extended"   # Partial --> Complete
    TRUNCATED = "Truncated" # Complete --> Partial


@dataclass
class PFAMHit:
    pfam_id: str
    complete: bool


@dataclass
class TranscriptAnnotation:
    transcript_id: str
    gene_id: str
    pfam_hits: tuple[PFAMHit]
    domain_annotation: DomainAnnotation
    domain_changes: set[DomainChange] = field(default_factory=set)

    def add_domain_change(self, domain_change: DomainChange):
        self.domain_changes.add(domain_change)


def get_gene_to_transcripts(gff_file) -> dict[str, list[str]]:
    gff_db = gffutils.create_db(gff_file, dbfn=":memory:", merge_strategy="merge")
    gene_to_transcripts = defaultdict(list)

    for gene in gff_db.features_of_type("gene"):
        transcript_ids = [
            child.id
            for child in gff_db.children(gene, level=1)
            if child.featuretype.lower() in {"mrna", "transcript"}
        ]

        if transcript_ids:
            gene_to_transcripts[gene.id] = transcript_ids
        else:
            gene_to_transcripts[gene.id] = [gene.id]

    return gene_to_transcripts


def get_pfam_hits(hmmscan_file) -> dict[str, list[PFAMHit]]:
    pfam_hits = defaultdict(list)
    with open(hmmscan_file, "r", encoding="utf-8") as infile:
        for line in infile:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split()
            pfam_id = fields[1]
            transcript_id = fields[3]
            hmm_length = int(fields[2])
            aligned_length = int(fields[17]) - int(fields[16]) + 1
            complete = aligned_length >= 0.8 * hmm_length
            hit = PFAMHit(pfam_id=pfam_id, complete=complete)
            pfam_hits[transcript_id].append(hit)
    return pfam_hits


def add_domain_changes(transcripts: list[TranscriptAnnotation],
                       reference: TranscriptAnnotation) -> None:
    ref_hits = {hit.pfam_id: hit.complete for hit in reference.pfam_hits}
    ref_pfams = set(ref_hits.keys())
    for transcript in transcripts:
        transcript_hits = {hit.pfam_id: hit.complete for hit in transcript.pfam_hits}
        transcript_pfams = set(transcript_hits.keys())
        added_pfams = transcript_pfams - ref_pfams
        lost_pfams = ref_pfams - transcript_pfams
        common_pfams = ref_pfams & transcript_pfams
        if added_pfams:
            transcript.add_domain_change(DomainChange.GAIN)
        if lost_pfams:
            transcript.add_domain_change(DomainChange.LOSS)
        for pfam_id in common_pfams:
            if ref_hits[pfam_id] and not transcript_hits[pfam_id]:
                transcript.add_domain_change(DomainChange.TRUNCATED)
            elif not ref_hits[pfam_id] and transcript_hits[pfam_id]:
                transcript.add_domain_change(DomainChange.EXTENDED)
        if not transcript.domain_changes:
            transcript.add_domain_change(DomainChange.UNCHANGED)


def get_annotated_transcripts(gene_to_transcripts, pfam_hits) -> list[TranscriptAnnotation]:
    transcripts: list[TranscriptAnnotation] = []
    for gene_id, transcript_ids in gene_to_transcripts.items():
        gene_transcripts = []
        for transcript_id in transcript_ids:
            # Hmmscan sometimes splits IDs by "-", god knows why
            try:
                hits = pfam_hits[transcript_id]
            except KeyError:
                transcript_id_shorter = transcript_id.split("-")[-1]
                hits = pfam_hits.get(transcript_id_shorter, [])
            domain_annotation = DomainAnnotation.NONE

            if hits:
                if any(hit.complete for hit in hits):
                    domain_annotation = DomainAnnotation.COMPLETE
                else:
                    domain_annotation = DomainAnnotation.PARTIAL

            transcript_annotation = TranscriptAnnotation(
                transcript_id=transcript_id,
                gene_id=gene_id,
                pfam_hits=tuple(hits),
                domain_annotation=domain_annotation
            )
            gene_transcripts.append(transcript_annotation)
        if len(gene_transcripts) > 1 and re.fullmatch(r"GML\d+_mRNA1", gene_transcripts[0].transcript_id):
            add_domain_changes(gene_transcripts[1:], gene_transcripts[0])
        transcripts.extend(gene_transcripts)
    return transcripts


def write_transcripts(transcripts, output_file) -> None:
    with open(output_file, "w", encoding="utf-8") as outfile:
        outfile.write("gene_id\ttranscript_id\tpfam_ids\tdomain_annotation\tdomain_changes\n")
        for transcript in transcripts:
            pfam_ids_str = ",".join(sorted([hit.pfam_id for hit in transcript.pfam_hits]))
            domain_changes_str = ",".join([change.value for change in transcript.domain_changes])
            domain_changes = domain_changes_str if domain_changes_str else None
            outfile.write(f"{transcript.gene_id}\t{transcript.transcript_id}\t{pfam_ids_str}\t{transcript.domain_annotation.value}\t{domain_changes}\n")


def main(hmmscan, gff, output):
    gene_to_transcripts = get_gene_to_transcripts(gff)
    pfam_hits = get_pfam_hits(hmmscan)
    transcripts = get_annotated_transcripts(gene_to_transcripts, pfam_hits)
    write_transcripts(transcripts, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--hmmscan", type=str,
                        help="Hmmscan table mapping to PFAM (domtblout format)")
    parser.add_argument("--gff", type=str,
                        help="GFF file with gene annotations")
    parser.add_argument("--output", type=str,
                        help="Output file with extracted gene functions")
    args = parser.parse_args()
    main(**vars(args))
