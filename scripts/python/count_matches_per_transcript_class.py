import argparse
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
    parser = argparse.ArgumentParser()
    parser.add_argument("gff3",  help="Reference GFF3 file with transcript class annotations")
    parser.add_argument("tracking", help="Gffcompare tracking file comparing predicted and reference transcripts")
    parser.add_argument("output", help="Output TSV file with counts of predicted transcripts per reference transcript class")
    args = parser.parse_args()

    mrnas_by_class = defaultdict(list)  # class -> [mrna_id, ...]
    matches_by_class = defaultdict(int)  # class -> count of correctly predicted transcripts

    with open(args.gff3) as s:
        for line in s:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                continue

            if parts[2] == "mRNA":
                attrs = parse_attributes(parts[8])
                transcript_id = attrs.get("ID")
                try:
                    transcript_class = attrs.get("TranscriptVariant")
                except KeyError:
                    raise KeyError(f"Missing TranscriptVariant attribute for transcript {transcript_id}")
                mrnas_by_class[transcript_class].append(transcript_id)

    with open(args.tracking) as s:
        for line in s:
            if not line.strip():
                continue
            parts = line.split("\t")

            ref_id = parts[2]
            if ref_id == "-":
                continue
            transcript_id = ref_id.split('|')[1]
            match_type = parts[3]

            found = False
            for transcript_class, ids in mrnas_by_class.items():
                if transcript_id in ids:
                    found = True
                    if match_type == "=":
                        matches_by_class[transcript_class] += 1
                        break
            if not found:
                raise ValueError(f"Reference transcript {transcript_id} not found in GFF3 file")

    with open(args.output, "w") as out:
        out.write("TranscriptVariant\tMatches\tTotal\n")
        for transcript_class, count in matches_by_class.items():
            total = len(mrnas_by_class[transcript_class])
            out.write(f"{transcript_class}\t{count}\t{total}\n")


if __name__ == "__main__":
    main()
