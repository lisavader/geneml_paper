from Bio import SeqIO
from pathlib import Path
import sys

fasta = sys.argv[1]
path = Path(fasta).with_suffix("")

records = list(SeqIO.parse(fasta, "fasta"))

with open(f"{path}_lengths.tsv", "w") as f:
    f.write("protein_id\tlength\n")
    for i, r in enumerate(records):
        if i % 100 == 0:    # sample every 100th protein
            f.write(f"{r.id}\t{len(r.seq)}\n")
