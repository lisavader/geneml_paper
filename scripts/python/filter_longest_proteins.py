from Bio import SeqIO
from pathlib import Path
import sys

THRESHOLD = 100

fasta = sys.argv[1]
path = Path(fasta).with_suffix("")

records = list(SeqIO.parse(fasta, "fasta"))
long  = [r for r in records if len(r.seq) >=  THRESHOLD]

SeqIO.write(long,  f"{path}_long.faa",  "fasta")
