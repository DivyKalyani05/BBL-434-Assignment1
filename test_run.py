"""
Automated test for plasmid redesign tool.
"""

from Bio import SeqIO

seq = str(SeqIO.read("pUC19_modified.fa", "fasta").seq)

# EcoRI must be removed (not listed in design file)
assert "GAATTC" not in seq, "EcoRI site still present"

print("All tests passed.")
