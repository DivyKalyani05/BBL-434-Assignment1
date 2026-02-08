import sys
# fasta reader
def read_fasta(file):
    seq = ""
    with open(file) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip().upper()
    return seq


# GC Skew based ORI Finder

def find_ori(sequence):
    skew = []
    s = 0
    for base in sequence:
        if base == 'G':
            s += 1
        elif base == 'C':
            s -= 1
        skew.append(s)
    return skew.index(min(skew))


# Read Design File

def read_design(file):
    enzymes = []
    markers = []
    with open(file) as f:
        for line in f:
            if not line.strip():
                continue
            part, name = [x.strip() for x in line.split(",")]
            if "site" in part.lower():
                enzymes.append(name)
            else:
                markers.append(name)
    return enzymes, markers

# Read markers.tab database

def read_marker_db(file):
    db = {}
    with open(file) as f:
        for line in f:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Skip table headers or separator lines
            if "|" in line or "---" in line:
                continue

            # Only process valid CSV lines
            parts = line.split(",")
            if len(parts) != 2:
                continue

            part, name = [x.strip() for x in parts]
            db[name] = part

    return db

# Restriction enzyme recognition sequences

RE_SITE = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "PstI": "CTGCAG",
    "SphI": "GCATGC",
    "SalI": "GTCGAC",
    "XbaI": "TCTAGA",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "SmaI": "CCCGGG",
    "NotI": "GCGGCCGC"
}

MARKER_SEQ = {
    "Ampicillin": "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTG",
    "Kanamycin": "ATGAGCCATATTCAACGGGAAACGTCTTGCTCGAGGCC",
    "Chloramphenicol": "ATGGAGAAAAAAATCACTGGATATACCACCGTTGATATATCC",
    "Blue_White_Selection": "ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGAC"
}

# Main

if len(sys.argv) != 4:
    print("Usage: python universal_plasmid_maker.py Input.fa Design.txt markers.tab")
    sys.exit(1)

genome_file = sys.argv[1]
design_file = sys.argv[2]
marker_file = sys.argv[3]

genome = read_fasta(genome_file)

# Find ORI
ori_index = find_ori(genome)
ori_seq = genome[ori_index-500 : ori_index+500]

# Read design file
enzymes, markers = read_design(design_file)

# Load marker database
marker_db = read_marker_db(marker_file)

# Build plasmid
plasmid = ""

# Add ORI backbone
plasmid += ori_seq

# Add restriction sites (MCS region)
for enzyme in enzymes:
    if enzyme in RE_SITE:
        plasmid += RE_SITE[enzyme]
    else:
        print(f"Warning: {enzyme} not found in enzyme database")

# Add marker genes
for marker in markers:
    if marker in MARKER_SEQ:
        plasmid += MARKER_SEQ[marker]
    else:
        print(f"Warning: {marker} not found in marker sequence database")

# Write Output
with open("Output.Fa", "w") as f:
    f.write(">Universal_Plasmid\n")
    for i in range(0, len(plasmid), 70):
        f.write(plasmid[i:i+70] + "\n")

print("Plasmid successfully generated → Output.Fa")
print("Final plasmid length:", len(plasmid))
