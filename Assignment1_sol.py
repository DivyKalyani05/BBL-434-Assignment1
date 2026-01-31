"""
Plasmid Redesign Tool

This script takes an already functional plasmid (FASTA),
detects and preserves its origin of replication (ORI),
and edits restriction sites based on a user-provided
design specification file.

The design file is treated as the FINAL specification
of what features the plasmid should contain.
"""

from Bio import SeqIO
from Bio.Seq import Seq

# ============================================================
# Restriction enzyme recognition sequences (from markers.tab)
# ============================================================
RESTRICTION_SITES = {
    "EcoRI":  "GAATTC",
    "BamHI":  "GGATCC",
    "HindIII":"AAGCTT",
    "PstI":   "CTGCAG",
    "KpnI":   "GGTACC",
    "SacI":   "GAGCTC",
    "SalI":   "GTCGAC",
    "XbaI":   "TCTAGA",
    "NotI":   "GCGGCCGC",
    "SmaI":   "CCCGGG"
}

# ============================================================
# Parse design file
# ============================================================
def parse_design_file(design_file):
    """
    Reads the design file and extracts:
    - Allowed restriction sites
    - Required markers (if any)
    - Required ORI type (declarative)

    The design file specifies what the FINAL plasmid
    should contain.
    """

    allowed_sites = set()
    required_markers = set()
    required_ori = None

    with open(design_file) as f:
        for line in f:
            if not line.strip():
                continue

            key, value = [x.strip() for x in line.split(",")]

            if key.endswith("_site"):
                allowed_sites.add(value)

            elif key.endswith("_gene"):
                required_markers.add(value)

            elif key.startswith("ori"):
                required_ori = value

    return allowed_sites, required_markers, required_ori

# ============================================================
# ORI detection (AT-rich / low-GC sliding window)
# ============================================================
def detect_ori(sequence, window=100):
    """
    Detects ORI as the most AT-rich (lowest GC%) region.

    This is a GC-content method, NOT GC skew.
    Appropriate for plasmids like pUC19.
    """

    seq = str(sequence).upper()
    min_gc = 1.0
    ori_start = 0

    for i in range(len(seq) - window):
        fragment = seq[i:i + window]
        gc = (fragment.count("G") + fragment.count("C")) / window

        if gc < min_gc:
            min_gc = gc
            ori_start = i

    return ori_start, ori_start + window

# ============================================================
# Remove disallowed restriction sites
# ============================================================
def remove_disallowed_sites(sequence, allowed_sites, protected_regions):
    """
    Removes restriction sites that are NOT listed
    in the design file. Protected regions (ORI)
    are never modified.
    """

    seq = str(sequence)

    for enzyme, site in RESTRICTION_SITES.items():

        if enzyme in allowed_sites:
            continue

        search_pos = 0
        while True:
            idx = seq.find(site, search_pos)
            if idx == -1:
                break

            overlaps_protected = False
            for start, end in protected_regions:
                if idx < end and idx + len(site) > start:
                    overlaps_protected = True
                    break

            if overlaps_protected:
                search_pos = idx + 1
                continue

            seq = seq[:idx] + seq[idx + len(site):]

    return Seq(seq)

# ============================================================
# Main workflow
# ============================================================
def build_plasmid(fasta_file, design_file, output_file):
    """
    Main pipeline for plasmid redesign.
    """

    record = SeqIO.read(fasta_file, "fasta")
    sequence = record.seq

    allowed_sites, required_markers, required_ori = parse_design_file(design_file)

    ori_start, ori_end = detect_ori(sequence)
    protected_regions = [(ori_start, ori_end)]

    modified_seq = remove_disallowed_sites(
        sequence,
        allowed_sites,
        protected_regions
    )

    record.seq = modified_seq
    record.description = (
        f"Modified plasmid | ORI:{ori_start}-{ori_end} | "
        f"Allowed sites: {','.join(sorted(allowed_sites))}"
    )

    SeqIO.write(record, output_file, "fasta")

# ============================================================
# Run
# ============================================================
if __name__ == "__main__":
    build_plasmid(
        "pUC19.fa",
        "Design_pUC19.txt",
        "pUC19_modified.fa"
    )
