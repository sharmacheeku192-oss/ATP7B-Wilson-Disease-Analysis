from pathlib import Path

# File locations
variant_file = Path("data/variants/variants.csv")
sequence_file = Path("data/reference_sequences/ATP7B_mRNA_NM_000053.4.fasta")

# ATP7B reference CDS coordinates
cds_start = 114
cds_end = 4511

# Genetic code
genetic_code = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",

    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",

    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

# Read variant data
with open(variant_file, "r") as file:
    lines = file.readlines()

data = lines[1].strip().split(",")

variant = data[0]
dna_change = data[1]
protein_change = data[2]

print("Variant:", variant)
print("DNA change:", dna_change)
print("Protein change:", protein_change)

# Parse DNA change
position = int(dna_change[2:6])
reference_base = dna_change[6]
alternate_base = dna_change[8]

print("Coding position:", position)
print("Reference base:", reference_base)
print("Alternate base:", alternate_base)

# Read FASTA sequence
with open(sequence_file, "r") as file:
    fasta_lines = file.readlines()

sequence = "".join(line.strip() for line in fasta_lines[1:])

# Extract CDS
cds_sequence = sequence[cds_start - 1:cds_end]

print("CDS length:", len(cds_sequence))

# Convert coding position to mRNA position
mrna_position = cds_start + position - 1

print("mRNA position:", mrna_position)

# Validate reference nucleotide
reference_nucleotide = sequence[mrna_position - 1]

print("Reference nucleotide in FASTA:", reference_nucleotide)

if reference_nucleotide == reference_base:
    print("Reference base validated: True")
else:
    print("Reference base validated: False")

# Find codon
codon_start = position - ((position - 1) % 3)

print("Codon starts at CDS position:", codon_start)

codon_start_index = codon_start - 1

reference_codon = cds_sequence[codon_start_index:codon_start_index + 3]

print("Reference codon:", reference_codon)

# Position within codon
codon_position = (position - 1) % 3

print("Position within codon:", codon_position + 1)

# Create mutated codon
mutated_codon = list(reference_codon)
mutated_codon[codon_position] = alternate_base
mutated_codon = "".join(mutated_codon)

print("Mutated codon:", mutated_codon)

# Translate codons
reference_amino_acid = genetic_code.get(reference_codon, "?")
mutated_amino_acid = genetic_code.get(mutated_codon, "?")

print("Reference amino acid:", reference_amino_acid)
print("Mutated amino acid:", mutated_amino_acid)

# Final interpretation
if reference_amino_acid != mutated_amino_acid:
    print("Amino acid change: Missense")
else:
    print("Amino acid change: Synonymous")
    
