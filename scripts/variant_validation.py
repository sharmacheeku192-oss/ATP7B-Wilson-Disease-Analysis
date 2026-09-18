from pathlib import Path

variant_file = Path("data/variants/variants.csv")
sequence_file = Path("data/reference_sequences/ATP7B_mRNA_NM_000053.4.fasta")

with open(variant_file, "r") as file:
    lines = file.readlines()

print(lines)
data = lines[1].strip().split(",")
print(data)
dna_change = data[1]

position = int(dna_change[2:6])
reference_base = dna_change[6]
alternate_base = dna_change[8]

print("DNA change:", dna_change)
print("Coding position:", position)
print("Reference base:", reference_base)
print("Alternate base:", alternate_base)
with open(sequence_file, "r") as file:
    fasta_lines = file.readlines()
codon_start = position - ((position - 1) % 3)

print("Codon starts at CDS position:", codon_start)
sequence = "".join(line.strip() for line in fasta_lines[1:])
cds_start = 114

mrna_position = cds_start + position - 1

print("mRNA position:", mrna_position)

reference_nucleotide = sequence[mrna_position - 1]

print("Reference nucleotide in FASTA:", reference_nucleotide)

if reference_nucleotide == reference_base:
    print("Reference base validated: True")
else:
    print("Reference base validated: False")
cds_start = 114
cds_end = 4511

cds_sequence = sequence[cds_start - 1:cds_end]

codon_start = position - ((position - 1) % 3)

codon_start_index = codon_start - 1

reference_codon = cds_sequence[codon_start_index:codon_start_index + 3]

print("CDS length:", len(cds_sequence))
print("Reference codon:", reference_codon)
mutant_codon = (
    reference_codon[:codon_position_in_codon]
    + alternate_base
    + reference_codon[codon_position_in_codon + 1:]
)
