from pathlib import Path

variant_file = Path("data/variants/variants.csv")

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
