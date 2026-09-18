from pathlib import Path

variant_file = Path("data/variants/variants.csv")

with open(variant_file, "r") as file:
    lines = file.readlines()

print(lines)
data = lines[1].strip().split(",")
print(data)
