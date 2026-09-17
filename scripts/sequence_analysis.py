from pathlib import Path

mRNA_file = Path("data/reference_sequences/ATP7B_mRNA_NM_000053.4.fasta")
protein_file = Path("data/reference_sequences/ATP7B_protein_NP_000044.2.fasta")


def read_fasta(filename):
    sequence = ""

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if not line.startswith(">"):
                sequence += line

    return sequence


mRNA = read_fasta(mRNA_file)
protein = read_fasta(protein_file)


print("ATP7B mRNA length:", len(mRNA), "nt")
print("ATP7B protein length:", len(protein), "aa")
print("ATP7B protein starts with methionine (M):", protein.startswith("M"))
