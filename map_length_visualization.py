import pandas as pd
import matplotlib.pyplot as plt

PAF_FILE = "alignment.paf" # input PAF
cols = [
    "query", "query_len", "query_start", "query_end",
    "strand", "target", "target_len", "target_start", "target_end",
    "matches", "aln_len", "mapq"
]
df = pd.read_csv(PAF_FILE, sep='\t', header=None, names=cols, usecols=range(12))
df = df[df["mapq"] >= 20] # alignments quality filtering

chromosomes = {
    "CM099860.1": "1",
    "CM099861.1": "2",
    "CM099862.1": "3",
    "CM099863.1": "4",
    "CM099864.1": "5",
    "CM099865.1": "6",
    "CM099866.1": "7",
    "CM099867.1": "8",
    "CM099868.1": "9",
    "CM099869.1": "10",
    "CM099870.1": "11",
    "CM099871.1": "12",
    "CM099872.1": "13",
    "CM099873.1": "14",
    "CM099874.1": "15",
    "CM099875.1": "16",
    "CM099876.1": "X",
    "CM099877.1": "Y",
    "CM099878.1": "MT"
}

df = df[df["target"].isin(chromosomes.keys())] # chromosomes filtration and renaming
df["chromosome"] = df["target"].map(chromosomes)
df["alignment_length"] = df["query_end"] - df["query_start"] # alignment length calculation
alignment_sums = df.groupby("chromosome")["alignment_length"].sum() # group alignments length per chromosome
chrom_order = [str(i) for i in range(1, 17)] + ["X", "Y", "MT"] # chromosome sorting
alignment_sums = alignment_sums.reindex(chrom_order, fill_value=0)

plt.figure(figsize=(10, 5))
alignment_sums.plot(kind="bar", color="darkorange")
plt.title("Total alignment length per chromosome")
plt.xlabel("Chromosome")
plt.ylabel("Total alignment length (bp)")
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig("scaffold_total_length_per_chromosome.png", dpi=300)
plt.close()

print("Result saved as: scaffold_total_length_per_chromosome.png")
