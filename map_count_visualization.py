import pandas as pd
import matplotlib.pyplot as plt

PAF_FILE = "alignment.paf"
cols = [
    "query", "query_len", "query_start", "query_end",
    "strand", "target", "target_len", "target_start", "target_end",
    "matches", "aln_len", "mapq"
]

data = pd.read_csv(PAF_FILE, sep='\t', header=None, names=cols, usecols=range(12))
data = data[data["mapq"] >= 20]

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

data = data[data["target"].isin(chromosomes.keys())]
data["chromosome"] = data["target"].map(chromosomes)
scaffold_counts = data.groupby("chromosome")["query"].nunique()
chrom_order = [str(i) for i in range(1, 17)] + ["X", "Y", "MT"]
scaffold_counts = scaffold_counts.reindex(chrom_order, fill_value=0)
plt.figure(figsize=(10, 5))
scaffold_counts.plot(kind="bar", color="steelblue")
plt.title("Number of scaffolds aligned to each chromosome")
plt.xlabel("Chromosome")
plt.ylabel("Number of scaffolds")
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig("scaffold_distribution_named.png", dpi=300)
plt.close()

print("Plot saved as: scaffold_distribution_named.png")
