import matplotlib.pyplot as plt
import seaborn as sns

busco_stats = {
    'Complete (Single-copy)': 9689,
    'Complete (Duplicated)': 194,
    'Fragmented': 1286,
    'Missing': 1108
}
color_list = ['#1f78b4', '#33a02c', '#ff7f00', '#e31a1c']
labels = list(busco_stats.keys())
sizes = list(busco_stats.values())

fig, ax = plt.subplots(figsize=(7, 7))
wedges, texts, autotexts = ax.pie(
    sizes,
    labels=None,
    autopct='%1.1f%%',
    startangle=140,
    colors=color_list,
    textprops=dict(color="black", fontsize=12),
    pctdistance=1.15
)
ax.legend(wedges, labels, title="BUSCO Categories", loc="lower left", bbox_to_anchor=(1, 0.5), fontsize=12)

plt.setp(autotexts, size=15, weight="demibold")
plt.title("BUSCO Assessment", fontsize=20, weight='bold', pad=20, x=0.85, y=1.05)
plt.tight_layout()
plt.show()
