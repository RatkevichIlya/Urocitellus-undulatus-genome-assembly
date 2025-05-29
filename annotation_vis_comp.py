import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

features=[
    'Genes',
    'mRNAs',
    'CDSs',
    'exons in CDSs',
    'introns\n in exon/CDSs'
]

undulatus_values=[27921, 27921, 27921, 199793, 171872]
parryii_values= [26655, 36308, 36389, 216202, 193983]

width=0.4
x = np.arange(5)
plt.figure(figsize=(12,6))
plt.bar(x, undulatus_values, width, log=True, color='#1e96fc')
plt.bar(x+0.4, parryii_values, width, log=True, color='#9a031e')
plt.xlabel("Feature", fontsize=14)
plt.xticks(x, features, fontsize=12)
plt.ylabel('Count (log scale)', fontsize=14)
plt.title('Genome Annotation Comparison', fontsize=16, fontweight='bold', pad=20)
plt.grid(axis='y', which='both', linestyle='--', linewidth=0.5)
plt.legend(["${U. undulatus}$", "${U. parryii}$"], fontsize=14)
plt.tight_layout()
plt.show()
