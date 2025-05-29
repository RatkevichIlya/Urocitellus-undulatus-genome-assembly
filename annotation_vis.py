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

isoform_values=[27921, 101849, 101849, 1137621, 1035772]
noiso_values=[27921, 27921, 27921, 199793, 171872]

width=0.4
x = np.arange(5)
plt.figure(figsize=(12,6))
plt.bar(x, isoform_values, width, log=True, color='#ff9505')
plt.bar(x+0.4, noiso_values, width, log=True, color='#1b998b')
plt.xlabel("Feature", fontsize=14)
plt.xticks(x, features, fontsize=12)
plt.ylabel('Count (log scale)', fontsize=14)
plt.title('Genome Annotation Statistics', fontsize=16, fontweight='bold', pad=20)
plt.grid(axis='y', which='both', linestyle='--', linewidth=0.5)
plt.legend(["Isoforms included", "Isoforms filtered"], fontsize=14)
plt.tight_layout()
plt.show()
