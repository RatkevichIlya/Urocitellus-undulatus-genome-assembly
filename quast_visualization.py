import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler


rows = ['Number of contigs', 'Total length', 'Largest scaffold', 'N50', 'Average length']
data =  {
    'MaSuRCA(Celera)':[28699, 23593.24609, 15046.24, 17796.0, 8206.06],
    'MaSuRCA(Flye)':[19190, 13042.82694, 5849.48, 8295.8, 6796.68],
    'MAECI-based':[16432, 7247.40478, 4511.75, 5114.3, 4410.54]
}

df = pd.DataFrame(data=data)
df.insert(0, 'Parameter', rows)
dfm = pd.melt(df, id_vars='Parameter', var_name="Assembler", value_name="Length")

sns.set_theme(style='white', font_scale=1.1)
g = sns.catplot(x='Parameter', y='Length', hue='Assembler', data=dfm, kind='bar',
                height=8, aspect=1, palette='bright', legend_out=False)
g.set_xlabels("Parameter", fontweight='demibold', fontsize=16)
g.set_ylabels("Length", fontweight='demibold', fontsize=16, labelpad=15)

plt.title('Assemblies statistics', fontsize=20, fontweight='bold', pad=20)
plt.xticks(rotation=30)
plt.tight_layout()
plt.show()
