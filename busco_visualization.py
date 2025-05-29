import re
import matplotlib.pyplot as plt

busco_files = {
    'MaSuRCA(Celera)': '/home/ilya/Desktop/DIPLOMA/Ground_squirrel/busco/miniprot/busco_masurca_celera/run_mammalia_odb12/short_summary.txt',
    'MaSuRCA(Flye)': '/home/ilya/Desktop/DIPLOMA/Ground_squirrel/busco/miniprot/busco_masurca_flye/run_mammalia_odb12/short_summary.txt',
    'MAECI': '/home/ilya/Desktop/DIPLOMA/Ground_squirrel/busco/miniprot/busco_maeci/run_mammalia_odb12/short_summary.txt'
}

busco_data = {}
for filename, filepath in busco_files.items():
    with open(filepath) as file:
        for line in file:
            if line.startswith('\tC:'):
                match = re.search(
                    r'C:(\d+\.\d+)%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%',
                    line
                )
                if match:
                    busco_data[filename] = {
                        'Complete': float(match.group(1)),
                        'Fragmented': float(match.group(4)),
                        'Missing': float(match.group(5))
                    }

complete, fragmented, missing =[], [], []
assemblies = list(busco_data.keys())[::-1]
for k in assemblies:
    complete.append(busco_data[k]['Complete'])
    fragmented.append(busco_data[k]['Fragmented'])
    missing.append(busco_data[k]['Missing'])

font = {'family' : 'normal',
        'size'   : 16}
plt.rc('font', **font)
plt.style.use("seaborn-v0_8-muted")
fig, ax = plt.subplots(figsize=(14, 6))
y = range(len(assemblies))
ax.barh(y, complete, label='Complete', color='green')
ax.barh(y, fragmented, left=complete, label='Fragmented', color='orange')
ax.barh(y, missing, left=[complete[i] + fragmented[i] for i in y], label='Missing', color='red')
ax.set_yticks(y)
ax.set_facecolor("#eeefff")
ax.set_yticklabels(['MAECI-based approach', 'MaSuRCA(Flye)', 'MaSuRCA(Celera)'])
ax.set_xlabel('BUSCO %', fontsize=16, labelpad=15)
ax.set_title('BUSCO Summary per Assembly', fontsize=20, fontweight='bold', pad=20)
ax.legend(bbox_to_anchor=(1,1))
plt.tight_layout()
plt.show()
