import matplotlib.pyplot as plt

repeat_class={'ARTEFACT':4,
              'DNA':312779,
              'LINE':795816,
              'LTR':643549,
              'RC (Helitron)':6571,
              'Retroposon':21808,
              'SINE':1457631,
              'Unknown':172939,
              'Low_complexity':201468,
              'Satellite':129209, 'Simple_repeat':1087882,
              'rRNA':746,
              'scRNA':53365,
              'snRNA':2414,
              'srpRNA':171,
              'tRNA':7560
              }

sorted_data = dict(sorted(repeat_class.items(), key=lambda item: item[1], reverse=True))
labels = list(sorted_data.keys())
values = list(sorted_data.values())

plt.figure(figsize=(14, 7))
bars = plt.barh(labels, values, color="c")
plt.xscale('log')  # логарифмическая шкала
plt.xlim(1, max(values) * 3)
plt.xlabel('Count (log scale)', fontweight='bold', fontsize=12)
plt.ylabel('Repeat class', fontweight='bold', fontsize=12)
plt.title('Repeat Annotation Summary', fontweight='bold', fontsize=14, pad=15)
plt.grid(axis='x', which='both', linestyle='--', linewidth=0.8)

for bar in bars:
    width = bar.get_width()
    plt.text(width * 1.05, bar.get_y() + bar.get_height() / 2,
             f'{width:,}', va='center', fontsize=10)

plt.tight_layout()
plt.show()
