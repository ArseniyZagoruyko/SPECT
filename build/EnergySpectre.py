import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("/home/zas/CERN/SPECT/build/deposits(1detector).txt")

event_id = data[:, 0]
energy_scatterer = data[:, 1]
energy_absorber = data[:, 2]

percent_scatterer = (np.sum(energy_scatterer <= 5) / len(energy_scatterer)) * 100
percent_absorber = (np.sum(energy_absorber <= 5) / len(energy_absorber)) * 100

bin_edges = np.arange(0, 200, 2)

plt.figure(figsize=(12, 6))

# Гистограмма для рассеивателя
plt.subplot(1, 2, 1)
plt.hist(energy_scatterer, bins=bin_edges, color='blue', alpha=0.7, edgecolor='black')
plt.title("Спектр потерь энергии в рассеивателе", fontsize=24)
plt.xlabel("Потери энергии (кэВ)",fontsize=24)
plt.ylabel("Количество событий",fontsize=24)
plt.grid(True, which="both", ls="--")  # Сетка 
plt.yscale('log')  # Логарифмическая шкала по Y
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

plt.axvline(x=5, color='red', linestyle='--', linewidth=1.5, label=f'{percent_scatterer:.2f}% ≤5 кэВ')
plt.text(5.2, 1e4, f'{percent_scatterer:.2f}%', color='red', fontsize=24)

# Гистограмма для поглотителя
plt.subplot(1, 2, 2)
plt.hist(energy_absorber, bins=bin_edges, color='green', alpha=0.7, edgecolor='black')
plt.title("Спектр потерь энергии в поглотителе",fontsize=24)
plt.xlabel("Потери энергии (кэВ)",fontsize=24)
plt.ylabel("Количество событий",fontsize=24)
plt.grid(True, which="both", ls="--")  
plt.yscale('log')  
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

plt.tight_layout()
plt.show()