import matplotlib.pyplot as plt
import numpy as np

# Данные
x = [0.01, 0.03, 0.05, 0.07, 0.10]  


sigma_without_cut = [9.0,9.8,10.3,10.7,12.2]
sigma_with_cut = [7.6, 8.4, 9.1, 9.5, 11.2]
# sigma_without_cut_without_noise = [1.61,1.78,1.87,2.04,2.6]

yerr1 = (sigma_without_cut)/(np.sqrt(2*1000))
yerr2 = (sigma_with_cut)/(np.sqrt(2*1000))
# yerr3 = (sigma_without_cut_without_noise)/(np.sqrt(2*1000))



plt.figure(figsize=(10, 6))
plt.errorbar(x, sigma_without_cut,yerr=yerr1, marker='o',capsize=4, label='Без отсечки')
plt.errorbar(x, sigma_with_cut,yerr=yerr2, marker='s',capsize=4, label='С отсечкой 5 кэВ')
# plt.errorbar(x, sigma_without_cut_without_noise,yerr=yerr2, marker='v',capsize=4, label='Без отсечки')


plt.xlabel('Значение параметра a',fontsize=24)
plt.ylabel('ПШПВ (мм)',fontsize=24)
plt.grid(True)
plt.legend(fontsize=24)
plt.xticks(x)  
plt.tight_layout()
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

plt.show()