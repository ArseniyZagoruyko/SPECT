import numpy as np
import matplotlib.pyplot as plt


def f(theta, E0):
    alpha = E0 / 511
    dE = 4 + 0.29 * np.sqrt(E0)
    return (180.0 / np.pi) * dE * (1 + alpha * (1 - np.cos(theta * np.pi / 180.0)))**2 / \
           (E0 * alpha * np.sin(theta * np.pi / 180.0))

E0 = 140  
theta_values = np.linspace(0.1, 180, 500) 


f_values = f(theta_values, E0)


plt.figure(figsize=(10, 6))
plt.plot(theta_values, f_values, label=f'E0 = {E0}')
plt.title('140 кэВ')
plt.xlabel('угол рассеяния θ')
plt.ylabel('угловая неопределенность')
plt.ylim(0,25)
plt.grid(True)
plt.legend()
plt.show()