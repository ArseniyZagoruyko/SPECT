import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit

file_path = 'photon_coordinates(0,0,0)rand.txt'
df = pd.read_csv(file_path, sep='\t', header=None, names=['Detector', 'EventID', 'X', 'Z'])

df_counts = df.groupby(['Detector', 'EventID', 'X', 'Z']).size().reset_index(name='Counts')

def gaussian_2d(xy, x0, y0, sigma_x, sigma_y, amplitude):
    x, y = xy
    return amplitude * np.exp(-(((x - x0) ** 2) / (2 * sigma_x ** 2) + ((y - y0) ** 2) / (2 * sigma_y ** 2)))

detectors = df_counts['Detector'].unique()

for detector in detectors:
    event_ids = df_counts[df_counts['Detector'] == detector]['EventID'].unique()

    for event_id in event_ids:
        subset = df_counts[(df_counts['Detector'] == detector) & (df_counts['EventID'] == event_id)]
        if subset.empty:
            continue
        
        x_data = subset['X'].values
        z_data = subset['Z'].values
        counts = subset['Counts'].values

        if len(x_data) < 5:
            print(f"Not enough data to fit Gaussian for Detector {detector}, EventID {event_id}")
            continue
        
        xy_data = np.vstack((x_data, z_data))
        
        initial_guess = (np.mean(x_data), np.mean(z_data), np.std(x_data), np.std(z_data), np.max(counts))
        
        try:
            popt, _ = curve_fit(gaussian_2d, xy_data, counts, p0=initial_guess)
            x0, z0, sigma_x, sigma_z, amplitude = popt
            fwhm_x = 2 * np.sqrt(2 * np.log(2)) * sigma_x
            fwhm_z = 2 * np.sqrt(2 * np.log(2)) * sigma_z   

            print(f"Detector {detector}, EventID {event_id}:")
            print(f"  Mean X: {x0:.2f}")
            print(f"  Mean Z: {z0:.2f}")
            print(f"  FWHM X: {fwhm_x:.2f}")
            print(f"  FWHM Z: {fwhm_z:.2f}")

        except (RuntimeError, TypeError) as e:
            print(f"Error fitting Gaussian for Detector {detector}, EventID {event_id}: {e}")
            continue

        x = np.arange(min(x_data), max(x_data) + 1, 1)
        z = np.arange(min(z_data), max(z_data) + 1, 1)
        X, Z = np.meshgrid(x, z)
        fitted_counts = gaussian_2d((X, Z), x0, z0, sigma_x, sigma_z, amplitude)
        
        plt.figure(figsize=(8, 6))
        ax = sns.heatmap(fitted_counts, cmap='viridis', cbar=True, 
                         cbar_kws={'label': 'Counts', 'shrink': 0.75})

        xticks = np.arange(0, len(x), 5)
        yticks = np.arange(0, len(z), 5)
        plt.xticks(ticks=xticks, labels=np.arange(min(x_data), max(x_data) + 1, 5), fontsize=24)
        plt.yticks(ticks=yticks, labels=np.arange(min(z_data), max(z_data) + 1, 5), fontsize=24)

        plt.title(f'Detector {detector}, EventID {event_id}', fontsize=24)
        plt.xlabel('X', fontsize=24)
        plt.ylabel('Z', fontsize=24)

        colorbar = ax.collections[0].colorbar
        colorbar.ax.tick_params(labelsize=24)  

        
        ax.tick_params(axis='both', which='major', labelsize=24)
        
        plt.show()