import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


base_path = '/home/zas/CERN/SPECT/build/'
file_count = 56  

dtype = np.dtype([
    ('EventID', 'int32'),  # Event ID (4 bytes)
    ('Detector', 'int32'), # Radiator ID (4 bytes)
    ('X', 'float64'),      # X coordinate (8 bytes)
    ('Z', 'float64')       # Z coordinate (8 bytes)
])


all_data = []


for i in range(file_count):
    file_path = f'{base_path}photon_coordinates_{i}.bin'
    try:
        data = np.fromfile(file_path, dtype=dtype)
        all_data.append(data)
    except FileNotFoundError:
        print(f"File {file_path} not found. Skipping...")
        continue


if not all_data:
    raise ValueError("No binary files were found or read successfully.")


combined_data = np.concatenate(all_data)

#  DataFrame
df = pd.DataFrame(combined_data)

# группировка данных 
df_counts = df.groupby(['Detector', 'EventID', 'X', 'Z']).size().reset_index(name='Counts')


def gaussian_2d(xy, x0, y0, sigma_x, sigma_y, amplitude):
    x, y = xy
    return amplitude * np.exp(-(((x - x0) ** 2) / (2 * sigma_x ** 2) + ((y - y0) ** 2) / (2 * sigma_y ** 2)))


event_ids = df_counts['EventID'].unique()

results = []

for event_id in event_ids:
    detectors_data = []
    
    for detector in df_counts['Detector'].unique():
        subset = df_counts[(df_counts['EventID'] == event_id) & (df_counts['Detector'] == detector)]
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

            print(f"Detector {detector}, EventID {event_id}:")
            print(f"  Mean X: {x0}")
            print(f"  Mean Z: {z0}")
            
            detectors_data.append((x0, z0, sigma_x, sigma_z))

        except (RuntimeError, TypeError) as e:
            print(f"Error fitting Gaussian for Detector {detector}, EventID {event_id}: {e}")
            continue

    if len(detectors_data) == 2:
        x1, z1, sigma_x1, sigma_z1 = detectors_data[0]
        x2, z2, sigma_x2, sigma_z2 = detectors_data[1]
        results.append([event_id, x1, "54", z1, x2, "84", z2])


output_file = 'calculated_coordinates.txt'
with open(output_file, 'w') as f:
    f.write("EventID\tX1\tY1\tZ1\tX2\tY2\tZ2\n")  
    for result in results:
        f.write("\t".join(map(str, result)) + "\n")

print(f"Results saved to {output_file}")