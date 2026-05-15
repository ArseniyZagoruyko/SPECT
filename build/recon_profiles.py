import os
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tqdm import tqdm

from recon_utils import (
    discover_pairs,
    iter_events_from_pair,
    find_theta,
    find_vector_P,
    find_intersection_with_XZ_plane,
    find_vector_h,
    find_vector_d,
    rotate_and_intersect,
)

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
MAX_LINES_PER_FILE = None

cell_size = 0.1
x_range = (-50, 50)
z_range = (-50, 50)

x_bins = np.arange(x_range[0] - cell_size / 2, x_range[1] + cell_size, cell_size)
z_bins = np.arange(z_range[0] - cell_size / 2, z_range[1] + cell_size, cell_size)

all_x_coords = []
all_z_coords = []

pairs = discover_pairs(DATA_DIR)
print("=" * 60)
print(f"Найдено пар файлов: {len(pairs)}")
print("=" * 60)

total_events = 0
used_events = 0

for seed, dep_path, coord_path in tqdm(pairs, desc="Файлы", unit="пара"):
    for ev, E_lost, E_remaines, Vx, Vy, Vz, X2, Y2, Z2 in iter_events_from_pair(
        dep_path, coord_path, max_lines=MAX_LINES_PER_FILE
    ):
        total_events += 1
        if E_lost > E_remaines:
            continue

        theta = find_theta(E_remaines, E_lost, smear=False)
        if theta is None:
            continue

        P = find_vector_P(Vx, Vy, Vz, X2, Y2, Z2)
        H = find_intersection_with_XZ_plane(Vx, Vy, Vz, P[0], P[1], P[2])
        if H is None:
            continue
        h = find_vector_h(Vx, Vy, Vz, H[0], 0, H[2])
        d = find_vector_d(h, theta)
        if d is None:
            continue
        V = np.array([Vx, Vy, Vz])
        pts = rotate_and_intersect(V, d, h)
        if pts.size == 0:
            continue
        all_x_coords.extend(pts[:, 0])
        all_z_coords.extend(pts[:, 2])
        used_events += 1

print(f"\nВсего прочитано совпадений: {total_events}")
print(f"Использовано для реконструкции: {used_events}")
print(f"Накоплено точек пересечения: {len(all_x_coords)}")

point_counts = Counter()
for x, z in zip(all_x_coords, all_z_coords):
    if not (np.isfinite(x) and np.isfinite(z)):
        continue
    point_counts[(round(x, 1), round(z, 1))] += 1

filtered_points = {p: c for p, c in point_counts.items() if c > 2}
filtered_x, filtered_z, filtered_counts = [], [], []
for (x, z), c in filtered_points.items():
    filtered_x.append(x)
    filtered_z.append(z)
    filtered_counts.append(c)

x_vals = filtered_x
z_vals = filtered_z
counts = filtered_counts

print(f"Уникальных точек после фильтрации (count > 2): {len(filtered_points)}")


def gaussian(x, A, mu, sigma):
    return A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))


x_bins = np.arange(x_range[0] - cell_size / 2, x_range[1] + cell_size, cell_size)
z_bins = np.arange(z_range[0] - cell_size / 2, z_range[1] + cell_size, cell_size)

x_profile, x_edges = np.histogram(x_vals, bins=x_bins, weights=counts)
z_profile, z_edges = np.histogram(z_vals, bins=z_bins, weights=counts)

x_centers = (x_edges[:-1] + x_edges[1:]) / 2
z_centers = (z_edges[:-1] + z_edges[1:]) / 2

total_w = np.sum(counts) if counts else 0.0
if total_w > 0:
    weighted_mean_x = np.sum(np.array(x_vals) * np.array(counts)) / total_w
    weighted_mean_z = np.sum(np.array(z_vals) * np.array(counts)) / total_w
    std_x = np.std(x_vals)
    std_z = np.std(z_vals)
else:
    weighted_mean_x = weighted_mean_z = 0.0
    std_x = std_z = 1.0

try:
    popt_x, _ = curve_fit(
        gaussian, x_centers, x_profile,
        p0=[np.max(x_profile) if x_profile.size else 1.0, weighted_mean_x, max(std_x, 1e-3)],
    )
except Exception as e:
    print(f"Предупреждение: фит X не сошёлся: {e}")
    popt_x = None

try:
    popt_z, _ = curve_fit(
        gaussian, z_centers, z_profile,
        p0=[np.max(z_profile) if z_profile.size else 1.0, weighted_mean_z, max(std_z, 1e-3)],
    )
except Exception as e:
    print(f"Предупреждение: фит Z не сошёлся: {e}")
    popt_z = None

peak_position_x = x_centers[np.argmax(x_profile)] if x_profile.size else 0.0
peak_position_z = z_centers[np.argmax(z_profile)] if z_profile.size else 0.0

print(f"\nМаксимальная интенсивность по X: {peak_position_x:.2f}")
if popt_x is not None:
    print(f"Гауссовский пик по X: A={popt_x[0]:.2f}, μ={popt_x[1]:.2f}, σ={popt_x[2]:.2f}")
print(f"\nМаксимальная интенсивность по Z: {peak_position_z:.2f}")
if popt_z is not None:
    print(f"Гауссовский пик по Z: A={popt_z[0]:.2f}, μ={popt_z[1]:.2f}, σ={popt_z[2]:.2f}")

plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(x_centers, x_profile, 'bo', label='Data')
if popt_x is not None:
    plt.plot(x_centers, gaussian(x_centers, *popt_x), 'r-', label='Gaussian fit')
plt.axvline(peak_position_x, color='g', linestyle='--', label='Peak position')
plt.xlabel('X', fontsize=14)
plt.ylabel('Counts', fontsize=14)
plt.title('X Profile')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(z_centers, z_profile, 'go', label='Data')
if popt_z is not None:
    plt.plot(z_centers, gaussian(z_centers, *popt_z), 'm-', label='Gaussian fit')
plt.axvline(peak_position_z, color='g', linestyle='--', label='Peak position')
plt.xlabel('Z', fontsize=14)
plt.ylabel('Counts', fontsize=14)
plt.title('Z Profile')
plt.legend()

plt.tight_layout()
plt.show()
