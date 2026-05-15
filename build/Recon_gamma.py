import os
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tqdm import tqdm
import imageio

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

# Лимит строк на один файл (None = без ограничений). Каждый файл моделирования
# содержит ~неcколько тысяч событий, всего ~1900 пар файлов.
MAX_LINES_PER_FILE = None

cell_size = 0.1
x_range = (-50, 50)
z_range = (-50, 50)

x_bins = np.arange(x_range[0] - cell_size / 2, x_range[1] + cell_size, cell_size)
z_bins = np.arange(z_range[0] - cell_size / 2, z_range[1] + cell_size, cell_size)

all_x_coords = []
all_z_coords = []

pairs = discover_pairs(DATA_DIR)
total_events = 0
used_events = 0

print("=" * 60)
print(f"Найдено пар файлов: {len(pairs)}")
print("=" * 60)

# Обходим пары файлов в внешнем прогресс-баре, события внутри файла — без баром,
# чтобы tqdm не мигал на каждом файле.
for seed, dep_path, coord_path in tqdm(pairs, desc="Файлы", unit="пара"):
    for ev, E_lost, E_remaines, Vx, Vy, Vz, X2, Y2, Z2 in iter_events_from_pair(
        dep_path, coord_path, max_lines=MAX_LINES_PER_FILE
    ):
        total_events += 1
        if E_lost > E_remaines:
            continue

        theta = find_theta(E_remaines, E_lost, smear=True)
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

point_counts = defaultdict(int)
for x, z in zip(all_x_coords, all_z_coords):
    if not (np.isfinite(x) and np.isfinite(z)):
        continue
    point_counts[(round(x, 1), round(z, 1))] += 1

filtered_points = {p: c for p, c in point_counts.items() if c > 10}

print(f"\nДИАГНОСТИКА ФИЛЬТРАЦИИ:")
print(f"Уникальных точек до фильтрации: {len(point_counts)}")
print(f"Точек после фильтрации (count > 10): {len(filtered_points)}")
if filtered_points:
    print(f"Максимальное количество пересечений: {max(filtered_points.values())}")
    print(f"Среднее количество пересечений: {np.mean(list(filtered_points.values())):.2f}")

x_vals = [p[0] for p in filtered_points.keys()]
z_vals = [p[1] for p in filtered_points.keys()]
counts = list(filtered_points.values())

x_range = (-20, 20)
z_range = (-20, 20)
x_bins = np.arange(x_range[0], x_range[1] + cell_size, cell_size)
z_bins = np.arange(z_range[0], z_range[1] + cell_size, cell_size)

heatmap, x_edges, z_edges = np.histogram2d(x_vals, z_vals, bins=[x_bins, z_bins], weights=counts)

print(f"\nДИАГНОСТИКА HEATMAP:")
print(f"Размер heatmap: {heatmap.shape}")
print(f"Максимальное значение: {np.max(heatmap)}")
print(f"Ненулевых элементов: {np.count_nonzero(heatmap)}")

x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
max_z_idx = np.argmax(np.max(heatmap, axis=0))
profile_x = heatmap[:, max_z_idx]

if np.sum(profile_x) > 0:
    x_mean = float(np.average(x_centers, weights=profile_x))
else:
    x_mean = float(np.mean(x_centers)) if np.size(x_centers) else 0.0
y_peak = float(np.max(profile_x)) if np.size(profile_x) else 0.0


def gaussian_fixed_mean_and_y(x, A, sigma):
    return y_peak + A * (np.exp(-((x - x_mean) ** 2) / (2 * sigma ** 2)) - 1.0)


if np.sum(profile_x) > 0:
    var_est_x = np.average((x_centers - x_mean) ** 2, weights=profile_x)
    sigma0_x = float(np.sqrt(max(var_est_x, 1e-6)))
else:
    sigma0_x = 2.0
A0_x = float(np.max(profile_x) - y_peak) if np.size(profile_x) else 1.0

fit_params_x = None
try:
    fit_params_x, _ = curve_fit(
        gaussian_fixed_mean_and_y,
        x_centers, profile_x,
        p0=[A0_x, sigma0_x],
        bounds=([0.0, 1e-6], [np.inf, np.inf]),
        maxfev=20000,
    )
except Exception as e:
    print(f"Предупреждение: не удалось выполнить гауссовый фит профиля X: {e}")

fig, (ax_img, ax_prof_x) = plt.subplots(1, 2, figsize=(14, 5))

im = ax_img.imshow(
    heatmap.T, origin='lower',
    extent=[x_bins[0], x_bins[-1], z_bins[0], z_bins[-1]],
    aspect='equal', cmap='hot',
)
cbar = fig.colorbar(im, ax=ax_img)
cbar.set_label('Интенсивность (а.е.)', fontsize=18)
cbar.ax.tick_params(labelsize=18)

ax_img.set_xlabel('X (мм)', fontsize=18)
ax_img.set_ylabel('Z (мм)', fontsize=18)
ax_img.tick_params(axis='both', which='major', labelsize=18)
ax_img.grid(True, alpha=0.3)

if np.max(heatmap) > 0:
    normalized = (heatmap.T / np.max(heatmap.T) * 255).astype('uint8')
else:
    normalized = np.zeros_like(heatmap.T, dtype='uint8')
imageio.imwrite('data_test.jpg', normalized)

ax_prof_x.plot(x_centers, profile_x, 'b-', linewidth=2, label='Данные (профиль X)')
if fit_params_x is not None:
    x_fit = np.linspace(x_centers[0], x_centers[-1], 1000)
    ax_prof_x.plot(x_fit, gaussian_fixed_mean_and_y(x_fit, *fit_params_x),
                   '--', color='red', linewidth=1.5, label='Фитирование')

ax_prof_x.set_xlabel('X мм', fontsize=18)
ax_prof_x.set_ylabel('Интенсивность (а.е.)', fontsize=18)
ax_prof_x.tick_params(axis='both', which='major', labelsize=18)
ax_prof_x.grid(True, alpha=0.3)
ax_prof_x.set_xlim(x_centers[0], x_centers[-1])
ax_prof_x.legend(fontsize=12, loc='upper right')

plt.tight_layout()
plt.show()


def calculate_fwhm(heatmap, edges, axis):
    if axis == 'x':
        idx = np.argmax(np.max(heatmap, axis=0))
        profile = heatmap[:, idx]
        bins = edges[0]
    elif axis == 'z':
        idx = np.argmax(np.max(heatmap, axis=1))
        profile = heatmap[idx, :]
        bins = edges[1]
    else:
        raise ValueError("axis must be 'x' or 'z'")
    if profile.size == 0:
        return None
    mx = np.max(profile)
    if mx <= 0:
        return None
    idxs = np.where(profile >= mx / 2)[0]
    if len(idxs) == 0:
        return None
    centers = 0.5 * (bins[:-1] + bins[1:])
    return float(centers[idxs[-1]] - centers[idxs[0]])


fwhm_x = calculate_fwhm(heatmap, [x_edges, z_edges], axis='x')
fwhm_z = calculate_fwhm(heatmap, [x_edges, z_edges], axis='z')

print(f"\nРЕЗУЛЬТАТЫ РЕКОНСТРУКЦИИ:")
if fwhm_x is not None and fwhm_z is not None:
    print(f"FWHM по оси X: {fwhm_x:.2f} мм")
    print(f"FWHM по оси Z: {fwhm_z:.2f} мм")
else:
    print("FWHM не может быть вычислен — недостаточно данных")
print("=" * 60)
