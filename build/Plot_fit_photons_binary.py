import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import struct
import glob
import os

# --- Чтение всех бинарных файлов ---
# Формат каждой записи: [eventID: int32 (4б), radiatorID: int32 (4б),
#                        X: int32 (4б) + 4б padding, Z: int32 (4б) + 4б padding]
# Итого 24 байта на запись (pair<int,int> записывается с sizeof(double)=8б на каждое поле)
RECORD_SIZE = 24  # байт

records = []

script_dir = os.path.dirname(os.path.abspath(__file__))
bin_files = sorted(glob.glob(os.path.join(script_dir, 'photon_coordinates_*.bin')),
                   key=lambda f: int(os.path.basename(f).split('_')[-1].split('.')[0]))

if not bin_files:
    print("Бинарные файлы photon_coordinates_*.bin не найдены в текущей директории.")
    exit(1)

print(f"Найдено файлов: {len(bin_files)}: {bin_files}")

for fname in bin_files:
    file_size = os.path.getsize(fname)
    n_records = file_size // RECORD_SIZE
    remainder = file_size % RECORD_SIZE
    if remainder != 0:
        print(f"Предупреждение: файл {fname} содержит неполную запись ({remainder} лишних байт), они будут пропущены.")
    with open(fname, 'rb') as f:
        for _ in range(n_records):
            raw = f.read(RECORD_SIZE)
            if len(raw) < RECORD_SIZE:
                break
            event_id   = struct.unpack('i', raw[0:4])[0]
            radiator_id = struct.unpack('i', raw[4:8])[0]
            x = struct.unpack('i', raw[8:12])[0]    # int32, следующие 4б — padding
            z = struct.unpack('i', raw[16:20])[0]   # int32, следующие 4б — padding
            records.append((radiator_id, event_id, x, z))

if not records:
    print("Данные не найдены в бинарных файлах.")
    exit(1)

arr = np.array(records)  # columns: [Detector, EventID, X, Z]
print(f"Всего записей: {len(arr)}")

# --- 2D Гауссиан ---
def gaussian_2d(xy, x0, z0, sigma_x, sigma_z, amplitude):
    x, z = xy
    return amplitude * np.exp(-(((x - x0) ** 2) / (2 * sigma_x ** 2) +
                                ((z - z0) ** 2) / (2 * sigma_z ** 2)))

detectors = np.unique(arr[:, 0])
event_ids  = np.unique(arr[:, 1])

for detector in detectors:
    for event_id in event_ids:
        mask = (arr[:, 0] == detector) & (arr[:, 1] == event_id)
        sub = arr[mask]
        if len(sub) == 0:
            continue

        x_vals = sub[:, 2].astype(float)
        z_vals = sub[:, 3].astype(float)

        # --- Подсчёт числа фотонов в каждой ячейке ---
        cell_size = 5  # мм
        x_min, x_max = int(x_vals.min()), int(x_vals.max())
        z_min, z_max = int(z_vals.min()), int(z_vals.max())
        x_bins = np.arange(x_min, x_max + cell_size, cell_size)
        z_bins = np.arange(z_min, z_max + cell_size, cell_size)

        counts_raw, z_edges, x_edges = np.histogram2d(z_vals, x_vals,
                                                       bins=[z_bins - cell_size/2,
                                                             x_bins - cell_size/2])

        x_centers = (x_edges[:-1] + x_edges[1:]) / 2  # центры по X
        z_centers = (z_edges[:-1] + z_edges[1:]) / 2  # центры по Z

        if counts_raw.size < 5:
            print(f"Недостаточно данных для фита: Detector {detector}, EventID {event_id}")
            continue

        # Плоские массивы для фита
        Xc, Zc = np.meshgrid(x_centers, z_centers)
        xy_flat = np.vstack([Xc.ravel(), Zc.ravel()])
        counts_flat = counts_raw.ravel()

        # Начальные приближения
        total = counts_flat.sum()
        x0_guess = np.average(Xc.ravel(), weights=counts_flat + 1e-10)
        z0_guess = np.average(Zc.ravel(), weights=counts_flat + 1e-10)
        sigma_guess = max(np.std(x_vals), cell_size)
        amp_guess = counts_flat.max()
        p0 = (x0_guess, z0_guess, sigma_guess, sigma_guess, amp_guess)

        try:
            popt, _ = curve_fit(gaussian_2d, xy_flat, counts_flat, p0=p0,
                                maxfev=5000,
                                bounds=([-np.inf, -np.inf, 0.1, 0.1, 0],
                                        [np.inf,  np.inf,  200, 200, np.inf]))
            x0_fit, z0_fit, sigma_x, sigma_z, amplitude = popt
            fwhm_x = 2 * np.sqrt(2 * np.log(2)) * abs(sigma_x)
            fwhm_z = 2 * np.sqrt(2 * np.log(2)) * abs(sigma_z)

            print(f"Detector {detector}, EventID {event_id}:")
            print(f"  Centre X: {x0_fit:.2f} mm,  FWHM X: {fwhm_x:.2f} mm")
            print(f"  Centre Z: {z0_fit:.2f} mm,  FWHM Z: {fwhm_z:.2f} mm")
            fit_ok = True
        except (RuntimeError, TypeError, ValueError) as e:
            print(f"Ошибка фита: Detector {detector}, EventID {event_id}: {e}")
            fit_ok = False

        # --- Мелкая сетка для сглаженного фита ---
        if fit_ok:
            fine = 0.5  # мм — шаг мелкой сетки
            x_fine = np.arange(x_min - cell_size, x_max + cell_size + fine, fine)
            z_fine = np.arange(z_min - cell_size, z_max + cell_size + fine, fine)
            Xf, Zf = np.meshgrid(x_fine, z_fine)
            fitted_fine = gaussian_2d((Xf, Zf), *popt)

        # ================================================================
        # РИСУЕМ: левый — сырые данные, правый — фит (если успешен)
        # ================================================================
        ncols = 2 if fit_ok else 1
        fig, axes = plt.subplots(1, ncols, figsize=(8 * ncols, 6))
        if ncols == 1:
            axes = [axes]

        # --- Левый график: сырые данные ---
        ax0 = axes[0]
        im0 = ax0.imshow(counts_raw,
                         extent=[x_bins[0] - cell_size/2, x_bins[-1] + cell_size/2,
                                 z_bins[0] - cell_size/2, z_bins[-1] + cell_size/2],
                         origin='lower', aspect='equal', cmap='viridis',
                         interpolation='nearest')
        cb0 = fig.colorbar(im0, ax=ax0, shrink=0.75)
        cb0.set_label('Отсчеты', fontsize=20)
        cb0.ax.tick_params(labelsize=18)
        ax0.set_title(f'Данные \n Детектор {detector}, Номер события {event_id}', fontsize=20)
        ax0.set_xlabel('X, мм', fontsize=20)
        ax0.set_ylabel('Z, мм', fontsize=20)
        ax0.tick_params(labelsize=18)

        # --- Правый график: сглаженный фит ---
        if fit_ok:
            ax1 = axes[1]
            im1 = ax1.imshow(fitted_fine,
                             extent=[x_fine[0], x_fine[-1], z_fine[0], z_fine[-1]],
                             origin='lower', aspect='equal', cmap='viridis',
                             interpolation='bilinear')
            cb1 = fig.colorbar(im1, ax=ax1, shrink=0.75)
            cb1.set_label('Отсчеты', fontsize=20)
            cb1.ax.tick_params(labelsize=18)
            ax1.set_title(f'Фитирование \n Детектор {detector}, Номер события {event_id}\n'
                          , fontsize=20)
            ax1.set_xlabel('X, мм', fontsize=20)
            ax1.set_ylabel('Z, мм', fontsize=20)
            ax1.tick_params(labelsize=18)


        plt.tight_layout()
        plt.show()
