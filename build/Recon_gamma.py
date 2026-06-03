import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tqdm import tqdm
import imageio

from recon_utils import (
    load_events,
    compute_thetas,
    cone_intersections_xz,
    HeatmapAccumulator,
)

# ---------------- Параметры ------------------------------------------------

CACHE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "events_raw.npy")
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

ENERGY_SIGMA_REL = 0.00          # статистическая часть, в сотых
SIGMA_ELEC_KEV = 0.0       # электронный (аддитивный) шум, σ в кэВ (~ 1 кэВ FWHM, клип ±1 кэВ)
E_LOST_MIN = 0.0              # отсечка снизу по ИЗМЕРЕННОЙ E_lost, кэВ
X_RANGE = (-50.0, 50.0)
Z_RANGE = (-50.0, 50.0)
CELL_SIZE = 0.1
STEPS = 360
RNG_SEED = 12345
FILTER_MIN_COUNT = 5             # порог count > N

# Область для финального построения и фитирования
X_VIEW = (-20.0, 20.0)
Z_VIEW = (-20.0, 20.0)

# ---------------- Загрузка событий -----------------------------------------

events = load_events(cache_path=CACHE_PATH, data_dir=DATA_DIR)
E_lost = events[:, 0]; E_rem = events[:, 1]
Vx = events[:, 2]; Vy = events[:, 3]; Vz = events[:, 4]
X2 = events[:, 5]; Y2 = events[:, 6]; Z2 = events[:, 7]

keep = E_lost <= E_rem
print(f"После отбора E_lost <= E_rem: {keep.sum()} / {len(events)}")
E_lost = E_lost[keep]; E_rem = E_rem[keep]
Vx = Vx[keep]; Vy = Vy[keep]; Vz = Vz[keep]
X2 = X2[keep]; Y2 = Y2[keep]; Z2 = Z2[keep]

rng = np.random.default_rng(RNG_SEED)
theta, theta_ok = compute_thetas(E_lost, E_rem, sigma_rel=ENERGY_SIGMA_REL,
                                 sigma_elec_keV=SIGMA_ELEC_KEV,
                                 e_lost_min=E_LOST_MIN, rng=rng)
print(f"Валидных θ (σ_rel={ENERGY_SIGMA_REL*100:.1f}%, σ_elec={SIGMA_ELEC_KEV} кэВ, "
      f"E_lost_m > {E_LOST_MIN} кэВ): {theta_ok.sum()} / {len(theta)}")
idx = np.where(theta_ok)[0]

# ---------------- Потоковая аккумуляция heatmap ----------------------------

acc = HeatmapAccumulator(X_RANGE, Z_RANGE, CELL_SIZE)

used = 0
for i in tqdm(idx, desc="События", unit="ev"):
    xs, zs = cone_intersections_xz(
        Vx[i], Vy[i], Vz[i],
        X2[i], Y2[i], Z2[i],
        theta[i], steps=STEPS,
    )
    if xs is None:
        continue
    acc.add(xs, zs)
    used += 1

print(f"\nИспользовано событий: {used}")
print(f"Heatmap shape: {acc.heatmap.shape}, ненулевых ячеек: {np.count_nonzero(acc.heatmap)}")
print(f"Max в heatmap: {acc.heatmap.max():.0f}")

# Порог
heatmap_full = acc.heatmap.copy()
heatmap_full[heatmap_full <= FILTER_MIN_COUNT] = 0

# Обрезка до области отображения
x_centers_full = acc.x_centers
z_centers_full = acc.z_centers
ix_mask = (x_centers_full >= X_VIEW[0]) & (x_centers_full <= X_VIEW[1])
iz_mask = (z_centers_full >= Z_VIEW[0]) & (z_centers_full <= Z_VIEW[1])

heatmap = heatmap_full[np.ix_(ix_mask, iz_mask)]
x_centers = x_centers_full[ix_mask]
z_centers = z_centers_full[iz_mask]

# Кромки для imshow
x_edges = np.concatenate([
    x_centers - CELL_SIZE / 2,
    [x_centers[-1] + CELL_SIZE / 2]
])
z_edges = np.concatenate([
    z_centers - CELL_SIZE / 2,
    [z_centers[-1] + CELL_SIZE / 2]
])

# ---------------- Профили X, Z через срез по точке максимума ---------------

max_z_idx = int(np.argmax(np.max(heatmap, axis=0)))
max_x_idx = int(np.argmax(np.max(heatmap, axis=1)))
profile_x = heatmap[:, max_z_idx]
profile_z = heatmap[max_x_idx, :]


# ---------------- Псевдо-Войт фит -------------------------------------------
# Параметризация с общим FWHM: оба компонента нормированы на 1 в вершине,
# имеют ОДИНАКОВЫЙ FWHM, eta - доля Лоренца (0 = чистый Гаусс).
LN2 = np.log(2.0)


def pseudo_voigt(x, A, x0, fwhm, eta, C):
    # Обе компоненты имеют один и тот же FWHM, нормированы на 1 в вершине.
    # Gauss:  exp(-4 ln2 · (x-x0)² / FWHM²)   -> при |x-x0|=FWHM/2 значение = 1/2.
    # Lorentz: 1 / (1 + (2(x-x0)/FWHM)²)
    half = fwhm / 2.0
    G = np.exp(-LN2 * ((x - x0) / half) ** 2)
    L = 1.0 / (1.0 + ((x - x0) / half) ** 2)
    return C + A * (eta * L + (1.0 - eta) * G)


def fit_pseudo_voigt(centers, profile):
    if profile.sum() <= 0:
        return None
    s = profile.sum()
    mu0 = float(np.sum(centers * profile) / s)
    var0 = float(np.sum(((centers - mu0) ** 2) * profile) / s)
    sigma0 = float(np.sqrt(max(var0, 1e-6)))
    A0 = float(profile.max() - profile.min())
    C0 = float(profile.min())
    # FWHM_gauss = 2*sqrt(2*ln2)*sigma ≈ 2.355*sigma — это начальное приближение
    fwhm0 = 2.355 * sigma0
    p0 = [A0, mu0, fwhm0, 0.5, C0]
    bounds = (
        [0.0, centers[0], 1e-3, 0.0, 0.0],
        [np.inf, centers[-1], (centers[-1] - centers[0]), 1.0, np.inf],
    )
    try:
        popt, _ = curve_fit(pseudo_voigt, centers, profile, p0=p0, bounds=bounds, maxfev=40000)
        return popt
    except Exception as e:
        print(f"Предупреждение: фит pseudo-Voigt не сошёлся: {e}")
        return None


pv_x = fit_pseudo_voigt(x_centers, profile_x)
pv_z = fit_pseudo_voigt(z_centers, profile_z)

if pv_x is not None:
    A, x0, fwhm_fit_x, eta_x, C = pv_x
    print(f"\npseudo-Voigt X: x0={x0:.2f} мм, FWHM={fwhm_fit_x:.2f} мм, eta={eta_x:.2f}, A={A:.0f}, C={C:.0f}")
if pv_z is not None:
    A, z0, fwhm_fit_z, eta_z, C = pv_z
    print(f"pseudo-Voigt Z: z0={z0:.2f} мм, FWHM={fwhm_fit_z:.2f} мм, eta={eta_z:.2f}, A={A:.0f}, C={C:.0f}")

# ---------------- Визуализация --------------------------------------------

fig, (ax_img, ax_prof) = plt.subplots(1, 2, figsize=(14, 5))

im = ax_img.imshow(
    heatmap.T, origin='lower',
    extent=[x_edges[0], x_edges[-1], z_edges[0], z_edges[-1]],
    aspect='equal', cmap='hot',
)
cbar = fig.colorbar(im, ax=ax_img)
cbar.set_label('Интенсивность (а.е.)', fontsize=18)
cbar.ax.tick_params(labelsize=18)
ax_img.set_xlabel('X (мм)', fontsize=18)
ax_img.set_ylabel('Z (мм)', fontsize=18)
ax_img.tick_params(axis='both', which='major', labelsize=18)
ax_img.grid(True, alpha=0.3)

if heatmap.max() > 0:
    norm_img = (heatmap.T / heatmap.max() * 255).astype('uint8')
else:
    norm_img = np.zeros_like(heatmap.T, dtype='uint8')
imageio.imwrite('data_test.jpg', norm_img)

ax_prof.plot(x_centers, profile_x, 'b-', linewidth=2, label='Профиль X')
if pv_x is not None:
    xf = np.linspace(x_centers[0], x_centers[-1], 1000)
    ax_prof.plot(xf, pseudo_voigt(xf, *pv_x), '--', color='red', linewidth=1.5,
                 label=f'псевдо-Войт фит (η={pv_x[3]:.2f})')
ax_prof.set_xlabel('X мм', fontsize=18)
ax_prof.set_ylabel('Интенсивность (а.е.)', fontsize=18)
ax_prof.tick_params(axis='both', which='major', labelsize=18)
ax_prof.grid(True, alpha=0.3)
ax_prof.set_xlim(x_centers[0], x_centers[-1])
ax_prof.legend(fontsize=12, loc='upper right')

plt.tight_layout()
plt.show()


def calculate_fwhm(profile, centers):
    """FWHM с вычетом фона (минимум профиля в окне) и поиском связной области
    вокруг точки максимума: half = base + (peak - base) / 2.
    """
    if profile.size == 0:
        return None
    base = float(profile.min())
    peak = float(profile.max())
    if peak - base <= 0:
        return None
    half = base + (peak - base) / 2
    imax = int(np.argmax(profile))
    above = profile >= half
    if not above[imax]:
        return None
    left = imax
    while left > 0 and above[left - 1]:
        left -= 1
    right = imax
    while right < len(profile) - 1 and above[right + 1]:
        right += 1
    return float(centers[right] - centers[left])


fwhm_x_num = calculate_fwhm(profile_x, x_centers)
fwhm_z_num = calculate_fwhm(profile_z, z_centers)

print(f"\nРЕЗУЛЬТАТЫ:")
print("Численный FWHM (по половине максимума, без модели):")
if fwhm_x_num is not None: print(f"  X: {fwhm_x_num:.2f} мм")
if fwhm_z_num is not None: print(f"  Z: {fwhm_z_num:.2f} мм")
print("pseudo-Voigt FWHM (из фита):")
if pv_x is not None: print(f"  X: {pv_x[2]:.2f} мм  (eta={pv_x[3]:.2f})")
if pv_z is not None: print(f"  Z: {pv_z[2]:.2f} мм  (eta={pv_z[3]:.2f})")
print("=" * 60)
