"""List-mode OSEM реконструкция точечного источника для Комптон-камеры.

В качестве начального приближения используется простая backprojection
(тот же алгоритм, что в Recon_gamma.py — конус, повёрнутый на 360°,
пересечённый с плоскостью y=0).

Системная матрица: A_ij = exp(-((α_ij - θ_i)^2) / (2 σ_θ^2)),
где α_ij — угол между осью конуса i (V_i → второй детектор) и вектором
V_i → пиксель j на плоскости y=0; θ_i — комптоновский угол события.

Сетка изображения здесь грубее, чем у backprojection (по умолчанию
0.5 мм против 0.1 мм), потому что OSEM требует хранения и пересчёта
A для каждого пикселя на каждой итерации. Для сравнения формы
ядра пика 0.5 мм достаточно.
"""

import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
from tqdm import tqdm

from recon_utils import (
    load_events,
    compute_thetas,
    cone_intersections_xz,
    HeatmapAccumulator,
)

# ============================================================================
# Параметры
# ============================================================================

CACHE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "events_raw.npy")
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

ENERGY_SIGMA_REL = 0.05          # статистическая часть, в сотых
SIGMA_ELEC_KEV = 0.42            # электронный (аддитивный) шум, σ в кэВ (~ 1 кэВ FWHM, клип ±1 кэВ)
E_LOST_MIN = 0.0                 # отсечка снизу по ИЗМЕРЕННОЙ E_lost, кэВ
RNG_SEED = 12345

# Сетка изображения OSEM
X_RANGE = (-20.0, 20.0)
Z_RANGE = (-20.0, 20.0)
CELL_SIZE = 0.5                  # шаг сетки OSEM, мм

# Параметры OSEM
N_SUBSETS = 4
N_ITERATIONS = 8
SIGMA_THETA_DEG = 1.5            # σ углового ядра системной матрицы, градусы
MAX_EVENTS = 150_000             # подвыборка событий (None — все); масштабирует время линейно
CHUNK_SIZE = 4000                # размер чанка событий при векторизации
USE_UNIFORM_SENSITIVITY = True   # True: S_j = const; False: data-driven S_j = Σ_i A_ij
                                  # Для list-mode при реконструкции точечного источника
                                  # data-driven sensitivity завышена в центре и подавляет пик.

# Фильтрация
FILTER_MIN_COUNT = 0             # порог для финального отображения (0 — без отсечки)
POST_SMOOTH_MM = 1.0             # σ гауссова пост-фильтра, мм (0 — выключить)
                                  # σ ≪ ожидаемой FWHM (~3–5 мм), чтобы не уширить пик,
                                  # но достаточно, чтобы погасить редкие «горячие» пиксели по краям.

# ============================================================================
# Загрузка событий (использует тот же кэш, что и Recon_gamma.py)
# ============================================================================

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
idx_all = np.where(theta_ok)[0]

if MAX_EVENTS is not None and len(idx_all) > MAX_EVENTS:
    idx_all = rng.choice(idx_all, size=MAX_EVENTS, replace=False)
    idx_all.sort()
print(f"OSEM: используется {len(idx_all)} событий")

Vx_e = Vx[idx_all].astype(np.float32)
Vy_e = Vy[idx_all].astype(np.float32)
Vz_e = Vz[idx_all].astype(np.float32)
X2_e = X2[idx_all].astype(np.float32)
Y2_e = Y2[idx_all].astype(np.float32)
Z2_e = Z2[idx_all].astype(np.float32)
theta_e = theta[idx_all].astype(np.float32)
N = len(idx_all)

# Ось конуса возможных позиций ИСТОЧНИКА: -d_out = -(V1 → V2).
# Источник находится на конусе с вершиной в V1, осью -d_out и полу-углом θ.
# (См. кинематику Комптона: θ — угол между ИСХОДНЫМ направлением фотона и рассеянным.)
Px = (X2_e - Vx_e)
Py = (Y2_e - Vy_e)
Pz = (Z2_e - Vz_e)
P_norm = np.sqrt(Px * Px + Py * Py + Pz * Pz)
nPx = (-Px / P_norm).astype(np.float32)
nPy = (-Py / P_norm).astype(np.float32)
nPz = (-Pz / P_norm).astype(np.float32)

# ============================================================================
# Сетка изображения
# ============================================================================

x_edges = np.arange(X_RANGE[0], X_RANGE[1] + CELL_SIZE, CELL_SIZE)
z_edges = np.arange(Z_RANGE[0], Z_RANGE[1] + CELL_SIZE, CELL_SIZE)
x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
nx, nz = len(x_centers), len(z_centers)
Npix = nx * nz
print(f"Сетка OSEM: {nx} x {nz} = {Npix} пикселей, шаг {CELL_SIZE} мм")

XX, ZZ = np.meshgrid(x_centers, z_centers, indexing='ij')
pix_x = XX.ravel().astype(np.float32)              # (Npix,)
pix_z = ZZ.ravel().astype(np.float32)

# ============================================================================
# Начальное приближение λ⁰ — backprojection через cone_intersections_xz
# ============================================================================

print("\nНачальное приближение (backprojection)...")
acc = HeatmapAccumulator(X_RANGE, Z_RANGE, CELL_SIZE)
for k in tqdm(range(N), desc="Backproj", unit="ev"):
    xs, zs = cone_intersections_xz(
        Vx_e[k], Vy_e[k], Vz_e[k],
        X2_e[k], Y2_e[k], Z2_e[k],
        float(theta_e[k]),
        steps=360,
    )
    if xs is None:
        continue
    acc.add(xs, zs)

lam0_raw = acc.heatmap.astype(np.float64)          # (nx, nz)
print(f"  max(λ⁰_raw) = {lam0_raw.max():.0f}, ненулевых = {np.count_nonzero(lam0_raw)}")

# Для list-mode MLEM с y_i ≡ 1 правильная нормировка такая, чтобы forward
# A·λ для типичного события был ~1. Эквивалентно: mean(λ) ≈ 1 / mean(A_chunk).
# Грубая прикидка mean(A) — отношение ширины ядра к телесному углу пикселей,
# но проще пронормировать так, чтобы mean(λ⁰) = 1. EM это поправит на масштабе.
if lam0_raw.sum() > 0:
    lam0 = lam0_raw / lam0_raw.mean()
else:
    lam0 = np.ones_like(lam0_raw)

# EM запрещает нули — они «замораживают» пиксель навсегда.
EPS = max(1e-6 * lam0.max(), 1e-12)
lam = np.maximum(lam0, EPS).copy()
print(f"  λ⁰: mean(λ)={lam.mean():.3f}, max(λ)={lam.max():.3f}")


# ============================================================================
# Системная матрица: гауссово ядро в угловом пространстве (чанками)
# ============================================================================

SIGMA_TH = np.radians(SIGMA_THETA_DEG)
TWO_SIG2 = 2.0 * SIGMA_TH * SIGMA_TH
THETA_RAD = np.radians(theta_e).astype(np.float32)  # (N,)


def system_matrix_chunk(event_idx):
    """Строит фрагмент A для chunk событий: shape (Ne, Npix).
    Использует cos между (пиксель−V) и нормированной осью конуса.
    """
    Ne = len(event_idx)

    Vx_b = Vx_e[event_idx, None]    # (Ne, 1)
    Vy_b = Vy_e[event_idx, None]
    Vz_b = Vz_e[event_idx, None]
    nPx_b = nPx[event_idx, None]
    nPy_b = nPy[event_idx, None]
    nPz_b = nPz[event_idx, None]
    th_b = THETA_RAD[event_idx, None]   # (Ne, 1)

    ux = pix_x[None, :] - Vx_b      # (Ne, Npix)
    uy = -Vy_b                       # (Ne, 1) — broadcast
    uz = pix_z[None, :] - Vz_b

    u_norm = np.sqrt(ux * ux + uy * uy + uz * uz)
    dot = ux * nPx_b + uy * nPy_b + uz * nPz_b
    cos_alpha = np.clip(dot / np.maximum(u_norm, 1e-12), -1.0, 1.0)
    alpha = np.arccos(cos_alpha)
    diff = alpha - th_b

    A = np.exp(-(diff * diff) / TWO_SIG2).astype(np.float32)
    return A


# ============================================================================
# Глобальная sensitivity S_j = Σ_i A_ij — нужна для финальной нормализации.
# Для OSEM с per-subset sensitivity её можно не считать; но полезно иметь
# для диагностики.
# ============================================================================

# В классическом OSEM используется per-subset sensitivity, и это уже учтено в
# цикле ниже. Однако подсчёт глобальной S дополнительно стоит ровно один проход
# по всем событиям — посчитаем для диагностики.
print("\nПодсчёт sensitivity (1 проход по событиям)...")
S_global = np.zeros(Npix, dtype=np.float64)
for cstart in tqdm(range(0, N, CHUNK_SIZE), desc="Sensitivity", unit="chunk"):
    cend = min(cstart + CHUNK_SIZE, N)
    A = system_matrix_chunk(np.arange(cstart, cend, dtype=np.int64))
    S_global += A.sum(axis=0)
print(f"  max(S) = {S_global.max():.2f}, min(S) = {S_global.min():.4f}")

# ============================================================================
# OSEM итерации
# ============================================================================

# Случайное разбиение событий на подмножества
perm = rng.permutation(N)
subsets = [np.sort(s) for s in np.array_split(perm, N_SUBSETS)]

print(f"\nOSEM: {N_ITERATIONS} итераций × {N_SUBSETS} подмножеств")

for it in range(N_ITERATIONS):
    for sub_i, sub in enumerate(subsets):
        S_sub = np.zeros(Npix, dtype=np.float64)
        update = np.zeros(Npix, dtype=np.float64)
        lam_flat = lam.ravel()

        for cstart in range(0, len(sub), CHUNK_SIZE):
            cend = min(cstart + CHUNK_SIZE, len(sub))
            ev_idx = sub[cstart:cend]
            A = system_matrix_chunk(ev_idx)            # (chunk, Npix)
            forward = A @ lam_flat                      # (chunk,)
            forward = np.maximum(forward, EPS)
            update += A.T @ (1.0 / forward)
            S_sub += A.sum(axis=0)

        if USE_UNIFORM_SENSITIVITY:
            # Усредняем S по полю — теряем правильную нормировку, но не
            # подавляем центральный пик (data-driven S завышена в центре,
            # потому что туда смотрит больше конусов).
            S_eff = np.full_like(S_sub, S_sub.mean())
        else:
            S_eff = np.maximum(S_sub, EPS)

        lam_flat = lam_flat * (update / S_eff)
        np.maximum(lam_flat, EPS, out=lam_flat)
        # Стабилизируем масштаб: каждый шаг возвращаем mean(λ) к 1.
        lam_flat *= (1.0 / lam_flat.mean())
        lam = lam_flat.reshape(nx, nz)

        print(f"  iter {it+1}/{N_ITERATIONS}  subset {sub_i+1}/{N_SUBSETS}  "
              f"max(λ)={lam.max():.3f}  mean(λ)={lam.mean():.3f}")

# ============================================================================
# Финальные карты и фит профилей (срез по точке максимума)
# ============================================================================

heatmap_bp = lam0.copy()
heatmap_osem = lam.copy()

if POST_SMOOTH_MM > 0:
    sigma_px = POST_SMOOTH_MM / CELL_SIZE
    heatmap_osem = gaussian_filter(heatmap_osem, sigma=sigma_px)
    print(f"Пост-фильтр: gaussian σ={POST_SMOOTH_MM} мм ({sigma_px:.2f} px)")

if FILTER_MIN_COUNT > 0:
    heatmap_bp[heatmap_bp <= FILTER_MIN_COUNT] = 0
    heatmap_osem[heatmap_osem <= FILTER_MIN_COUNT] = 0


# Куда «протыкать» срез для профиля. Для точечного источника в (0,0) это
# (0,0); ближайший к нулю пиксель — индекс ix0, iz0.
SLICE_THROUGH = (0.0, 0.0)


def slice_profiles(heatmap, x_c, z_c, through=SLICE_THROUGH):
    ix0 = int(np.argmin(np.abs(x_c - through[0])))
    iz0 = int(np.argmin(np.abs(z_c - through[1])))
    return heatmap[:, iz0], heatmap[ix0, :]


prof_x_bp, prof_z_bp = slice_profiles(heatmap_bp, x_centers, z_centers)
prof_x_os, prof_z_os = slice_profiles(heatmap_osem, x_centers, z_centers)


LN2 = np.log(2.0)


def pseudo_voigt(x, A, x0, fwhm, eta, C):
    half = fwhm / 2.0
    G = np.exp(-LN2 * ((x - x0) / half) ** 2)
    L = 1.0 / (1.0 + ((x - x0) / half) ** 2)
    return C + A * (eta * L + (1.0 - eta) * G)


def fit_pv(centers, profile):
    if profile.sum() <= 0:
        return None
    s = profile.sum()
    mu0 = float(np.sum(centers * profile) / s)
    var0 = float(np.sum(((centers - mu0) ** 2) * profile) / s)
    sigma0 = float(np.sqrt(max(var0, 1e-6)))
    A0 = float(profile.max() - profile.min())
    C0 = float(profile.min())
    rng_c = float(centers[-1] - centers[0])
    pad = 0.05 * rng_c
    # Если weighted-mean выскочил за границы (типично для плоского профиля) —
    # стартуем из центра, а bounds расширяем на 5%.
    mu0 = float(np.clip(mu0, centers[0] - pad + 1e-6, centers[-1] + pad - 1e-6))
    p0 = [max(A0, 1e-6), mu0, 2.355 * sigma0, 0.5, C0]
    bounds = (
        [0.0, centers[0] - pad, 1e-3, 0.0, 0.0],
        [np.inf, centers[-1] + pad, 2.0 * rng_c, 1.0, np.inf],
    )
    try:
        popt, _ = curve_fit(pseudo_voigt, centers, profile, p0=p0, bounds=bounds, maxfev=40000)
        return popt
    except Exception as e:
        print(f"  pseudo-Voigt не сошёлся: {e}")
        return None


def numeric_fwhm(profile, centers):
    """FWHM с вычетом фона (минимум профиля) и поиском связной области
    вокруг точки максимума, чтобы шум на других концах не расширял оценку."""
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
    # связная область вокруг imax
    left = imax
    while left > 0 and above[left - 1]:
        left -= 1
    right = imax
    while right < len(profile) - 1 and above[right + 1]:
        right += 1
    return float(centers[right] - centers[left])


pv_bp_x = fit_pv(x_centers, prof_x_bp)
pv_bp_z = fit_pv(z_centers, prof_z_bp)
pv_os_x = fit_pv(x_centers, prof_x_os)
pv_os_z = fit_pv(z_centers, prof_z_os)

print(f"\n{'='*60}")
print("ИТОГ (численный FWHM / pseudo-Voigt FWHM, η):")
def line(tag, prof, c, pv):
    nf = numeric_fwhm(prof, c)
    if pv is not None:
        print(f"  {tag}: {nf:.2f} мм / {pv[2]:.2f} мм, η={pv[3]:.2f}")
    else:
        print(f"  {tag}: {nf:.2f} мм / —")
line("Backprojection X", prof_x_bp, x_centers, pv_bp_x)
line("Backprojection Z", prof_z_bp, z_centers, pv_bp_z)
line("OSEM           X", prof_x_os, x_centers, pv_os_x)
line("OSEM           Z", prof_z_os, z_centers, pv_os_z)
print("=" * 60)

# ============================================================================
# Визуализация: 2x2 — backprojection и OSEM, картина + профиль X
# ============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 11))

def show_image(ax, hm, title):
    im = ax.imshow(
        hm.T, origin='lower',
        extent=[x_edges[0], x_edges[-1], z_edges[0], z_edges[-1]],
        aspect='equal', cmap='hot',
    )
    cb = fig.colorbar(im, ax=ax)
    cb.set_label('Интенсивность (а.е.)', fontsize=14)
    cb.ax.tick_params(labelsize=12)
    ax.set_xlabel('X (мм)', fontsize=14)
    ax.set_ylabel('Z (мм)', fontsize=14)
    ax.set_title(title, fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, alpha=0.3)


def show_profile(ax, profile, centers, pv, title):
    ax.plot(centers, profile, 'b-', linewidth=2, label='Профиль X')
    if pv is not None:
        xf = np.linspace(centers[0], centers[-1], 1000)
        ax.plot(xf, pseudo_voigt(xf, *pv), '--', color='red', linewidth=1.5,
                label=f'pV (FWHM={pv[2]:.2f} мм, η={pv[3]:.2f})')
    ax.set_xlabel('X (мм)', fontsize=14)
    ax.set_ylabel('Интенсивность (а.е.)', fontsize=14)
    ax.set_title(title, fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(centers[0], centers[-1])
    ax.legend(fontsize=10, loc='upper right')


show_image(axes[0, 0], heatmap_bp,
           f"Backprojection (init λ⁰)  σE/E={ENERGY_SIGMA_REL*100:.1f}%")
show_profile(axes[0, 1], prof_x_bp, x_centers, pv_bp_x,
             "Профиль X — backprojection")
osem_title = f"OSEM ({N_ITERATIONS}×{N_SUBSETS}, σθ={SIGMA_THETA_DEG}°)"
if POST_SMOOTH_MM > 0:
    osem_title += f", post σ={POST_SMOOTH_MM} мм"
show_image(axes[1, 0], heatmap_osem, osem_title)
show_profile(axes[1, 1], prof_x_os, x_centers, pv_os_x,
             "Профиль X — OSEM")

plt.tight_layout()
plt.show()
