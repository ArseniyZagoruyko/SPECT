"""Прогон recon_OSEM по нескольким значениям σ_E/E и двум отсечкам E_lost.

Аналог recon_gamma_sweep.py, но с OSEM-реконструкцией. Для каждой пары
(E_LOST_MIN, σ_E/E):
  • строится backprojection (начальное приближение);
  • запускается OSEM (фикс параметры N_SUBSETS × N_ITERATIONS);
  • применяется gaussian post-filter;
  • профили + pseudo-Voigt fit;
  • PNG + строка FWHM в fwhm_summary.txt.

Структура повторяет recon_OSEM.py, чтобы результаты совпадали при тех же
параметрах. Идентичность с одиночным запуском обеспечивается тем, что rng
переинициализируется внутри цикла с одинаковым RNG_SEED.
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
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

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE_PATH = os.path.join(HERE, "events_raw.npy")
DATA_DIR = os.path.join(HERE, "data")

SIGMA_LIST = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07]
E_LOST_MIN_LIST = [0.0, 5.0]
SIGMA_ELEC_KEV = 0.42

RNG_SEED = 12345

# Сетка OSEM
X_RANGE = (-20.0, 20.0)
Z_RANGE = (-20.0, 20.0)
CELL_SIZE = 0.5

# OSEM
N_SUBSETS = 4
N_ITERATIONS = 8
SIGMA_THETA_DEG = 1.5
MAX_EVENTS = 150_000
CHUNK_SIZE = 4000
USE_UNIFORM_SENSITIVITY = True
POST_SMOOTH_MM = 1.0
FILTER_MIN_COUNT = 0
SLICE_THROUGH = (0.0, 0.0)

# ============================================================================
# Загрузка событий (один раз)
# ============================================================================

events = load_events(cache_path=CACHE_PATH, data_dir=DATA_DIR)
E_lost_all = events[:, 0]; E_rem_all = events[:, 1]
Vx_all = events[:, 2]; Vy_all = events[:, 3]; Vz_all = events[:, 4]
X2_all = events[:, 5]; Y2_all = events[:, 6]; Z2_all = events[:, 7]

keep = E_lost_all <= E_rem_all
print(f"После отбора E_lost <= E_rem: {keep.sum()} / {len(events)}")
E_lost_all = E_lost_all[keep]; E_rem_all = E_rem_all[keep]
Vx_all = Vx_all[keep]; Vy_all = Vy_all[keep]; Vz_all = Vz_all[keep]
X2_all = X2_all[keep]; Y2_all = Y2_all[keep]; Z2_all = Z2_all[keep]

# ============================================================================
# Псевдо-Войт фит и численный FWHM
# ============================================================================

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


def slice_profiles(heatmap, x_c, z_c, through=SLICE_THROUGH):
    ix0 = int(np.argmin(np.abs(x_c - through[0])))
    iz0 = int(np.argmin(np.abs(z_c - through[1])))
    return heatmap[:, iz0], heatmap[ix0, :]


# ============================================================================
# Сетка изображения (фиксирована для всего свипа)
# ============================================================================

x_edges = np.arange(X_RANGE[0], X_RANGE[1] + CELL_SIZE, CELL_SIZE)
z_edges = np.arange(Z_RANGE[0], Z_RANGE[1] + CELL_SIZE, CELL_SIZE)
x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
nx, nz = len(x_centers), len(z_centers)
Npix = nx * nz
print(f"Сетка OSEM: {nx} x {nz} = {Npix} пикселей, шаг {CELL_SIZE} мм")

XX, ZZ = np.meshgrid(x_centers, z_centers, indexing='ij')
pix_x = XX.ravel().astype(np.float32)
pix_z = ZZ.ravel().astype(np.float32)


# ============================================================================
# Цикл по (E_LOST_MIN, σ_rel)
# ============================================================================

for E_LOST_MIN in E_LOST_MIN_LIST:
    OUT_DIR = os.path.join(
        HERE,
        f"sweep_osem_Elost{E_LOST_MIN:g}keV_elec{SIGMA_ELEC_KEV:.2f}keV",
    )
    SUMMARY_PATH = os.path.join(OUT_DIR, "fwhm_summary.txt")
    os.makedirs(OUT_DIR, exist_ok=True)

    print("\n" + "#" * 70)
    print(f"# E_LOST_MIN = {E_LOST_MIN} кэВ   →   {OUT_DIR}")
    print("#" * 70)

    with open(SUMMARY_PATH, "w") as fsum:
        fsum.write(
            "# σ_E/E[%]  "
            "BP_num_X[mm]  BP_pV_X[mm]  eta_BP_X  "
            "OS_num_X[mm]  OS_pV_X[mm]  eta_OS_X  "
            "BP_num_Z[mm]  BP_pV_Z[mm]  eta_BP_Z  "
            "OS_num_Z[mm]  OS_pV_Z[mm]  eta_OS_Z\n"
        )

    for sigma_rel in SIGMA_LIST:
        print("\n" + "=" * 70)
        print(f"σ_E/E = {sigma_rel*100:.0f}%   (E_LOST_MIN={E_LOST_MIN} кэВ)")
        print("=" * 70)

        rng = np.random.default_rng(RNG_SEED)
        theta, theta_ok = compute_thetas(
            E_lost_all, E_rem_all,
            sigma_rel=sigma_rel,
            sigma_elec_keV=SIGMA_ELEC_KEV,
            e_lost_min=E_LOST_MIN, rng=rng,
        )
        print(f"Валидных θ: {theta_ok.sum()} / {len(theta)}")
        idx_all = np.where(theta_ok)[0]

        if MAX_EVENTS is not None and len(idx_all) > MAX_EVENTS:
            idx_all = rng.choice(idx_all, size=MAX_EVENTS, replace=False)
            idx_all.sort()
        print(f"OSEM: используется {len(idx_all)} событий")

        Vx_e = Vx_all[idx_all].astype(np.float32)
        Vy_e = Vy_all[idx_all].astype(np.float32)
        Vz_e = Vz_all[idx_all].astype(np.float32)
        X2_e = X2_all[idx_all].astype(np.float32)
        Y2_e = Y2_all[idx_all].astype(np.float32)
        Z2_e = Z2_all[idx_all].astype(np.float32)
        theta_e = theta[idx_all].astype(np.float32)
        N = len(idx_all)

        # Ось конуса возможных позиций источника
        Px = (X2_e - Vx_e); Py = (Y2_e - Vy_e); Pz = (Z2_e - Vz_e)
        P_norm = np.sqrt(Px * Px + Py * Py + Pz * Pz)
        nPx = (-Px / P_norm).astype(np.float32)
        nPy = (-Py / P_norm).astype(np.float32)
        nPz = (-Pz / P_norm).astype(np.float32)

        # ----- Backprojection (начальное приближение) -----
        print("Backprojection...")
        acc = HeatmapAccumulator(X_RANGE, Z_RANGE, CELL_SIZE)
        for k in tqdm(range(N), desc=f"BP σ={sigma_rel*100:.0f}%, Emin={E_LOST_MIN}",
                      unit="ev"):
            xs, zs = cone_intersections_xz(
                Vx_e[k], Vy_e[k], Vz_e[k],
                X2_e[k], Y2_e[k], Z2_e[k],
                float(theta_e[k]),
                steps=360,
            )
            if xs is None:
                continue
            acc.add(xs, zs)
        lam0_raw = acc.heatmap.astype(np.float64)
        print(f"  max(λ⁰_raw) = {lam0_raw.max():.0f}, ненулевых = {np.count_nonzero(lam0_raw)}")

        if lam0_raw.sum() > 0:
            lam0 = lam0_raw / lam0_raw.mean()
        else:
            lam0 = np.ones_like(lam0_raw)
        EPS = max(1e-6 * lam0.max(), 1e-12)
        lam = np.maximum(lam0, EPS).copy()

        # ----- Системная матрица (внутри замыкания, чтобы видеть локальные nPx/Vx_e/...)
        SIGMA_TH = np.radians(SIGMA_THETA_DEG)
        TWO_SIG2 = 2.0 * SIGMA_TH * SIGMA_TH
        THETA_RAD = np.radians(theta_e).astype(np.float32)

        def system_matrix_chunk(event_idx):
            Vx_b = Vx_e[event_idx, None]; Vy_b = Vy_e[event_idx, None]; Vz_b = Vz_e[event_idx, None]
            nPx_b = nPx[event_idx, None]; nPy_b = nPy[event_idx, None]; nPz_b = nPz[event_idx, None]
            th_b = THETA_RAD[event_idx, None]
            ux = pix_x[None, :] - Vx_b
            uy = -Vy_b
            uz = pix_z[None, :] - Vz_b
            u_norm = np.sqrt(ux * ux + uy * uy + uz * uz)
            dot = ux * nPx_b + uy * nPy_b + uz * nPz_b
            cos_alpha = np.clip(dot / np.maximum(u_norm, 1e-12), -1.0, 1.0)
            alpha = np.arccos(cos_alpha)
            diff = alpha - th_b
            return np.exp(-(diff * diff) / TWO_SIG2).astype(np.float32)

        # ----- OSEM -----
        perm = rng.permutation(N)
        subsets = [np.sort(s) for s in np.array_split(perm, N_SUBSETS)]

        print(f"OSEM: {N_ITERATIONS} итераций × {N_SUBSETS} подмножеств")
        for it in range(N_ITERATIONS):
            for sub_i, sub in enumerate(subsets):
                S_sub = np.zeros(Npix, dtype=np.float64)
                update = np.zeros(Npix, dtype=np.float64)
                lam_flat = lam.ravel()
                for cstart in range(0, len(sub), CHUNK_SIZE):
                    cend = min(cstart + CHUNK_SIZE, len(sub))
                    ev_idx = sub[cstart:cend]
                    A = system_matrix_chunk(ev_idx)
                    forward = A @ lam_flat
                    forward = np.maximum(forward, EPS)
                    update += A.T @ (1.0 / forward)
                    S_sub += A.sum(axis=0)
                if USE_UNIFORM_SENSITIVITY:
                    S_eff = np.full_like(S_sub, S_sub.mean())
                else:
                    S_eff = np.maximum(S_sub, EPS)
                lam_flat = lam_flat * (update / S_eff)
                np.maximum(lam_flat, EPS, out=lam_flat)
                lam_flat *= (1.0 / lam_flat.mean())
                lam = lam_flat.reshape(nx, nz)

        # ----- Финальная обработка -----
        heatmap_bp = lam0.copy()
        heatmap_osem = lam.copy()
        if POST_SMOOTH_MM > 0:
            sigma_px = POST_SMOOTH_MM / CELL_SIZE
            heatmap_osem = gaussian_filter(heatmap_osem, sigma=sigma_px)
        if FILTER_MIN_COUNT > 0:
            heatmap_bp[heatmap_bp <= FILTER_MIN_COUNT] = 0
            heatmap_osem[heatmap_osem <= FILTER_MIN_COUNT] = 0

        prof_x_bp, prof_z_bp = slice_profiles(heatmap_bp, x_centers, z_centers)
        prof_x_os, prof_z_os = slice_profiles(heatmap_osem, x_centers, z_centers)

        pv_bp_x = fit_pv(x_centers, prof_x_bp)
        pv_bp_z = fit_pv(z_centers, prof_z_bp)
        pv_os_x = fit_pv(x_centers, prof_x_os)
        pv_os_z = fit_pv(z_centers, prof_z_os)

        def vals(pv):
            return (pv[2], pv[3]) if pv is not None else (float("nan"), float("nan"))

        bp_x_num = numeric_fwhm(prof_x_bp, x_centers)
        bp_z_num = numeric_fwhm(prof_z_bp, z_centers)
        os_x_num = numeric_fwhm(prof_x_os, x_centers)
        os_z_num = numeric_fwhm(prof_z_os, z_centers)
        bp_x_pv, bp_x_eta = vals(pv_bp_x)
        bp_z_pv, bp_z_eta = vals(pv_bp_z)
        os_x_pv, os_x_eta = vals(pv_os_x)
        os_z_pv, os_z_eta = vals(pv_os_z)

        print(f"  BP   X: num={bp_x_num:.2f}, pV={bp_x_pv:.2f}, η={bp_x_eta:.2f}")
        print(f"  OSEM X: num={os_x_num:.2f}, pV={os_x_pv:.2f}, η={os_x_eta:.2f}")
        print(f"  BP   Z: num={bp_z_num:.2f}, pV={bp_z_pv:.2f}, η={bp_z_eta:.2f}")
        print(f"  OSEM Z: num={os_z_num:.2f}, pV={os_z_pv:.2f}, η={os_z_eta:.2f}")

        with open(SUMMARY_PATH, "a") as fsum:
            fsum.write(
                f"{sigma_rel*100:5.1f}  "
                f"{bp_x_num:8.2f}  {bp_x_pv:8.2f}  {bp_x_eta:5.2f}  "
                f"{os_x_num:8.2f}  {os_x_pv:8.2f}  {os_x_eta:5.2f}  "
                f"{bp_z_num:8.2f}  {bp_z_pv:8.2f}  {bp_z_eta:5.2f}  "
                f"{os_z_num:8.2f}  {os_z_pv:8.2f}  {os_z_eta:5.2f}\n"
            )

        # ----- Картинка 2x2 -----
        fig, axes = plt.subplots(2, 2, figsize=(14, 11))

        def show_image(ax, hm, title):
            im = ax.imshow(
                hm.T, origin='lower',
                extent=[x_edges[0], x_edges[-1], z_edges[0], z_edges[-1]],
                aspect='equal', cmap='hot',
            )
            cb = fig.colorbar(im, ax=ax)
            cb.set_label('Интенсивность (а.е.)', fontsize=12)
            ax.set_xlabel('X (мм)', fontsize=12)
            ax.set_ylabel('Z (мм)', fontsize=12)
            ax.set_title(title, fontsize=12)
            ax.grid(True, alpha=0.3)

        def show_profile(ax, profile, centers, pv, title):
            ax.plot(centers, profile, 'b-', linewidth=2, label='Профиль X')
            if pv is not None:
                xf = np.linspace(centers[0], centers[-1], 1000)
                ax.plot(xf, pseudo_voigt(xf, *pv), '--', color='red', linewidth=1.5,
                        label=f'pV (FWHM={pv[2]:.2f} мм, η={pv[3]:.2f})')
            ax.set_xlabel('X (мм)', fontsize=12)
            ax.set_ylabel('Интенсивность (а.е.)', fontsize=12)
            ax.set_title(title, fontsize=12)
            ax.grid(True, alpha=0.3)
            ax.set_xlim(centers[0], centers[-1])
            ax.legend(fontsize=10, loc='upper right')

        show_image(axes[0, 0], heatmap_bp,
                   f"BP   σE/E={sigma_rel*100:.0f}%, E_lost>{E_LOST_MIN} кэВ")
        show_profile(axes[0, 1], prof_x_bp, x_centers, pv_bp_x, "BP — профиль X")
        osem_title = f"OSEM ({N_ITERATIONS}×{N_SUBSETS}, σθ={SIGMA_THETA_DEG}°)"
        if POST_SMOOTH_MM > 0:
            osem_title += f", post σ={POST_SMOOTH_MM} мм"
        show_image(axes[1, 0], heatmap_osem, osem_title)
        show_profile(axes[1, 1], prof_x_os, x_centers, pv_os_x, "OSEM — профиль X")

        plt.tight_layout()
        out_png = os.path.join(OUT_DIR, f"recon_smear{int(round(sigma_rel*100)):02d}.png")
        fig.savefig(out_png, dpi=130)
        plt.close(fig)
        print(f"  сохранено: {out_png}")

    print(f"\nГотово для E_LOST_MIN={E_LOST_MIN} кэВ. Свод: {SUMMARY_PATH}")

print("\nВсе OSEM прогоны завершены.")
