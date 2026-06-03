"""Прогон Recon_gamma по нескольким значениям энергетического смеаринга.

Использует тот же кэш events_raw.npy, что и Recon_gamma.py. Cmearинг
применяется при вычислении θ из RAW энергий, поэтому пересобирать
кэш не надо.

Для каждого σ_E/E из SIGMA_LIST:
  • строится heatmap (потоковый аккумулятор);
  • из среза по точке максимума фитится pseudo-Voigt;
  • картинка + профиль X сохраняются в sweep_smear/recon_smearXX.png;
  • строка с FWHM (численный + из фита, X и Z) пишется в fwhm_summary.txt.
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")  # без интерактивного окна
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tqdm import tqdm

from recon_utils import (
    load_events,
    compute_thetas,
    cone_intersections_xz,
    HeatmapAccumulator,
)

# ---------------- Параметры ------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE_PATH = os.path.join(HERE, "events_raw.npy")
DATA_DIR = os.path.join(HERE, "data")
E_LOST_MIN_LIST = [0.0, 5.0]   # пройдёмся по обеим отсечкам в одном запуске
SIGMA_ELEC_KEV = 0.42          # электронный шум, σ в кэВ (~ 1 кэВ FWHM, клип ±1 кэВ)
SIGMA_LIST = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07]

X_RANGE = (-50.0, 50.0)
Z_RANGE = (-50.0, 50.0)
CELL_SIZE = 0.1
STEPS = 360
RNG_SEED = 12345
FILTER_MIN_COUNT = 5

X_VIEW = (-20.0, 20.0)
Z_VIEW = (-20.0, 20.0)

# ---------------- Загрузка событий (один раз) ------------------------------

events = load_events(cache_path=CACHE_PATH, data_dir=DATA_DIR)
E_lost_all = events[:, 0]; E_rem_all = events[:, 1]
Vx_all = events[:, 2]; Vy_all = events[:, 3]; Vz_all = events[:, 4]
X2_all = events[:, 5]; Y2_all = events[:, 6]; Z2_all = events[:, 7]

keep = E_lost_all <= E_rem_all
print(f"После отбора E_lost <= E_rem: {keep.sum()} / {len(events)}")
E_lost_all = E_lost_all[keep]; E_rem_all = E_rem_all[keep]
Vx_all = Vx_all[keep]; Vy_all = Vy_all[keep]; Vz_all = Vz_all[keep]
X2_all = X2_all[keep]; Y2_all = Y2_all[keep]; Z2_all = Z2_all[keep]

# ---------------- Псевдо-Войт фит -------------------------------------------

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
        print(f"  pseudo-Voigt не сошёлся: {e}")
        return None


def fwhm_numeric(profile, centers):
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


# ---------------- Цикл по отсечкам и σ -------------------------------------

for E_LOST_MIN in E_LOST_MIN_LIST:
    OUT_DIR = os.path.join(
        HERE,
        f"sweep_smear_Elost{E_LOST_MIN:g}keV_elec{SIGMA_ELEC_KEV:.2f}keV",
    )
    SUMMARY_PATH = os.path.join(OUT_DIR, "fwhm_summary.txt")
    os.makedirs(OUT_DIR, exist_ok=True)

    print("\n" + "#" * 70)
    print(f"# E_LOST_MIN = {E_LOST_MIN} кэВ   →   {OUT_DIR}")
    print("#" * 70)

    with open(SUMMARY_PATH, "w") as fsum:
        fsum.write("# σ_E/E[%]  FWHM_num_X[mm]  FWHM_num_Z[mm]  FWHM_pV_X[mm]  eta_X  FWHM_pV_Z[mm]  eta_Z\n")

    for sigma_rel in SIGMA_LIST:
        print("\n" + "=" * 70)
        print(f"σ_E/E = {sigma_rel*100:.0f}%   (E_LOST_MIN={E_LOST_MIN} кэВ)")
        print("=" * 70)

        rng = np.random.default_rng(RNG_SEED)
        theta, theta_ok = compute_thetas(E_lost_all, E_rem_all, sigma_rel=sigma_rel,
                                         sigma_elec_keV=SIGMA_ELEC_KEV,
                                         e_lost_min=E_LOST_MIN, rng=rng)
        print(f"Валидных θ (σ_elec={SIGMA_ELEC_KEV} кэВ, E_lost_m > {E_LOST_MIN} кэВ): "
              f"{theta_ok.sum()} / {len(theta)}")
        idx = np.where(theta_ok)[0]

        acc = HeatmapAccumulator(X_RANGE, Z_RANGE, CELL_SIZE)
        used = 0
        for i in tqdm(idx, desc=f"σ={sigma_rel*100:.0f}%, Emin={E_LOST_MIN}", unit="ev"):
            xs, zs = cone_intersections_xz(
                Vx_all[i], Vy_all[i], Vz_all[i],
                X2_all[i], Y2_all[i], Z2_all[i],
                theta[i], steps=STEPS,
            )
            if xs is None:
                continue
            acc.add(xs, zs)
            used += 1
        print(f"Использовано событий: {used}, max(heatmap)={acc.heatmap.max():.0f}")

        heatmap_full = acc.heatmap.copy()
        heatmap_full[heatmap_full <= FILTER_MIN_COUNT] = 0

        x_centers_full = acc.x_centers
        z_centers_full = acc.z_centers
        ix_mask = (x_centers_full >= X_VIEW[0]) & (x_centers_full <= X_VIEW[1])
        iz_mask = (z_centers_full >= Z_VIEW[0]) & (z_centers_full <= Z_VIEW[1])

        heatmap = heatmap_full[np.ix_(ix_mask, iz_mask)]
        x_centers = x_centers_full[ix_mask]
        z_centers = z_centers_full[iz_mask]
        x_edges = np.concatenate([x_centers - CELL_SIZE / 2, [x_centers[-1] + CELL_SIZE / 2]])
        z_edges = np.concatenate([z_centers - CELL_SIZE / 2, [z_centers[-1] + CELL_SIZE / 2]])

        max_z_idx = int(np.argmax(np.max(heatmap, axis=0)))
        max_x_idx = int(np.argmax(np.max(heatmap, axis=1)))
        profile_x = heatmap[:, max_z_idx]
        profile_z = heatmap[max_x_idx, :]

        pv_x = fit_pv(x_centers, profile_x)
        pv_z = fit_pv(z_centers, profile_z)

        fwhm_num_x = fwhm_numeric(profile_x, x_centers)
        fwhm_num_z = fwhm_numeric(profile_z, z_centers)

        fwhm_pv_x = pv_x[2] if pv_x is not None else float("nan")
        eta_x = pv_x[3] if pv_x is not None else float("nan")
        fwhm_pv_z = pv_z[2] if pv_z is not None else float("nan")
        eta_z = pv_z[3] if pv_z is not None else float("nan")

        print(f"  FWHM численный:   X={fwhm_num_x:.2f} мм, Z={fwhm_num_z:.2f} мм")
        print(f"  FWHM pseudo-Voigt: X={fwhm_pv_x:.2f} мм (η={eta_x:.2f}), Z={fwhm_pv_z:.2f} мм (η={eta_z:.2f})")

        with open(SUMMARY_PATH, "a") as fsum:
            fsum.write(
                f"{sigma_rel*100:5.1f}  "
                f"{fwhm_num_x:8.2f}  {fwhm_num_z:8.2f}  "
                f"{fwhm_pv_x:8.2f}  {eta_x:5.2f}  "
                f"{fwhm_pv_z:8.2f}  {eta_z:5.2f}\n"
            )

        # ----- картинка -----
        fig, (ax_img, ax_prof) = plt.subplots(1, 2, figsize=(14, 5))

        im = ax_img.imshow(
            heatmap.T, origin='lower',
            extent=[x_edges[0], x_edges[-1], z_edges[0], z_edges[-1]],
            aspect='equal', cmap='hot',
        )
        cbar = fig.colorbar(im, ax=ax_img)
        cbar.set_label('Интенсивность (а.е.)', fontsize=14)
        cbar.ax.tick_params(labelsize=12)
        ax_img.set_xlabel('X (мм)', fontsize=14)
        ax_img.set_ylabel('Z (мм)', fontsize=14)
        ax_img.set_title(
            f"σE/E = {sigma_rel*100:.0f}%, σ_elec = {SIGMA_ELEC_KEV} кэВ, E_lost > {E_LOST_MIN} кэВ",
            fontsize=12,
        )
        ax_img.tick_params(axis='both', which='major', labelsize=12)
        ax_img.grid(True, alpha=0.3)

        ax_prof.plot(x_centers, profile_x, 'b-', linewidth=2, label='Профиль X')
        if pv_x is not None:
            xf = np.linspace(x_centers[0], x_centers[-1], 1000)
            ax_prof.plot(xf, pseudo_voigt(xf, *pv_x), '--', color='red', linewidth=1.5,
                         label=f"pV: FWHM={pv_x[2]:.2f} мм, η={pv_x[3]:.2f}")
        ax_prof.set_xlabel('X (мм)', fontsize=14)
        ax_prof.set_ylabel('Интенсивность (а.е.)', fontsize=14)
        ax_prof.set_title(f"σE/E = {sigma_rel*100:.0f}%   FWHM_num(X) = {fwhm_num_x:.2f} мм", fontsize=14)
        ax_prof.tick_params(axis='both', which='major', labelsize=12)
        ax_prof.grid(True, alpha=0.3)
        ax_prof.set_xlim(x_centers[0], x_centers[-1])
        ax_prof.legend(fontsize=11, loc='upper right')

        plt.tight_layout()
        out_png = os.path.join(OUT_DIR, f"recon_smear{int(round(sigma_rel*100)):02d}.png")
        fig.savefig(out_png, dpi=130)
        plt.close(fig)
        print(f"  сохранено: {out_png}")

    print(f"\nГотово для E_LOST_MIN={E_LOST_MIN} кэВ. Свод: {SUMMARY_PATH}")

print("\nВсе прогоны завершены.")
