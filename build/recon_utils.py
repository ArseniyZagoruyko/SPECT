"""Общие функции для обработки данных моделирования Комптон-камеры.

Подгружает все пары файлов deposits_FirstDet_<seed>.txt /
GammaCoordinates_FirstDet_<seed>.txt из каталога data/ и итеративно
возвращает события. Каждый seed соответствует отдельному запуску на
узле суперкомпьютера и нумерация eventID независима, поэтому склейка
через dict[eventID] = ... приведёт к молчаливой потере событий
(один и тот же eventID встречается во многих файлах). Здесь
используется построчная итерация в пределах одной пары, что делает
коллизии eventID между файлами невозможными.

Энергетический смеаринг (моделирование ухудшения энергетического
разрешения детекторов) применяется на этапе реконструкции, а не при
кэшировании сырых данных в .npy — это позволяет переиспользовать дамп
для разных значений σ.
"""

import os
import re
from glob import glob

import numpy as np
from tqdm import tqdm

DATA_DIR_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
CACHE_PATH_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "events_raw.npy")

_SEED_RE = re.compile(r"deposits_FirstDet_(\d+)\.txt$")


# =============================================================================
# Чтение данных
# =============================================================================

def discover_pairs(data_dir=DATA_DIR_DEFAULT):
    deposit_files = sorted(glob(os.path.join(data_dir, "deposits_FirstDet_*.txt")))
    pairs = []
    missing = 0
    for dep_path in deposit_files:
        m = _SEED_RE.search(os.path.basename(dep_path))
        if not m:
            continue
        seed = m.group(1)
        coord_path = os.path.join(data_dir, f"GammaCoordinates_FirstDet_{seed}.txt")
        if not os.path.isfile(coord_path):
            missing += 1
            continue
        pairs.append((seed, dep_path, coord_path))
    if missing:
        print(f"[recon_utils] Пропущено {missing} файл(ов) deposits без парного GammaCoordinates")
    return pairs


def iter_events_from_pair(dep_path, coord_path, max_lines=None):
    """Парная построчная итерация для одной пары файлов.
    Возвращает кортежи (eventID, E_lost, E_remaines, Vx, Vy, Vz, X2, Y2, Z2).
    """
    coords = {}
    with open(coord_path, "r") as fc:
        for i, line in enumerate(fc):
            if max_lines is not None and i >= max_lines:
                break
            parts = line.strip().split()
            if len(parts) < 7:
                continue
            try:
                ev = int(parts[0])
                vals = tuple(map(float, parts[1:7]))
                coords[ev] = vals
            except ValueError:
                continue

    with open(dep_path, "r") as fd:
        for i, line in enumerate(fd):
            if max_lines is not None and i >= max_lines:
                break
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                ev = int(parts[0])
                E_lost = float(parts[1])
                E_remaines = float(parts[2])
            except ValueError:
                continue
            if ev not in coords:
                continue
            Vx, Vy, Vz, X2, Y2, Z2 = coords[ev]
            yield ev, E_lost, E_remaines, Vx, Vy, Vz, X2, Y2, Z2


# =============================================================================
# Кэш в .npy
# =============================================================================

def build_cache(data_dir=DATA_DIR_DEFAULT, cache_path=CACHE_PATH_DEFAULT, dtype=np.float32):
    """Один проход по всем txt-парам, сохранение совпадений в один .npy.

    Сохраняются ИСХОДНЫЕ значения энергий, без смеаринга — смеаринг применяется
    при реконструкции. Колонки:
        [E_lost, E_remaines, Vx, Vy, Vz, X2, Y2, Z2]
    eventID и seed не сохраняются — для обратного проецирования они не нужны,
    и опускание экономит ~25% памяти/диска.

    Для 3.8M событий .npy займёт ~120 МБ (float32) — это умещается в ОЗУ.
    """
    pairs = discover_pairs(data_dir)
    print(f"[build_cache] Пар файлов: {len(pairs)}")

    chunks = []
    total = 0
    for seed, dep_path, coord_path in tqdm(pairs, desc="Чтение", unit="пара"):
        buf = []
        for ev, E_lost, E_rem, Vx, Vy, Vz, X2, Y2, Z2 in iter_events_from_pair(dep_path, coord_path):
            buf.append((E_lost, E_rem, Vx, Vy, Vz, X2, Y2, Z2))
        if buf:
            chunks.append(np.asarray(buf, dtype=dtype))
            total += len(buf)

    if not chunks:
        raise RuntimeError("Не найдено ни одного события")
    arr = np.concatenate(chunks, axis=0)
    np.save(cache_path, arr)
    print(f"[build_cache] Сохранено {total} событий -> {cache_path} ({arr.nbytes/1e6:.1f} МБ)")
    return arr


def load_events(cache_path=CACHE_PATH_DEFAULT, data_dir=DATA_DIR_DEFAULT, rebuild=False):
    """Возвращает массив событий формы (N, 8). При отсутствии кэша строит его."""
    if rebuild or not os.path.isfile(cache_path):
        return build_cache(data_dir=data_dir, cache_path=cache_path)
    arr = np.load(cache_path)
    print(f"[load_events] Загружено {len(arr)} событий из {cache_path}")
    return arr


# =============================================================================
# Физика — векторизованная проекция конуса
# =============================================================================

def compute_thetas(E_lost, E_rem, sigma_rel=0.05, sigma_elec_keV=0.0,
                   elec_clip_keV=1.0, e_lost_min=0.5, rng=None):
    """Векторизованный расчёт θ для массива событий.

    sigma_rel       — статистическая (мультипликативная) часть разрешения,
                      σ_stat = sigma_rel · E. Например, 0.05 для 5%.
                      Применяется как мультипликативный гаусс с σ = sigma_rel · E
                      на обе энергии независимо.
    sigma_elec_keV  — σ электронного (аддитивного) шума в кэВ.
                      Не зависит от E. Реализован как N(0, sigma_elec_keV),
                      потом КЛИПИРУЕТСЯ в [-elec_clip_keV, +elec_clip_keV] —
                      это срезает хвосты гаусса (физически: реальный baseline
                      не уходит произвольно далеко, его ограничивают
                      разрядность АЦП и pile-up rejection).
    elec_clip_keV   — граница клиппинга шумовой добавки, кэВ.
    Установи sigma_rel=0 и sigma_elec_keV=0, чтобы отключить смеаринг.
    Возвращает (theta_deg, mask): theta только для тех событий, где mask=True.
    """
    if rng is None:
        rng = np.random.default_rng()
    E_lost = np.asarray(E_lost, dtype=np.float64)
    E_rem = np.asarray(E_rem, dtype=np.float64)

    if sigma_rel > 0:
        E_lost_m = rng.normal(E_lost, sigma_rel * E_lost)
        E_rem_m = rng.normal(E_rem, sigma_rel * E_rem)
    else:
        E_lost_m = E_lost.copy()
        E_rem_m = E_rem.copy()

    if sigma_elec_keV > 0:
        delta_lost = rng.normal(0.0, sigma_elec_keV, size=E_lost.shape)
        delta_rem = rng.normal(0.0, sigma_elec_keV, size=E_rem.shape)
        np.clip(delta_lost, -elec_clip_keV, elec_clip_keV, out=delta_lost)
        np.clip(delta_rem, -elec_clip_keV, elec_clip_keV, out=delta_rem)
        E_lost_m = E_lost_m + delta_lost
        E_rem_m = E_rem_m + delta_rem

    E_init = E_lost_m + E_rem_m
    denom = E_init * (E_init - E_lost_m)
    with np.errstate(divide='ignore', invalid='ignore'):
        cos_theta = 1.0 - (511.0 * E_lost_m) / denom
    # Отсечка применяется к ИЗМЕРЕННОЙ энергии: реальный триггер/дискриминатор
    # видит только смеаренный сигнал. Дополнительно требуем E_rem_m > 0 — события,
    # где электронный шум «утопил» baseline ниже нуля, в формулу Комптона ставить
    # нельзя (нефизичный знаменатель). При sigma_rel=sigma_elec=0 это сводится
    # к старому поведению (E_lost_m == E_lost).
    mask = (
        (E_lost_m > e_lost_min) &
        (E_rem_m > 0) &
        np.isfinite(cos_theta) &
        (cos_theta >= -1.0) & (cos_theta <= 1.0)
    )
    theta = np.full_like(cos_theta, np.nan)
    theta[mask] = np.degrees(np.arccos(cos_theta[mask]))
    return theta, mask


def cone_intersections_xz(Vx, Vy, Vz, X2, Y2, Z2, theta_deg, steps=360):
    """Векторизованная версия rotate_and_intersect для ОДНОГО события.

    Сначала строится h (проекция вектора V→второй детектор на плоскость XZ
    через продолжение P до XZ), затем перпендикуляр в XZ-плоскости, отклонение
    d на угол θ от h, и потом d вращается вокруг h на 0..steps-1 градусов.
    Каждое повёрнутое направление продолжается до плоскости y=0; возвращаются
    (x, z) точек пересечения.

    Возвращает (x_arr, z_arr) или (None, None) если событие отбраковано
    (вырожденная геометрия).
    """
    Px = X2 - Vx
    Py = Y2 - Vy
    Pz = Z2 - Vz
    if abs(Py) < 1e-10:
        return None, None
    t0 = -Vy / Py
    Hx = Vx + t0 * Px
    Hz = Vz + t0 * Pz

    hx = Hx - Vx
    hy = -Vy
    hz = Hz - Vz
    h = np.array([hx, hy, hz], dtype=np.float64)
    hn = np.linalg.norm(h)
    if hn < 1e-10:
        return None, None
    h /= hn

    perp = np.array([-h[2], 0.0, h[0]], dtype=np.float64)
    pn = np.linalg.norm(perp)
    if pn < 1e-10:
        return None, None
    perp /= pn

    th = np.radians(theta_deg)
    d = h * np.cos(th) + perp * np.sin(th)

    angles = np.radians(np.arange(steps, dtype=np.float64))
    ca = np.cos(angles)
    sa = np.sin(angles)

    # Поворот d вокруг оси h на каждый угол (формула Родрига).
    # d, h фиксированы; (h·d) — скаляр, (h×d) — постоянный вектор.
    h_dot_d = float(np.dot(h, d))
    h_cross_d = np.cross(h, d)
    # d_rot(a) = d*cos a + (h × d)*sin a + h*(h·d)*(1 - cos a)
    drx = d[0] * ca + h_cross_d[0] * sa + h[0] * h_dot_d * (1 - ca)
    dry = d[1] * ca + h_cross_d[1] * sa + h[1] * h_dot_d * (1 - ca)
    drz = d[2] * ca + h_cross_d[2] * sa + h[2] * h_dot_d * (1 - ca)

    # Пересечение луча (V + t*d_rot) с плоскостью y=0: t = -Vy / dry
    valid = np.abs(dry) > 1e-10
    if not np.any(valid):
        return None, None
    t = -Vy / dry[valid]
    x_out = Vx + t * drx[valid]
    z_out = Vz + t * drz[valid]
    finite = np.isfinite(x_out) & np.isfinite(z_out)
    return x_out[finite], z_out[finite]


# =============================================================================
# Потоковый аккумулятор heatmap
# =============================================================================

class HeatmapAccumulator:
    """Аккумулирует точки пересечения сразу в 2D-гистограмму,
    не храня сами точки. Память O(nx*nz), постоянная.
    """
    def __init__(self, x_range, z_range, cell_size):
        self.x_edges = np.arange(x_range[0], x_range[1] + cell_size, cell_size)
        self.z_edges = np.arange(z_range[0], z_range[1] + cell_size, cell_size)
        self.nx = len(self.x_edges) - 1
        self.nz = len(self.z_edges) - 1
        self.heatmap = np.zeros((self.nx, self.nz), dtype=np.float64)
        self.x_min = self.x_edges[0]
        self.z_min = self.z_edges[0]
        self.inv_cell = 1.0 / cell_size

    def add(self, x_arr, z_arr):
        if x_arr is None or len(x_arr) == 0:
            return
        ix = np.floor((x_arr - self.x_min) * self.inv_cell).astype(np.int64)
        iz = np.floor((z_arr - self.z_min) * self.inv_cell).astype(np.int64)
        mask = (ix >= 0) & (ix < self.nx) & (iz >= 0) & (iz < self.nz)
        if not np.any(mask):
            return
        np.add.at(self.heatmap, (ix[mask], iz[mask]), 1.0)

    @property
    def x_centers(self):
        return 0.5 * (self.x_edges[:-1] + self.x_edges[1:])

    @property
    def z_centers(self):
        return 0.5 * (self.z_edges[:-1] + self.z_edges[1:])
