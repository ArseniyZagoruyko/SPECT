"""Общие функции для обработки данных моделирования Комптон-камеры.

Подгружает все пары файлов deposits_FirstDet_<seed>.txt /
GammaCoordinates_FirstDet_<seed>.txt из каталога data/ и итеративно
возвращает события. Каждый seed соответствует отдельному запуску на
узле суперкомпьютера и нумерация eventID независима, поэтому склейка
через dict[eventID] = ... приведёт к молчаливой потере событий
(один и тот же eventID встречается во многих файлах). Здесь
используется построчная итерация в пределах одной пары, что делает
коллизии eventID между файлами невозможными.
"""

import os
import re
from glob import glob

import numpy as np

DATA_DIR_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

_SEED_RE = re.compile(r"deposits_FirstDet_(\d+)\.txt$")


def discover_pairs(data_dir=DATA_DIR_DEFAULT):
    """Возвращает список (seed, deposits_path, coords_path) для всех файлов в data/.

    Пары формируются по seed в имени файла. Если для deposits нет парного
    GammaCoordinates с тем же seed — такая пара пропускается с предупреждением.
    """
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
    """Парная построчная итерация: для каждой строки deposits ищет такую же
    eventID в coordinates ВНУТРИ ОДНОГО файла.

    В пределах одного файла моделирования eventID уникален. Возвращает кортежи
    (seed_unused, eventID, E_lost, E_remaines, Vx, Vy, Vz, X2, Y2, Z2).

    Если в файле deposits eventID встречается, а в coordinates — нет (записи
    могут расходиться, если триггер записи в SteppingAction различался),
    такая строка пропускается.
    """
    # Coordinates строятся в dict в пределах одного файла — это безопасно,
    # потому что внутри файла eventID уникальны.
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


def iter_all_events(data_dir=DATA_DIR_DEFAULT, max_lines_per_file=None):
    """Итератор по событиям из всех пар файлов в data_dir.

    Коллизий eventID между файлами не возникает — каждое событие из каждого
    файла отдаётся независимо. Возвращает те же поля, что iter_events_from_pair,
    дополненные seed на первой позиции:
    (seed, eventID, E_lost, E_remaines, Vx, Vy, Vz, X2, Y2, Z2).
    """
    pairs = discover_pairs(data_dir)
    print(f"[recon_utils] Найдено пар файлов: {len(pairs)} (каталог {data_dir})")
    for seed, dep_path, coord_path in pairs:
        for tup in iter_events_from_pair(dep_path, coord_path, max_lines=max_lines_per_file):
            yield (seed,) + tup


# ---- Физика конусов Комптона ----------------------------------------------

def find_theta(E_remaines, E_lost, smear=True, e_lost_min=0.5):
    if E_lost <= e_lost_min:
        return None
    if smear:
        E_lost_m = np.random.normal(E_lost, E_lost * 0.05)
        E_rem_m = np.random.normal(E_remaines, E_remaines * 0.05)
    else:
        E_lost_m = E_lost
        E_rem_m = E_remaines
    E_init = E_lost_m + E_rem_m
    denom = E_init * (E_init - E_lost_m)
    if denom == 0:
        return None
    cos_theta = 1 - (511 * E_lost_m) / denom
    if not -1 <= cos_theta <= 1:
        return None
    return float(np.degrees(np.arccos(cos_theta)))


def find_vector_P(Vx, Vy, Vz, X2, Y2, Z2):
    return np.array([X2 - Vx, Y2 - Vy, Z2 - Vz])


def find_intersection_with_XZ_plane(Vx, Vy, Vz, Px, Py, Pz):
    if abs(Py) < 1e-10:
        return None
    t = -Vy / Py
    return np.array([Vx + t * Px, 0.0, Vz + t * Pz])


def find_vector_h(Vx, Vy, Vz, Hx, Hy, Hz):
    return np.array([Hx - Vx, Hy - Vy, Hz - Vz])


def find_vector_d(h, theta):
    norm = np.linalg.norm(h)
    if norm < 1e-10:
        return None
    h_norm = h / norm
    hx, hy, hz = h_norm
    perp = np.array([-hz, 0.0, hx])
    pn = np.linalg.norm(perp)
    if pn < 1e-10:
        return None
    perp_norm = perp / pn
    return h_norm * np.cos(np.radians(theta)) + perp_norm * np.sin(np.radians(theta))


def rotation_matrix(axis, theta_deg):
    norm = np.linalg.norm(axis)
    if norm < 1e-10:
        return np.eye(3)
    axis = axis / norm
    a = np.cos(np.radians(theta_deg / 2))
    b, c, d = -axis * np.sin(np.radians(theta_deg / 2))
    return np.array([
        [a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
        [2 * (b * c + a * d), a * a + c * c - b * b - d * d, 2 * (c * d - a * b)],
        [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - b * b - c * c],
    ])


def rotate_and_intersect(V, d, h, steps=360):
    pts = []
    for angle in range(steps):
        R = rotation_matrix(h, angle)
        dr = R @ d
        if abs(dr[1]) < 1e-10:
            continue
        t = -V[1] / dr[1]
        p = V + t * dr
        if np.isfinite(p).all():
            pts.append(p)
    return np.array(pts) if pts else np.empty((0, 3))
