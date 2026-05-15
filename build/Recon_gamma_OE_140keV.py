#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import math
import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
from tqdm import tqdm
from scipy import stats

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input files
DEPOSITS_FILE = "/home/zas/CERN/SPECT/build/deposits(1detector).txt"
COORDS_FILE = "/home/zas/CERN/SPECT/build/GammaCoordinates(1detector).txt"

# Physical constants
ELECTRON_MASS_KEV = 511.0
E0_KEV = 140.0  # Initial gamma energy

# Reconstruction grid
CELL_SIZE = 0.1  # mm
X_RANGE = (-50.0, 50.0)  # mm
Z_RANGE = (-50.0, 50.0)  # mm
GRID_Y = 0.0  # mm, reconstruction plane

# Energy cuts
MIN_E_LOST = 0.001  # keV
MAX_E_LOST = 140.0  # keV

# OE-specific parameters (KEY DIFFERENCES FROM BACKPROJECTION)
SAMPLES_PER_EVENT = 50  # Number of probabilistic samples per event (increased for ideal system)
ENERGY_RESOLUTION = 0.1  # 0.1% energy resolution - nearly ideal system
POSITION_RESOLUTION = 0.01  # mm, detector position resolution (σ) - nearly ideal
KDE_BANDWIDTH = 0.1  # mm, Gaussian kernel bandwidth - very tight for maximum sharpness
KDE_CUTOFF_SIGMA = 1.0  # Cutoff kernel at this many sigmas (for speed)
USE_DISCRETE_MODE = False  # TRUE OE without KDE smoothing (like backprojection but with sampling)
USE_ENERGY_WEIGHTING = True  # Weight by photopeak proximity
PHOTOPEAK_WINDOW = 2.0  # keV, half-width around 140 keV for weighting

# Processing limits
MAX_EVENTS_PROCESS = 10000  # None = all events, set to int for testing

# Output
SAVE_PREFIX = "oe_true_140keV"

# =============================================================================
# DATA LOADING
# =============================================================================

def load_deposits(path):
    """Load energy deposits: eventID, E_lost, E_remains"""
    data = {}
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            try:
                event_id = int(parts[0])
                e_lost = float(parts[1])
                e_rem = float(parts[2])
                data[event_id] = (e_lost, e_rem)
            except ValueError:
                continue
    return data


def load_coordinates(path):
    """Load interaction coordinates: eventID, Vx, Vy, Vz, X2, Y2, Z2"""
    data = {}
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 7:
                continue
            try:
                event_id = int(parts[0])
                vals = tuple(float(x) for x in parts[1:7])
                data[event_id] = vals
            except ValueError:
                continue
    return data


# =============================================================================
# COMPTON KINEMATICS WITH UNCERTAINTY
# =============================================================================

def compton_angle_with_uncertainty(e_lost, e_rem, n_samples=50):
    """
    Calculate Compton angle distribution accounting for energy resolution.

    Returns:
        theta_samples: array of sampled angles (degrees)
        weights: corresponding probability weights
    """
    # Convert FWHM to sigma
    sigma_rel = ENERGY_RESOLUTION / 2.355

    # Sample energies from Gaussian distributions
    e_lost_samples = np.random.normal(e_lost, e_lost * sigma_rel, n_samples)
    e_rem_samples = np.random.normal(e_rem, e_rem * sigma_rel, n_samples)

    # Keep only physical values
    mask = (e_lost_samples > 0) & (e_rem_samples > 0)
    e_lost_samples = e_lost_samples[mask]
    e_rem_samples = e_rem_samples[mask]

    if len(e_lost_samples) == 0:
        return None, None

    # Calculate angles
    E_initial = e_lost_samples + e_rem_samples
    mask2 = (E_initial > e_lost_samples) & (E_initial > 0)

    if np.sum(mask2) == 0:
        return None, None

    e_lost_valid = e_lost_samples[mask2]
    E_initial_valid = E_initial[mask2]

    cos_theta = 1.0 - (ELECTRON_MASS_KEV * e_lost_valid) / (E_initial_valid * (E_initial_valid - e_lost_valid))
    mask3 = (cos_theta >= -1.0) & (cos_theta <= 1.0)

    if np.sum(mask3) == 0:
        return None, None

    theta_rad = np.arccos(np.clip(cos_theta[mask3], -1.0, 1.0))
    theta_deg = np.degrees(theta_rad)

    # Uniform weights for now (could use Klein-Nishina here)
    weights = np.ones(len(theta_deg)) / len(theta_deg)

    return theta_deg, weights


# =============================================================================
# ORIGIN ENSEMBLE CORE: PROBABILISTIC SAMPLING
# =============================================================================

def sample_source_positions(V, X2, e_lost, e_rem, n_samples=SAMPLES_PER_EVENT):
    """
    TRUE ORIGIN ENSEMBLE METHOD:

    For a given event, generate an ensemble of probable source positions
    by sampling the measurement uncertainties.

    This is NOT backprojection! We:
    1. Sample angle θ from energy uncertainty distribution
    2. Sample azimuth φ uniformly on cone
    3. Sample positions V and X2 from detector position uncertainty
    4. For each sample, compute one potential source position
    5. Return cloud of positions with weights

    Returns:
        source_samples: array of (x, z) source position samples
        weights: probability weight for each sample
    """
    # Sample angles accounting for energy uncertainty
    theta_samples, theta_weights = compton_angle_with_uncertainty(
        e_lost, e_rem, n_samples=max(10, n_samples // 10)
    )

    if theta_samples is None:
        return None, None

    # Sample position uncertainties
    rng = np.random.default_rng()

    source_samples = []
    weights = []

    for _ in range(n_samples):
        # Sample one angle from distribution
        theta_idx = rng.choice(len(theta_samples), p=theta_weights)
        theta_deg = theta_samples[theta_idx]
        theta_rad = math.radians(theta_deg)

        # Sample position uncertainties (Gaussian around measured positions)
        if POSITION_RESOLUTION > 0:
            V_sample = V + rng.normal(0, POSITION_RESOLUTION, 3)
            X2_sample = X2 + rng.normal(0, POSITION_RESOLUTION, 3)
        else:
            V_sample = V.copy()
            X2_sample = X2.copy()

        # Scattered direction
        d_vec = X2_sample - V_sample
        d_norm = np.linalg.norm(d_vec)
        if d_norm < 1e-9:
            continue
        d_hat = d_vec / d_norm

        # Sample azimuth uniformly
        phi = rng.uniform(0, 2 * math.pi)

        # Create orthonormal basis for cone
        if abs(d_hat[0]) < 0.9:
            ref = np.array([1.0, 0.0, 0.0])
        else:
            ref = np.array([0.0, 1.0, 0.0])

        b1 = np.cross(d_hat, ref)
        b1_norm = np.linalg.norm(b1)
        if b1_norm < 1e-12:
            continue
        b1 = b1 / b1_norm
        b2 = np.cross(d_hat, b1)

        # Incident direction on cone
        u_hat = (math.cos(theta_rad) * d_hat +
                 math.sin(theta_rad) * (math.cos(phi) * b1 + math.sin(phi) * b2))

        # Backproject to find source position at Y=GRID_Y
        if abs(u_hat[1]) < 1e-9:
            continue

        t = (V_sample[1] - GRID_Y) / u_hat[1]
        if t <= 0:
            continue

        source_pos = V_sample - t * u_hat

        # Store sample
        source_samples.append((source_pos[0], source_pos[2]))
        weights.append(1.0)  # Equal weight for now

    if len(source_samples) == 0:
        return None, None

    return np.array(source_samples), np.array(weights) / np.sum(weights)


def event_quality_weight(e_lost, e_rem):
    """
    Weight event by proximity to photopeak and energy balance.
    This is a soft, probabilistic weight (unlike backprojection's hard cuts).
    """
    if not USE_ENERGY_WEIGHTING:
        return 1.0

    E_total = e_lost + e_rem

    # Gaussian weight around 140 keV photopeak
    sigma_window = PHOTOPEAK_WINDOW / 2.355
    w_energy = math.exp(-0.5 * ((E_total - E0_KEV) / sigma_window) ** 2)

    # Prefer events with reasonable energy split
    ratio = min(e_lost, e_rem) / max(e_lost, e_rem)
    w_balance = math.sqrt(ratio)

    return w_energy * w_balance


# =============================================================================
# KERNEL DENSITY ESTIMATION
# =============================================================================

def accumulate_kde(grid_x, grid_z, samples, weights, bandwidth=KDE_BANDWIDTH):
    """
    Accumulate probability density.

    Two modes:
    1. DISCRETE MODE (USE_DISCRETE_MODE=True): Like backprojection - count hits in bins
    2. KDE MODE (USE_DISCRETE_MODE=False): Gaussian smoothing

    FULLY VECTORIZED for speed.
    """
    density = np.zeros_like(grid_x, dtype=np.float64)

    if len(samples) == 0:
        return density

    # Convert samples to array
    samples_arr = np.array(samples)  # Shape: (N, 2)
    weights_arr = np.array(weights)  # Shape: (N,)

    # Grid parameters
    dx_grid = grid_x[0, 1] - grid_x[0, 0] if grid_x.shape[1] > 1 else 1.0
    dz_grid = grid_z[1, 0] - grid_z[0, 0] if grid_z.shape[0] > 1 else 1.0
    x_min, x_max = grid_x[0, 0], grid_x[0, -1]
    z_min, z_max = grid_z[0, 0], grid_z[-1, 0]

    # Filter samples outside grid bounds
    margin = bandwidth * KDE_CUTOFF_SIGMA if not USE_DISCRETE_MODE else 0
    mask = ((samples_arr[:, 0] >= x_min - margin) &
            (samples_arr[:, 0] <= x_max + margin) &
            (samples_arr[:, 1] >= z_min - margin) &
            (samples_arr[:, 1] <= z_max + margin))

    samples_arr = samples_arr[mask]
    weights_arr = weights_arr[mask]

    if len(samples_arr) == 0:
        return density

    if USE_DISCRETE_MODE:
        # DISCRETE MODE: Histogram binning (like backprojection)
        # Find which bin each sample belongs to
        i_indices = np.floor((samples_arr[:, 0] - x_min) / dx_grid).astype(int)
        j_indices = np.floor((samples_arr[:, 1] - z_min) / dz_grid).astype(int)

        # Clip to grid bounds
        i_indices = np.clip(i_indices, 0, grid_x.shape[1] - 1)
        j_indices = np.clip(j_indices, 0, grid_z.shape[0] - 1)

        # Accumulate weights in bins
        for j, i, w in zip(j_indices, i_indices, weights_arr):
            density[j, i] += w
    else:
        # KDE MODE: Gaussian smoothing
        cutoff_radius = KDE_CUTOFF_SIGMA * bandwidth
        inv_bandwidth_sq = 1.0 / (bandwidth**2)
        norm_const = 1.0 / (2 * math.pi * bandwidth**2)

        # Vectorized processing
        grid_3d = np.stack([grid_x, grid_z], axis=-1)  # Shape: (M, N, 2)
        samples_3d = samples_arr[:, np.newaxis, np.newaxis, :]  # Shape: (K, 1, 1, 2)

        # Compute all distances at once
        diff = grid_3d - samples_3d  # Shape: (K, M, N, 2)
        dist_sq = np.sum(diff**2, axis=-1)  # Shape: (K, M, N)

        # Apply cutoff mask
        cutoff_mask = dist_sq <= (cutoff_radius**2)

        # Compute kernels with cutoff
        kernels = np.zeros_like(dist_sq)
        kernels[cutoff_mask] = np.exp(-0.5 * dist_sq[cutoff_mask] * inv_bandwidth_sq)

        # Apply weights and normalization
        kernels *= (weights_arr[:, np.newaxis, np.newaxis] * norm_const)

        # Sum all kernels
        density = np.sum(kernels, axis=0)

    return density


# =============================================================================
# MAIN RECONSTRUCTION
# =============================================================================

def reconstruct_oe():
    """
    Perform TRUE Origin Ensemble reconstruction.
    """
    print("=" * 70)
    print("TRUE ORIGIN ENSEMBLE (OE) RECONSTRUCTION FOR 140 keV")
    print("=" * 70)
    mode_str = "DISCRETE (no smoothing)" if USE_DISCRETE_MODE else "KDE (Gaussian smoothing)"
    print(f"Method: Probabilistic sampling - {mode_str}")
    print(f"Samples per event: {SAMPLES_PER_EVENT}")
    print(f"Energy resolution: {ENERGY_RESOLUTION * 100:.3f}%")
    print(f"Position resolution: {POSITION_RESOLUTION:.2f} mm")
    if not USE_DISCRETE_MODE:
        print(f"KDE bandwidth: {KDE_BANDWIDTH:.2f} mm")
    print("=" * 70)

    # Load data
    deposits = load_deposits(DEPOSITS_FILE)
    coords = load_coordinates(COORDS_FILE)

    common_ids = sorted(set(deposits.keys()) & set(coords.keys()))
    print(f"\nEvents in deposits: {len(deposits)}")
    print(f"Events in coordinates: {len(coords)}")
    print(f"Common events: {len(common_ids)}")

    # Create grid
    bins_x = np.arange(X_RANGE[0], X_RANGE[1] + CELL_SIZE, CELL_SIZE)
    bins_z = np.arange(Z_RANGE[0], Z_RANGE[1] + CELL_SIZE, CELL_SIZE)
    centers_x = (bins_x[:-1] + bins_x[1:]) / 2.0
    centers_z = (bins_z[:-1] + bins_z[1:]) / 2.0
    grid_x, grid_z = np.meshgrid(centers_x, centers_z, indexing='xy')

    # Initialize density map
    density_map = np.zeros_like(grid_x, dtype=np.float64)
    total_weight = 0.0

    # Process events
    events_processed = 0
    events_used = 0

    if MAX_EVENTS_PROCESS is not None:
        common_ids = common_ids[:MAX_EVENTS_PROCESS]

    print(f"\nProcessing {len(common_ids)} events...")

    for eid in tqdm(common_ids, desc="OE reconstruction"):
        e_lost, e_rem = deposits[eid]

        # Filter events
        if e_lost < MIN_E_LOST or e_lost > MAX_E_LOST:
            continue
        if e_lost > e_rem:
            continue

        Vx, Vy, Vz, X2, Y2, Z2 = coords[eid]
        V = np.array([Vx, Vy, Vz])
        P2 = np.array([X2, Y2, Z2])

        # Generate origin ensemble (probabilistic samples)
        samples, sample_weights = sample_source_positions(
            V, P2, e_lost, e_rem, n_samples=SAMPLES_PER_EVENT
        )

        if samples is None:
            continue

        # Event quality weight
        event_weight = event_quality_weight(e_lost, e_rem)

        # Accumulate KDE
        event_density = accumulate_kde(
            grid_x, grid_z, samples,
            sample_weights * event_weight,
            bandwidth=KDE_BANDWIDTH
        )

        density_map += event_density
        total_weight += event_weight
        events_used += 1
        events_processed += 1

    print(f"\nEvents processed: {events_processed}")
    print(f"Events with contribution: {events_used}")

    # Normalize
    if total_weight > 0:
        density_map /= total_weight

    return density_map, grid_x, grid_z, bins_x, bins_z


# =============================================================================
# METRICS AND VISUALIZATION
# =============================================================================

def calculate_metrics(recon_map, grid_x, grid_z, bins_x, bins_z):
    """Calculate reconstruction quality metrics."""
    # Peak
    max_idx = np.unravel_index(np.argmax(recon_map), recon_map.shape)
    peak_x = grid_x[max_idx]
    peak_z = grid_z[max_idx]
    peak_val = float(recon_map[max_idx])

    # Profiles
    x_profile = recon_map.sum(axis=0)
    z_profile = recon_map.sum(axis=1)

    def fwhm_from_profile(profile, edges):
        if profile.size == 0 or np.max(profile) <= 0:
            return None
        half = np.max(profile) * 0.5
        idx = np.where(profile >= half)[0]
        if idx.size == 0:
            return None
        left = edges[idx[0]]
        right = edges[idx[-1] + 1] if (idx[-1] + 1) < len(edges) else edges[-1]
        return float(right - left)

    fwhm_x = fwhm_from_profile(x_profile, bins_x)
    fwhm_z = fwhm_from_profile(z_profile, bins_z)

    return {
        'peak_position': (float(peak_x), float(peak_z)),
        'peak_value': peak_val,
        'fwhm_x': fwhm_x,
        'fwhm_z': fwhm_z,
    }


def save_results(recon_map, grid_x, grid_z, metrics):
    """Save reconstruction results."""
    extent = [grid_x.min(), grid_x.max(), grid_z.min(), grid_z.max()]

    # Plot
    plt.figure(figsize=(8, 7))
    plt.imshow(recon_map.T, origin='lower', extent=extent,
               aspect='equal', cmap='hot', interpolation='bilinear')
    plt.colorbar(label='Probability Density')
    plt.xlabel('X (mm)', fontsize=14)
    plt.ylabel('Z (mm)', fontsize=14)
    plt.title('True OE Reconstruction (140 keV)', fontsize=16)
    plt.tight_layout()
    plt.savefig(f'{SAVE_PREFIX}_raw.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Save JPG
    if np.max(recon_map) > 0:
        norm = (recon_map.T / np.max(recon_map.T) * 255).astype(np.uint8)
    else:
        norm = np.zeros_like(recon_map.T, dtype=np.uint8)
    imageio.imwrite(f'{SAVE_PREFIX}_final.jpg', norm)

    # Print results
    print("\n" + "=" * 70)
    print("RECONSTRUCTION RESULTS")
    print("=" * 70)
    print(f"Peak position: X={metrics['peak_position'][0]:.2f} mm, "
          f"Z={metrics['peak_position'][1]:.2f} mm")
    print(f"Peak value: {metrics['peak_value']:.6f}")
    if metrics['fwhm_x'] and metrics['fwhm_z']:
        print(f"FWHM X: {metrics['fwhm_x']:.2f} mm")
        print(f"FWHM Z: {metrics['fwhm_z']:.2f} mm")
    print("=" * 70)
    print(f"Images saved: {SAVE_PREFIX}_raw.png, {SAVE_PREFIX}_final.jpg")
    print("=" * 70)


# =============================================================================
# MAIN
# =============================================================================

def main():
    density_map, grid_x, grid_z, bins_x, bins_z = reconstruct_oe()
    metrics = calculate_metrics(density_map, grid_x, grid_z, bins_x, bins_z)
    save_results(density_map, grid_x, grid_z, metrics)


if __name__ == "__main__":
    main()
