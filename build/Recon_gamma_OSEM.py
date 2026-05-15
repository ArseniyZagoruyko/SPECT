#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""


OSEM is an iterative statistical reconstruction method that:
1. Starts with initial estimate (backprojection)
2. Forward projects current estimate to predict measurements
3. Compares with actual data
4. Backprojects the ratio to update estimate
5. Repeats for N iterations

Formula: λ^(n+1) = λ^(n) × (backproj(measured/forward_proj(λ^(n))) / sensitivity)


"""

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from tqdm import tqdm
import imageio
import math

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input files
DEPOSITS_FILE = "/home/zas/CERN/SPECT/build/deposits(1detector).txt"
COORDS_FILE = "/home/zas/CERN/SPECT/build/GammaCoordinates(1detector).txt"

# Reconstruction grid
CELL_SIZE = 0.1  # mm
X_RANGE = (-50.0, 50.0)  # mm
Z_RANGE = (-50.0, 50.0)  # mm
GRID_Y = 0.0  # mm, reconstruction plane

# Energy cuts
MIN_E_LOST = 0.01  # keV
ELECTRON_MASS_KEV = 511.0  # keV

# OSEM parameters
N_ITERATIONS = 3  # Number of iterations (fewer is often better)
N_SUBSETS = 5  # Number of ordered subsets (OSEM speedup)
INITIAL_VALUE = 0.1  # Initial uniform value for image
CONE_ROTATION_STEPS = 360  # Angular sampling for cone intersection
USE_BACKPROJ_INIT = True  # Use backprojection for initialization
RELAXATION_PARAM = 1.0  # Relaxation parameter (1.0 = standard OSEM)
MAX_IMAGE_VALUE = 10.0  # Maximum allowed pixel value (prevent overflow)
SMOOTHING_SIGMA = 0.0  # Small Gaussian smoothing after each iteration (0 = disabled)

# Output
SAVE_PREFIX = "osem_140keV"
MAX_EVENTS = 10000  # None = all events

# =============================================================================
# DATA LOADING
# =============================================================================

def load_deposits(filename, max_events=None):
    """Load energy deposits"""
    deposits = {}
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            if max_events and i >= max_events:
                break
            parts = line.split()
            eventID = int(parts[0])
            E_lost = float(parts[1])
            E_remains = float(parts[2])
            deposits[eventID] = (E_lost, E_remains)
    return deposits

def load_coordinates(filename, max_events=None):
    """Load interaction coordinates"""
    coordinates = {}
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            if max_events and i >= max_events:
                break
            parts = line.strip().split()
            if len(parts) < 7:
                continue
            try:
                eventID = int(parts[0])
                Vx, Vy, Vz = map(float, parts[1:4])
                X2, Y2, Z2 = map(float, parts[4:7])
                coordinates[eventID] = (Vx, Vy, Vz, X2, Y2, Z2)
            except ValueError:
                continue
    return coordinates

# =============================================================================
# GEOMETRY UTILITIES (from backprojection)
# =============================================================================

def compute_theta(E_remains, E_lost):
    """Compute Compton scatter angle from energies"""
    if E_lost <= MIN_E_LOST or E_remains <= MIN_E_LOST:
        return None

    E_initial = E_lost + E_remains
    cos_theta = 1 - (ELECTRON_MASS_KEV * E_lost) / (E_initial * (E_initial - E_lost))

    if not -1 <= cos_theta <= 1:
        return None

    theta = np.degrees(np.arccos(cos_theta))
    return theta

def rotation_matrix(axis, theta_deg):
    """Rodrigues rotation formula"""
    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-10:
        return np.eye(3)
    axis = axis / axis_norm
    theta = np.radians(theta_deg)
    a = np.cos(theta / 2)
    b, c, d = -axis * np.sin(theta / 2)
    return np.array([
        [a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
        [2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
        [2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c]
    ])

def get_cone_intersections(V, X2, theta, steps=CONE_ROTATION_STEPS):
    """
    Get intersection points of Compton cone with Y=0 plane.
    Returns list of (x, z) coordinates.
    """
    # Scattered direction
    P = X2 - V

    # Intersection of scattered ray with Y=0
    if abs(P[1]) < 1e-10:
        return []

    t = -V[1] / P[1]
    H = V + t * P

    # Axis of cone rotation
    h = H - V
    h_norm = np.linalg.norm(h)
    if h_norm < 1e-10:
        return []
    h_unit = h / h_norm

    # Perpendicular vector in XZ plane
    perp = np.array([-h_unit[2], 0, h_unit[0]])
    perp_norm = np.linalg.norm(perp)
    if perp_norm < 1e-10:
        return []
    perp_unit = perp / perp_norm

    # Initial direction on cone
    theta_rad = np.radians(theta)
    d = h_unit * np.cos(theta_rad) + perp_unit * np.sin(theta_rad)

    # Rotate and find intersections
    intersections = []
    for angle in range(steps):
        R = rotation_matrix(h, angle)
        d_rotated = R @ d

        if abs(d_rotated[1]) < 1e-10:
            continue

        t_int = -V[1] / d_rotated[1]
        if t_int <= 0:
            continue

        point = V + t_int * d_rotated
        if np.isfinite(point).all():
            intersections.append((point[0], point[2]))

    return intersections

# =============================================================================
# BACKPROJECTION (for initialization and backward step)
# =============================================================================

def backproject(events_data, grid_shape, x_bins, z_bins, subset_indices=None, desc="Backprojection"):
    """
    Perform backprojection for given events.

    events_data: list of (eventID, theta, V, X2, weight)
    subset_indices: if not None, only use these event indices
    Returns: 2D array of backprojected image

    CRITICAL: Normalize weight by number of valid intersection points per event
    to prevent ring artifacts and ensure proper weighting.
    """
    backproj_image = np.zeros(grid_shape, dtype=np.float64)

    if subset_indices is None:
        subset_indices = range(len(events_data))

    for idx in tqdm(subset_indices, desc=desc, leave=False):
        eventID, theta, V, X2, weight = events_data[idx]

        # Get cone intersections
        intersections = get_cone_intersections(V, X2, theta)

        if len(intersections) == 0:
            continue

        # Count valid intersection points for normalization
        valid_points = []
        for x, z in intersections:
            if not (np.isfinite(x) and np.isfinite(z)):
                continue

            i = int((x - x_bins[0]) / CELL_SIZE)
            j = int((z - z_bins[0]) / CELL_SIZE)

            if 0 <= i < grid_shape[1] and 0 <= j < grid_shape[0]:
                valid_points.append((i, j))

        # Normalize weight by number of valid points
        if len(valid_points) > 0:
            normalized_weight = weight / len(valid_points)

            # Add normalized weight to all valid intersection points
            for i, j in valid_points:
                backproj_image[j, i] += normalized_weight

    return backproj_image

# =============================================================================
# FORWARD PROJECTION
# =============================================================================

def forward_project_event(image, x_bins, z_bins, theta, V, X2):
    """
    Forward project image for one event.
    Returns expected count for this event based on current image estimate.

    CRITICAL: Must normalize by number of intersection points to get proper
    expected value for each event. Otherwise creates ring artifacts.
    """
    intersections = get_cone_intersections(V, X2, theta)

    if len(intersections) == 0:
        return 0.0

    expected_count = 0.0
    n_valid_points = 0

    for x, z in intersections:
        if not (np.isfinite(x) and np.isfinite(z)):
            continue

        # Find grid indices
        i = int((x - x_bins[0]) / CELL_SIZE)
        j = int((z - z_bins[0]) / CELL_SIZE)

        if 0 <= i < image.shape[1] and 0 <= j < image.shape[0]:
            expected_count += image[j, i]
            n_valid_points += 1

    # MUST normalize by number of points - average intensity along cone
    if n_valid_points > 0:
        return expected_count / n_valid_points
    else:
        return 0.0


# =============================================================================
# OSEM RECONSTRUCTION
# =============================================================================

def osem_reconstruct(events_data, grid_shape, x_bins, z_bins,
                     n_iterations=N_ITERATIONS, n_subsets=N_SUBSETS):
    """
    Perform OSEM reconstruction.

    OSEM update formula:
    λ^(n+1) = λ^(n) × (backproj(data/forward_proj(λ^(n))) / sensitivity)
    """
    print("=" * 70)
    print("OSEM RECONSTRUCTION FOR COMPTON CAMERA")
    print("=" * 70)
    print(f"Grid size: {grid_shape}")
    print(f"Total events: {len(events_data)}")
    print(f"Iterations: {n_iterations}")
    print(f"Subsets: {n_subsets}")
    print(f"Relaxation parameter: {RELAXATION_PARAM}")
    print("=" * 70)

    # Initialize image
    if USE_BACKPROJ_INIT:
        print("\nInitializing with backprojection...")
        initial_image = backproject(events_data, grid_shape, x_bins, z_bins,
                                   desc="Initial backprojection")
        initial_image[initial_image < 0] = 0
        if np.max(initial_image) > 0:
            initial_image = initial_image / np.max(initial_image) * INITIAL_VALUE
        else:
            initial_image = np.ones(grid_shape) * INITIAL_VALUE
    else:
        initial_image = np.ones(grid_shape) * INITIAL_VALUE

    current_image = initial_image.copy()

    # Compute sensitivity image (normalization)
    print("Computing sensitivity image...")
    sens_events = [(eid, th, v, x2, 1.0) for eid, th, v, x2, _ in events_data]
    sensitivity = backproject(sens_events, grid_shape, x_bins, z_bins,
                             desc="Sensitivity computation")
    sensitivity[sensitivity < 1e-6] = 1.0  # Avoid division by zero

    # Divide events into ordered subsets
    n_events = len(events_data)
    subset_size = n_events // n_subsets
    subsets = []
    for s in range(n_subsets):
        start = s * subset_size
        end = (s + 1) * subset_size if s < n_subsets - 1 else n_events
        subsets.append(list(range(start, end)))

    # OSEM iterations
    iteration_images = [initial_image.copy()]

    print(f"\nStarting OSEM iterations...")
    for iteration in range(n_iterations):
        for subset_idx, subset_indices in enumerate(subsets):
            # Forward projection for this subset
            ratios = []
            pbar_desc = f"Iter {iteration+1}/{n_iterations} Sub {subset_idx+1}/{n_subsets} - Forward"

            for idx in tqdm(subset_indices, desc=pbar_desc, leave=False):
                eventID, theta, V, X2, measured = events_data[idx]

                expected = forward_project_event(current_image, x_bins, z_bins,
                                                theta, V, X2)

                # Compute ratio (measured / expected)
                if expected > 1e-6:
                    ratio = measured / expected
                else:
                    ratio = 1.0

                ratios.append((eventID, theta, V, X2, ratio))

            # Backproject ratios
            backproj_desc = f"Iter {iteration+1}/{n_iterations} Sub {subset_idx+1}/{n_subsets} - Backproj"
            correction = backproject(ratios, grid_shape, x_bins, z_bins,
                                   range(len(ratios)), desc=backproj_desc)

            # Multiplicative update with relaxation
            correction = correction / sensitivity * n_subsets  # Normalize by subset
            update_factor = np.power(correction, RELAXATION_PARAM)
            current_image = current_image * update_factor

            # Enforce non-negativity and maximum value
            current_image[current_image < 0] = 0
            current_image[current_image > MAX_IMAGE_VALUE] = MAX_IMAGE_VALUE

            # Optional: light smoothing to prevent noise amplification
            if SMOOTHING_SIGMA > 0:
                from scipy import ndimage
                current_image = ndimage.gaussian_filter(current_image, sigma=SMOOTHING_SIGMA)

        iteration_images.append(current_image.copy())

        # Diagnostics (print every few iterations)
        if (iteration + 1) % 2 == 0 or iteration == n_iterations - 1:
            peak_val = np.max(current_image)
            peak_idx = np.unravel_index(np.argmax(current_image), current_image.shape)
            x_centers = (x_bins[:-1] + x_bins[1:]) / 2
            z_centers = (z_bins[:-1] + z_bins[1:]) / 2
            peak_x = x_centers[peak_idx[1]]
            peak_z = z_centers[peak_idx[0]]
            print(f"Iteration {iteration+1}: Peak at ({peak_x:.2f}, {peak_z:.2f}) = {peak_val:.6f}")

    return current_image, iteration_images

# =============================================================================
# METRICS
# =============================================================================

def calculate_fwhm(image, x_bins, z_bins):
    """Calculate FWHM along X and Z axes"""
    # Profiles
    x_profile = np.sum(image, axis=0)
    z_profile = np.sum(image, axis=1)

    def fwhm_from_profile(profile, bins):
        if len(profile) == 0 or np.max(profile) <= 0:
            return None
        half_max = np.max(profile) / 2
        indices = np.where(profile >= half_max)[0]
        if len(indices) == 0:
            return None
        fwhm = bins[indices[-1]] - bins[indices[0]]
        return fwhm

    fwhm_x = fwhm_from_profile(x_profile, x_bins[:-1])
    fwhm_z = fwhm_from_profile(z_profile, z_bins[:-1])

    return fwhm_x, fwhm_z

# =============================================================================
# VISUALIZATION
# =============================================================================

def visualize_results(final_image, iteration_images, x_bins, z_bins):
    """Visualize reconstruction results"""
    extent = [x_bins[0], x_bins[-1], z_bins[0], z_bins[-1]]

    # Plot final result
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Initial (backprojection)
    im1 = axes[0].imshow(iteration_images[0].T, origin='lower',
                         extent=extent, aspect='equal', cmap='hot')
    axes[0].set_title('Initial (Backprojection)', fontsize=14)
    axes[0].set_xlabel('X (mm)', fontsize=12)
    axes[0].set_ylabel('Z (mm)', fontsize=12)
    plt.colorbar(im1, ax=axes[0], label='Intensity')

    # Final OSEM result
    im2 = axes[1].imshow(final_image.T, origin='lower',
                         extent=extent, aspect='equal', cmap='hot')
    axes[1].set_title(f'OSEM ({N_ITERATIONS} iterations)', fontsize=14)
    axes[1].set_xlabel('X (mm)', fontsize=12)
    axes[1].set_ylabel('Z (mm)', fontsize=12)
    plt.colorbar(im2, ax=axes[1], label='Intensity')

    plt.tight_layout()
    plt.savefig(f'{SAVE_PREFIX}_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Save final as JPG
    if np.max(final_image) > 0:
        normalized = (final_image.T / np.max(final_image.T) * 255).astype('uint8')
    else:
        normalized = np.zeros_like(final_image.T, dtype='uint8')
    imageio.imwrite(f'{SAVE_PREFIX}_final.jpg', normalized)

    # Plot convergence (FWHM vs iteration)
    fwhm_x_history = []
    fwhm_z_history = []
    for img in iteration_images:
        fx, fz = calculate_fwhm(img, x_bins, z_bins)
        fwhm_x_history.append(fx if fx else np.nan)
        fwhm_z_history.append(fz if fz else np.nan)

    fig, ax = plt.subplots(figsize=(10, 6))
    iterations = range(len(fwhm_x_history))
    ax.plot(iterations, fwhm_x_history, 'o-', label='FWHM X', linewidth=2)
    ax.plot(iterations, fwhm_z_history, 's-', label='FWHM Z', linewidth=2)
    ax.set_xlabel('Iteration', fontsize=14)
    ax.set_ylabel('FWHM (mm)', fontsize=14)
    ax.set_title('OSEM Convergence', fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{SAVE_PREFIX}_convergence.png', dpi=300, bbox_inches='tight')
    plt.close()

# =============================================================================
# MAIN
# =============================================================================

def main():
    print("Loading data...")
    deposits = load_deposits(DEPOSITS_FILE, MAX_EVENTS)
    coordinates = load_coordinates(COORDS_FILE, MAX_EVENTS)

    common_ids = sorted(set(deposits.keys()) & set(coordinates.keys()))
    print(f"Total events: {len(common_ids)}")

    # Prepare event data
    print("Preparing events...")
    events_data = []
    for eventID in tqdm(common_ids, desc="Processing events"):
        E_lost, E_remains = deposits[eventID]

        # Filter
        if E_lost > E_remains:
            continue

        theta = compute_theta(E_remains, E_lost)
        if theta is None:
            continue

        Vx, Vy, Vz, X2, Y2, Z2 = coordinates[eventID]
        V = np.array([Vx, Vy, Vz])
        X2_pos = np.array([X2, Y2, Z2])

        # Measured data (for this event, we count it as "1" detection)
        measured = 1.0

        events_data.append((eventID, theta, V, X2_pos, measured))

    print(f"Valid events for reconstruction: {len(events_data)}")

    # Setup grid
    x_bins = np.arange(X_RANGE[0], X_RANGE[1] + CELL_SIZE, CELL_SIZE)
    z_bins = np.arange(Z_RANGE[0], Z_RANGE[1] + CELL_SIZE, CELL_SIZE)
    grid_shape = (len(z_bins) - 1, len(x_bins) - 1)

    # Run OSEM
    final_image, iteration_images = osem_reconstruct(
        events_data, grid_shape, x_bins, z_bins,
        n_iterations=N_ITERATIONS, n_subsets=N_SUBSETS
    )

    # Calculate metrics
    print("\n" + "=" * 70)
    print("FINAL RESULTS")
    print("=" * 70)

    # Peak position
    peak_idx = np.unravel_index(np.argmax(final_image), final_image.shape)
    x_centers = (x_bins[:-1] + x_bins[1:]) / 2
    z_centers = (z_bins[:-1] + z_bins[1:]) / 2
    peak_x = x_centers[peak_idx[1]]
    peak_z = z_centers[peak_idx[0]]
    peak_val = final_image[peak_idx]

    print(f"Peak position: X={peak_x:.2f} mm, Z={peak_z:.2f} mm")
    print(f"Peak value: {peak_val:.6f}")

    # FWHM
    fwhm_x, fwhm_z = calculate_fwhm(final_image, x_bins, z_bins)
    if fwhm_x and fwhm_z:
        print(f"FWHM X: {fwhm_x:.2f} mm")
        print(f"FWHM Z: {fwhm_z:.2f} mm")

    # Compare with initial
    init_fwhm_x, init_fwhm_z = calculate_fwhm(iteration_images[0], x_bins, z_bins)
    if init_fwhm_x and fwhm_x:
        improvement_x = (init_fwhm_x - fwhm_x) / init_fwhm_x * 100
        improvement_z = (init_fwhm_z - fwhm_z) / init_fwhm_z * 100
        print(f"\nImprovement vs. backprojection:")
        print(f"  FWHM X: {improvement_x:.1f}% better")
        print(f"  FWHM Z: {improvement_z:.1f}% better")

    print("=" * 70)

    # Visualize
    visualize_results(final_image, iteration_images, x_bins, z_bins)
    print(f"\nImages saved: {SAVE_PREFIX}_*.png/jpg")

if __name__ == "__main__":
    main()
