
import numpy as np
import matplotlib
matplotlib.use('Agg')  # headless by default; remove if you want interactive windows
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
from shapely.geometry import LineString, Point

# Detector noise simulation
def simulate_detector_noise(E_lost, energy_resolution=True, noise_sigma_base=1.0, noise_sigma_factor=0.05):
    """
    Simulate detector energy resolution noise.
    
    Parameters:
    - E_lost: deposited energy (keV)
    - energy_resolution: if True, apply noise model dE = ±(1 + 0.05*E_lost)
    - noise_sigma_base: base noise level (keV)
    - noise_sigma_factor: energy-dependent factor
    
    Returns:
    - E_lost_noisy: energy with detector noise
    """
    if not energy_resolution:
        return E_lost
    
    # Energy resolution model: dE = ±(sigma_base + sigma_factor * E_lost)
    sigma = noise_sigma_base + noise_sigma_factor * E_lost
    
    # Add Gaussian noise
    noise = np.random.normal(0, sigma)
    E_lost_noisy = E_lost + noise
    
    # Ensure positive energy
    return max(0.1, E_lost_noisy)  # minimum 0.1 keV



def read_deposits(filename):
    """Read energy deposits: eventID R1 R2 (keV)."""
    deposits = {}
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                eventID = int(parts[0])
                R1 = float(parts[1])
                R2 = float(parts[2])
                deposits[eventID] = (R1, R2)
    return deposits


def read_gamma_coordinates(filename):
    """Read gamma interaction coordinates: eventID Rx1 Ry1 Rz1 Rx2 Ry2 Rz2."""
    coords = {}
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 7:
                eventID = int(parts[0])
                Rx1, Ry1, Rz1 = map(float, parts[1:4])
                Rx2, Ry2, Rz2 = map(float, parts[4:7])
                coords[eventID] = (Rx1, Ry1, Rz1, Rx2, Ry2, Rz2)
    return coords


def calculate_compton_angle(E_lost, E_remaining):
    """Compute Compton scatter angle theta (rad) with quality filters."""
    if E_lost <= 1.5 or E_lost >= 25.0:
        return None
    E_total = E_lost + E_remaining
    if E_total <= 0:
        return None
    cos_theta = 1 - (511.0 * E_lost) / (E_total * (E_total - E_lost))
    if not -1 <= cos_theta <= 1:
        return None
    theta = np.arccos(cos_theta)
    theta_deg = np.degrees(theta)
    if theta_deg < 30 or theta_deg > 150:
        return None
    return theta


def get_cone_ellipse_points_direct(V1, V2, theta, plane_y=0, num_points=180,
                                   x_range=(-50, 50), z_range=(-50, 50)):
    """Direct construction of cone-plane intersection ellipse (x,z points)."""
    V1 = np.array(V1, dtype=np.float64)
    V2 = np.array(V2, dtype=np.float64)
    P = V2 - V1
    P_norm = np.linalg.norm(P)
    if P_norm < 1e-10:
        return []
    if abs(P[1]) < 1e-10:
        return []
    t = (plane_y - V1[1]) / P[1]
    H = V1 + t * P
    h = H - V1
    h_norm = np.linalg.norm(h)
    if h_norm < 1e-10:
        return []
    h_unit = h / h_norm
    perp = np.array([-h_unit[2], 0, h_unit[0]])
    perp_norm = np.linalg.norm(perp)
    if perp_norm < 1e-10:
        perp = np.array([h_unit[0], 0, h_unit[2]])
        perp_norm = np.linalg.norm(perp)
        if perp_norm < 1e-10:
            return []
    perp_unit = perp / perp_norm
    perp2 = np.cross(h_unit, perp_unit)
    perp2_norm = np.linalg.norm(perp2)
    if perp2_norm < 1e-10:
        return []
    perp2_unit = perp2 / perp2_norm
    intersections = []
    phi_values = np.linspace(0, 2 * np.pi, num_points)
    for phi in phi_values:
        d_rot = h_unit * np.cos(theta) + (perp_unit * np.cos(phi) + perp2_unit * np.sin(phi)) * np.sin(theta)
        dn = np.linalg.norm(d_rot)
        if dn < 1e-10:
            continue
        d_rot = d_rot / dn
        if abs(d_rot[1]) < 1e-10:
            continue
        t_int = (plane_y - V1[1]) / d_rot[1]
        if t_int <= 0:
            continue
        pt = V1 + t_int * d_rot
        if (x_range[0] <= pt[0] <= x_range[1]) and (z_range[0] <= pt[2] <= z_range[1]):
            intersections.append((pt[0], pt[2]))
    return intersections


# -----------------------------
# Core: adaptive consensus intersections
# -----------------------------

def find_adaptive_consensus(ellipses, image_size=200, voxel_size=0.5,
                           x_range=(-50, 50), z_range=(-50, 50),
                           min_ellipses=3, coarse_step=8, detail_radius=16,
                           coarse_buffer_factor=2.0, fine_buffer_factor=0.8):
    """
    Two-stage search: coarse hotspot detection, then fine voting only nearby.
    Returns the intensity map (image_size x image_size).
    """
    x_min, x_max = x_range
    z_min, z_max = z_range

    # Coarse stage
    coarse_size = max(1, image_size // coarse_step)
    coarse_map = np.zeros((coarse_size, coarse_size), dtype=np.float32)

    buffer_distance = voxel_size * coarse_buffer_factor
    ellipse_polygons = []
    for ellipse in ellipses:
        if len(ellipse) < 3:
            ellipse_polygons.append(None)
            continue
        try:
            line = LineString(ellipse)
            buffered = line.buffer(buffer_distance)
            ellipse_polygons.append(buffered)
        except Exception:
            ellipse_polygons.append(None)

    for i in range(coarse_size):
        for k in range(coarse_size):
            x = x_min + voxel_size/2 + i * coarse_step * voxel_size
            z = z_min + voxel_size/2 + k * coarse_step * voxel_size
            p = Point(x, z)
            count = 0
            for poly in ellipse_polygons:
                if poly is not None and poly.contains(p):
                    count += 1
            coarse_map[i, k] = count

    # Hotspots threshold (fraction of max)
    if np.max(coarse_map) <= 0:
        return coarse_map.repeat(coarse_step, axis=0)[:image_size, :].repeat(coarse_step, axis=1)[:,:image_size]
    hotspots = np.where(coarse_map >= 0.4 * np.max(coarse_map))

    # Fine stage around hotspots
    detailed_map = np.zeros((image_size, image_size), dtype=np.float32)
    fine_buffer = voxel_size * fine_buffer_factor
    fine_polygons = []
    for ellipse in ellipses:
        if len(ellipse) < 3:
            fine_polygons.append(None)
            continue
        try:
            line = LineString(ellipse)
            buffered = line.buffer(fine_buffer)
            fine_polygons.append(buffered)
        except Exception:
            fine_polygons.append(None)

    for hot_i, hot_k in tqdm(list(zip(hotspots[0], hotspots[1])), desc='Refine'):
        center_i = hot_i * coarse_step
        center_k = hot_k * coarse_step
        for di in range(-detail_radius, detail_radius + 1):
            for dk in range(-detail_radius, detail_radius + 1):
                i = center_i + di
                k = center_k + dk
                if not (0 <= i < image_size and 0 <= k < image_size):
                    continue
                x = x_min + voxel_size/2 + i * voxel_size
                z = z_min + voxel_size/2 + k * voxel_size
                p = Point(x, z)
                cnt = 0
                for poly in fine_polygons:
                    if poly is not None and poly.contains(p):
                        cnt += 1
                if cnt >= min_ellipses:
                    w = cnt * cnt
                    if cnt >= 5:
                        w *= 1.5
                    detailed_map[i, k] = max(detailed_map[i, k], w)

    return detailed_map


# -----------------------------
# Reconstruction pipeline and simple visualization
# -----------------------------

def reconstruct(deposits_file, coordinates_file,
                image_size=200, voxel_size=0.5,
                plane_y=0.0, num_cone_points=180,
                x_range=(-50, 50), z_range=(-50, 50),
                min_ellipses=3, coarse_step=8, detail_radius=16,
                filter_sigma=0.7, max_events=None,
                apply_detector_noise=True, noise_sigma_base=1.0, noise_sigma_factor=0.05):
    """Build ellipses and run adaptive consensus reconstruction."""
    deposits = read_deposits(deposits_file)
    coords = read_gamma_coordinates(coordinates_file)
    common = list(set(deposits.keys()) & set(coords.keys()))
    if max_events is not None and max_events < len(common):
        common = common[:max_events]

    ellipses = []
    noise_stats = {'original_count': 0, 'noise_rejected': 0}
    
    for eventID in tqdm(common, desc='Ellipses'):
        R1_orig, R2_orig = deposits[eventID]
        
        # Apply detector noise simulation
        if apply_detector_noise:
            R1 = simulate_detector_noise(R1_orig, True, noise_sigma_base, noise_sigma_factor)
            R2 = simulate_detector_noise(R2_orig, True, noise_sigma_base, noise_sigma_factor)
            noise_stats['original_count'] += 1
        else:
            R1, R2 = R1_orig, R2_orig
        
        Rx1, Ry1, Rz1, Rx2, Ry2, Rz2 = coords[eventID]
        theta = calculate_compton_angle(R1, R2)
        if theta is None:
            if apply_detector_noise:
                noise_stats['noise_rejected'] += 1
            continue
        V1 = (Rx1, Ry1, Rz1)
        V2 = (Rx2, Ry2, Rz2)
        pts = get_cone_ellipse_points_direct(
            V1, V2, theta, plane_y,
            num_points=num_cone_points,
            x_range=x_range, z_range=z_range
        )
        if len(pts) > 3:
            ellipses.append(pts)

    if not ellipses:
        return None
    
    # Print noise statistics if noise was applied
    if apply_detector_noise and noise_stats['original_count'] > 0:
        rejection_rate = 100 * noise_stats['noise_rejected'] / noise_stats['original_count']
        print(f"Detector noise applied: {noise_stats['noise_rejected']}/{noise_stats['original_count']} events rejected ({rejection_rate:.1f}%)")
        print(f"Noise model: dE = ±({noise_sigma_base:.1f} + {noise_sigma_factor:.2f}*E_lost) keV")

    result_map = find_adaptive_consensus(
        ellipses, image_size, voxel_size, x_range, z_range,
        min_ellipses=min_ellipses, coarse_step=coarse_step,
        detail_radius=detail_radius
    )

    if filter_sigma and filter_sigma > 0:
        result_map = gaussian_filter(result_map, sigma=filter_sigma)

    return result_map


def gaussian_func(x, amplitude, center, sigma, baseline):
    """Gaussian function for fitting."""
    return amplitude * np.exp(-((x - center) ** 2) / (2 * sigma ** 2)) + baseline

def _compute_fwhm_1d(profile, spacing):
    """Compute FWHM of a 1D profile using linear interpolation."""
    if profile is None or len(profile) == 0:
        return np.nan
    y = np.asarray(profile, dtype=float)
    y = np.maximum(y, 0.0)
    peak = np.max(y)
    if peak <= 0:
        return np.nan
    half = peak / 2.0
    idx_max = int(np.argmax(y))
    # left half-maximum crossing
    i_left = idx_max
    while i_left > 0 and y[i_left] > half:
        i_left -= 1
    if y[i_left] == y[i_left+1]:
        x_left = i_left * spacing
    else:
        frac = (half - y[i_left]) / (y[i_left+1] - y[i_left])
        x_left = (i_left + frac) * spacing
    # right half-maximum crossing
    i_right = idx_max
    n = len(y)
    while i_right < n-1 and y[i_right] > half:
        i_right += 1
    if y[i_right] == y[i_right-1]:
        x_right = i_right * spacing
    else:
        frac = (half - y[i_right]) / (y[i_right-1] - y[i_right])
        x_right = (i_right - frac) * spacing
    return max(0.0, x_right - x_left)

def fit_gaussian_to_profile(x_coords, profile, fix_center_to_max=True):
    """Fit Gaussian to profile with center fixed to maximum position."""
    try:
        # Find maximum position and validate data
        peak_idx = np.argmax(profile)
        amplitude_guess = np.max(profile)
        center_fixed = x_coords[peak_idx]
        baseline_guess = np.min(profile)
        
        # Better sigma estimation from FWHM approximation
        # Estimate FWHM first, then convert to sigma
        half_max = amplitude_guess / 2.0 + baseline_guess
        indices_above_half = np.where(profile > half_max)[0]
        if len(indices_above_half) > 0:
            fwhm_estimate = (indices_above_half[-1] - indices_above_half[0]) * abs(x_coords[1] - x_coords[0])
            sigma_guess = max(0.5, fwhm_estimate / 2.355)  # FWHM = 2.355 * sigma
        else:
            sigma_guess = 2.0
        
        # Ensure reasonable bounds
        sigma_guess = max(0.1, min(sigma_guess, 20.0))
        
        if fix_center_to_max:
            # Define constrained Gaussian function with fixed center AND amplitude
            amplitude_fixed = amplitude_guess  # Fix amplitude to data peak
            
            def gaussian_fixed_center_and_amplitude(x, sigma, baseline):
                return amplitude_fixed * np.exp(-((x - center_fixed) ** 2) / (2 * sigma ** 2)) + baseline
            
            initial_guess = [sigma_guess, baseline_guess]
            
            # Set bounds to ensure reasonable parameters
            bounds = ([0.1, -np.inf],          # lower bounds: sigma > 0.1, baseline unbounded
                     [20.0, amplitude_fixed])  # upper bounds: sigma < 20, baseline < peak
            
            # Fit with fixed center and amplitude
            popt, pcov = curve_fit(gaussian_fixed_center_and_amplitude, x_coords, profile, 
                                 p0=initial_guess, bounds=bounds, maxfev=5000)
            sigma, baseline = popt
            amplitude = amplitude_fixed  # Amplitude is fixed
            center = center_fixed        # Center is fixed
            
            # Generate fit curve
            y_pred = gaussian_fixed_center_and_amplitude(x_coords, sigma, baseline)
            
        else:
            # Original free-center fitting
            initial_guess = [amplitude_guess, center_fixed, sigma_guess, baseline_guess]
            popt, pcov = curve_fit(gaussian_func, x_coords, profile, p0=initial_guess, maxfev=5000)
            amplitude, center, sigma, baseline = popt
            y_pred = gaussian_func(x_coords, *popt)
        
        # Calculate FWHM from fitted sigma: FWHM = 2.355 * sigma
        fwhm_fit = 2.355 * abs(sigma)
        
        # Calculate R-squared (y_pred already calculated above)
        ss_res = np.sum((profile - y_pred) ** 2)
        ss_tot = np.sum((profile - np.mean(profile)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        return {
            'amplitude': amplitude,
            'center': center,
            'sigma': abs(sigma),
            'baseline': baseline,
            'fwhm_fit': fwhm_fit,
            'r_squared': r_squared,
            'fit_curve': y_pred
        }
    except Exception as e:
        # Print debug info for troubleshooting
        print(f"Gaussian fit failed: {e}")
        print(f"Peak amplitude: {np.max(profile):.1f}, Profile range: {np.min(profile):.1f} to {np.max(profile):.1f}")
        return None


def simple_visualize(result_map, voxel_size=0.5,
                      x_range=(-50, 50), z_range=(-50, 50),
                      output_file='reconstruction.png'):
    """Simple figure: left = reconstruction heatmap; right = X-axis profile only."""
    if result_map is None:
        print('No result to visualize')
        return
    x_min, x_max = x_range
    z_min, z_max = z_range
    max_idx = np.unravel_index(np.argmax(result_map), result_map.shape)
    max_x = x_min + voxel_size/2 + max_idx[0] * voxel_size
    max_z = z_min + voxel_size/2 + max_idx[1] * voxel_size

    # Profile along X at Z = max_z (this is result_map[:, max_idx[1]])
    x_profile = result_map[:, max_idx[1]]
    x_coords = np.linspace(x_min + voxel_size/2, x_max - voxel_size/2, len(x_profile))
    
    fwhm_x = _compute_fwhm_1d(x_profile, voxel_size)
    
    # Gaussian fitting with center fixed to maximum
    fit_result = fit_gaussian_to_profile(x_coords, x_profile, fix_center_to_max=True)

    # Console output
    print(f"Found position (x,z): ({max_x:.3f} mm, {max_z:.3f} mm)")
    print(f"FWHM (direct): {fwhm_x:.3f} mm")
    
    if fit_result is not None:
        print(f"Gaussian fit:")
        print(f"  Center: {fit_result['center']:.3f} mm")
        print(f"  Sigma: {fit_result['sigma']:.3f} mm")
        print(f"  FWHM (fit): {fit_result['fwhm_fit']:.3f} mm")
        print(f"  R²: {fit_result['r_squared']:.3f}")
    else:
        print("Gaussian fit failed")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: heatmap (no markers, no titles)
    im = axes[0].imshow(result_map.T, extent=[x_min, x_max, z_min, z_max],
                        origin='lower', cmap='hot', aspect='equal')
    axes[0].set_xlabel('X (mm)')
    axes[0].set_ylabel('Z (mm)')
    axes[0].grid(True, alpha=0.2)
    plt.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04, label='Intensity')

    # Right: X profile with Gaussian fit overlay
    axes[1].plot(x_coords, x_profile, 'b-', linewidth=2, label='Data')
    
    # Add Gaussian fit if available
    if fit_result is not None:
        axes[1].plot(x_coords, fit_result['fit_curve'], 'r--', linewidth=1.5, label='Gaussian fit')
        axes[1].legend(fontsize=10)
    
    axes[1].set_xlabel('X (mm)')
    axes[1].set_ylabel('Intensity')
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved: {output_file}")
    # No plt.show() for headless


# -----------------------------
# CLI entry point (adjust parameters here)
# -----------------------------
if __name__ == '__main__':
    # User-adjustable parameters
    deposits_file = 'build/deposits_FirstDet.txt'
    coordinates_file = 'build/GammaCoordinates_FirstDet.txt'

    # Fast test parameters for noise validation
    image_size = 200   
    voxel_size = 0.5   # back to 0.5 for speed
    plane_y = 0.0
    num_cone_points = 180
    x_range = (-50, 50)
    z_range = (-50, 50)

    # Original working parameters
    min_ellipses = 3
    coarse_step = 8    # original values
    detail_radius = 16 
    filter_sigma = 0.7
    max_events = 5000  # test with more events to see noise effect
    
    # Detector noise simulation parameters
    apply_detector_noise = True    # Enable/disable noise simulation
    noise_sigma_base = 1.0         # Base noise: ±1 keV
    noise_sigma_factor = 0.05      # Energy-dependent: ±0.05*E_lost keV

    result = reconstruct(
        deposits_file, coordinates_file,
        image_size=image_size, voxel_size=voxel_size,
        plane_y=plane_y, num_cone_points=num_cone_points,
        x_range=x_range, z_range=z_range,
        min_ellipses=min_ellipses,
        coarse_step=coarse_step, detail_radius=detail_radius,
        filter_sigma=filter_sigma, max_events=max_events,
        apply_detector_noise=apply_detector_noise,
        noise_sigma_base=noise_sigma_base, noise_sigma_factor=noise_sigma_factor,
    )

    simple_visualize(result, voxel_size=voxel_size,
                     x_range=x_range, z_range=z_range,
                     output_file='reconstruction.png')
