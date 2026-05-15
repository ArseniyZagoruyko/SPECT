import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from tqdm import tqdm
import imageio
from scipy import stats, ndimage
from scipy.optimize import minimize
import warnings
import time
warnings.filterwarnings('ignore')

# Try to import GPU libraries
try:
    import cupy as cp
    GPU_AVAILABLE = True
    print("GPU (CuPy) доступен")
except ImportError:
    GPU_AVAILABLE = False
    print("GPU недоступен, используем CPU")

try:
    from numba import jit, cuda
    NUMBA_AVAILABLE = True
    print("Numba доступен для JIT компиляции")
except ImportError:
    NUMBA_AVAILABLE = False
    print("Numba недоступен")

class EnhancedOriginEnsembleReconstructor:
    """
    Enhanced Origin Ensemble (OE) reconstruction method for Compton camera
    with GPU acceleration, advanced noise handling, and iterative improvements.
    """

    def __init__(self, cell_size=0.1, x_range=(-50, 50), z_range=(-50, 50),
                 use_gpu=False, use_iterative=True):
        self.cell_size = cell_size
        self.x_range = x_range
        self.z_range = z_range
        self.use_gpu = use_gpu and GPU_AVAILABLE
        self.use_iterative = use_iterative

        # Grid setup
        self.x_bins = np.arange(x_range[0], x_range[1] + cell_size, cell_size)
        self.z_bins = np.arange(z_range[0], z_range[1] + cell_size, cell_size)
        self.grid_x, self.grid_z = np.meshgrid(
            (self.x_bins[:-1] + self.x_bins[1:]) / 2,
            (self.z_bins[:-1] + self.z_bins[1:]) / 2
        )

        # Enhanced reconstruction parameters
        self.energy_resolution = 0.05  # 5% energy resolution
        self.position_uncertainty = 0.5  # mm
        self.min_energy_threshold = 10  # keV minimum energy
        self.max_energy_threshold = 800  # keV maximum energy

        # Advanced filtering parameters
        self.adaptive_threshold = True
        self.use_bilateral_filter = True
        self.use_morphological_ops = True

        # Iterative reconstruction parameters
        self.max_iterations = 5
        self.convergence_threshold = 1e-4
        self.regularization_lambda = 0.01

        # Noise robustness parameters
        self.outlier_threshold = 3.0
        self.robust_estimation = True
        self.confidence_weighting = True

        if self.use_gpu:
            self._setup_gpu()

    def _setup_gpu(self):
        """Setup GPU arrays and functions"""
        self.grid_x_gpu = cp.asarray(self.grid_x)
        self.grid_z_gpu = cp.asarray(self.grid_z)
        print(f"GPU память: {cp.cuda.MemoryPool().used_bytes() / 1024**2:.1f} MB")

    def read_deposits(self, filename, max_lines=10000):
        """Read energy deposit data with validation"""
        deposits = {}
        invalid_count = 0

        with open(filename, 'r') as file:
            for i, line in enumerate(file):
                if i >= max_lines:
                    break
                try:
                    parts = line.split()
                    eventID = int(parts[0])
                    E_lost = float(parts[1])
                    E_remaines = float(parts[2])

                    # Enhanced validation
                    if (E_lost > 0 and E_remaines > 0 and
                        E_lost < self.max_energy_threshold and
                        E_remaines < self.max_energy_threshold):
                        deposits[eventID] = (E_lost, E_remaines)
                    else:
                        invalid_count += 1

                except (ValueError, IndexError):
                    invalid_count += 1
                    continue

        print(f"Загружено {len(deposits)} валидных событий, пропущено {invalid_count} невалидных")
        return deposits

    def read_gamma_coordinates(self, filename, max_lines=10000):
        """Read gamma-ray coordinate data with validation"""
        coordinates = {}
        invalid_count = 0

        with open(filename, 'r') as file:
            for i, line in enumerate(file):
                if i >= max_lines:
                    break
                parts = line.strip().split()
                if len(parts) < 7:
                    invalid_count += 1
                    continue
                try:
                    eventID = int(parts[0])
                    coords = list(map(float, parts[1:7]))

                    # Validate coordinates are finite
                    if all(np.isfinite(coords)):
                        coordinates[eventID] = coords
                    else:
                        invalid_count += 1

                except ValueError:
                    invalid_count += 1
                    continue

        print(f"Загружено {len(coordinates)} валидных координат, пропущено {invalid_count} невалидных")
        return coordinates

    def calculate_scattering_angle_robust(self, E_remaines, E_lost, n_samples=100):
        """
        Robust scattering angle calculation with outlier handling
        """
        # Generate samples with energy resolution
        E_lost_samples = np.random.normal(E_lost, E_lost * self.energy_resolution, n_samples)
        E_remaines_samples = np.random.normal(E_remaines, E_remaines * self.energy_resolution, n_samples)

        # Physical constraints
        valid_mask = ((E_lost_samples > 0) & (E_remaines_samples > 0) &
                     (E_lost_samples < self.max_energy_threshold) &
                     (E_remaines_samples < self.max_energy_threshold))

        if np.sum(valid_mask) < 10:  # Need at least 10 valid samples
            return None, None, 0.0

        E_lost_valid = E_lost_samples[valid_mask]
        E_remaines_valid = E_remaines_samples[valid_mask]
        E_initial = E_lost_valid + E_remaines_valid

        # Compton scattering formula
        cos_theta = 1 - (511 * E_lost_valid) / (E_initial * (E_initial - E_lost_valid))

        # Physical constraint: -1 <= cos(theta) <= 1
        angle_mask = (cos_theta >= -1.001) & (cos_theta <= 1.001)  # Small tolerance for numerical errors
        cos_theta = np.clip(cos_theta[angle_mask], -1, 1)

        if len(cos_theta) < 5:
            return None, None, 0.0

        angles = np.degrees(np.arccos(cos_theta))

        # Robust statistics - remove outliers
        if self.robust_estimation and len(angles) > 10:
            q75, q25 = np.percentile(angles, [75, 25])
            iqr = q75 - q25
            lower_bound = q25 - 1.5 * iqr
            upper_bound = q75 + 1.5 * iqr
            angles = angles[(angles >= lower_bound) & (angles <= upper_bound)]

        if len(angles) == 0:
            return None, None, 0.0

        # Calculate confidence based on sample consistency
        confidence = 1.0 / (1.0 + np.std(angles) / max(np.mean(angles), 1.0))

        return np.mean(angles), max(np.std(angles), 1.0), confidence

    def calculate_likelihood_map_gpu(self, V, direction, theta_mean, theta_std, weight=1.0):
        """GPU-accelerated likelihood map calculation"""
        if not self.use_gpu:
            return self.calculate_likelihood_map_cpu(V, direction, theta_mean, theta_std, weight)

        # Convert to GPU arrays
        V_gpu = cp.asarray(V)
        direction_gpu = cp.asarray(direction)
        direction_norm_gpu = direction_gpu / cp.linalg.norm(direction_gpu)

        # Create grid points
        grid_points_gpu = cp.stack([self.grid_x_gpu, cp.zeros_like(self.grid_x_gpu), self.grid_z_gpu], axis=2)

        # Vector calculations
        source_vectors_gpu = grid_points_gpu - V_gpu.reshape(1, 1, 3)
        source_distances_gpu = cp.linalg.norm(source_vectors_gpu, axis=2)

        # Avoid division by zero
        valid_mask_gpu = source_distances_gpu > 1e-10
        source_vectors_norm_gpu = cp.zeros_like(source_vectors_gpu)

        # Normalize valid vectors
        valid_indices = cp.where(valid_mask_gpu)
        if len(valid_indices[0]) > 0:
            source_vectors_norm_gpu[valid_indices] = (source_vectors_gpu[valid_indices] /
                                                     source_distances_gpu[valid_indices, cp.newaxis])

        # Calculate angles
        cos_angles_gpu = cp.sum(source_vectors_norm_gpu * direction_norm_gpu, axis=2)
        cos_angles_gpu = cp.clip(cos_angles_gpu, -1, 1)
        calculated_angles_gpu = cp.degrees(cp.arccos(cos_angles_gpu))

        # Calculate likelihood
        if theta_std > 0:
            likelihood_map_gpu = cp.exp(-0.5 * ((calculated_angles_gpu - theta_mean) / theta_std) ** 2)
        else:
            likelihood_map_gpu = cp.where(cp.abs(calculated_angles_gpu - theta_mean) < 1.0, 1.0, 0.0)

        # Apply weights
        distance_weights_gpu = 1.0 / (source_distances_gpu**2 + 1e-6)
        likelihood_map_gpu = likelihood_map_gpu * distance_weights_gpu * weight * valid_mask_gpu

        # Convert back to CPU
        return cp.asnumpy(likelihood_map_gpu)

    def calculate_likelihood_map_cpu(self, V, direction, theta_mean, theta_std, weight=1.0):
        """CPU version of likelihood map calculation"""
        direction_norm = direction / np.linalg.norm(direction)
        grid_points = np.stack([self.grid_x, np.zeros_like(self.grid_x), self.grid_z], axis=2)

        source_vectors = grid_points - V.reshape(1, 1, 3)
        source_distances = np.linalg.norm(source_vectors, axis=2)

        valid_mask = source_distances > 1e-10
        source_vectors_norm = np.zeros_like(source_vectors)
        source_vectors_norm[valid_mask] = source_vectors[valid_mask] / source_distances[valid_mask, np.newaxis]

        cos_angles = np.sum(source_vectors_norm * direction_norm, axis=2)
        cos_angles = np.clip(cos_angles, -1, 1)
        calculated_angles = np.degrees(np.arccos(cos_angles))

        if theta_std > 0:
            likelihood_map = np.exp(-0.5 * ((calculated_angles - theta_mean) / theta_std) ** 2)
        else:
            likelihood_map = np.where(np.abs(calculated_angles - theta_mean) < 1.0, 1.0, 0.0)

        distance_weights = 1.0 / (source_distances**2 + 1e-6)
        likelihood_map = likelihood_map * distance_weights * weight * valid_mask

        return likelihood_map

    def estimate_event_quality_advanced(self, E_lost, E_remaines, V, X2_pos, confidence):
        """Advanced event quality estimation"""
        # Energy quality
        total_energy = E_lost + E_remaines
        energy_quality = np.exp(-abs(total_energy - 511) / 100)  # Prefer 511 keV

        # Geometry quality
        scatter_distance = np.linalg.norm(np.array(X2_pos) - np.array(V))
        geometry_quality = 1.0 / (1.0 + scatter_distance / 100.0)

        # Energy ratio quality (prefer balanced splits)
        energy_ratio = min(E_lost, E_remaines) / max(E_lost, E_remaines)
        ratio_quality = energy_ratio ** 0.5

        # Angular quality (prefer forward scattering for better resolution)
        scatter_angle = np.arccos(np.clip(1 - (511 * E_lost) / (total_energy * (total_energy - E_lost)), -1, 1))
        angle_quality = np.exp(-scatter_angle / (np.pi/2))

        # Confidence from angle calculation
        confidence_quality = confidence if self.confidence_weighting else 1.0

        return energy_quality * geometry_quality * ratio_quality * angle_quality * confidence_quality

    def iterative_reconstruction(self, initial_map, events_data, max_iterations=5):
        """
        Iterative reconstruction using Maximum Likelihood Expectation Maximization (MLEM)
        """
        if not self.use_iterative:
            return initial_map

        print("Запуск итеративной реконструкции MLEM...")

        current_map = initial_map.copy()
        current_map[current_map <= 0] = 1e-10  # Avoid zeros

        for iteration in range(max_iterations):
            print(f"  Итерация {iteration + 1}/{max_iterations}")

            # Forward projection
            forward_proj = np.zeros_like(current_map)
            back_proj = np.zeros_like(current_map)

            for event_data in tqdm(events_data[:100], desc=f"Итерация {iteration+1}", leave=False):  # Limit for speed
                V, direction, theta_mean, theta_std, weight = event_data

                # Calculate likelihood map for this event
                likelihood = self.calculate_likelihood_map_cpu(V, direction, theta_mean, theta_std, 1.0)

                # Forward projection
                expected_count = np.sum(current_map * likelihood)
                if expected_count > 1e-10:
                    # Back projection with correction
                    correction = weight / expected_count
                    back_proj += likelihood * correction

            # MLEM update
            old_map = current_map.copy()
            current_map = current_map * back_proj / (np.sum(back_proj) + 1e-10)

            # Apply regularization
            if self.regularization_lambda > 0:
                smoothed = ndimage.gaussian_filter(current_map, sigma=1.0)
                current_map = (1 - self.regularization_lambda) * current_map + self.regularization_lambda * smoothed

            # Check convergence
            change = np.sum(np.abs(current_map - old_map)) / np.sum(old_map)
            if change < self.convergence_threshold:
                print(f"  Сходимость достигнута на итерации {iteration + 1}")
                break

        return current_map

    def advanced_filtering(self, reconstruction_map):
        """
        Advanced filtering with multiple techniques
        """
        filtered_map = reconstruction_map.copy()

        # 1. Adaptive thresholding
        if self.adaptive_threshold:
            # Calculate local thresholds
            local_mean = ndimage.uniform_filter(filtered_map, size=5)
            local_std = ndimage.uniform_filter(filtered_map**2, size=5) - local_mean**2
            local_std = np.sqrt(np.maximum(local_std, 0))
            threshold_map = local_mean + 2 * local_std
            filtered_map = np.where(filtered_map > threshold_map, filtered_map, 0)

        # 2. Bilateral filtering for edge preservation
        if self.use_bilateral_filter:
            try:
                from skimage.restoration import denoise_bilateral
                filtered_map = denoise_bilateral(filtered_map, sigma_color=0.1, sigma_spatial=1.0)
            except ImportError:
                # Fallback to Gaussian filter
                filtered_map = ndimage.gaussian_filter(filtered_map, sigma=0.8)

        # 3. Morphological operations
        if self.use_morphological_ops:
            from scipy.ndimage import binary_opening, binary_closing
            # Create binary mask of significant regions
            mask = filtered_map > 0.1 * np.max(filtered_map)
            # Clean up small artifacts
            mask = binary_opening(mask, iterations=1)
            mask = binary_closing(mask, iterations=1)
            filtered_map = filtered_map * mask

        # 4. Final smoothing
        filtered_map = ndimage.gaussian_filter(filtered_map, sigma=0.5)

        return filtered_map

    def reconstruct_enhanced_oe(self, deposits_file, coordinates_file, max_events=10000):
        """
        Enhanced OE reconstruction with all advanced features
        """
        print("=" * 70)
        print("ENHANCED ORIGIN ENSEMBLE (OE) RECONSTRUCTION")
        print("=" * 70)

        start_time = time.time()

        # Load and validate data
        deposits = self.read_deposits(deposits_file, max_events)
        coordinates = self.read_gamma_coordinates(coordinates_file, max_events)

        common_events = list(set(deposits.keys()) & set(coordinates.keys()))
        print(f"Общие события: {len(common_events)}")

        # Pre-filter events
        valid_events = []
        events_data = []  # For iterative reconstruction

        for eventID in common_events:
            E_lost, E_remaines = deposits[eventID]
            if (E_lost >= self.min_energy_threshold and
                E_remaines >= self.min_energy_threshold and
                E_lost <= E_remaines):
                valid_events.append(eventID)

        print(f"Валидные события: {len(valid_events)}")

        # Initialize reconstruction
        reconstruction_map = np.zeros_like(self.grid_x)
        total_weight = 0.0
        processed_events = 0

        print("\nОбработка событий...")

        # Main reconstruction loop
        batch_size = 50
        for i in tqdm(range(0, len(valid_events), batch_size), desc="Enhanced OE"):
            batch_events = valid_events[i:i + batch_size]

            for eventID in batch_events:
                E_lost, E_remaines = deposits[eventID]
                Vx, Vy, Vz, X2, Y2, Z2 = coordinates[eventID]

                # Robust angle calculation
                theta_mean, theta_std, confidence = self.calculate_scattering_angle_robust(E_remaines, E_lost)
                if theta_mean is None:
                    continue

                if theta_std < 1.0:
                    theta_std = 2.0

                # Advanced quality estimation
                V = np.array([Vx, Vy, Vz])
                X2_pos = [X2, Y2, Z2]
                event_weight = self.estimate_event_quality_advanced(E_lost, E_remaines, V, X2_pos, confidence)

                direction = np.array([X2 - Vx, Y2 - Vy, Z2 - Vz])

                # Calculate likelihood map (GPU or CPU)
                if self.use_gpu:
                    event_likelihood = self.calculate_likelihood_map_gpu(
                        V, direction, theta_mean, theta_std, event_weight
                    )
                else:
                    event_likelihood = self.calculate_likelihood_map_cpu(
                        V, direction, theta_mean, theta_std, event_weight
                    )

                reconstruction_map += event_likelihood
                total_weight += event_weight
                processed_events += 1

                # Store for iterative reconstruction
                if self.use_iterative and len(events_data) < 1000:  # Limit for memory
                    events_data.append((V, direction, theta_mean, theta_std, event_weight))

        print(f"\nОбработано событий: {processed_events}")
        print(f"Общий вес: {total_weight:.2f}")

        # Normalize
        if total_weight > 0:
            reconstruction_map /= total_weight

        # Iterative reconstruction
        if self.use_iterative and len(events_data) > 0:
            reconstruction_map = self.iterative_reconstruction(reconstruction_map, events_data)

        # Advanced filtering
        filtered_map = self.advanced_filtering(reconstruction_map)

        processing_time = time.time() - start_time
        print(f"Время обработки: {processing_time:.2f} секунд")

        return reconstruction_map, filtered_map

    def comprehensive_analysis(self, raw_map, filtered_map):
        """
        Comprehensive analysis of reconstruction results
        """
        def analyze_map(rmap, name):
            # Basic metrics
            max_idx = np.unravel_index(np.argmax(rmap), rmap.shape)
            peak_pos = (self.grid_x[max_idx], self.grid_z[max_idx])
            peak_val = rmap[max_idx]

            # FWHM calculation
            x_profile = np.sum(rmap, axis=0)
            z_profile = np.sum(rmap, axis=1)

            def calc_fwhm(profile, bins):
                if np.max(profile) <= 0:
                    return None
                half_max = np.max(profile) / 2
                indices = np.where(profile >= half_max)[0]
                if len(indices) == 0:
                    return None
                return bins[indices[-1]] - bins[indices[0]]

            fwhm_x = calc_fwhm(x_profile, self.x_bins[:-1])
            fwhm_z = calc_fwhm(z_profile, self.z_bins[:-1])

            # Advanced metrics
            # Signal-to-noise ratio
            signal_region = rmap > 0.5 * peak_val
            noise_region = rmap < 0.1 * peak_val

            if np.sum(noise_region) > 0:
                snr = np.mean(rmap[signal_region]) / (np.std(rmap[noise_region]) + 1e-10)
            else:
                snr = float('inf')

            # Localization accuracy (assuming source at origin)
            localization_error = np.sqrt(peak_pos[0]**2 + peak_pos[1]**2)

            # Reconstruction efficiency
            efficiency = np.sum(rmap[signal_region]) / np.sum(rmap)

            return {
                'name': name,
                'peak_position': peak_pos,
                'peak_value': peak_val,
                'fwhm_x': fwhm_x,
                'fwhm_z': fwhm_z,
                'snr': snr,
                'localization_error': localization_error,
                'efficiency': efficiency
            }

        raw_results = analyze_map(raw_map, "Raw")
        filtered_results = analyze_map(filtered_map, "Filtered")

        return raw_results, filtered_results

def main():
    """Enhanced main function with comprehensive testing"""

    print("Инициализация Enhanced OE Reconstructor...")

    # Initialize with GPU if available
    reconstructor = EnhancedOriginEnsembleReconstructor(
        cell_size=0.1,
        x_range=(-50, 50),
        z_range=(-50, 50),
        use_gpu=GPU_AVAILABLE,
        use_iterative=True
    )

    # File paths
    deposits_file = '/home/zas/CERN/SPECT/build/deposits_FirstDet.txt'
    coordinates_file = '/home/zas/CERN/SPECT/build/GammaCoordinates_FirstDet.txt'

    # Perform enhanced reconstruction
    raw_map, filtered_map = reconstructor.reconstruct_enhanced_oe(
        deposits_file, coordinates_file, max_events=10000
    )

    # Comprehensive analysis
    raw_results, filtered_results = reconstructor.comprehensive_analysis(raw_map, filtered_map)

    # Display results
    print("\n" + "=" * 70)
    print("РЕЗУЛЬТАТЫ ENHANCED OE РЕКОНСТРУКЦИИ:")
    print("=" * 70)

    for results in [raw_results, filtered_results]:
        print(f"\n{results['name']} реконструкция:")
        print(f"  Пик: X={results['peak_position'][0]:.3f}, Z={results['peak_position'][1]:.3f} мм")
        print(f"  Значение пика: {results['peak_value']:.6f}")
        if results['fwhm_x'] and results['fwhm_z']:
            print(f"  FWHM: X={results['fwhm_x']:.2f}, Z={results['fwhm_z']:.2f} мм")
        print(f"  SNR: {results['snr']:.2f}")
        print(f"  Ошибка локализации: {results['localization_error']:.3f} мм")
        print(f"  Эффективность: {results['efficiency']:.3f}")

    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # Raw reconstruction
    im1 = axes[0].imshow(
        raw_map.T, origin='lower',
        extent=[reconstructor.x_range[0], reconstructor.x_range[1],
                reconstructor.z_range[0], reconstructor.z_range[1]],
        aspect='equal', cmap='hot', interpolation='bilinear'
    )
    axes[0].set_title('Enhanced OE (Raw)', fontsize=14)
    axes[0].set_xlabel('X (mm)', fontsize=12)
    axes[0].set_ylabel('Z (mm)', fontsize=12)
    plt.colorbar(im1, ax=axes[0])

    # Filtered reconstruction
    im2 = axes[1].imshow(
        filtered_map.T, origin='lower',
        extent=[reconstructor.x_range[0], reconstructor.x_range[1],
                reconstructor.z_range[0], reconstructor.z_range[1]],
        aspect='equal', cmap='hot', interpolation='bilinear'
    )
    axes[1].set_title('Enhanced OE (Filtered)', fontsize=14)
    axes[1].set_xlabel('X (mm)', fontsize=12)
    axes[1].set_ylabel('Z (mm)', fontsize=12)
    plt.colorbar(im2, ax=axes[1])

    plt.tight_layout()
    plt.savefig('enhanced_oe_reconstruction.png', dpi=300, bbox_inches='tight')
    plt.show()

    # Save final result
    if np.max(filtered_map) > 0:
        normalized_map = (filtered_map.T / np.max(filtered_map.T) * 255).astype('uint8')
    else:
        normalized_map = np.zeros_like(filtered_map.T, dtype='uint8')

    imageio.imwrite('enhanced_oe_final.jpg', normalized_map)

    print("=" * 70)
    return raw_map, filtered_map, raw_results, filtered_results

if __name__ == "__main__":
    raw_map, filtered_map, raw_results, filtered_results = main()
