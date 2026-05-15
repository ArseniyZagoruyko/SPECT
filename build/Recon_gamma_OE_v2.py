import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from tqdm import tqdm
import imageio
from scipy import stats
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

class OriginEnsembleReconstructor:
    """
    Origin Ensemble (OE) reconstruction method for Compton camera.

    This method focuses on statistical reconstruction of source positions
    by analyzing the ensemble of possible gamma-ray origins, incorporating
    measurement uncertainties and using Bayesian-like statistical approaches.
    """

    def __init__(self, cell_size=0.1, x_range=(-50, 50), z_range=(-50, 50)):
        self.cell_size = cell_size
        self.x_range = x_range
        self.z_range = z_range

        # Grid setup
        self.x_bins = np.arange(x_range[0], x_range[1] + cell_size, cell_size)
        self.z_bins = np.arange(z_range[0], z_range[1] + cell_size, cell_size)
        self.grid_x, self.grid_z = np.meshgrid(
            (self.x_bins[:-1] + self.x_bins[1:]) / 2,
            (self.z_bins[:-1] + self.z_bins[1:]) / 2
        )

        # Reconstruction parameters
        self.energy_resolution = 0.05  # 5% energy resolution
        self.position_uncertainty = 0.5  # mm
        self.min_energy_threshold = 10  # keV minimum energy

    def read_deposits(self, filename, max_lines=10000):
        """Read energy deposit data"""
        deposits = {}
        with open(filename, 'r') as file:
            for i, line in enumerate(file):
                if i >= max_lines:
                    break
                parts = line.split()
                eventID = int(parts[0])
                E_lost = float(parts[1])
                E_remaines = float(parts[2])
                deposits[eventID] = (E_lost, E_remaines)
        return deposits

    def read_gamma_coordinates(self, filename, max_lines=10000):
        """Read gamma-ray coordinate data"""
        coordinates = {}
        with open(filename, 'r') as file:
            for i, line in enumerate(file):
                if i >= max_lines:
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

    def calculate_scattering_angle_with_uncertainty(self, E_remaines, E_lost, n_samples=50):
        """
        Calculate scattering angle with Monte Carlo uncertainty propagation (optimized)
        """
        # Apply energy resolution (vectorized)
        E_lost_samples = np.random.normal(E_lost, E_lost * self.energy_resolution, n_samples)
        E_remaines_samples = np.random.normal(E_remaines, E_remaines * self.energy_resolution, n_samples)

        # Vectorized calculations
        valid_mask = (E_lost_samples > 0) & (E_remaines_samples > 0)
        E_lost_valid = E_lost_samples[valid_mask]
        E_remaines_valid = E_remaines_samples[valid_mask]

        if len(E_lost_valid) == 0:
            return None, None

        E_initial = E_lost_valid + E_remaines_valid
        valid_mask2 = E_initial > 0

        if np.sum(valid_mask2) == 0:
            return None, None

        E_lost_final = E_lost_valid[valid_mask2]
        E_initial_final = E_initial[valid_mask2]

        cos_theta = 1 - (511 * E_lost_final) / (E_initial_final * (E_initial_final - E_lost_final))
        angle_mask = (cos_theta >= -1) & (cos_theta <= 1)

        if np.sum(angle_mask) == 0:
            return None, None

        valid_cos_theta = cos_theta[angle_mask]
        valid_angles = np.degrees(np.arccos(valid_cos_theta))

        return np.mean(valid_angles), np.std(valid_angles)

    def calculate_likelihood_map(self, V, direction, theta_mean, theta_std, weight=1.0):
        """
        Calculate likelihood map for a single event using vectorized approach
        """
        # Normalize direction vector
        direction_norm = direction / np.linalg.norm(direction)

        # Create grid points (vectorized)
        grid_points = np.stack([self.grid_x, np.zeros_like(self.grid_x), self.grid_z], axis=2)

        # Vector from scatter point to all potential sources (vectorized)
        source_vectors = grid_points - V.reshape(1, 1, 3)
        source_distances = np.linalg.norm(source_vectors, axis=2)

        # Avoid division by zero
        valid_mask = source_distances > 1e-10
        source_vectors_norm = np.zeros_like(source_vectors)
        source_vectors_norm[valid_mask] = source_vectors[valid_mask] / source_distances[valid_mask, np.newaxis]

        # Calculate angles (vectorized)
        cos_angles = np.sum(source_vectors_norm * direction_norm, axis=2)
        cos_angles = np.clip(cos_angles, -1, 1)
        calculated_angles = np.degrees(np.arccos(cos_angles))

        # Calculate likelihood
        if theta_std > 0:
            likelihood_map = stats.norm.pdf(calculated_angles, theta_mean, theta_std)
        else:
            likelihood_map = np.where(np.abs(calculated_angles - theta_mean) < 1.0, 1.0, 0.0)

        # Apply distance weighting and validity mask
        distance_weights = 1.0 / (source_distances**2 + 1e-6)
        likelihood_map = likelihood_map * distance_weights * weight * valid_mask

        return likelihood_map

    def estimate_event_quality(self, E_lost, E_remaines, V, X2_pos):
        """
        Estimate the quality/reliability of an event for weighting
        """
        # Energy quality factor
        total_energy = E_lost + E_remaines
        energy_quality = 1.0 / (1.0 + abs(total_energy - 511) / 511)  # Prefer 511 keV events

        # Geometry quality factor
        scatter_distance = np.linalg.norm(np.array(X2_pos) - np.array(V))
        geometry_quality = 1.0 / (1.0 + scatter_distance / 50.0)  # Prefer closer interactions

        # Energy ratio quality
        energy_ratio = min(E_lost, E_remaines) / max(E_lost, E_remaines)
        ratio_quality = energy_ratio  # Prefer more balanced energy splits

        return energy_quality * geometry_quality * ratio_quality

    def reconstruct_oe(self, deposits_file, coordinates_file, max_events=5000):
        """
        Main OE reconstruction method (optimized)
        """
        print("=" * 60)
        print("ORIGIN ENSEMBLE (OE) RECONSTRUCTION")
        print("=" * 60)

        # Load data
        deposits = self.read_deposits(deposits_file, max_events)
        coordinates = self.read_gamma_coordinates(coordinates_file, max_events)

        common_events = list(set(deposits.keys()) & set(coordinates.keys()))
        print(f"События в deposits: {len(deposits)}")
        print(f"События в coordinates: {len(coordinates)}")
        print(f"Общие события: {len(common_events)}")

        # Pre-filter events for better performance
        valid_events = []
        for eventID in common_events:
            E_lost, E_remaines = deposits[eventID]
            if (E_lost >= self.min_energy_threshold and
                E_remaines >= self.min_energy_threshold and
                E_lost <= E_remaines):
                valid_events.append(eventID)

        print(f"Валидные события после фильтрации: {len(valid_events)}")

        # Initialize reconstruction map
        reconstruction_map = np.zeros_like(self.grid_x)
        total_weight = 0.0
        processed_events = 0

        print("\nОбработка событий для OE реконструкции...")

        # Process in batches for better memory management
        batch_size = 100
        for i in tqdm(range(0, len(valid_events), batch_size), desc="OE реконструкция"):
            batch_events = valid_events[i:i + batch_size]

            for eventID in batch_events:
                E_lost, E_remaines = deposits[eventID]
                Vx, Vy, Vz, X2, Y2, Z2 = coordinates[eventID]

                # Calculate scattering angle with uncertainty
                theta_mean, theta_std = self.calculate_scattering_angle_with_uncertainty(E_remaines, E_lost)
                if theta_mean is None:
                    continue

                # Set minimum uncertainty if too small
                if theta_std is None or theta_std < 1.0:
                    theta_std = 2.0

                # Calculate event weight based on quality
                V = np.array([Vx, Vy, Vz])
                X2_pos = [X2, Y2, Z2]
                event_weight = self.estimate_event_quality(E_lost, E_remaines, V, X2_pos)

                # Calculate direction of scattered gamma
                direction = np.array([X2 - Vx, Y2 - Vy, Z2 - Vz])

                # Calculate likelihood map for this event
                event_likelihood = self.calculate_likelihood_map(
                    V, direction, theta_mean, theta_std, event_weight
                )

                # Add to reconstruction
                reconstruction_map += event_likelihood
                total_weight += event_weight
                processed_events += 1

        print(f"\nОбработано событий: {processed_events}")
        print(f"Общий вес: {total_weight:.2f}")

        # Normalize reconstruction map
        if total_weight > 0:
            reconstruction_map /= total_weight

        return reconstruction_map

    def apply_statistical_filters(self, reconstruction_map, percentile_threshold=95):
        """
        Apply statistical filtering to remove noise and enhance signal
        """
        # Calculate threshold based on percentile
        threshold = np.percentile(reconstruction_map.flatten(), percentile_threshold)

        # Apply threshold
        filtered_map = reconstruction_map.copy()
        filtered_map[filtered_map < threshold] = 0

        # Apply Gaussian smoothing for noise reduction
        from scipy import ndimage
        smoothed_map = ndimage.gaussian_filter(filtered_map, sigma=1.0)

        return smoothed_map

    def calculate_reconstruction_metrics(self, reconstruction_map):
        """
        Calculate quality metrics for the reconstruction
        """
        # Find peak location
        max_idx = np.unravel_index(np.argmax(reconstruction_map), reconstruction_map.shape)
        peak_x = self.grid_x[max_idx]
        peak_z = self.grid_z[max_idx]
        peak_value = reconstruction_map[max_idx]

        # Calculate FWHM
        x_profile = np.sum(reconstruction_map, axis=0)
        z_profile = np.sum(reconstruction_map, axis=1)

        def calculate_fwhm_profile(profile, bins):
            if np.max(profile) <= 0:
                return None
            half_max = np.max(profile) / 2
            indices = np.where(profile >= half_max)[0]
            if len(indices) == 0:
                return None
            return bins[indices[-1]] - bins[indices[0]]

        fwhm_x = calculate_fwhm_profile(x_profile, self.x_bins[:-1])
        fwhm_z = calculate_fwhm_profile(z_profile, self.z_bins[:-1])

        # Calculate contrast-to-noise ratio
        signal_region = reconstruction_map[reconstruction_map > 0.8 * peak_value]
        noise_region = reconstruction_map[reconstruction_map < 0.1 * peak_value]

        if len(noise_region) > 0 and np.std(noise_region) > 0:
            cnr = (np.mean(signal_region) - np.mean(noise_region)) / np.std(noise_region)
        else:
            cnr = float('inf')

        return {
            'peak_position': (peak_x, peak_z),
            'peak_value': peak_value,
            'fwhm_x': fwhm_x,
            'fwhm_z': fwhm_z,
            'contrast_to_noise_ratio': cnr
        }

    def visualize_results(self, reconstruction_map, filtered_map=None, save_prefix="oe_reconstruction"):
        """
        Visualize reconstruction results
        """
        fig, axes = plt.subplots(1, 2 if filtered_map is not None else 1, figsize=(15, 6))
        if filtered_map is None:
            axes = [axes]

        # Original reconstruction
        im1 = axes[0].imshow(
            reconstruction_map.T,
            origin='lower',
            extent=[self.x_range[0], self.x_range[1], self.z_range[0], self.z_range[1]],
            aspect='equal',
            cmap='hot',
            interpolation='bilinear'
        )
        axes[0].set_xlabel('X (mm)', fontsize=14)
        axes[0].set_ylabel('Z (mm)', fontsize=14)
        axes[0].set_title('OE Reconstruction (Raw)', fontsize=16)
        plt.colorbar(im1, ax=axes[0], label='Likelihood Density')

        # Filtered reconstruction
        if filtered_map is not None:
            im2 = axes[1].imshow(
                filtered_map.T,
                origin='lower',
                extent=[self.x_range[0], self.x_range[1], self.z_range[0], self.z_range[1]],
                aspect='equal',
                cmap='hot',
                interpolation='bilinear'
            )
            axes[1].set_xlabel('X (mm)', fontsize=14)
            axes[1].set_ylabel('Z (mm)', fontsize=14)
            axes[1].set_title('OE Reconstruction (Filtered)', fontsize=16)
            plt.colorbar(im2, ax=axes[1], label='Likelihood Density')

        plt.tight_layout()
        plt.savefig(f'{save_prefix}_comparison.png', dpi=300, bbox_inches='tight')
        plt.show()

        # Save final image
        final_map = filtered_map if filtered_map is not None else reconstruction_map
        if np.max(final_map) > 0:
            normalized_map = (final_map.T / np.max(final_map.T) * 255).astype('uint8')
        else:
            normalized_map = np.zeros_like(final_map.T, dtype='uint8')

        imageio.imwrite(f'{save_prefix}_final.jpg', normalized_map)

        return final_map


def main():
    """
    Main execution function
    """
    # Initialize reconstructor
    reconstructor = OriginEnsembleReconstructor(
        cell_size=0.1,
        x_range=(-50, 50),
        z_range=(-50, 50)
    )

    # File paths
    deposits_file = '/home/zas/CERN/SPECT/build/deposits_FirstDet(1source).txt'
    coordinates_file = '/home/zas/CERN/SPECT/build/GammaCoordinates_FirstDet(1source).txt'

    # Perform OE reconstruction
    reconstruction_map = reconstructor.reconstruct_oe(deposits_file, coordinates_file, max_events=5000)

    # Apply statistical filtering
    filtered_map = reconstructor.apply_statistical_filters(reconstruction_map, percentile_threshold=90)

    # Calculate metrics
    raw_metrics = reconstructor.calculate_reconstruction_metrics(reconstruction_map)
    filtered_metrics = reconstructor.calculate_reconstruction_metrics(filtered_map)

    # Display results
    print("\n" + "=" * 60)
    print("РЕЗУЛЬТАТЫ OE РЕКОНСТРУКЦИИ:")
    print("=" * 60)

    print("\nСырая реконструкция:")
    print(f"Пик в позиции: X={raw_metrics['peak_position'][0]:.2f}, Z={raw_metrics['peak_position'][1]:.2f}")
    print(f"Значение пика: {raw_metrics['peak_value']:.6f}")
    if raw_metrics['fwhm_x'] and raw_metrics['fwhm_z']:
        print(f"FWHM X: {raw_metrics['fwhm_x']:.2f} мм")
        print(f"FWHM Z: {raw_metrics['fwhm_z']:.2f} мм")
    print(f"Контраст к шуму: {raw_metrics['contrast_to_noise_ratio']:.2f}")

    print("\nФильтрованная реконструкция:")
    print(f"Пик в позиции: X={filtered_metrics['peak_position'][0]:.2f}, Z={filtered_metrics['peak_position'][1]:.2f}")
    print(f"Значение пика: {filtered_metrics['peak_value']:.6f}")
    if filtered_metrics['fwhm_x'] and filtered_metrics['fwhm_z']:
        print(f"FWHM X: {filtered_metrics['fwhm_x']:.2f} мм")
        print(f"FWHM Z: {filtered_metrics['fwhm_z']:.2f} мм")
    print(f"Контраст к шуму: {filtered_metrics['contrast_to_noise_ratio']:.2f}")

    print("=" * 60)

    # Visualize results
    final_map = reconstructor.visualize_results(reconstruction_map, filtered_map)

    return reconstruction_map, filtered_map, raw_metrics, filtered_metrics


if __name__ == "__main__":
    reconstruction_map, filtered_map, raw_metrics, filtered_metrics = main()
