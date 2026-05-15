import numpy as np
import os
from scipy.optimize import curve_fit
from scipy.ndimage import center_of_mass

mec2 = 511.0  
E0 = 140.0    
GRID_SIZE = 50
RANGE = (-25, 25)

def load_data(deposits_file, coords_file, max_lines=10000):
    events = {}
    with open(deposits_file, 'r') as f:
        for i, line in enumerate(f):
            if i >= max_lines: break
            event_id, e1, e2 = map(float, line.split())
            events[event_id] = {'E1': e1, 'E2': e2}
    
    with open(coords_file, 'r') as f:
        for i, line in enumerate(f):
            if i >= max_lines: break
            parts = list(map(float, line.split()))
            event_id = parts[0]
            if event_id in events:
                events[event_id].update({
                    'x1': parts[1], 'y1': parts[2], 'z1': parts[3],
                    'x2': parts[4], 'y2': parts[5], 'z2': parts[6]
                })
    
    print(f"Loaded {len(events)} events")
    return events

def calculate_theta(e1, e2, e0=E0):
    if e1 <= 0 or e2 <= 0 or (e1 + e2) > (e0 + 1):
        return None
    
    e_scattered = e0 - e1
    cos_theta = 1 - mec2 * (1/e_scattered - 1/e0)
    
    if not (-1 <= cos_theta <= 1):
        return None
        
    theta = np.arccos(cos_theta)
    return theta if theta <= np.pi/2 else None

def determine_primary_axis(events, tolerance=1.0):
    coords = {
        'x': [e['x1'] for e in events.values()],
        'y': [e['y1'] for e in events.values()],
        'z': [e['z1'] for e in events.values()]
    }
    
    stds = {axis: np.std(vals) for axis, vals in coords.items()}
    fixed_axes = [axis for axis, std in stds.items() if std < tolerance]
    
    if len(fixed_axes) == 1:
        return fixed_axes[0]
    raise ValueError("Cannot determine primary axis")

def backprojection_3d(events):
    primary_axis = determine_primary_axis(events)
    print(f"Primary axis: {primary_axis.upper()}")

    voxel_size = (RANGE[1] - RANGE[0]) / GRID_SIZE
    x = np.linspace(RANGE[0] + voxel_size/2, RANGE[1] - voxel_size/2, GRID_SIZE)
    y = np.linspace(RANGE[0] + voxel_size/2, RANGE[1] - voxel_size/2, GRID_SIZE)
    z = np.linspace(RANGE[0] + voxel_size/2, RANGE[1] - voxel_size/2, GRID_SIZE)
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    image = np.zeros((GRID_SIZE, GRID_SIZE, GRID_SIZE), dtype=np.int32)
    valid_events = 0

    for event_id, data in events.items():
        if not all(key in data for key in ['E1','E2','x1','y1','z1','x2','y2','z2']):
            continue
            
        e1, e2 = data['E1'], data['E2']
        x1, y1, z1 = data['x1'], data['y1'], data['z1']
        x2, y2, z2 = data['x2'], data['y2'], data['z2']
        
        theta = calculate_theta(e1, e2)
        if theta is None: continue
        
        if primary_axis == 'x':
            point_scatter = np.array([x1, y1, z1])
            point_absorber = np.array([x2, y2, z2])
            points_grid = np.stack([X, Y, Z], axis=-1)
        elif primary_axis == 'y':
            point_scatter = np.array([y1, x1, z1])
            point_absorber = np.array([y2, x2, z2])
            points_grid = np.stack([Y, X, Z], axis=-1)
        else:
            point_scatter = np.array([z1, x1, y1])
            point_absorber = np.array([z2, x2, y2])
            points_grid = np.stack([Z, X, Y], axis=-1)
        
        direction = point_absorber - point_scatter
        length = np.linalg.norm(direction)
        if length < 1e-6: continue
            
        u = direction / length
        
        points = points_grid - point_scatter
        cross = np.cross(points, u)
        dist_to_axis = np.linalg.norm(cross, axis=-1)
        dist_along_axis = np.dot(points, u)
        
        r = np.abs(dist_along_axis) * np.tan(theta)
        max_dist = np.linalg.norm(point_scatter)
        
        mask = (
            (np.abs(dist_to_axis - r) < 0.2) & 
            (np.abs(dist_along_axis) <= max_dist) & 
            (X >= RANGE[0]) & (X <= RANGE[1]) & 
            (Y >= RANGE[0]) & (Y <= RANGE[1]) & 
            (Z >= RANGE[0]) & (Z <= RANGE[1])
        )
        
        image[mask] += 1
        valid_events += 1
        
    print(f"Processed events: {valid_events}/{len(events)}")
    return image

def save_volume_data(data, filename="/home/zas/CERN/SPECT/build/voxels_data.txt"):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    data_3d = np.transpose(data, (0, 1, 2)) 
    np.savetxt(filename, data_3d.ravel(), fmt='%d')
    print(f"Data saved to {filename}")

def find_max_voxel_coordinates(image):

    max_index = np.unravel_index(np.argmax(image), image.shape)
    max_value = image[max_index]
    voxel_size = (RANGE[1] - RANGE[0]) / GRID_SIZE
    x = RANGE[0] + voxel_size * (max_index[0] + 0.5)
    y = RANGE[0] + voxel_size * (max_index[1] + 0.5)
    z = RANGE[0] + voxel_size * (max_index[2] + 0.5)
    print(f"Voxel with maximum intersections: ({x:.2f}, {y:.2f}, {z:.2f}), Count: {max_value}")
    return (x, y, z), max_value

def gaussian_3d(coords, amplitude, x0, y0, z0, sigma_x, sigma_y, sigma_z):

    x, y, z = coords
    return amplitude * np.exp(-((x - x0) ** 2 / (2 * sigma_x ** 2) +
                                 (y - y0) ** 2 / (2 * sigma_y ** 2) +
                                 (z - z0) ** 2 / (2 * sigma_z ** 2)))

def fit_gaussian_and_calculate_fwhm(image, threshold=0.1):

    voxel_size = (RANGE[1] - RANGE[0]) / GRID_SIZE
    x = np.linspace(RANGE[0] + voxel_size/2, RANGE[1] - voxel_size/2, GRID_SIZE)
    y = np.linspace(RANGE[0] + voxel_size/2, RANGE[1] - voxel_size/2, GRID_SIZE)
    z = np.linspace(RANGE[0] + voxel_size/2, RANGE[1] - voxel_size/2, GRID_SIZE)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    threshold_value = threshold * np.max(image)
    mask = image > threshold_value
    filtered_data = image[mask].ravel()
    coords = np.vstack((X[mask].ravel(), Y[mask].ravel(), Z[mask].ravel()))

    max_index = np.unravel_index(np.argmax(image), image.shape)
    initial_guess = (np.max(filtered_data), 
                     x[max_index[0]], y[max_index[1]], z[max_index[2]], 
                     1.0, 1.0, 1.0)


    popt, _ = curve_fit(gaussian_3d, coords, filtered_data, p0=initial_guess)

    fwhm_x = 2 * np.sqrt(2 * np.log(2)) * popt[4]
    fwhm_y = 2 * np.sqrt(2 * np.log(2)) * popt[5]
    fwhm_z = 2 * np.sqrt(2 * np.log(2)) * popt[6]

    print(f"Fitted Gaussian parameters: Amplitude={popt[0]:.2f}, "
          f"Center=({popt[1]:.2f}, {popt[2]:.2f}, {popt[3]:.2f}), "
          f"Sigma=({popt[4]:.2f}, {popt[5]:.2f}, {popt[6]:.2f})")
    print(f"FWHM: X={fwhm_x:.2f}, Y={fwhm_y:.2f}, Z={fwhm_z:.2f}")

    return popt, (fwhm_x, fwhm_y, fwhm_z)

if __name__ == "__main__":

    deposits_x = '/home/zas/CERN/SPECT/build/deposits_SecondDet1.txt'
    coords_x = '/home/zas/CERN/SPECT/build/GammaCoordinates_SecondDet1.txt'
    deposits_y = '/home/zas/CERN/SPECT/build/deposits_FirstDet1.txt'
    coords_y = '/home/zas/CERN/SPECT/build/GammaCoordinates_FirstDet1.txt'


    events_x = load_data(deposits_x, coords_x)
    events_y = load_data(deposits_y, coords_y)
    
    image_x = backprojection_3d(events_x)
    image_y = backprojection_3d(events_y)
    

    # combined_image = image_x + image_y

    combined_image = image_y
    save_volume_data(combined_image)
    
    max_coords, max_value = find_max_voxel_coordinates(combined_image)

    fit_gaussian_and_calculate_fwhm(combined_image)