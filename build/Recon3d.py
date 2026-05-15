import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

mec2 = 511.0  
E0 = 140.0    

def load_data(deposits_file, coords_file, max_lines=10000):
    events = {}
    with open(deposits_file, 'r') as f:
        for i, line in enumerate(f):
            if i >= max_lines:
                break
            event_id, e1, e2 = map(float, line.split())
            events[event_id] = {'E1': e1, 'E2': e2}
    
    with open(coords_file, 'r') as f:
        for i, line in enumerate(f):
            if i >= max_lines:
                break
            parts = list(map(float, line.split()))
            event_id = parts[0]
            if event_id in events:
                events[event_id].update({
                    'x1': parts[1], 
                    'y1': parts[2], 
                    'z1': parts[3],
                    'x2': parts[4], 
                    'y2': parts[5], 
                    'z2': parts[6]
                })
    
    print(f"Loaded {len(events)} events with energies and coordinates")
    return events

def calculate_theta(e1, e2, e0=E0):
    if e1 <= 0 or e2 <= 0 or (e1 + e2) > (e0 + 1):
        return None
    
    e_scattered = e0 - e1
    cos_theta = 1 - mec2 * (1 / e_scattered - 1 / e0)
    
    if not (-1 <= cos_theta <= 1):
        return None
        
    theta = np.arccos(cos_theta)
    if theta > np.pi / 2: 
        return None

    return theta

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
    else:
        raise ValueError("Невозможно определить ось детекторов. Проверьте входные данные.")

def backprojection_3d(events, grid_size=50, range=(-25, 25)):
    primary_axis = determine_primary_axis(events)
    print(f"Основная ось детекторов: {primary_axis.upper()}")

    voxel_size = (range[1] - range[0]) / grid_size
    x = np.linspace(range[0] + voxel_size/2, range[1] - voxel_size/2, grid_size)
    y = np.linspace(range[0] + voxel_size/2, range[1] - voxel_size/2, grid_size)
    z = np.linspace(range[0] + voxel_size/2, range[1] - voxel_size/2, grid_size)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    image = np.zeros((grid_size, grid_size, grid_size))
    valid_events = 0

    for event_id, data in events.items():
        if not all(key in data for key in ['E1', 'E2', 'x1', 'y1', 'z1', 'x2', 'y2', 'z2']):
            continue
            
        e1, e2 = data['E1'], data['E2']
        x1, y1, z1 = data['x1'], data['y1'], data['z1']
        x2, y2, z2 = data['x2'], data['y2'], data['z2']
        
        theta = calculate_theta(e1, e2)
        if theta is None:
            continue
        
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
        if length < 1e-6:
            continue
            
        u = direction / length
        
        if primary_axis == 'x':
            points = np.stack([X, Y, Z], axis=-1) - point_scatter
        else:
            points = np.stack([Y, X, Z], axis=-1) - point_scatter
        
        cross = np.cross(points, u)
        dist_to_axis = np.linalg.norm(cross, axis=-1)
        dist_along_axis = np.dot(points, u)
        
        r = np.abs(dist_along_axis) * np.tan(theta)
        max_dist = np.linalg.norm(point_scatter)
        
        mask = (np.abs(dist_to_axis - r) < 0.2) & \
               (np.abs(dist_along_axis) <= max_dist) & \
               (X >= range[0]) & (X <= range[1]) & \
               (Y >= range[0]) & (Y <= range[1]) & \
               (Z >= range[0]) & (Z <= range[1])
        
        image += mask.astype(int)
        valid_events += 1
        
    if image.max() > 0:
        image /= image.max()
    
    print(f"Обработано событий: {valid_events}/{len(events)}")
    max_pos = np.unravel_index(np.argmax(image), image.shape)
    print(f"Максимум интенсивности: X={x[max_pos[0]]:.1f}, Y={y[max_pos[1]]:.1f}, Z={z[max_pos[2]]:.1f}")
    
    return x, y, z, image

def calculate_fwhm(profile):
    half_max = profile.max() / 2
    indices = np.where(profile >= half_max)[0]
    return indices[-1] - indices[0] if len(indices) > 1 else 0

def combine_images(images):
    combined = np.sum(images, axis=0)
    if combined.max() > 0:
        combined /= combined.max()
    return combined



# def visualize_combined_results(x, y, z, images, primary_axes):
#     # Объединяем изображения
#     combined_image = combine_images(images)
    
#     # Визуализация
#     fig = plt.figure(figsize=(18, 12))
    
#     # 3D визуализация
#     ax3d = fig.add_subplot(221, projection='3d')
#     threshold = 0.1
#     x_idx, y_idx, z_idx = np.where(combined_image >= threshold)
    
#     scatter = ax3d.scatter(
#         x[x_idx], y[y_idx], z[z_idx],
#         c=combined_image[x_idx, y_idx, z_idx],
#         cmap='hot', alpha=0.6, s=20
#     )
    
#     ax_range = [-25, 25]
#     ax3d.set_box_aspect([1, 1, 1])
#     ax3d.set_xlim(ax_range)
#     ax3d.set_ylim(ax_range)
#     ax3d.set_zlim(ax_range)
#     ax3d.set_xlabel('X (мм)')
#     ax3d.set_ylabel('Y (мм)')
#     ax3d.set_zlabel('Z (мм)')
#     plt.colorbar(scatter, ax=ax3d, label='Интенсивность')
    
#     # 2D срезы
#     slices = {
#         'xy': (combined_image[:, :, int(len(z)/2)], [x.min(), x.max()], [y.min(), y.max()], 'XY (Z=0)'),
#         'xz': (combined_image[:, int(len(y)/2), :], [x.min(), x.max()], [z.min(), z.max()], 'XZ (Y=0)'), 
#         'yz': (combined_image[int(len(x)/2), :, :], [y.min(), y.max()], [z.min(), z.max()], 'YZ (X=0)')
#     }
    
#     positions = [222, 223, 224]
#     for pos, (slice_key, (slice_data, x_ext, y_ext, title)) in zip(positions, slices.items()):
#         ax = fig.add_subplot(pos)
#         im = ax.imshow(slice_data.T, extent=[x_ext[0], x_ext[1], y_ext[0], y_ext[1]], 
#                       cmap='hot', origin='lower')
#         ax.set_title(title)
#         ax.set_xlabel('X (мм)' if 'x' in slice_key else 'Y (мм)')
#         ax.set_ylabel('Y (мм)' if 'y' in slice_key else 'Z (мм)')
#         plt.colorbar(im, ax=ax, label='Интенсивность')
    
#     plt.suptitle("Комбинированная реконструкция", y=1.02)
#     plt.tight_layout()
#     plt.show()


def get_intersection_mask(images, threshold=0.1):

    masks = [img > threshold for img in images]
    intersection_mask = np.ones_like(masks[0], dtype=bool)
    for mask in masks:
        intersection_mask &= mask
    return intersection_mask



def save_voxels_for_mayavi_txt(x, y, z, intersection_image, filename="voxel_data.txt"):
    """
    Сохраняет воксели в текстовом формате, подходящем для загрузки и визуализации.
    
    Параметры:
    - x, y, z: координатные массивы (1D)
    - intersection_image: 3D массив интенсивностей
    - filename: имя файла для сохранения (.txt)
    """
    grid_size_x = len(x)
    grid_size_y = len(y)
    grid_size_z = len(z)
    
    with open(filename, 'w') as f:
        # Записываем размерность сетки в начало файла
        f.write(f"{grid_size_x} {grid_size_y} {grid_size_z}\n")
        
        # Записываем координаты вокселей и их интенсивности
        for i in range(grid_size_x):
            for j in range(grid_size_y):
                for k in range(grid_size_z):
                    intensity = intersection_image[i, j, k]
                    if intensity > 0:  # Сохраняем только ненулевые воксели
                        f.write(f"{x[i]:.2f} {y[j]:.2f} {z[k]:.2f} {intensity:.6f}\n")
    
    print(f"Сохранено вокселей в файл {filename}")


def visualize_intersection_results(x, y, z, image_x, image_y):
    threshold = 0.1
    mask_x = image_x > threshold
    mask_y = image_y > threshold
    intersection_mask = np.logical_and(mask_x, mask_y)
    
    intersection_image = np.zeros_like(image_x)
    intersection_image[intersection_mask] = (image_x[intersection_mask] + image_y[intersection_mask]) / 2
    
    if intersection_image.max() > 0:
        intersection_image /= intersection_image.max()
    
    save_voxels_for_mayavi_txt(x, y, z, intersection_image, filename="voxel_data.txt")
    
    #  точка максимальной интенсивности
    max_pos = np.unravel_index(np.argmax(intersection_image), intersection_image.shape)
    max_coords = (x[max_pos[0]], y[max_pos[1]], z[max_pos[2]])
    max_intensity = intersection_image[max_pos]
    
    print(f"\nАнализ пересечения данных:")
    print(f"Максимальная интенсивность: {max_intensity:.3f} в точке ({max_coords[0]:.1f}, {max_coords[1]:.1f}, {max_coords[2]:.1f})")
    
    fwhm_x = calculate_fwhm(intersection_image[:, max_pos[1], max_pos[2]])
    fwhm_y = calculate_fwhm(intersection_image[max_pos[0], :, max_pos[2]])
    fwhm_z = calculate_fwhm(intersection_image[max_pos[0], max_pos[1], :])
    
    # Перевод из индексов в миллиметры
    voxel_size = x[1] - x[0]
    fwhm_x_mm = fwhm_x * voxel_size
    fwhm_y_mm = fwhm_y * voxel_size
    fwhm_z_mm = fwhm_z * voxel_size
    
    print(f"FWHM X: {fwhm_x_mm:.1f} мм")
    print(f"FWHM Y: {fwhm_y_mm:.1f} мм")
    print(f"FWHM Z: {fwhm_z_mm:.1f} мм")
    
    # Визуализация
    fig = plt.figure(figsize=(18, 12))
    
    # 3D визуализация
    ax3d = fig.add_subplot(221, projection='3d')
    display_threshold = 0.1
    x_idx, y_idx, z_idx = np.where(np.logical_and(intersection_image >= display_threshold, intersection_mask))
    
    scatter = ax3d.scatter(
        x[x_idx], y[y_idx], z[z_idx],
        c=intersection_image[x_idx, y_idx, z_idx],
        cmap='hot', alpha=0.6, s=20
    )
    
    ax3d.scatter([max_coords[0]], [max_coords[1]], [max_coords[2]], 
                c='blue', s=100, label='Максимум')
    
    ax_range = [-25, 25]
    ax3d.set_box_aspect([1, 1, 1])
    ax3d.set_xlim(ax_range)
    ax3d.set_ylim(ax_range)
    ax3d.set_zlim(ax_range)
    ax3d.set_xlabel('X (мм)')
    ax3d.set_ylabel('Y (мм)')
    ax3d.set_zlabel('Z (мм)')
    plt.colorbar(scatter, ax=ax3d, label='Норм. интенсивность')
    ax3d.legend()
    
    # 2D срезы с FWHM
    slices = {
        'xy': (intersection_image[:, :, max_pos[2]], [x.min(), x.max()], [y.min(), y.max()], 'XY'),
        'xz': (intersection_image[:, max_pos[1], :], [x.min(), x.max()], [z.min(), z.max()], 'XZ'), 
        'yz': (intersection_image[max_pos[0], :, :], [y.min(), y.max()], [z.min(), z.max()], 'YZ')
    }
    
    positions = [222, 223, 224]
    for pos, (slice_key, (slice_data, x_ext, y_ext, title)) in zip(positions, slices.items()):
        ax = fig.add_subplot(pos)
        
        # Отмечаем FWHM на срезах
        if slice_key == 'xy':
            profile_x = slice_data[:, max_pos[1]]
            profile_y = slice_data[max_pos[0], :]
        elif slice_key == 'xz':
            profile_x = slice_data[:, max_pos[2]]
            profile_z = slice_data[max_pos[0], :]
        else:  # 'yz'
            profile_y = slice_data[:, max_pos[2]]
            profile_z = slice_data[max_pos[1], :]
        
        im = ax.imshow(slice_data.T, extent=[x_ext[0], x_ext[1], y_ext[0], y_ext[1]], 
                      cmap='hot', origin='lower')
        ax.set_title(f"{title} (FWHM X:{fwhm_x_mm:.1f}mm, Y:{fwhm_y_mm:.1f}mm, Z:{fwhm_z_mm:.1f}mm)")
        ax.set_xlabel('X (мм)' if 'x' in slice_key else 'Y (мм)')
        ax.set_ylabel('Y (мм)' if 'y' in slice_key else 'Z (мм)')
        plt.colorbar(im, ax=ax, label='Норм. интенсивность')
    
    plt.suptitle(f"Пересечение данных | Макс: {max_intensity:.2f} @ ({max_coords[0]:.1f},{max_coords[1]:.1f},{max_coords[2]:.1f})", y=1.02)
    plt.tight_layout()
    plt.show()


deposits_x = '/home/zas/CERN/SPECT/build/deposits_SecondDet.txt'
coords_x = '/home/zas/CERN/SPECT/build/GammaCoordinates_SecondDet.txt'
deposits_y = '/home/zas/CERN/SPECT/build/deposits_FirstDet.txt'
coords_y = '/home/zas/CERN/SPECT/build/GammaCoordinates_FirstDet.txt'

events_x = load_data(deposits_x, coords_x)
events_y = load_data(deposits_y, coords_y)

x_x, y_x, z_x, image_x = backprojection_3d(events_x)
x_y, y_y, z_y, image_y = backprojection_3d(events_y)

# visualize_combined_results(x_x, y_x, z_x, [image_x, image_y], ['x', 'y'])

visualize_intersection_results(x_x, y_x, z_x, image_x, image_y)