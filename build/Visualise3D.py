from mayavi import mlab
import numpy as np
from traits.api import HasTraits, Bool, Instance, on_trait_change, Range
from traitsui.api import View, Item, Group
from mayavi.core.ui.api import MayaviScene, SceneEditor, MlabSceneModel
import matplotlib.pyplot as plt

# Границы осей
x_min = -500.0
x_max = 500.0
y_min = -2000.0
y_max = 2000.0
z_min = -1000.0
z_max = 1400.0

# Загрузка данных
data = np.loadtxt("/home/zas/CERN/SPECT/build/voxels_data(1detY).txt")
data = data.reshape((50, 50, 50))


class SliceVisualizer(HasTraits):
    show_slice_x = Bool(False)
    show_slice_y = Bool(False)
    show_slice_z = Bool(False)
    show_volume = Bool(True)

    # Настройка контраста
    vmin = Range(low=float(np.min(data)), high=float(np.max(data)), value=float(np.min(data)))
    vmax = Range(low=float(np.min(data)), high=float(np.max(data)), value=float(np.max(data)))

    # Настройка прозрачности
    opacity = Range(low=0.0, high=1.0, value=0.2)

    # Индексы срезов
    Xl = Range(low=0, high=data.shape[0] - 1, value=0)
    Yl = Range(low=0, high=data.shape[1] - 1, value=0)
    Zl = Range(low=0, high=data.shape[2] - 1, value=0)

    # Сцена Mayavi
    scene = Instance(MlabSceneModel, ())

    def __init__(self):
        super(SliceVisualizer, self).__init__()
        self.src = None
        self.contour = None
        self.data = data
        self.axas = None

        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.z_min = z_min
        self.z_max = z_max

        self.max = np.max(data)

    def plot_slice(self, axis, index):

        if axis not in ['x', 'y', 'z']:
            raise ValueError("Параметр axis должен быть 'x', 'y' или 'z'.")
    
        if axis == 'x':
            if index < 0 or index >= self.data.shape[0]:
                raise ValueError(f"Индекс {index} выходит за границы данных по оси X (0-{self.data.shape[0]-1}).")
            slice_data = self.data[index, :, :]
            title = f"Срез по X "
            xlabel = "Y"
            ylabel = "Z"
        elif axis == 'y':
            if index < 0 or index >= self.data.shape[1]:
                raise ValueError(f"Индекс {index} выходит за границы данных по оси Y (0-{self.data.shape[1]-1}).")
            slice_data = self.data[:, index, :]
            title = f"Срез по Y "
            xlabel = "X"
            ylabel = "Z"
        elif axis == 'z':
            if index < 0 or index >= self.data.shape[2]:
                raise ValueError(f"Индекс {index} выходит за границы данных по оси Z (0-{self.data.shape[2]-1}).")
            slice_data = self.data[:, :, index]
            title = f"Срез по Z"
            xlabel = "X"
            ylabel = "Y"

        shape = slice_data.shape
        extent_x = [-shape[1] // 2, shape[1] // 2, -shape[0] // 2, shape[0] // 2]
        

        plt.figure(figsize=(8, 6))
        plt.imshow(slice_data, cmap='viridis', origin='lower', extent=extent_x )
        cbar = plt.colorbar(label='Количество пересечений')
        cbar.ax.tick_params(labelsize=24)
        cbar.set_label('Количество пересечений', fontsize=24)   
        plt.title(title, fontsize=24)
        plt.xlabel(xlabel, fontsize=24)
        plt.ylabel(ylabel, fontsize=24)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.show()
        
        max_intensity = np.max(slice_data)
        half_max = max_intensity / 2   

        central_row_profile = slice_data[slice_data.shape[0] // 2, :]  
        above_half_max_row = central_row_profile >= half_max
        edges_row = np.where(np.diff(above_half_max_row.astype(int)) != 0)[0]  

        if len(edges_row) >= 2:
            fwhm_row = edges_row[-1] - edges_row[0]  
            print(f"FWHM по горизонтальной оси для среза по {axis} (индекс={index}): {fwhm_row} пикселей")
        
        central_col_profile = slice_data[:, slice_data.shape[1] // 2] 
        above_half_max_col = central_col_profile >= half_max
        edges_col = np.where(np.diff(above_half_max_col.astype(int)) != 0)[0]  

        if len(edges_col) >= 2:
            fwhm_col = edges_col[-1] - edges_col[0]  
            print(f"FWHM по вертикальной оси для среза по {axis} (индекс={index}): {fwhm_col} пикселей")
        

    @on_trait_change('scene.activated')
    def init_scene(self):
        if self.data is None:
            raise ValueError("Данные не загружены.")
        
        if self.src is None:
            self.src = mlab.pipeline.scalar_field(self.data)

            # Создание объемного представления
            self.contour = mlab.pipeline.volume(self.src, vmin=self.vmin, vmax=self.vmax)
            self.contour._volume_property.scalar_opacity_unit_distance = self.opacity

            # Добавление плоскостей срезов
            self.slice_x = mlab.pipeline.image_plane_widget(self.src, plane_orientation='x_axes', slice_index=int(self.Xl))
            self.slice_y = mlab.pipeline.image_plane_widget(self.src, plane_orientation='y_axes', slice_index=int(self.Yl))
            self.slice_z = mlab.pipeline.image_plane_widget(self.src, plane_orientation='z_axes', slice_index=int(self.Zl))
            self.slice_x.visible = self.show_slice_x
            self.slice_y.visible = self.show_slice_y
            self.slice_z.visible = self.show_slice_z

            # Настройка вида
            mlab.view(azimuth=0, elevation=90, distance=200, focalpoint=(0, 30, 20), figure=self.scene.mayavi_scene)

            # Обновление осей
            self.axas = mlab.axes(xlabel='X [m]', ylabel='Y [m]', zlabel='Z [m]')
            self.axas.axes.use_ranges = False
            self.axas.axes.ranges = [self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max]

            self.colorbar = mlab.colorbar(title=r'Intensety', orientation='horizontal', label_fmt='%.2f')

    @on_trait_change('Xl')
    def SceneX(self):
        if self.slice_x:
            self.slice_x.ipw.slice_index = int(self.Xl)
            self.scene.mayavi_scene.render()

    @on_trait_change('Yl')
    def SceneY(self):
        if self.slice_y:
            self.slice_y.ipw.slice_index = int(self.Yl)
            self.scene.mayavi_scene.render()

    @on_trait_change('Zl')
    def SceneZ(self):
        if self.slice_z:
            self.slice_z.ipw.slice_index = int(self.Zl)
            self.scene.mayavi_scene.render()

    @on_trait_change('show_slice_x, show_slice_y, show_slice_z')
    def update_slices(self):
        if hasattr(self, 'slice_x'):
            self.slice_x.visible = self.show_slice_x
        if hasattr(self, 'slice_y'):
            self.slice_y.visible = self.show_slice_y
        if hasattr(self, 'slice_z'):
            self.slice_z.visible = self.show_slice_z

    @on_trait_change('vmin, vmax, opacity')
    def update_volume_properties(self):
        if self.contour:
            azimuth, elevation, distance, focalpoint = mlab.view()
            del self.contour.actors
            self.contour = mlab.pipeline.volume(self.src, vmin=self.vmin, vmax=self.vmax)
            mlab.view(azimuth, elevation, distance, focalpoint)
            self.colorbar.data_range = (self.vmin, self.vmax)
            self.contour._volume_property.scalar_opacity_unit_distance = self.opacity
            self.contour.scene.render()

    @on_trait_change('show_volume')
    def update_volume(self):
        if self.contour:
            self.contour.visible = self.show_volume

    view = View(
        Item('scene', editor=SceneEditor(scene_class=MayaviScene), height=400, width=400, show_label=False),
        Group(
            Item('show_volume', label='Показать объем'),
            Item('show_slice_x', label='Показать срез X'),
            Item('show_slice_y', label='Показать срез Y'),
            Item('show_slice_z', label='Показать срез Z'),
            Item('vmin', label='Контраст: мин'),
            Item('vmax', label='Контраст: макс'),
            Item('opacity', label='Прозрачность'),
            Item('Xl', label='Срез X'),
            Item('Yl', label='Срез Y'),
            Item('Zl', label='Срез Z'),
        ),
        resizable=True
    )


if __name__ == '__main__':
    visualizer = SliceVisualizer()
    visualizer.configure_traits()
    visualizer.plot_slice('x', 25)  
    visualizer.plot_slice('y', 25) 
    visualizer.plot_slice('z', 25)  