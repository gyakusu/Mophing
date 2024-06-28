#In[0]:
# !maturin develop --release

#In[0]:
import sys
import numpy as np
import pyvista as pv
from PyQt5 import QtCore, QtWidgets
from pyqtgraph import opengl as gl
from pyqtgraph import mkColor
from matplotlib import pyplot as plt

from morphing.morphing import CageParameter, Brg

def map_stress_to_color(stress_values, cmap='jet'):
    norm = plt.Normalize(vmin=stress_values.min(), vmax=stress_values.max())
    return plt.get_cmap(cmap)(norm(stress_values))

#In[0]:
# hoge = Brg('../data/FEM/FEM2_solved.vtu', '../data/FEM/index2.xml', CageParameter(2.35e-3, 2.85e-3, 0.93e-3, 2.10e-3, 0.2, 0.840e-3, 1.70e-3, 1.20e-3, 2.45e-3, 0.2))
# index, face = hoge.extract_surface_as_list()

# np.array(index)[:5], np.array(face)[:5]

#In[0]:

#In[0]:
def file2mesh_item(filename, whiteness=0.8, opacity=0.3):
    mesh = pv.read(filename)
    points = np.array(mesh.points)
    faces = np.array(mesh.faces).reshape((-1, 4))[:, 1:]
    colors = np.ones((points.shape[0], 4)) * whiteness
    colors[:, 3] = opacity
    mesh_item = gl.GLMeshItem(vertexes=points, faces=faces, vertexColors=colors, glOptions='translucent')
    
    edges = mesh.extract_feature_edges()
    edge_coords = [edges.points[edge] for edge in edges.lines.reshape(-1, 3)[:, 1:]]
    edge_items  = [gl.GLLinePlotItem(pos=edge_coord, color=(0.0, 0.0, 0.0, opacity), antialias=True) for edge_coord in edge_coords]
    edge_items.append(mesh_item)
    
    return edge_items

def sphire_item(radius=0.76e-3, center=[0.0, 2.645e-3, 2.0e-3]):
    md = gl.MeshData.sphere(rows=20, cols=20, radius=radius)
    sphere_item = gl.GLMeshItem(meshdata=md, color=(0.9, 0.9, 0.9, 0.9), glOptions='translucent')
    sphere_item.translate(*center)
    return sphere_item

def file2arrow_item(filename, scale=2e-5):
    mesh = pv.read(filename)
    points = np.array(mesh.points)
    velocity = np.array(mesh['v'])
    vnorm = np.linalg.norm(velocity, axis=1)
    colors = np.repeat(map_stress_to_color(vnorm), 2, axis=0)

    vel2pos = velocity * scale
    vectors = np.stack([points, points + vel2pos], axis=1).reshape((-1, 3))
    arrow_item = gl.GLLinePlotItem(pos=vectors, color=colors, width=2.0, antialias=True, mode='lines')
    return arrow_item, [vel2pos, colors]

#In[0]:

#In[0]:

#In[0]:
class ControlGroup(QtWidgets.QWidget):
    def __init__(self, label, min_max, value_label, parent=None):
        super().__init__(parent)
        
        # Create the control bar
        self.control_bar = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.control_bar.setMinimum(min_max[0])
        self.control_bar.setMaximum(min_max[1])
        self.control_bar.setValue((min_max[0] + min_max[1]) // 2)
        self.control_bar.sliderReleased.connect(self.update_value)
        self.control_bar.setFixedSize(200, 30)
        self.control_bar.update()

        # Create the label
        self.label = QtWidgets.QLabel(self)
        self.label.setFixedSize(100, 30)
        self.label.setText(label)

        # Create the value label
        self.vl = value_label
        self.value_label = QtWidgets.QLabel(self)
        self.value_label.setFixedSize(100, 30)
        self.value_label.setText(self.vl(self.control_bar.value()))

        # Arrange the widgets in a layout
        layout = QtWidgets.QHBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.control_bar)
        layout.addWidget(self.value_label)
        self.setLayout(layout)

    def update_value(self):
        self.value_label.setText(self.vl(self.control_bar.value()))

    def get_value(self):
        return self.control_bar.value()

    def on_value_changed(self):
        self.value_changed.emit()
        
class Window(QtWidgets.QMainWindow):
    def __init__(self, control_groups, mesh_items, fem, fem_params, cfd, cfd_params):
        super().__init__()

        self.control_groups = control_groups
        self.fem = fem
        self.cfd = cfd
        self.cfd_params = cfd_params
        self.fem_params = fem_params
        
        self.total_label = QtWidgets.QLabel(self)
        self.total_label.setFixedSize(400, 30)
        self.total_label.setText(f'fluid torque:\t0 [Nmm]\nmises stress:\t0 [MPa]\n')
        
        # Create a 3D graphics widget
        self.view = gl.GLViewWidget(self)
        self.view.setFixedSize(800, 800)
        self.view.setCameraPosition(distance=1e-2, azimuth=30, elevation=-60)
        self.view.setBackgroundColor(mkColor(255, 220, 255))
        
        # メッシュアイテムの作成
        verts  = np.array(self.fem.get_points_as_list())[self.fem_params[0]]
        faces  = self.fem_params[1]
        vertex_stress_values = self.fem_params[2][self.fem_params[0]]
        colors = map_stress_to_color(vertex_stress_values)
        mesh_item = gl.GLMeshItem(vertexes=verts, faces=faces, vertexColors=colors, drawEdges=True)
        
        self.view.addItem(mesh_item)
        
        for mesh_item in mesh_items:
            self.view.addItem(mesh_item)
        
        # Create a splitter to hold the control bar and the view
        splitter = QtWidgets.QSplitter(self)
        splitter.addWidget(self.view)
        
        layout = QtWidgets.QVBoxLayout()  # Change here
        layout.addWidget(self.total_label)
        
        # Create the deformation toggle button
        h_layout = QtWidgets.QHBoxLayout()
        
        self.deformation_button = QtWidgets.QCheckBox("変形", self)
        self.deformation_button.setCheckable(True)
        self.deformation_button.setChecked(False)
        self.deformation_button.clicked.connect(self.update_value)
        
        # Create the fluid toggle button
        self.fluid_button = QtWidgets.QCheckBox("流体", self)
        self.fluid_button.setCheckable(True)
        self.fluid_button.setChecked(True)
        self.fluid_button.clicked.connect(self.update_value)

        # Add the deformation button and indicator lamp to the layout
        h_layout.addWidget(self.deformation_button)
        h_layout.addWidget(self.fluid_button)
        layout.addLayout(h_layout)
        
        for control_group in self.control_groups:
            layout.addWidget(control_group)
            control_group.control_bar.sliderReleased.connect(self.update_value)
        
        control_widget = QtWidgets.QWidget()
        control_widget.setLayout(layout)
        splitter.addWidget(control_widget)
            
        self.setCentralWidget(splitter)

    def update_value(self):
        cp = [control_group.get_value() for control_group in self.control_groups]
        for i in [0,1,2,3,5,6,7,8]:
            cp[i] *= 1e-6
        cp[4] *= 1e-2 * (cp[1] - cp[0])
        cp[9] *= 1e-2 * (cp[8] - cp[3])
        
        self.fem.reload_cage(CageParameter(*cp), 20, 0)
        self.cfd.reload_cage(CageParameter(*cp), 20, 10)
        
        verts  = np.array(self.fem.get_points_as_list())[self.fem_params[0]]
        faces  = self.fem_params[1]
        vertex_stress_values = self.fem_params[2][self.fem_params[0]]
        colors = map_stress_to_color(vertex_stress_values)
        
        if self.deformation_button.isChecked():
            deformation = self.fem_params[3][self.fem_params[0]]
            self.view.items[0].setMeshData(vertexes=verts+deformation, faces=faces, vertexColors=colors)
        else:
            self.view.items[0].setMeshData(vertexes=verts, faces=faces, vertexColors=colors)        
                
        cfd_points = np.array(self.cfd.get_points_as_list())
        
        if self.fluid_button.isChecked():
            new_vectors = np.stack([cfd_points, cfd_points + self.cfd_params[0]], axis=1).reshape((-1, 3))
        else:
            new_vectors = np.stack([cfd_points, cfd_points + 0], axis=1).reshape((-1, 3))
        
        self.view.items[-1] = gl.GLLinePlotItem(pos=new_vectors, color=self.cfd_params[1], width=2.0, antialias=True, mode='lines')
        
        total = sum([control_group.get_value() for control_group in self.control_groups])
        self.total_label.setText(f'fluid torque:\t{total} [Nmm]\nmises stress:\t{total*2} [MPa]\n')
        
if __name__ == '__main__':
    
    app = QtWidgets.QApplication(sys.argv)
    
    labels = ['内径', '外径', '底面高さ', '肩高さ', '底面取り', 'ポケット径', '肩中心高さ', '肩径', '爪高さ', '爪面取り']
    means = [2.35e-3, 2.85e-3, 0.93e-3, 2.10e-3, 0.2, 0.840e-3, 1.70e-3, 1.20e-3, 2.45e-3, 0.2]
    diffs = [0.10e-3, 0.10e-3, 0.10e-3, 0.10e-3, 0.1, 0.040e-3, 0.02e-3, 0.02e-3, 0.10e-3, 0.1]
    mm0 = lambda i: [int((means[i] - diffs[i])*1e6), int((means[i] + diffs[i])*1e6)]
    mm1 = lambda i: [int((means[i] - diffs[i])*1e2), int((means[i] + diffs[i])*1e2)]
    vl0 = lambda v: f'{v * 1e-3:.4g} [mm]'
    vl1 = lambda v: f'{v * 1e-0:.4g} [%]'
    min_maxs     = [mm0, mm0, mm0, mm0, mm1, mm0, mm0, mm0, mm0, mm1]
    value_labels = [vl0, vl0, vl0, vl0, vl1, vl0, vl0, vl0, vl0, vl1]
    control_groups = [ControlGroup(labels[i], min_maxs[i](i), value_labels[i]) for i in range(10)]
    
    cage = CageParameter(*means)
    fem = Brg('data/FEM/FEM2_solved.vtu', "data/FEM/index2.xml", cage)
    cfd = Brg('data/CFD/star_master4_rotated.vtu', "data/CFD/index4.xml", cage)
    index, face = fem.extract_surface_as_list()
    von_mises_stress = pv.read('data/FEM/FEM2_solved.vtu')['von_mises_stress']
    deformation = pv.read('data/FEM/FEM2_solved.vtu')['deformation']
    fem_params = [np.array(param) for param in [index, face, von_mises_stress, deformation]]
    
    sphire = sphire_item()
    outer  = file2mesh_item('data/A_Cut_Scaled.vtk', opacity=0.7)
    inner  = file2mesh_item('data/B_Cut_Scaled.vtk')
    arrow, cfd_params = file2arrow_item('data/CFD/star_master4_rotated.vtu')
    items = []
    items.extend([sphire])
    items.extend(outer)
    items.extend(inner)
    items.extend([arrow])
    
    window = Window(control_groups, items, fem, fem_params, cfd, cfd_params)
    window.show()
    sys.exit(app.exec_())
    
