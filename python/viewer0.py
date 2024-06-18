import sys
from PyQt5 import QtCore, QtWidgets
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
import pyvista as pv
from morphing.morphing import CageParameter, Brg

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
    def __init__(self, control_groups, brg, mesh):
        super().__init__()

        self.brg = brg
        self.mesh = mesh
        
        self.total_label = QtWidgets.QLabel(self)
        self.total_label.setFixedSize(150, 30)
        self.total_label.setText(f'Torque: 0 [Nmm]')
        
        # Create a 3D graphics widget
        self.view = gl.GLViewWidget(self)
        self.view.setFixedSize(800, 800)
        self.view.setCameraPosition(distance=1e-2, azimuth=30, elevation=-60)
        self.view.setBackgroundColor(pg.mkColor(255, 255, 255))
        
        # メッシュアイテムの作成
        surface = self.mesh.extract_geometry()
        verts = np.array(surface.points)
        faces = np.array(surface.faces).reshape((-1, 4))[:, 1:]
        colors = np.array([(i / len(faces), 0, 0, 0.5) for i in range(len(faces))])       
        
        mesh_item = gl.GLMeshItem(vertexes=verts, faces=faces, vertexColors=colors)        
        self.view.addItem(mesh_item)
        
        self.control_groups = control_groups
        
        # Create a splitter to hold the control bar and the view
        splitter = QtWidgets.QSplitter(self)
        splitter.addWidget(self.view)
        splitter.addWidget(self.total_label)
        
        layout = QtWidgets.QVBoxLayout()  # Change here
        for control_group in self.control_groups:
            layout.addWidget(control_group)
        
        control_widget = QtWidgets.QWidget()
        control_widget.setLayout(layout)
        splitter.addWidget(control_widget)
                
        # Set the splitter as the central widget
        self.setCentralWidget(splitter)
        
        for control_group in self.control_groups:
            control_group.control_bar.sliderReleased.connect(self.update_value)

    def update_value(self):
        cp = [control_group.get_value() for control_group in self.control_groups]
        for i in [0,1,2,3,5,6,7,8]:
            cp[i] *= 1e-6
        cp[4] *= 1e-2 * (cp[1] - cp[0])
        cp[9] *= 1e-2 * (cp[8] - cp[3])
        
        cage = CageParameter(*cp)    
        brg.reload_cage(cage)
        brg.linspace_all()
        
        for _ in range(50):
            brg.smooth_face()
            
        for _ in range(10):
            brg.smooth_ball()
            brg.smooth_inner()
        
        verts = np.array(brg.get_points_as_list())
        self.mesh.points = verts
        surface = self.mesh.extract_geometry()
        verts = np.array(surface.points)
        faces = np.array(surface.faces).reshape((-1, 4))[:, 1:]
        colors = np.array([(i / len(faces), 0, 0, 0.5) for i in range(len(faces))])       
        
        self.view.items[0].setMeshData(vertexes=verts, faces=faces, vertexColors=colors)
               
        total = sum([control_group.get_value() for control_group in self.control_groups])
        self.total_label.setText(f'Torque: {total} [Nmm]')
        
if __name__ == '__main__':
    
    labels = ['内径', '外径', '底面高さ', '肩高さ', '面取り', 'ポケット径', '肩中心高さ', '肩径', '爪高さ', '爪から面までの距離']
    means = [2.35e-3, 2.85e-3, 0.93e-3, 2.10e-3, 0.2, 0.825e-3, 1.70e-3, 1.20e-3, 2.45e-3, 0.2]
    diffs = [0.10e-3, 0.10e-3, 0.10e-3, 0.10e-3, 0.1, 0.001e-3, 0.02e-3, 0.02e-3, 0.10e-3, 0.1]
    mm0 = lambda i: [int((means[i] - diffs[i])*1e6), int((means[i] + diffs[i])*1e6)]
    mm1 = lambda i: [int((means[i] - diffs[i])*1e2), int((means[i] + diffs[i])*1e2)]
    vl0 = lambda v: f'{v * 1e-3:.4g} [mm]'
    vl1 = lambda v: f'{v * 1e-2:.4g} [-]'
    min_maxs     = [mm0, mm0, mm0, mm0, mm1, mm0, mm0, mm0, mm0, mm1]
    value_labels = [vl0, vl0, vl0, vl0, vl1, vl0, vl0, vl0, vl0, vl1]
    
    cage = CageParameter(means[0], means[1], means[2], means[3], means[4], means[5], means[6], means[7], means[8], means[9])
    brg = Brg("data/FEM/FEM2.vtu", "data/FEM/index2.xml", cage)
    filename = 'data/FEM/FEM2.vtu'
    mesh = pv.read(filename)
        
    app = QtWidgets.QApplication(sys.argv)
    control_groups = [ControlGroup(labels[i], min_maxs[i](i), value_labels[i]) for i in range(10)]
    window = Window(control_groups, brg, mesh)
    window.show()
    sys.exit(app.exec_())
    