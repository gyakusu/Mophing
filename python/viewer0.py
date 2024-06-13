import sys
from PyQt5 import QtCore, QtWidgets
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
import pyvista as pv

class Window(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()

        filename = 'data/FEM/FEM2.vtu'
        FEM2 = pv.read(filename)

        # Create a control bar
        self.control_bar = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
        self.control_bar.setMinimum(0)
        self.control_bar.setMaximum(100)
        self.control_bar.setValue(50)
        self.control_bar.sliderReleased.connect(self.update_color)  # Change here
        self.control_bar.setFixedSize(200, 30)
        
        # Create a label to display the description
        self.label = QtWidgets.QLabel(self)
        self.label.setText("保持器高さ")
        self.label.setFixedSize(100, 30)
        
        # Create a 3D graphics widget
        self.view = gl.GLViewWidget(self)
        self.view.setFixedSize(800, 800)

        # Create a label to display the value indicator
        self.value_label = QtWidgets.QLabel(self)
        self.value_label.setFixedSize(50, 30)
        self.value_label.setText(f'{50} [mm]')
        
        # Create a splitter to hold the control bar and the view
        splitter = QtWidgets.QSplitter(self)
        splitter.addWidget(self.view)
        splitter.addWidget(self.label)
        splitter.addWidget(self.control_bar)
        splitter.addWidget(self.value_label)
        
        # Set the splitter as the central widget
        self.setCentralWidget(splitter)        

        # Set the view to show the entire scene
        self.view.setCameraPosition(distance=1e-2)
        
        # Set the background color to white
        self.view.setBackgroundColor(pg.mkColor(255, 255, 255))

        # 頂点の座標
        surface = FEM2.extract_geometry()
        self.verts = np.array(surface.points)
        self.faces = np.array(surface.faces).reshape((-1, 4))[:, 1:]

        # Calculate the color for each face based on its index
        self.num_faces = len(self.faces)
        self.colors = np.zeros((self.num_faces, 4))
        for i in range(self.num_faces):
            color = (i / self.num_faces, 0, 0, 0.5)
            self.colors[i] = color

        # メッシュアイテムの作成
        mesh = gl.GLMeshItem(vertexes=self.verts, faces=self.faces, vertexColors=self.colors)        
        self.view.addItem(mesh)

    def update_color(self):
        value = self.control_bar.value()
        
        # Convert the value to a color
        for i in range(self.num_faces):
            color = (i / self.num_faces, 0, value / 100, 0.5)
            self.colors[i] = color
        
        self.view.items[0].setMeshData(vertexes=self.verts, faces=self.faces, vertexColors=self.colors)
        
        # Update the value label when the value changes
        self.value_label.setText(f'{value} [mm]')
        
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = Window()
    window.show()
    sys.exit(app.exec_())