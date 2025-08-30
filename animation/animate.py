import sys
import numpy as np
from vispy import app, scene, color
from pathlib import Path
from time import sleep
import numpy as np
from pathlib import Path
import struct

def list_bin_files(folder):
    path = Path(folder)
    return sorted(path.glob("*.bin"))

def load_particles(filename):
    with open(filename, 'rb') as f:
        # Read header
        header_fmt = '<4s I Q Q I'  # little-endian: magic[4], version, step, count, flags
        header_size = struct.calcsize(header_fmt)
        header = f.read(header_size)
        magic, version, step, count, flags = struct.unpack(header_fmt, header)

        if magic != b'PART':
            raise ValueError("Invalid magic number. Not a particle file.")
        
        # Define particle struct: x, y, z (float32), type (uint32), id (uint32), state (uint8)
        particle_dtype = np.dtype([
            ('x', 'f4'),
            ('y', 'f4'),
            ('z', 'f4'),
            ('type', 'u4'),
            ('id', 'u4'),
            ('state', 'u1')
        ])


        # Read the particle data
        particles = np.fromfile(f, dtype=particle_dtype, count=count)

        return particles, step

def add_orientation_axes(parent):
    axes = []

    # Length of the arrows
    L = 1.0

    # X-axis (red)
    axis_x = scene.visuals.Arrow(pos=np.array([[0, 0, 0], [L, 0, 0]]),
                                 color='red', width=3, arrow_size=10, parent=parent)
    axes.append(axis_x)

    # Y-axis (green)
    axis_y = scene.visuals.Arrow(pos=np.array([[0, 0, 0], [0, L, 0]]),
                                 color='green', width=3, arrow_size=10, parent=parent)
    axes.append(axis_y)

    # Z-axis (blue)
    axis_z = scene.visuals.Arrow(pos=np.array([[0, 0, 0], [0, 0, L]]),
                                 color='blue', width=3, arrow_size=10, parent=parent)
    axes.append(axis_z)

    return axes

class ParticleVisualizer(scene.SceneCanvas):
    def __init__(self, folder, delay=0.05):
        super().__init__(keys='interactive', size=(800, 600), title='Particle Viewer', show=True)

        # Set up a view and camera
        self.unfreeze()  # Allow dynamic attributes
        self.view = self.central_widget.add_view()
        self.view.camera = 'turntable'

        # self.axes = add_orientation_axes(self.view.scene)        

        # Set background and create a marker visual
        self.bgcolor = '#000000'
        self.scatter = scene.visuals.Markers(parent=self.view.scene)

        # Load data
        self.files = list_bin_files(folder)
        self.current = 0
        self.delay = delay

        # Start timer
        self.timer = app.Timer(interval=delay, connect=self.update_frame, start=True)

    def update_frame(self, event):
        if self.current >= len(self.files):
            self.current = 0  # Loop

        file = self.files[self.current]
        particles, step = load_particles(file)

        positions = np.vstack((particles['x'], particles['y'], particles['z'])).T
        colors = self.color_by_type(particles['type'])

        self.scatter.set_data(positions, face_color=colors, size=10)
        self.current += 1

    def color_by_type(self, types):
        cmap = color.get_colormap('viridis')
        ptp = np.ptp(types)
        norm = (types - types.min()) / (ptp if ptp > 0 else 1)
        return cmap.map(norm)




if __name__ == "__main__": 
    vis = ParticleVisualizer("../output/halleys_comet")
    app.run()

    # print(load_particles(files[0]))

   