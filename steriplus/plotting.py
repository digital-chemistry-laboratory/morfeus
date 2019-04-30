"""This module contains functions related to 3D plotting.

Functions:
    ax_3D: Yields axis object for 3D plotting.
    coordinate_axis: Yield axis object for 3D plotting with coordinate axes.
    set_axes_equal: Sets equal perspective for 3D plot.
"""

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from contextlib import contextmanager
from matplotlib.colors import hex2color
import numpy as np
import vpython as vp
from steriplus.data import jmol_colors, atomic_symbols
from steriplus.helpers import convert_element_ids
import math

class MoleculeScene:
    def __init__(self, elements, coordinates, radii, indices=[]):
        elements = convert_element_ids(elements)
        symbols = [atomic_symbols[element] for element in elements]
        hex_colors = [jmol_colors[i] for i in elements]
        colors = [hex2color(color) for color in hex_colors]
    
        if not indices:
            indices = range(1, len(coordinates) + 1)
        
        # Create scene
        scene = vp.canvas(background=vp.color.white)

        # Center camera on geometric center
        center = np.mean(np.array(coordinates), axis=0)
        scene.center = vp.vector(*center)

        spheres = []
        labels_symbols = []
        labels_numbers = []
    
        for coordinate, radius, color, index, symbol in zip(coordinates, radii, colors, indices, symbols):
            pos = vp.vector(*coordinate)
            sphere = vp.sphere(pos=pos, radius=radius, color=vp.vector(*color))
            spheres.append(sphere)
    
            label_symbol = vp.label(pos=pos, text=symbol, opacity=0, box=False, visible=False)
            labels_symbols.append(label_symbol)
            
            label_number = vp.label(pos=pos, text=str(index), opacity=0, box=False, visible=False)
            labels_numbers.append(label_number)
    
        slider_size = vp.slider(pos=scene.title_anchor, bind=self._scale_size, min=0, max=1,value=1)
        scene.append_to_title('    ')
        vp.wtext(pos=scene.title_anchor, text="Size")

        scene.append_to_title('\n')

        slider_sphere_opacity = vp.slider(pos=scene.title_anchor, bind=self._change_sphere_opacity, min=0, max=1,value=1)
        scene.append_to_title('    ')
        vp.wtext(pos=scene.title_anchor, text="Atom opacity")

        checkbox_symbols = vp.checkbox(pos=scene.caption_anchor, bind=lambda x: self._switch_visibility(labels_symbols), text='Display symbols')
        scene.append_to_caption('    ')
        checkbox_numbers = vp.checkbox(pos=scene.caption_anchor, bind=lambda x: self._switch_visibility(labels_numbers), text='Display numbers')
    
        self.scene = scene
        self.spheres = spheres
        self.radii = radii
        self.labels_symbols = labels_symbols
        self.labels_numbers = labels_numbers
        self.slider_size = slider_size
        self.slider_sphere_opacity = slider_sphere_opacity
        self.checkbox_symbols = checkbox_symbols
        self.checkbox_numbers = checkbox_numbers
        self.arrows = []
        self.arrow_labels = []
        self.points = []
        self.triangles = []
    
    def add_arrow(self, start, stop, length, text=""):
        direction = np.array(stop) - np.array(start)
        color = vp.vector(0.12156862745098039, 0.4666666666666667, 0.7058823529411765)
        arrow = vp.arrow(pos=vp.vector(*start), axis=vp.vector(*direction), shaftwidth=0.1, length=length, color=color)
        self.arrows.append(arrow)

        if text:
            arrow_label = vp.label(pos=vp.vector(*stop), text=text, yoffset=10, opacity=0, line=False, box=False, color=vp.color.red)
            self.arrow_labels.append(arrow_label)
    
    def add_cone(self, start, normal, angle, length):
        r = math.tan(angle) * length
        axis = vp.vector(*(-normal))
        pos = vp.vector(*(start + normal * length))
        color = vp.vector(0.12156862745098039, 0.4666666666666667, 0.7058823529411765)
        self.cone = vp.cone(pos=pos, axis=axis, length=length, radius=r, color=color, opacity=0.15)
    
    def add_points(self, points, color):
        color = vp.vector(*hex2color(color))
        points = vp.points(pos=[vp.vector(*point) for point in points], color=color, radius=2)
        self.points.append(points)

    def set_scale(self, value):
        self.slider_size.value = value
        self._scale_size(self.spheres, self.radii, value)

    @staticmethod
    def _switch_visibility(objects):
        for obj in objects:
            obj.visible = not obj.visible

    def _change_sphere_opacity(self):
        for sphere in self.spheres:
            sphere.opacity = self.slider_sphere_opacity.value

    def _scale_size(self):
        for sphere, radius in zip(self.spheres, self.radii):
            sphere.radius = radius * self.slider_size.value

# Code taken from stackexchange
def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.25 * np.max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

@contextmanager
def ax_3D():
    """Context manager for providing 3D axes with Matplotlib.

    Yields:
        ax (object): Axes for Matplotlib plotting
    """
    try:
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.set_aspect('equal')
        yield ax
    finally:
        set_axes_equal(ax)

@contextmanager
def coordinate_axes():
    """Context manager for providing 3D axes with xyz Cartesian vectors
    with Matplotlib.

    Yields:
        ax (object): Axes for Matplotlib plotting
    """
    try:
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        x, y, z = np.zeros((3,3))
        u, v, w = np.array([[1,0,0],[0,1,0],[0,0,1]])

        # Draw arrows for the Cartesian axes
        ax.quiver(x,y,z,u,v,w,arrow_length_ratio=0.1)

        # Print out the text labels for origin, x, y and z
        ax.text(0, 0, 0, "O")
        ax.text(1, 0, 0, "X")
        ax.text(0, 1, 0, "Y")
        ax.text(0, 0, 1, "Z")

        yield ax
    finally:
        set_axes_equal(ax)
