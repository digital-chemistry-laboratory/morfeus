"""Support for 3D plotting of steric descriptors.

Classes:
    MoleculeScene: Draw molecules and visualize steric descriptors.
"""
import math

from matplotlib.colors import hex2color
import numpy as np
import vpython as vp

from steriplus.data import atomic_symbols, jmol_colors
from steriplus.helpers import convert_elements

class MoleculeScene:
    """Draw molecules and visualize steric descriptors using vpython.

    Args:
        coordinates (list): Coordinates (Å)
        elements (list): Elements as atomic symbols or numbers
        indices (list): Atomic indices (starting from 1)
        radii (list): Sphere radii (Å)
    
    Attributes:
        arrow_labels (list): Arrow labels
        arrows (list): Arrows objects
        checkbox_numbers (obj): Checkbox controlling atomic number visibility. 
        checkbox_symbols (obj): Checkbox controlling atomic symbol visilibity.
        labels_numbers (list): Atomic number labels
        labels_symbols (list): Atomic symbol labels 
        points (list): Points objects
        radii (list): Sphere radii (Å)
        scene (obj): Main scene for display
        slider_opacity (obj): Slider controlling sphere opacity
        slider_size (obj): Slider controlling sphere size
        spheres (list): Sphere objects
    """
    def __init__(self, elements, coordinates, radii, indices=[]):
        # Set up atomic numbers, symbols and colors
        elements = convert_elements(elements)
        symbols = [atomic_symbols[element] for element in elements]
        hex_colors = [jmol_colors[i] for i in elements]
        colors = [hex2color(color) for color in hex_colors]

        # Create indices if these are not supplied
        if not indices:
            indices = range(1, len(coordinates) + 1)
        
        # Create scene
        scene = vp.canvas(background=vp.color.white)

        # Center camera on geometric center
        center = np.mean(np.array(coordinates), axis=0)
        scene.center = vp.vector(*center)

        # Draw atomic spheres and set up labels and numbers
        spheres = []
        labels_symbols = []
        labels_numbers = []
    
        for coordinate, radius, color, index, symbol in \
                zip(coordinates, radii, colors, indices, symbols):
            pos = vp.vector(*coordinate)
            sphere = vp.sphere(pos=pos, radius=radius, color=vp.vector(*color))
            spheres.append(sphere)
    
            label_symbol = vp.label(pos=pos, text=symbol, opacity=0, box=False,
                                    visible=False)
            labels_symbols.append(label_symbol)
            
            label_number = vp.label(pos=pos, text=str(index), opacity=0, 
                                    box=False, visible=False)
            labels_numbers.append(label_number)

        # Draw sliders and checkboxes
        slider_size = vp.slider(pos=scene.title_anchor, bind=self._scale_size,
                                min=0, max=1,value=1)
        scene.append_to_title('    ')
        vp.wtext(pos=scene.title_anchor, text="Size")

        scene.append_to_title('\n')

        slider_opacity = vp.slider(pos=scene.title_anchor,
                                   bind=self._change_opacity,
                                   min=0, max=1,value=1)
        scene.append_to_title('    ')
        vp.wtext(pos=scene.title_anchor, text="Atom opacity")

        checkbox_symbols = vp.checkbox(
            pos=scene.caption_anchor,
            bind=lambda x: self._switch_visibility(labels_symbols),
            text='Display symbols')
        
        scene.append_to_caption('    ')

        checkbox_numbers = vp.checkbox(
            pos=scene.caption_anchor,
            bind=lambda x: self._switch_visibility(labels_numbers),
            text='Display numbers')
    
        # Set up attributes
        self.arrow_labels = []
        self.arrows = []
        self.checkbox_numbers = checkbox_numbers
        self.checkbox_symbols = checkbox_symbols
        self.labels_numbers = labels_numbers
        self.labels_symbols = labels_symbols
        self.points = []
        self.radii = radii
        self.scene = scene
        self.slider_size = slider_size
        self.slider_opacity = slider_opacity
        self.spheres = spheres
    
    def add_arrow(self, start, stop, length, text=""):
        """Adds an arrow with optional text.

        Args:
            length (float): Length of arrow (Å).
            start (list): Coordinates for start of arrow (Å)
            stop (list): Coordinates for end of arrow (Å)
            text (str): Text to display at end of arrow
        """
        # Draw arrow
        direction = np.array(stop) - np.array(start)
        color = vp.vector(0.12156862745098039, 0.4666666666666667, 0.7058823529411765)
        arrow = vp.arrow(pos=vp.vector(*start), axis=vp.vector(*direction), shaftwidth=0.1, length=length, color=color)
        self.arrows.append(arrow)

        # Draw text if given
        if text:
            arrow_label = vp.label(pos=vp.vector(*stop), text=text, yoffset=10, opacity=0, line=False, box=False, color=vp.color.red)
            self.arrow_labels.append(arrow_label)
    
    def add_cone(self, start, normal, angle, length):
        """Adds cone with apex at starting point.

        Args:
            angle (float): Cone angle (deg)
            length (float): Cone length (Å)
            normal (list): Direction vector of cone (Å)
            start (list): Starting coordinates for cone (Å)
        """
        r = math.tan(angle) * length
        axis = vp.vector(*(-normal))
        pos = vp.vector(*(start + normal * length))
        color = vp.vector(0.12156862745098039, 0.4666666666666667,
                          0.7058823529411765)
        self.cone = vp.cone(pos=pos, axis=axis, length=length,
                            radius=r, color=color, opacity=0.15)
    
    def add_points(self, points, color):
        """Add points of certain color.

        Args:
            color (str): Color in hexademical code
            points (list): Coordinates of points (Å)
        """
        color = vp.vector(*hex2color(color))
        points = vp.points(pos=[vp.vector(*point) for point in points], color=color, radius=2)
        self.points.append(points)

    def set_scale(self, value):
        """Set the scale of the atomic spheres.
        
        Args:
            value (float): Atom size factor (between 0 and 1) 
        """
        self.slider_size.value = value
        self._scale_size()
    
    def _change_opacity(self):
        """Change the sphere opacity according to the slider."""
        for sphere in self.spheres:
            sphere.opacity = self.slider_opacity.value

    def _scale_size(self):
        """Scale the sphere sizes according to the slider."""
        for sphere, radius in zip(self.spheres, self.radii):
            sphere.radius = radius * self.slider_size.value

    @staticmethod
    def _switch_visibility(objects):
        """Switch visibility of supplied objects. Helper function for slider.

        Args:
            objects (list): Objects to switch
        """
        for obj in objects:
            obj.visible = not obj.visible