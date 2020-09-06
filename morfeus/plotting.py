"""Plotting functions."""

import vtk
import pyvista as pv
import numpy as np

def Arrow_3D(start=[0, 0, 0], direction=[1, 0, 0], length=1, shaft_radius=0.05, shaft_resolution=20, tip_length=0.25, tip_radius=0.1, tip_resolution=20):
    start = np.array(start)
    direction = np.array(direction) / np.linalg.norm(direction)
    cylinder_length = length - tip_length
    cylinder_center = start + length * direction / 2
    cyl = pv.Cylinder(center=cylinder_center, direction=direction, radius=shaft_radius, height=cylinder_length, resolution=shaft_resolution)
    
    cone_center = start + (cylinder_length + tip_length / 2) * direction
    cone = Cone_3D(center=cone_center, direction=direction, radius=tip_radius, height=tip_length, resolution=tip_resolution)
    arrow = pv.MultiBlock()
    arrow.append(cyl)
    arrow.append(cone)
    
    return arrow

def Cone_3D(center=(0., 0., 0.), direction=(1., 0., 0.), height=1.0, radius=None,
         capping=True, angle=None, resolution=6):
    """Create a cone
    Parameters
    ----------
    center : np.ndarray or list
        Center in [x, y, z]. middle of the axis of the cone.
    direction : np.ndarray or list
        direction vector in [x, y, z]. orientation vector of the cone.
    height : float
        height along the cone in its specified direction.
    radius : float
        base radius of the cone
    capping : bool
        Turn on/off whether to cap the base of the cone with a polygon.
    angle : float
        The angle degrees between the axis of the cone and a generatrix.
    resolution : int
        number of facets used to represent the cone
    """
    src = vtk.vtkConeSource()
    src.SetCapping(capping)
    src.SetDirection(direction)
    src.SetCenter(center)
    src.SetHeight(height)
    if angle and radius:
        raise Exception ("Both radius and angle specified. They are mutually exclusive.")
    elif angle and not radius:
        src.SetAngle(angle)
    elif not angle and radius:
        src.SetRadius(radius)
    elif not angle and not radius:
        src.SetRadius(0.5)
    src.SetResolution(resolution)
    src.Update()
    return pv.wrap(src.GetOutput())