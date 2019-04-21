"""This module contains functions related to 3D plotting.

Functions:
    ax_3D: Yields axis object for 3D plotting.
    coordinate_axis: Yield axis object for 3D plotting with coordinate axes.
    set_axes_equal: Sets equal perspective for 3D plot.
"""

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from contextlib import contextmanager
import numpy as np

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
