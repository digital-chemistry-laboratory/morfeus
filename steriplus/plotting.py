import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from contextlib import contextmanager

@contextmanager
def ax_3D():
    """Context manager for providing 3D axes with Matplotlib.

    Yields:
        ax (object): Axes for Matplotlib plotting
    """
    try:
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        yield ax
    finally:
        pass

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
        pass
