import numpy as np
from scipy.spatial.transform import Rotation

def rotate_coordinates(coordinates, vector, axis):
    """Rotates coordinates with the rotation that aligns vector with axis.

    Args:
        axis (ndarray)                  :   Vector to align vector to
        coordinates (ndarray)           :   Coordinates to rotate
        vector (ndarray)                :   Vector to align

    Returns:
        rotated_coordinates (ndarray)   :   Rotated coordinates.
    """
    # Get rotation quaternion that overlays vector with x-asis
    real = np.dot(vector, axis).reshape(1) + 1
    #  Handle case of antiparallel vectors
    if real < 1e-6:
        w = np.cross(vector, np.array([0, 0, 1]))
        if np.linalg.norm(w) < 1e-6:
            w = np.cross(vector, np.array([1, 0, 0]))
    else:
        w = np.cross(vector, axis)

    q = np.concatenate((w, real))
    q = q / np.linalg.norm(q)

    # Rotate atomic coordinates
    rot = Rotation.from_quat(q)
    rotated_coordinates = rot.apply(coordinates)
    return rotated_coordinates
