"""Bite angle code."""

from __future__ import annotations

from typing import Optional

import numpy as np

from morfeus.typing import ArrayLike1D, ArrayLike2D


class BiteAngle:
    """Performs and stores the results of a bite angle calculation.

    The bite angle can become 'inverted' if it goes beyond 180°. This cannot be detected
    automatically, so the user has to supply a reference vector that should point in the
    direction of the ligand 'center of mass'.

    Args:
        coordinates: Coordinates (Å)
        metal_index: Index of metal atom
        ligand_index_1: Index of first ligand atom
        ligand_index_2: Index of second ligand atom

    Attributes:
        angle: Bite angle
        inverted: Whether angle is 'inverted'
    """

    angle: float
    inverted: bool

    def __init__(
        self,
        coordinates: ArrayLike2D,
        metal_index: int,
        ligand_index_1: int,
        ligand_index_2: int,
        ref_vector: Optional[ArrayLike1D] = None,
    ) -> None:
        # Construct vectors
        v_1 = coordinates[ligand_index_1 - 1] - coordinates[metal_index - 1]
        v_1 /= np.linalg.norm(v_1)
        v_2 = coordinates[ligand_index_2 - 1] - coordinates[metal_index - 1]
        v_2 /= np.linalg.norm(v_2)

        # Calculate angle between vectors
        angle_rad = np.arctan2(np.linalg.norm(np.cross(v_1, v_2)), np.dot(v_1, v_2))
        angle = np.rad2deg(angle_rad)

        # Check if angle should be inverted
        inverted = False
        if ref_vector is not None:
            ref_vector_norm = ref_vector / np.linalg.norm(ref_vector)
            if np.dot(v_1 + v_2, ref_vector_norm) < 0:
                inverted = True
                angle = 360 - angle

        self.angle = angle
        self.inverted = inverted
