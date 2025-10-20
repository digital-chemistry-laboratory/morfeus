"""Bite angle code."""

from __future__ import annotations

from collections.abc import Iterable
import functools
from typing import Any, cast

import numpy as np

from morfeus.io import read_geometry
from morfeus.typing import Array1DFloat, Array2DFloat, ArrayLike1D, ArrayLike2D


class BiteAngle:
    """Performs and stores the results of a bite angle calculation.

    The bite angle can become 'inverted' if it goes beyond 180°. This cannot be detected
    automatically, so the user has to supply a reference vector that should point in the
    'direction' of the ligand. As a convenience, the vector can be constructed
    automatically from the metal atom to the geometric mean of a set of 'ref_atoms'.

    Args:
        coordinates: Coordinates (Å)
        metal_index: Index of metal atom (1-indexed)
        ligand_index_1: Index of first ligand atom (1-indexed)
        ligand_index_2: Index of second ligand atom (1-indexed)
        ref_vector: Reference vector for determining inverted bite angle
        ref_atoms: Reference atoms for determnining inverted bite angle (1-indexed)

    Attributes:
        angle: Bite angle (degrees)
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
        ref_atoms: Iterable[int] | None = None,
        ref_vector: ArrayLike1D | None = None,
    ) -> None:
        if 0 in {
            metal_index,
            ligand_index_1,
            ligand_index_2,
            *(ref_atoms or ()),
        }:
            raise IndexError("Atom indices should not be 0 (1-indexed).")

        # Check keywords
        if ref_atoms is not None and ref_vector is not None:
            raise ValueError("ref_atoms and ref_vector cannot be set at the same time.")

        coordinates: Array2DFloat = np.asarray(coordinates)

        # Construct vectors
        v_1 = coordinates[ligand_index_1 - 1] - coordinates[metal_index - 1]
        v_1_norm = v_1 / np.linalg.norm(v_1)
        v_2 = coordinates[ligand_index_2 - 1] - coordinates[metal_index - 1]
        v_2_norm = v_2 / np.linalg.norm(v_2)

        # Calculate angle between vectors
        angle_rad = np.arctan2(
            np.linalg.norm(np.cross(v_1_norm, v_2_norm)), np.dot(v_1_norm, v_2_norm)
        )
        angle = np.rad2deg(angle_rad)

        # Check if angle should be inverted
        if ref_atoms is not None:
            ref_vector = (
                np.mean(coordinates[[i - 1 for i in ref_atoms]], axis=0)
                - coordinates[metal_index - 1]
            )
        ref_vector = cast(Array1DFloat, ref_vector)
        inverted = False
        if ref_vector is not None:
            ref_vector /= np.linalg.norm(ref_vector)
            v_mean = (v_1 + v_2) / 2
            v_mean /= np.linalg.norm(v_mean)
            if np.dot(v_mean, ref_vector) < 0:
                inverted = True
                angle = 360 - angle

        self.angle = angle
        self.inverted = inverted


def cli(file: str) -> Any:
    """CLI for bite angle.

    Args:
        file: Geometry file

    Returns:
        Partially instantiated class
    """
    _, coordinates = read_geometry(file)
    return functools.partial(BiteAngle, coordinates)
