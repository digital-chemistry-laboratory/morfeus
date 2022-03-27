"""Plotting functions."""

from __future__ import annotations

from collections.abc import Sequence
import typing

import numpy as np

from morfeus.typing import Array1DFloat
from morfeus.utils import Import, requires_dependency

if typing.TYPE_CHECKING:
    import pyvista as pv
    import vtk


@requires_dependency([Import(module="pyvista", alias="pv")], globals())
def get_drawing_arrow(
    start: Sequence[float] | None = None,
    direction: Sequence[float] | None = None,
    length: float = 1,
    shaft_radius: float = 0.05,
    shaft_resolution: int = 20,
    tip_length: float = 0.25,
    tip_radius: float = 0.1,
    tip_resolution: int = 20,
) -> "pv.MultiBlock":
    """Creates PyVista 3D arrow from cone and cylider.

    Args:
        start: Starting point (Å)
        direction: Direction vector (Å)
        length: Length (Å)
        shaft_radius: Shaft radius (Å)
        shaft_resolution: Shaft resolution
        tip_length: Tip length (Å)
        tip_radius:: Tip radius (Å)
        tip_resolution: Tip resoluation

    Returns:
        arrow: 3D arrow
    """
    # Set start and direction
    if start is None:
        start = [0, 0, 0]
    if direction is None:
        direction = [1, 0, 0]
    start: Array1DFloat = np.array(start)
    direction = np.array(direction) / np.linalg.norm(direction)

    # Create cylinder
    cylinder_length = length - tip_length
    cylinder_center = start + length * direction / 2
    cyl = pv.Cylinder(
        center=cylinder_center,
        direction=direction,
        radius=shaft_radius,
        height=cylinder_length,
        resolution=shaft_resolution,
    )

    # Create cone
    cone_center = start + (cylinder_length + tip_length / 2) * direction
    cone = get_drawing_cone(
        center=cone_center,
        direction=direction,
        radius=tip_radius,
        height=tip_length,
        resolution=tip_resolution,
    )

    # Assemble to arrow
    arrow = pv.MultiBlock()
    arrow.append(cyl)
    arrow.append(cone)

    return arrow


@requires_dependency([Import(module="pyvista", alias="pv"), Import("vtk")], globals())
def get_drawing_cone(
    center: Sequence[float] | None = None,
    direction: Sequence[float] | None = None,
    height: float = 1.0,
    radius: float | None = None,
    capping: bool = True,
    angle: float | None = None,
    resolution: int = 6,
) -> "pv.PolyData":
    """Create a cone.

    Copy from the PyVista code.

    Args:
        center: Center in [x, y, z]. middle of the axis of the cone.
        direction: Direction vector in [x, y, z]. orientation vector of the cone.
        height: Height along the cone in its specified direction.
        radius: Base radius of the cone
        capping: Turn on/off whether to cap the base of the cone with a polygon.
        angle: The angle degrees between the axis of the cone and a generatrix.
        resolution: Number of facets used to represent the cone

    Returns:
        cone: 3D cone

    Raises:
        Exception: When both radius and angle are specified.
    """
    if center is None:
        center = [0.0, 0.0, 0.0]
    if direction is None:
        direction = [1.0, 0.0, 0.0]

    src = vtk.vtkConeSource()
    src.SetCapping(capping)
    src.SetDirection(direction)
    src.SetCenter(center)
    src.SetHeight(height)
    if angle and radius:
        raise Exception("Both radius and angle specified. They are mutually exclusive.")
    elif angle and not radius:
        src.SetAngle(angle)
    elif not angle and radius:
        src.SetRadius(radius)
    elif not angle and not radius:
        src.SetRadius(0.5)
    src.SetResolution(resolution)
    src.Update()
    cone = pv.wrap(src.GetOutput())

    return cone
