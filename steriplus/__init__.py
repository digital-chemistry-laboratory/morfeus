# __init__.py
"""Calculates steric descriptors for molecules.

Modules:
    data: Data dictionaries.
    geometry: Helper functions related to geometry.
    io: Geometry file parsers.
    plotting: 3D plotting functions.
    steriplus: Classes for steric descriptor calculations.

Scripts:
    script_buried_volume: Calculate buried volume.
    script_cone_angle: Calculate cone angle.
    script_sasa: Calculate solvent accessible surface area.
    script_steriplus: Calculate Sterimol parameters.
"""

# Version of the Steriplus package
__version__ = "0.3.0"

from steriplus.steriplus import Sterimol, BuriedVolume, SASA, ConeAngle
from steriplus.io import read_gjf, read_xyz
