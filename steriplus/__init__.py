# __init__.py
"""Calculates steric descriptors for molecules.

Modules:
    data: Data dictionaries.
    geometry: Geometry functions and classes.
    helpers: General helper functions.
    io: File parser functions and classes.
    plotting: 3D plotting classes.
    steriplus: Classes for steric descriptor calculations.

Scripts:
    script_buried_volume: Calculate buried volume.
    script_cone_angle: Calculate cone angle.
    script_dispersion: Calculate dispersion descriptor.
    script_sasa: Calculate solvent accessible surface area.
    script_steriplus: Calculate Sterimol parameters.
"""

# Version of the Steriplus package
__version__ = "0.3.0"

from steriplus.steriplus import BuriedVolume, ConeAngle, Dispersion
from steriplus.steriplus import SASA, Sterimol
from steriplus.io import read_gjf, read_xyz
