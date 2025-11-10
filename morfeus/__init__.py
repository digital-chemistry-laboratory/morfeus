"""Calculates steric descriptors for molecules.

Modules:
    bite_angle: Bite angle code
    buried_volume: Buried volume code
    calculators: Internal calculators
    cone_angle: Cone angle code
    conformer: Conformer tools
    d3_data: Reference data for D3 calculations
    data: Module related to data
    geometry: Help classes and functions related to geometry
    io: Input and output
    local_force: Local force constant code
    plotting: Plotting functions
    pyramidalization: Pyramidalization code
    qc: Interface to quantum-chemical programs
    sasa: Solvent accessible surface area code
    solid_angle: Solid angle code
    sterimol: Sterimol code
    typing: Typing code for arrays
    utils: Helper functions.
    visible_volume: Visible volume code
    xtb: xtb interface code

Scripts:
    morfeus: Command line interface
"""

from importlib import metadata

from morfeus.bite_angle import BiteAngle
from morfeus.buried_volume import BuriedVolume
from morfeus.cone_angle import ConeAngle
from morfeus.dispersion import Dispersion
from morfeus.io import read_geometry, read_gjf, read_xyz, write_xyz
from morfeus.local_force import LocalForce
from morfeus.pyramidalization import Pyramidalization
from morfeus.sasa import SASA
from morfeus.solid_angle import SolidAngle
from morfeus.sterimol import Sterimol
from morfeus.visible_volume import VisibleVolume
from morfeus.xtb import XTB

__all__ = [
    "read_geometry",
    "read_gjf",
    "read_xyz",
    "write_xyz",
    "BiteAngle",
    "BuriedVolume",
    "ConeAngle",
    "Dispersion",
    "LocalForce",
    "Pyramidalization",
    "SASA",
    "SolidAngle",
    "Sterimol",
    "VisibleVolume",
    "XTB",
]

# Version of the morfeus package
__version__ = metadata.version("morfeus-ml")
