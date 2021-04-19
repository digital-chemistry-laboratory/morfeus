"""Typing code for arrays.

Borrowed from gpaw/typing.py
"""

from typing import Any

import numpy as np

try:
    from numpy.typing import ArrayLike
except ImportError:
    ArrayLike = Any

ArrayLike1D = ArrayLike
ArrayLike2D = ArrayLike
ArrayLike3D = ArrayLike

ArrayND = np.ndarray
Array1D = ArrayND
Array2D = ArrayND
Array3D = ArrayND
