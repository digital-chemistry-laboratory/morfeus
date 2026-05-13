"""Typing code for arrays.

Borrowed from gpaw/typing.py
"""

from __future__ import annotations

from typing import Any

from numpy import integer, signedinteger
from numpy.typing import ArrayLike, NDArray

IntLike = int | integer | signedinteger

ArrayLike1D = ArrayLike
ArrayLike2D = ArrayLike
ArrayLike3D = ArrayLike

ArrayND = NDArray
Array1DAny = Any
# TODO: Remove these when https://github.com/python/mypy/issues/11347 is resolved
Array1DFloat = Any  # NDArray[np.float_]
Array2DFloat = Any  # NDArray[np.float_]
Array3DFloat = Any  # NDArray[np.float_]
Array1DInt = Any  # NDArray[np.int_]
Array2DInt = Any  # NDArray[np.int_]
Array3DInt = Any  # NDArray[np.int_]
Array1DBool = Any  # NDArray[np.bool_]
Array1DStr = Any  # NDArray[np.str_]
