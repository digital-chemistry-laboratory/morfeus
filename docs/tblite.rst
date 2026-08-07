======
TBLite
======

The :py:class:`TBLite <morfeus.tblite.TBLite>` class provides a native Python
interface to the tblite__ library for xTB single-point properties.

This is an optional backend for the subset of quantities that tblite exposes
directly. It is not a drop-in replacement for :py:class:`XTB <morfeus.xtb.XTB>`,
which still uses the ``xtb`` command line program for descriptors that depend on
specialized workflows such as Fukui functions, IP/EA corrections, fractional
occupation density, CM5 charges, or detailed solvation terms.

******
Module
******

The native backend can calculate total energies, Mulliken charges, bond orders,
molecular dipoles, and frontier orbital energies when those quantities are
available from the tblite result container.

.. code-block:: python
  :caption: Example

  >>> from morfeus import TBLite, read_xyz
  >>> elements, coordinates = read_xyz("ammonia.xyz")
  >>> tblite = TBLite(elements, coordinates)
  >>> tblite.get_energy()
  -4.425...
  >>> tblite.get_charges()
  {1: -0.434..., 2: 0.144..., 3: 0.144..., 4: 0.144...}
  >>> tblite.get_bond_order(1, 2)
  0.97...

The method can be selected with ``method=1``, ``method=2``, ``method="GFN1-xTB"``,
``method="GFN2-xTB"``, or ``method="IPEA1-xTB"``. Molecular charge and unpaired
electrons are set with ``charge=<int>`` and ``n_unpaired=<int>``.

Implicit solvation can be requested with ``solvent=<str>`` and
``solvation_model="alpb"`` or ``solvation_model="gbsa"``. Only the total
single-point properties exposed by tblite are available; the solvation-specific
descriptors from :py:class:`XTB <morfeus.xtb.XTB>` remain tied to the ``xtb``
command line backend.

.. __: https://tblite.readthedocs.io
