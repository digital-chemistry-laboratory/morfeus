========
Multiwfn
========

Morfeus provides an interface to the `Multiwfn program`_ for wavefunction analysis as described by Tian Lu :footcite:`lu_2024` :footcite:`lu_chen_2012`.

Multiwfn itself must be installed separately and be available as ``Multiwfn`` in
your shell environment. Refer to the `Multiwfn manual`_ for installation instructions.

The package `Pexpect`_ is also required to run any function in this module.

Morfeus controls Multiwfn's interactive menu system with scripted command
sequences. Results are parsed into Python dictionaries indexed by atom number
(1-based indexing) or descriptor names.

For Multiwfn menu-level details, see the `Multiwfn manual`_.

.. warning::

    This module has only been tested with molden files derived from xTB, PySCF and ORCA.
    The `molden2aim`_ utility is recommended to standardize wavefunction files.

.. __:

******
Module
******

The :py:class:`Multiwfn <morfeus.multiwfn.Multiwfn>` class allows calculation
of a variety of descriptors and properties from wavefunction (``.molden``/``.wfn``)
files.
Additionally, grid files (``.cub`` or ``.grd``) can be generated or integrated.

Several descriptors are spin-dependent and are handled by
``has_spin=<bool|None>``, detected automatically if not provided.

All available options can be easily retrieved with the ``list_options()`` method.

Citations used by the selected analyses can be collected directly:

.. code-block:: python
  :caption: Example

  >>> from morfeus import Multiwfn
  >>> mwfn = Multiwfn("molden.input")
  >>> mwfn.get_charges("hirshfeld")
  >>> mwfn.get_citations()

All available options can be accessesed with the ``list_options()`` method:

.. code-block:: python
  :caption: Example

  >>> options = mwfn.list_options()
  >>> options.keys()
  ['bond_order', 'charges', 'descriptors', 'descriptors_fast', 'grid_quality', 'surface']


#######################
Charges and bond orders
#######################

Several charge and bond order models are available and can be calculated as:

.. code-block:: python
  :caption: Example
  
  >>> mwfn = Multiwfn("molden.input")
  >>> charges = mwfn.get_charges(model="adch")
  >>> bond_orders = mwfn.get_bond_order(model="mayer")
  >>> charges[1]
  0.123
  >>> bond_orders[(1, 2)]
  0.456


#######################
Fuzzy-space descriptors
#######################

Atomic descriptors can be calculated for a variety of real-space functions.
The function names can be listed with
``mwfn.list_options()["descriptors"]``. A sublist excluding costly functions is
given as ``mwfn.list_options()["descriptors_fast"]``.

The descriptors can be calculated with:

.. code-block:: python
  :caption: Example

  >>> mwfn = Multiwfn("molden.input")
  >>> rho = mwfn.get_descriptor("rho")

Caution: Some functions may not be available or meaningful for some
wavefunction types or spin states.

###################
Surface descriptors
###################

Several global and atomic descriptors can be derived from the quantitative
analysis of molecular surfaces. The available surfaces can be listed with
``mwfn.list_options()["surface"]``. 
The atomic descriptors are only provided for atoms with significant contributions.
For example, the electrostatic potential (ESP) can be calculated with:

.. code-block:: python
  :caption: Example
  
  >>> mwfn = Multiwfn("molden.input")
  >>> esp_surface = mwfn.get_surface("esp")
  >>> esp_surface.keys()
  dict_keys(['atomic', 'global'])
  >>> esp_surface["global"].get("minimal_value")
  -58.90822
  >>> atom_3 = esp_surface["atomic"][3]
  >>> atom_3.get("area_positive"), atom_3.get("area_negative")
  (8.82908, -2.14958)


#######################
Conceptual DFT analyses
#######################

These analyses require a closed-shell wavefunction, so set ``has_spin=False``
or let spin be auto-detected for a closed-shell input.

.. code-block:: python
  :caption: Example

  >>> mwfn = Multiwfn("molden.input", run_path="mwfn_run", has_spin=False)
  >>> fukui = mwfn.get_fukui()
  >>> superdeloc = mwfn.get_superdelocalizabilities()


##########
Grid files
##########

Grid files can be generated for all descriptors with ``get_grid()``
(the grid quality significantly affects the computational cost and accuracy).
Cube files can directly be integrated per atom with ``grid_to_descriptors()``.

.. code-block:: python
  :caption: Example

  >>> mwfn = Multiwfn("molden.input", run_path="mwfn_run")
  >>> cube_path = mwfn.get_grid("rho", "low", grid_file_name="rho.cub")
  >>> cube_path.name
  'rho.cub'
  >>> integrated = mwfn.grid_to_descriptors(cube_path)


**********
Background
**********

#######################
Electrostatic potential
#######################

One commonly used descriptor is the electrostatic potential (ESP).
The molecular ESP at position :math:`\mathbf{r}` is defined as

.. math::

  V(\mathbf{r}) = \sum_A \frac{Z_A}{|\mathbf{R}_A - \mathbf{r}|}
  - \int \frac{\rho(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|} d\mathbf{r}'

where :math:`Z_A` and :math:`\mathbf{R}_A` are nuclear charges/positions and
:math:`\rho(\mathbf{r}')` is the electron density. 
Mapping :math:`V(\mathbf{r})` onto a molecular surface is commonly used to 
identify electron-rich (:math:`V < 0`) and electron-poor (:math:`V > 0`) regions.
:footcite:`murray_politzer_2011`

In this module, ``get_surface("esp")`` provides global and atom-resolved
surface statistics, while ``get_descriptor("esp_total")`` and related
ESP-derived functions provide fuzzy-space atomic integrals.

.. footbibliography::

.. _Multiwfn program: http://sobereva.com/multiwfn/
.. _Multiwfn manual: http://sobereva.com/multiwfn/Multiwfn_manual.html
.. _Pexpect: https://pexpect.readthedocs.io/en/stable/
.. _molden2aim: https://github.com/zorkzou/Molden2AIM
