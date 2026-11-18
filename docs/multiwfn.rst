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

    Tests cover Molden files from xTB, PySCF and ORCA, with additional WFN and
    WFX coverage.
    The `molden2aim`_ utility is recommended to standardize wavefunction files.

******
Module
******

The :py:class:`Multiwfn <morfeus.multiwfn.Multiwfn>` class allows calculation
of a variety of descriptors and properties from wavefunction (``.molden``/``.wfn``/``.wfx``)
files.
Additionally, grid files (``.cub`` or ``.grd``) can be generated or integrated.

Several descriptors are spin-dependent and are handled by
``has_spin=<bool|None>``, inferred from the loaded wavefunction if not provided.

Citations used by the selected analyses can be collected directly:

.. code-block:: python
  :caption: Example

  >>> from morfeus import Multiwfn
  >>> mwfn = Multiwfn("molden.input")
  >>> mwfn.get_charges("hirshfeld")
  >>> mwfn.get_citations()

All available options can be accessed with the ``list_options()`` method:

.. code-block:: python
  :caption: Example

  >>> options = mwfn.list_options()
  >>> options.keys()
  ['bond_order', 'charges', 'descriptors', 'descriptors_fast', 'grid_quality', 'interaction', 'surface']


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
wavefunction types or spin states. Atomic integration supports Becke, Hirshfeld,
and MBIS partitioning through ``get_descriptor(..., partitioning=...)``.
Hirshfeld-I requires a separate atomic-reference preparation workflow and is not
supported by this interface. Localization and fuzzy/Laplacian bond orders use
Becke partitioning independently of previous calculations.

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
  (3.28428, 0.29337)


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

  >>> with Multiwfn("molden.input", persistent=True) as mwfn:
  ...     cube_path = mwfn.get_grid("rho", "low", output_path="results/rho.cub")
  ...     integrated = mwfn.grid_to_descriptors(cube_path)
  >>> cube_path.name
  'rho.cub'

Exported cube files are retained at the requested destination. Intermediate
calculation files are temporary by default.

For quantitative total-density analysis from xTB, prefer the Molden wavefunction
and ``Multiwfn.get_descriptor("rho")`` or ``get_grid("rho", ...)``. 
Check electron and spin integrals before using exported cubes quantitatively.


#########################
NCI and IGMH interactions
#########################

Noncovalent interaction (NCI) analysis :footcite:`johnson_revealing_2010` and the
independent gradient model based on Hirshfeld partition (IGMH)
:footcite:`lu_chen_igmh_2022` characterize interactions from the electron density.
``get_interaction()`` exports these fields as cube files. IGMH also supports
fragment analysis and atomic or atom-pair contributions.

.. code-block:: python
  :caption: Example

  >>> with Multiwfn("complex.wfx") as mwfn:
  ...     nci = mwfn.get_interaction("nci", output_dir="results/nci")
  ...     igmh = mwfn.get_interaction(
  ...         "igmh", output_dir="results/igmh", fragments=[[1, 2, 3], [4, 5, 6]],
  ...     )
  >>> inter_cube = igmh["grids"]["delta_g_inter"]

These descriptors describe interaction patterns, not interaction energies.
For details, see :py:meth:`get_interaction <morfeus.multiwfn.Multiwfn.get_interaction>`.


#################################
Orbital records and frontier gaps
#################################

Orbital records and HOMO–LUMO gaps are available in eV. Gaps are calculated
within each spin channel:

.. code-block:: python
  :caption: Example

  >>> with Multiwfn("molden.input") as mwfn:
  ...     orbitals = mwfn.get_orbitals()
  ...     frontiers = mwfn.get_frontier_gaps()
  >>> gap = frontiers["restricted"]["gap_ev"]

The signed gap, :math:`\varepsilon_\mathrm{LUMO} - \varepsilon_\mathrm{HOMO}`,
can be negative for non-Aufbau occupations. These are orbital-energy differences,
not excitation energies.


####################
Persistent execution
####################

By default, each calculation starts a fresh process. To opt into process reuse:

.. code-block:: python
  :caption: Example

  >>> with Multiwfn("molden.input", persistent=True) as mwfn:
  ...     charges = mwfn.get_charges("hirshfeld")
  ...     moments = mwfn.get_electric_moments()

The loaded wavefunction is reused until the context exits or ``close()`` is called.


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
