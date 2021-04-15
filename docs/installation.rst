============
Installation
============

ᴍᴏʀғᴇᴜs can be installed using pip and works with Python >= 3.7.

************
Dependencies
************
Core dependencies:

* numpy_
* scipy_

Optional depedencies:

* dftd4_
* matplotlib_
* openbabel_
* pymeshfix_
* pyvista_
* pyvistaqt_
* qcengine_
* rdkit_
* spyrmsd_
* vtk_
* xtb-python_

The optional dependencies are used for 3D visualization and in some of the
Dispersion descriptor calculations. *xtb* is used for calculation of electronic
properties.

***
pip
***

Clone the repository from GitHub. Install core functionality with pip:

.. code-block:: console

  pip install .

Install with extra graphics functionality:

.. code-block:: console

  pip install .[extras]

.. _dftd4: https://github.com/dftd4/dftd4
.. _matplotlib: https://matplotlib.org
.. _numpy: https://numpy.org
.. _openbabel: http://openbabel.org/
.. _pymeshfix: https://github.com/pyvista/pymeshfix
.. _pyvista: https://github.com/pyvista/pyvista
.. _pyvistaqt: https://github.com/pyvista/pyvistaqt
.. _qcengine: https://github.com/MolSSI/QCEngine
.. _rdkit: https://www.rdkit.org
.. _scipy: https://github.com/pyvista/pyvistaqt
.. _spyrmsd: https://github.com/RMeli/spyrmsd
.. _vtk:  https://vtk.org
.. _xtb-python: https://github.com/grimme-lab/xtb-python
