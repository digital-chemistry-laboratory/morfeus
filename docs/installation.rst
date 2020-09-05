============
Installation
============

ᴍᴏʀғᴇᴜs can be installed either with pip Python 3.6 and later.

************
Dependencies
************
Core dependencies:

* numpy_
* scipy_

Optional depedencies:

* matplotlib_
* pyvista_
* vtk_

The optional dependencies are used for 3D visualization and in some of the
Dispersion descriptor calculations.

***
pip
***

Clone the repository from GitHub. Install core functionality with pip:

.. code-block:: console

  pip install .

Install with extra graphics functionality:

.. code-block:: console

  pip install .[extras]

.. _matplotlib: https://pypi.org/project/matplotlib/
.. _numpy: https://pypi.org/project/numpy/
.. _pyvista: https://pypi.org/project/pyvista/
.. _scipy: https://pypi.org/project/scipy/
.. _vtk: https://pypi.org/project/vtk/
