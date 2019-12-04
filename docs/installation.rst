============
Installation
============

morfeus can be installed either with pip or conda for Python 3.

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

The optional dependencies are used for 3D visualization and in the Dispersion
descriptor calculations.

***
pip
***

Installation of core functionality with pip

.. code-block:: console

  pip install morfeus

Installation with extra graphics functionality

.. code-block:: console

  pip install morfeus[extras]


*****
conda
*****

Conda currently only supports installation with extra graphics functionality.

.. code-block:: console

  conda install -c conda-forge morfeus


.. _matplotlib: https://pypi.org/project/matplotlib/
.. _numpy: https://pypi.org/project/numpy/
.. _pyvista: https://pypi.org/project/pyvista/
.. _scipy: https://pypi.org/project/scipy/
.. _vtk: https://pypi.org/project/vtk/
