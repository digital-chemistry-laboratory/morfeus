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

* matplotlib_
* pymeshfix_
* pyvista_
* pyvistaqt_
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

.. _matplotlib: https://pypi.org/project/matplotlib/
.. _numpy: https://pypi.org/project/numpy/
.. _pymeshfix: https://pypi.org/project/pymeshfix/
.. _pyvista: https://pypi.org/project/pyvista/
.. _pyvistaqt: https://pypi.org/project/pyvistaqt/
.. _scipy: https://pypi.org/project/scipy/
.. _vtk:  https://pypi.org/project/vtk/
.. _xtb-python: https://github.com/grimme-lab/xtb-python
