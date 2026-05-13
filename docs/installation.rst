============
Installation
============

ᴍᴏʀғᴇᴜs can be installed using pip or conda and works with Python >= 3.10.

.. tab:: pip

  .. code-block:: shell

    $ pip install morfeus-ml

.. tab:: conda

  .. code-block:: shell

    $ conda install -c conda-forge morfeus-ml

************
Dependencies
************

It is recommended that dependencies are installed with conda if possible as the
performance of the linear algebra backends is generally better than with pip.

In some instances, such as when using `xtb`, it can help to switch the
`BLAS implementation`_ used by conda. On Linux or Windows, MKL often works
better than OpenBLAS, and it can be used with:

.. code-block:: shell

  conda install "libblas=*=*mkl"

Core dependencies:

* fire_
* numpy_
* packaging_
* scipy_

Optional depedencies:

* cclib_
* dftd4_
* geometric_
* libconeangle_
* matplotlib_
* openbabel_
* pymeshfix_
* pyvista_
* pyvistaqt_
* qcengine_
* rdkit_
* spyrmsd_
* vtk_
* xtb_
* xtb-python_

.. _cclib: https://cclib.github.io/
.. _dftd4: https://github.com/dftd4/dftd4
.. _fire: https://github.com/google/python-fire
.. _geometric: https://github.com/leeping/geomeTRIC
.. _libconeangle: https://github.com/digital-chemistry-laboratory/libconeangle
.. _matplotlib: https://matplotlib.org
.. _numpy: https://numpy.org
.. _openbabel: http://openbabel.org/
.. _packaging: https://github.com/pypa/packaging
.. _pymeshfix: https://github.com/pyvista/pymeshfix
.. _pyvista: https://github.com/pyvista/pyvista
.. _pyvistaqt: https://github.com/pyvista/pyvistaqt
.. _qcengine: https://github.com/MolSSI/QCEngine
.. _rdkit: https://www.rdkit.org
.. _scipy: https://github.com/pyvista/pyvistaqt
.. _spyrmsd: https://github.com/RMeli/spyrmsd
.. _vtk:  https://vtk.org
.. _xtb: https://xtb-docs.readthedocs.io/en/latest/
.. _xtb-python: https://github.com/grimme-lab/xtb-python
.. _BLAS implementation: https://conda-forge.org/docs/maintainer/knowledge_base.html#switching-blas-implementation
