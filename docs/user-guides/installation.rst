Installation
============

Dragonfly is distributed on PyPI as ``dragonfly-spi`` and the source repository is hosted at ``github.com/SerialXI/dragonfly``. The Python package you import is named ``dragonfly``.

Installing the Python package builds Cython extensions, so a working compiler toolchain and the external scientific libraries used by Dragonfly must be available at install time.

Prerequisites
-------------

System dependencies:

* **MPI** (OpenMPI is recommended; the current build uses its ``mpicc --showme`` flags)
* **GSL** (GNU Scientific Library)
* **HDF5**
* **OpenMP** (usually included with compiler)

Python dependencies:

* Python 3.8+
* numpy
* scipy
* h5py
* mpi4py
* cython
* scikit-learn
* matplotlib and PyQt5 for the GUIs

We recommend using a virtual environment and the environment modules provided by your cluster.

Install from PyPI
-----------------

Install the published package with ``pip``::

    python -m pip install mpi4py
    python -m pip install dragonfly-spi

This installs the ``dragonfly`` Python package and its command-line tools.

On HPC systems, load the MPI module you want to use before installation so ``mpi4py`` and the compiled extensions link against the same MPI library::

    module load mpi/openmpi-x.x.x
    python -m pip install mpi4py
    python -m pip install dragonfly-spi

.. note::

    The build process expects tools such as ``mpicc``, ``gsl-config``, and ``h5cc`` to be available on ``PATH``.


Install from Source
-------------------

Clone the repository::

    git clone https://github.com/SerialXI/dragonfly.git
    cd dragonfly

Install Python dependencies::

    python -m pip install -r requirements.txt

.. note::

    On managed systems, load the required compiler, MPI, HDF5, and GSL modules
    before running the install commands.

Install the package in editable mode::

    python -m pip install -e .

If your environment has trouble with isolated builds, use::

    python -m pip install -e . --no-build-isolation

Verify Installation
-------------------

Check that the package imports::

    python -c "import dragonfly; print(dragonfly.__file__)"

Check that MPI is available::

    mpirun --version
    python -c "from mpi4py import MPI; print(MPI.Get_library_version())"

You can also verify that the installed command-line entry points resolve::

    dragonfly.init -h
    dragonfly.emc -h

The second command also checks that the reconstruction extension and ``mpi4py`` can be imported. On a machine with a display, ``dragonfly.autoplot -h`` and ``dragonfly.frameviewer -h`` check the GUI imports without opening a window.

Uninstall
---------

Uninstall the PyPI or editable install with::

    python -m pip uninstall dragonfly-spi

Troubleshooting
---------------

**Build cannot find MPI, HDF5, or GSL**

    Make sure the corresponding development packages or environment modules are loaded and that ``mpicc``, ``h5cc``, and ``gsl-config`` are on ``PATH``.

**mpi4py installation fails or imports the wrong MPI**

    Compare ``MPI.Get_library_version()`` with your selected MPI module. If necessary, rebuild ``mpi4py`` from source against that module::

        module load mpi/openmpi-x.x.x
        MPICC=mpicc python -m pip install --force-reinstall --no-cache-dir --no-binary=mpi4py mpi4py

    If you also change the MPI library used to compile Dragonfly, rebuild or reinstall Dragonfly against the same environment. Do not mix a conda MPI runtime with a different cluster MPI runtime in the same process.

**Import errors after installation**

    Try reinstalling with a fresh build::

        python -m pip uninstall dragonfly-spi
        python -m pip install -e . --no-build-isolation

**Missing GSL library**

    Install GSL development libraries (e.g., ``libgsl-dev`` on Debian/Ubuntu).

Next Steps
----------

* :doc:`simulation-quickstart` - Run your first simulation
* :doc:`config` - Configure classification or reconstruction using your own data
* :doc:`experimental-quickstart` - A 3D reconstruction tutorial using experimental data
