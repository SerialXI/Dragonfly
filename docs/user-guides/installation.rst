Installation
============

Dragonfly is distributed on PyPI as ``dragonfly-spi`` and the source repository
is hosted at ``github.com/SerialXI/dragonfly``. The Python package you import is
named ``dragonfly``.

Installing the Python package builds Cython extensions, so a working compiler
toolchain and the external scientific libraries used by Dragonfly must be
available at install time.

Prerequisites
-------------

System dependencies:

* **MPI** (OpenMPI or MPICH)
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

Install from PyPI
-----------------

Install the published package with ``pip``::

    python -m pip install dragonfly-spi

This installs the ``dragonfly`` Python package and its command-line tools.

On HPC systems, load the MPI module you want to use before installation so
``mpi4py`` and the compiled extensions link against the same MPI library::

    module load mpi/openmpi-x.x.x
    python -m pip install dragonfly-spi

.. note::

    The build process expects tools such as ``mpicc``, ``gsl-config``, and
    ``h5cc`` to be available on ``PATH``.

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

You can also verify that the installed command-line entry points resolve::

    dragonfly.init -h

Uninstall
---------

Uninstall the PyPI or editable install with::

    python -m pip uninstall dragonfly-spi

Troubleshooting
---------------

**Build cannot find MPI, HDF5, or GSL**

    Make sure the corresponding development packages or environment modules are
    loaded and that ``mpicc``, ``h5cc``, and ``gsl-config`` are on ``PATH``.

**mpi4py installation fails or imports the wrong MPI**

    Reinstall after loading the MPI module you actually intend to use::

        module load mpi/openmpi-x.x.x
        python -m pip install --force-reinstall mpi4py dragonfly-spi

**Import errors after installation**

    Try reinstalling with a fresh build::

        python -m pip uninstall dragonfly-spi
        python -m pip install -e . --no-build-isolation

**Missing GSL library**

    Install GSL development libraries (e.g., ``libgsl-dev`` on Debian/Ubuntu).

Next Steps
----------

* :doc:`simulation-quickstart` - Run your first simulation
* :doc:`experimental-quickstart` - Process publicly available experimental data
