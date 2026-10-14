Data Stream Simulator
=====================

The data stream simulator generates synthetic single-particle diffraction data
from PDB files, including realistic noise and photon statistics.

sim_setup
---------

Convenience utility that runs the simulator pipeline sequentially:

.. py:function:: dragonfly.utils.sim_setup.main()
   :noindex:

    Runs through simulation utilities to generate data.

    This is a convenience wrapper that sequentially executes:
    * :func:`make_densities <dragonfly.utils.make_densities.make_dens>`
    * :func:`make_intensities <dragonfly.utils.make_intensities.make_intens>`
    * :func:`make_detector <dragonfly.utils.make_detector.make_detector>`
    * :func:`make_data <dragonfly.utils.make_data.make_data_cdef>`

    Usage::

        dragonfly.utils.sim_setup
        dragonfly.utils.sim_setup -c config.ini
        dragonfly.utils.sim_setup -y

make_densities
--------------

Creates electron density map from PDB file using atomic scattering factors.

.. py:function:: dragonfly.utils.make_densities.make_dens(config_fname, yes=False, verbose=False)
   :noindex:


    Generate density map from parameters in config file.

    :param config_fname: Path to configuration file.
    :type config_fname: str
    :param yes: Skip confirmation prompts.
    :type yes: bool
    :param verbose: Enable verbose logging.
    :type verbose: bool
    :returns: Outputs binary density file specified in config.
    :rtype: None

    **Usage**::

        dragonfly.utils.make_densities
        dragonfly.utils.make_densities -c config.ini

    **Output**: Binary file containing 3D cubic electron density distribution
    at the appropriate voxel size.

    A low-pass filter eliminates artifacts from electron density
    discretization on a cubic grid.

make_intensities
----------------

Creates 3D intensity map from electron density with low-pass filtering.

.. py:function:: dragonfly.utils.make_intensities.make_intens(config_fname, yes=False, verbose=False)
   :noindex:


    Generate intensity volume from config file parameters.

    :param config_fname: Path to configuration file.
    :type config_fname: str
    :param yes: Skip confirmation prompts.
    :type yes: bool
    :param verbose: Enable verbose logging.
    :type verbose: bool
    :returns: Outputs binary intensity file specified in config.
    :rtype: None

    **Usage**::

        dragonfly.utils.make_intensities
        dragonfly.utils.make_intensities -c config.ini

    **Output**: Binary file containing 3D intensity distribution in reciprocal space.

make_detector
-------------

Generates detector geometry file with pixel information.

.. py:function:: dragonfly.utils.make_detector.make_detector(config_fname, yes=False, verbose=False)
   :noindex:


    Generate detector file from parameters in config file.

    :param config_fname: Path to configuration file.
    :type config_fname: str
    :param yes: Skip confirmation prompts.
    :type yes: bool
    :param verbose: Enable verbose logging.
    :type verbose: bool
    :returns: Outputs detector file specified in config.
    :rtype: None

    **Usage**::

        dragonfly.utils.make_detector
        dragonfly.utils.make_detector -c config.ini

    **Output**: HDF5 detector geometry file (conventionally ``.h5``) containing:

    * Reciprocal-space pixel coordinates
    * Solid-angle and polarization correction factors
    * Mask information
    * Detector distance in detector-pixel units
    * Ewald sphere radius

    See :doc:`detector-format` for the dataset layout, units, and mask values.

make_data
---------

Simulates photon diffraction patterns from 3D intensity volume.

.. py:function:: dragonfly.utils.make_data.make_data_cdef(config_fname, yes=False, verbose=False)


    Generate simulated photon data from configuration file.

    :param config_fname: Path to configuration file.
    :type config_fname: str
    :param yes: Skip confirmation prompts.
    :type yes: bool
    :param verbose: Enable verbose logging.
    :type verbose: bool
    :returns: Outputs photon data file specified in config.
    :rtype: None

    **Usage**::

        dragonfly.utils.make_data
        dragonfly.utils.make_data -c config.ini

    **Output**: Binary file (``.emc`` or ``.h5``) containing sparse photon data.

Data Format
-----------

The simulator outputs photon data in sparse format with three arrays:

#. **Single-photon locations** (variable):
   - Pixel indices where exactly one photon was detected

#. **Multi-photon locations** (variable):
   - Pixel indices for multi-photon events

#. **Multi-photon counts** (variable):
   - Photon count at each multi-photon pixel

This sparse format is efficient for high-resolution SPI where most
pixels receive zero or one photon.

See :doc:`emc-format` for the binary block layout, historical format diagram,
header extensions, and HDF5 representation.

Pipeline Flow
-------------

::

    PDB file
       |
       v
    make_densities --> densityMap.bin
       |
       v
    make_intensities --> intensities.bin
       |
       v
    make_detector --> det_geom.h5
       |
       v
    make_data --> photons.emc
