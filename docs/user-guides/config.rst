Configuration Guide
===================

Dragonfly uses configuration files (typically ``config.ini``) that organize parameters
by workflow section. The first few sections correspond to the data stream simulator workflow.

If processing experimental data with a known detector file, only the ``[emc]`` section is needed.

Geometry Parameters
-------------------

The ``[parameters]`` section defines the experimental geometry:

.. list-table::
   :header-rows: 1
   :widths: 20 60

   * - Parameter
     - Description
   * - ``detd``
     - Detector distance (millimeters)
   * - ``lambda``
     - Photon wavelength (Angstrom)
   * - ``detsize``
     - Linear detector size (pixels, integer)
   * - ``pixsize``
     - Pixel size (millimeters)
   * - ``stoprad``
     - Beamstop radius (pixels, float)
   * - ``polarization``
     - Polarization correction ('x', 'y', or 'none')

make_densities Section
----------------------

Creates electron density map from PDB file.

.. list-table::
   :header-rows: 1
   :widths: 25 55

   * - Parameter
     - Description
   * - ``pdb_code``
     - PDB code to fetch (alternative to ``in_pdb_file``)
   * - ``in_pdb_file``
     - Path to custom PDB file
   * - ``scatt_dir``
     - Directory with atomic scattering factors (Henke table)
   * - ``out_density_file``
     - Output electron density file

make_intensities Section
------------------------

Creates 3D intensity map from electron density with low-pass filtering.

.. list-table::
   :header-rows: 1
   :widths: 25 55

   * - Parameter
     - Description
   * - ``in_density_file``
     - Input density file (e.g., ``make_densities:::out_density_file``)
   * - ``out_intensity_file``
     - Output 3D intensity file

make_detector Section
---------------------

Generates detector geometry file.

.. list-table::
   :header-rows: 1
   :widths: 25 55

   * - Parameter
     - Description
   * - ``out_detector_file``
     - Output detector geometry file

make_data Section
-----------------

Simulates photon diffraction patterns.

.. list-table::
   :header-rows: 1
   :widths: 25 55

   * - Parameter
     - Description
   * - ``num_data``
     - Number of data frames to simulate
   * - ``fluence``
     - Incident photon fluence (photons/um^2)
   * - ``in_detector_file``
     - Detector geometry file
   * - ``in_intensity_file``
     - Input 3D intensity volume
   * - ``out_photons_file``
     - Output photon data file

EMC Section
-----------

Parameters for the EMC reconstruction.

.. list-table::
   :header-rows: 1
   :widths: 25 55

   * - Parameter
     - Description
   * - ``in_photons_file``
     - Input photon data file
   * - ``in_detector_file``
     - Detector geometry file
   * - ``num_div``
     - Quaternion sampling refinement (``10(5n^3 + n)`` rotations)
   * - ``beta``
     - Initial annealing value. Use ``auto`` to estimate per-frame starting values automatically.
   * - ``beta_factor``
     - Additional multiplicative factor applied to the starting ``beta`` values
   * - ``beta_schedule``
     - Increase the effective ``beta`` by multiplying by <first> every <second> iterations
   * - ``output_folder``
     - Directory for output files
   * - ``log_file``
     - Reconstruction log file
   * - ``need_scaling``
     - Fluence scaling (0=off, 1=on)

Configuration Chaining
----------------------

Use ``:::`` syntax to link outputs between sections::

    in_detector_file = make_detector:::out_detector_file
    in_photons_file = make_data:::out_photons_file

This chains the output of one module as input to another.
