Configuration Guide
===================

Dragonfly uses INI files (typically ``config.ini``) to describe input data, geometry, classification or reconstruction settings, and output locations.  The same configuration is used by ``dragonfly.emc``, ``dragonfly.frameviewer``, and ``dragonfly.autoplot``.

For existing photon data and an HDF5 detector file, only **[emc]** is needed.  The **[parameters]** and **[make_detector]** sections are needed if you generate geometry from experimental parameters. The other **[make_*]** sections are for simulation, not prerequisites for analyzing an experimental dataset.

The examples below are starting points, not settings specific to the original Dragonfly paper or prescriptions for all experiments. Defaults listed in the tables are the parser defaults when an option is omitted; the config copied by ``dragonfly.init`` may explicitly choose different values.

Starting Configurations
-----------------------

2D Classification
~~~~~~~~~~~~~~~~~

For exploratory analysis of serial diffraction patterns, learn a set of 2D classes without requiring a 3D model. Given matching photon and detector files::

    [emc]
    in_photons_file = data/photons.emc
    in_detector_file = data/detector.h5
    recon_type = 2d
    num_modes = 50
    num_rot = 180
    friedel_sym = 0
    need_scaling = 1
    beta = auto
    beta_schedule = 2.0 10
    output_folder = results_2d/
    log_file = EMC_2d.log

Save this as ``config_2d.ini`` and run::

    dragonfly.emc -c config_2d.ini -t 4 5
    dragonfly.autoplot -c config_2d.ini

``num_rot`` is the number of in-plane angular samples; here 180 samples cover 360 degrees, so their spacing is two degrees. Increase sampling if finer features require it. 

.. note::

   If ``friedel_sym = 1``, the range of angular sampling is only 180 degrees. Then ``num_rot = 180`` corresponds to a one degree spacing.

``num_modes`` is the number of learned classes. Choose a smaller ``num_modes`` especially if you have fewer than around 500 patterns per mode. With ``num_rot = 1``, no in-plane rotation search is performed, which can be useful when laboratory-frame orientation is meaningful. In practice, we have found more than 100 classes to not be practical, although your mileage may vary.

.. note::

   Choosing ``num_modes = 2*n^2`` like 32, 50 and 98 results in a neat rendering in the autoplot GUI without empty slots. This is, of course, a purely aesthetic condition!
   
Scaling allows patterns to be correctly compared while accounting for exposure variation. Class separation can still reflect backgrounds or artifacts; inspect member frames and occupancies as described in :doc:`autoplot`.

3D Reconstruction
~~~~~~~~~~~~~~~~~

For patterns expected to come from a common 3D intensity distribution::

    [emc]
    in_photons_file = data/photons.emc
    in_detector_file = data/detector.h5
    recon_type = 3d
    num_div = 6
    need_scaling = 1
    beta = auto
    beta_schedule = 2.0 10
    output_folder = results_3d/
    log_file = EMC_3d.log

``num_div`` controls 3D orientation sampling and is not the 2D rotation count.  See :doc:`../background/emc-algorithm` for the distinction and :doc:`experimental-quickstart` for an experimental 3D tutorial.

Use different output folders and log files for independent runs. A run without ``-r`` starts again and will overwrite earlier results. The positional CLI argument is the number of iterations to run, including when resuming; it is not an option inside ``[emc]``.

.. _orientation-refinement:

3D Orientation Refinement
~~~~~~~~~~~~~~~~~~~~~~~~~

Coarse-to-fine orientation refinement is a common way to reach high angular sensitivity in 3D without performing a global search on an increasingly dense grid. Two values of ``num_div`` specify the **fine sampling first, coarse sampling second**. For example, after a global reconstruction with ``num_div = 6``, the refinement setting is::

    num_div = 12 6

The first value sets the new orientation grid; the second identifies the coarse grid used for the earlier assignments. Fine orientations are mapped to their nearest coarse orientation (a Voronoi cell). Each frame is then compared only with fine orientations belonging to its allowed coarse cells, rather than with the entire fine grid. This reduces the cost of the search, but also means that an incorrect coarse assignment can constrain the subsequent refinement.

In that workflow, refinement continues from a coarse reconstruction's model, scales, and orientation information. The previous ``log_file`` and ``output_folder`` identify the checkpoint when resuming; the photon-file order, geometry, mode count, and symmetry must remain consistent with the earlier assignments. The second ``num_div`` value is the checkpoint's sampling, not an iteration count or a second simultaneous global search.

With ``save_prob = 1`` in the coarse run, its output retains multiple plausible coarse orientations per frame under ``probabilities/place``. Without saved probabilities, the standalone implementation falls back to each frame's most likely orientation. Saving probabilities before the final coarse iteration preserves more of the assignment uncertainty for refinement.

The standalone ``probs_file`` option can pin the coarse restrictions to a particular checkpoint, for example::

    num_div = 12 6
    probs_file = data/output_100.h5

Here ``output_100.h5`` must come from the coarse grid of 6 subdivisions. Without an explicit ``probs_file``, the standalone resume path uses the last iteration output. Keeping the coarse checkpoint separately and referring to it explicitly avoids confusing fine-grid probabilities with coarse-grid indices when continuing a refinement stage. ``probs_file`` is not currently read by the Python setup.

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
     - One integer for a square detector, or two integers for its x/y sizes,
       in pixels. Photon array ordering must match this geometry.
   * - ``pixsize``
     - Pixel size (millimeters)
   * - ``stoprad``
     - Beamstop radius (pixels, float)
   * - ``polarization``
     - ``x`` or ``y`` for linear polarization. ``none`` uses the generator's
       average of the two polarization factors, not a unit correction.
   * - ``ewald_rad``
     - Ewald radius in model voxels. Defaults to ``detd / pixsize`` if omitted.
       It sets reciprocal-space voxel size; see :doc:`detector-format`.

These configuration distances are in millimeters; the detector file stores ``detd`` in detector-pixel units. Do not interchange the two representations.

make_densities Section
----------------------

Creates an electron density map from a PDB file. This section is only needed for the simulation pipeline.

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
     - Directory containing ``henke_data.h5`` and ``atom_mass.txt``.
       The installed template uses ``aux/``.
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
     - Output HDF5 detector geometry file (use a ``.h5`` name).
   * - ``in_mask_file``
     - Optional raw binary ``uint8`` mask, one value per pixel in detector
       order. Values are ``0`` good, ``1`` non-orientational, ``2`` bad.
   * - ``center``
     - Optional beam-center pixel indices: one integer for a square detector
       or two for a rectangular detector. Defaults to the geometric center.

Without an input mask, the generator marks the beamstop as bad and detector corners outside its inscribed circle as non-orientational. With a supplied mask, it also marks otherwise good corner pixels as non-orientational.  See :doc:`detector-format` for mask meanings and the generated datasets.

For many experimental datasets, you will not use this section, but rather create a detector file using an available detector geometry file in another format and a bad pixel mask. See :doc:`../api/detector` for details about how to do this.

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
   * - ``mean_count``
     - Alternative target mean photons per pattern. Specify either a positive
       fluence or a positive mean count.
   * - ``do_gamma``
     - ``0`` (default) or ``1`` to simulate variable incident intensity.
   * - ``in_detector_file``
     - Detector geometry file
   * - ``in_intensity_file``
     - Input 3D intensity volume
   * - ``out_photons_file``
     - Output photon data file. ``.emc`` selects binary output; other
       extensions select sparse HDF5 output in the simulator.
   * - ``background_file``
     - Optional HDF5 file with a root-level ``background`` array in detector
       pixel order, used when simulating background.

See :doc:`simulator` for the simulation commands.

EMC Section
-----------

The **[emc]** section controls both 2D classification and 3D reconstruction.  All input pixel arrays must match the associated detector's count and order.

Inputs and Outputs
~~~~~~~~~~~~~~~~~~

Ordinary filename options are resolved relative to the configuration file.  For photon and detector lists, **use absolute paths for their entries** so the same files work with both reconstruction and GUI readers. Keep list order fixed when reusing frame-wise scales, blacklists, or external metadata.

.. list-table::
   :header-rows: 1
   :widths: 25 55

   * - Parameter
     - Description
   * - ``in_photons_file``
     - One photon data file. Use this or ``in_photons_list``, not both.
   * - ``in_photons_list``
     - | Text file containing one photon-file path per line. 
       | Frames are concatenated in list order.
   * - ``in_detector_file``
     - | One HDF5 detector geometry file, shared by the input photon files.
       | Use this or ``in_detector_list``, not both.
   * - ``in_detector_list``
     - | Text file with one detector-file path per photon file, in matching order. 
       | A single shared detector may instead use ``in_detector_file``.
   * - ``output_folder``
     - | Default ``data/``. Output HDF5 files are written here for each iteration
       | The folder is created if needed.
   * - ``log_file``
     - | Reconstruction log. Default ``EMC.log``.
       | Create its parent folder if you choose a subdirectory.
   * - ``start_model_file``
     - | Optional starting model; an HDF5 file uses its ``intens`` dataset.
       | Without a file, the model is randomly initialized.
   * - ``blacklist_file``
     - | Optional text file with one number per frame ``0`` to keep or ``1`` 
       | to exclude. See :doc:`autoplot` for generating from class selection.
   * - ``selection``
     - | Optional ``odd_only`` or ``even_only`` to split the otherwise usable 
       | frames into alternating subsets for independent analyses.

Sampling and Symmetry
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 55

   * - Parameter
     - Description
   * - ``recon_type``
     - ``3d`` (default) for volume reconstruction or ``2d`` for classification.
   * - ``num_modes``
     - Number of learned models/classes. Default ``1``; choose a positive
       integer.
   * - ``num_div``
     - | Positive subdivision count for 3D orientation sampling, required for a 3D run.
       | One value selects a global grid; two values mean fine/coarse refinement.
       | See :ref:`orientation-refinement` for the workflow and implementation status.
       | Not used to set 2D sampling.
   * - ``num_rot``
     - | Positive number of in-plane angular samples, required for a 2D run.
       | Covers 360 degrees, or 180 degrees with Friedel symmetry enabled.
   * - ``friedel_sym``
     - | Default ``0``. Set ``1`` to impose centrosymmetry on 2D classes and restrict 
       | sampling to 180 degrees. Enable only when justified; finite Ewald curvature can 
       | break centrosymmetry in a measured pattern even for real-valued objects.
   * - ``point_group``
     - 3D point group: ``1`` (default), ``A5`` (icosahedral), or ``S4`` (octahedral).
   * - ``axial_sym``
     - | Integer rotational symmetry order about the model's z axis.
       | Default ``1`` means no additional axial symmetry.

Annealing and Model Updates
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 55

   * - Parameter
     - Description
   * - ``beta``
     - | Initial inverse-temperature value. Default ``auto`` estimates a starting 
       | value per frame; a positive number sets the same start for every frame. 
       | Smaller values broaden assignment probabilities.
   * - ``beta_factor``
     - Additional positive multiplier on starting beta values. Default ``1``.
   * - ``beta_schedule``
     - | Two values: multiplier and period in iterations. Default ``1.0 100`` leaves 
       | beta unchanged. For ``2.0 10``, beta doubles at iterations 11, 21, etc. 
       | For each frame, beta is capped at one.
   * - ``alpha``
     - | Fraction of the previous model retained in an update. 
       | Default ``0`` uses the new estimate directly. 
       | Values between zero and one add inertia.

Intensity Scaling and Background
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 55

   * - Parameter
     - Description
   * - ``need_scaling``
     - | ``0`` (default) or ``1`` to account for per-frame intensity variation
       | due to incident fluence fluctuations.
   * - ``update_scale``
     - | ``1`` (default) updates scale estimates; ``0`` keeps their supplied or 
       | initialized values. Requires scaling to be enabled to be relevant.
   * - ``scale_file``
     - | Optional HDF5 file with a one-dimensional ``scale`` dataset in frame order, 
       | used to initialize scales. Omitted scales start at one.
   * - ``background_file``
     - | Optional HDF5 file with a ``background`` dataset, one value per pixel in 
       | detector order. Use this or ``background_list``, not both.
   * - ``background_list``
     - Optional list of background files with one entry per supplied detector entry.
   * - ``bgscale_file``
     - | Optional HDF5 file with a ``scale`` dataset giving per-frame background 
       | multipliers. Note that this is an input and is not updated.

The background-aware photon likelihood is used for sparse data with ``need_scaling = 1``. Merely storing ``background`` inside a detector file does not enable it: configure ``background_file`` explicitly, which may point to that same detector file.

Noise Model and Saved Probabilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sparse photon data and dense integer binary EMC use a Poisson model. Dense double binary EMC uses the Gaussian path: set a positive ``gaussian_sigma`` (standard deviation, in the units of the input intensities) and use ``need_scaling = 0``. Its likelihood does not fit per-frame intensity scales.  The default ``gaussian_sigma = 0`` is not valid for Gaussian data.

Set ``save_prob = 1`` to save sparse joint model/orientation probabilities in ``probabilities/place`` and ``probabilities/prob``. The default is ``0`` because these arrays can be large. For multi-mode runs, per-frame ``occupancies`` are saved even without ``save_prob``; they summarize class probabilities after summing over orientations. See :doc:`../api/recon` for other output datasets.

Stochastic Data Subsets
~~~~~~~~~~~~~~~~~~~~~~~

For large datasets, optional stochastic updates use a subset of frames per iteration. Keep the default full-data settings for a first run:

.. list-table::
   :header-rows: 1
   :widths: 25 55

   * - Parameter
     - Description
   * - ``data_fraction``
     - | Fraction of usable frames selected, greater than zero and at most one.
       | Default ``1`` uses all usable frames.
   * - ``data_fraction_period``
     - Positive number of iterations to keep a subset before refreshing it. Default ``1``.
   * - ``data_fraction_until``
     - | Last iteration using stochastic subsets; later iterations use all usable frames. 
       | Default ``0`` means no cutoff. Useful for doing a last 1-2 iterations for classification.
   * - ``coverage_bias``
     - | Positive values favor frames not recently selected. 
       | Default ``0`` disables this bias; values must be greater than ``-1``.

Metrics from different subsets are not directly comparable as though they were computed on a fixed dataset. For frame classification tasks, set ``data_fraction_until`` to 1 or 2 less than the targeted number of iterations to finally classify all frames with the final model.

Configuration Chaining
--------------------------------

Use ``:::`` syntax to link outputs between sections::

    in_detector_file = make_detector:::out_detector_file
    in_photons_file = make_data:::out_photons_file

This refers to the other option's value; it does not run that module.  When replacing the simulation template with experimental inputs, remove the old linked options rather than adding conflicting input file/list settings.

.. raw:: html

   <style>
   .rst-content table .line-block {
       margin-bottom: 0;
   }
   </style>
