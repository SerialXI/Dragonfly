Detector File Format
====================

A detector file describes how each pixel samples reciprocal space, its
intensity correction factor, and whether it should be used in reconstruction.
It contains geometry, not diffraction patterns; photon counts are stored
separately in the :doc:`emc-format`.

The current Python package reads and writes **HDF5** detector files,
conventionally named ``detector.h5``. The historical ASCII format is described
below for reference, but is no longer accepted by ``Detector`` or ``CDetector``.
See :doc:`../api/detector` for the Python API.

HDF5 Datasets
-------------

All datasets are at the root of the HDF5 file. Let ``N`` be the number of
detector pixels. The current writer produces the following datasets:

.. list-table:: Detector datasets
   :header-rows: 1
   :widths: 20 15 15 50

   * - Dataset
     - Shape
     - Data type
     - Meaning
   * - ``qx``
     - ``(N,)``
     - ``float64``
     - x coordinate in reciprocal-space voxel units.
   * - ``qy``
     - ``(N,)``
     - ``float64``
     - y coordinate in reciprocal-space voxel units.
   * - ``qz``
     - ``(N,)``
     - ``float64``
     - z coordinate in reciprocal-space voxel units.
   * - ``corr``
     - ``(N,)``
     - ``float64``
     - Multiplicative intensity correction for each pixel.
   * - ``mask``
     - ``(N,)``
     - ``uint8``
     - Pixel classification: ``0``, ``1``, or ``2`` (see below).
   * - ``detd``
     - Scalar
     - ``float64``
     - Detector distance in detector-pixel units.
   * - ``ewald_rad``
     - Scalar
     - ``float64``
     - Ewald sphere radius in reciprocal-space voxel units.
   * - ``background``
     - ``(N,)``
     - ``float64``
     - Optional per-pixel background intensity.

The number of pixels is inferred from the length of ``qx``; there is no
required ``num_pix`` dataset. All five required per-pixel arrays must be
one-dimensional and have the same length. The file does not require separate
pixel-size, beam-center, or assembled-image-shape datasets.

``Detector.write()`` always writes HDF5, regardless of the filename extension.
Use ``.h5`` to make the format clear. Readers identify HDF5 from the contents,
not the filename.

Pixel Ordering
--------------

Entry ``i`` in each per-pixel array describes pixel ``i`` in the photon data.
Indices are zero-based. The detector and photon files must have the same pixel
count **and ordering**, including masked pixels. Changing the detector order
without changing the photon data associates counts with the wrong geometry.

For a rectangular detector, ``make_detector`` flattens arrays of shape
``(dets_x, dets_y)`` in NumPy C order: the second index varies fastest. Custom
detectors can use any ordering, provided the photon data use that same ordering.
Assembled image coordinates are derived from the geometry when loading a
``Detector``; they are not the photon file's pixel indices.

Coordinates and Units
---------------------

``qx``, ``qy``, and ``qz`` give scattering-vector coordinates relative to the
origin of the reciprocal-space model, measured in **model voxels**, not inverse
Angstroms. ``ewald_rad`` uses the same voxel units.

The Ewald radius therefore sets the **reciprocal-space voxel size**. For
elastic scattering at wavelength :math:`\lambda`, the physical Ewald sphere
radius is :math:`1/\lambda`. Writing :math:`R` for ``ewald_rad``, that radius
measured in voxels, one voxel corresponds to

.. math::

   \Delta q = \frac{1}{\lambda R}.

If the wavelength is in Angstroms, :math:`\Delta q` is in inverse Angstroms.
Multiply the stored ``qx``, ``qy``, and ``qz`` by :math:`\Delta q` to obtain
physical reciprocal-space coordinates. At fixed wavelength, increasing
``ewald_rad`` decreases the voxel size, giving finer reciprocal-space sampling.

Under the alternative scattering-vector convention where the Ewald sphere
radius is :math:`2\pi/\lambda`, the voxel size is instead

.. math::

   \Delta q = \frac{2\pi}{\lambda R}.

Use the same convention consistently when converting voxel coordinates to
physical units.

The stored ``detd`` is not the millimeter value from ``config.ini``.
``make_detector`` converts it using::

    detd_in_file = detd_in_config / pixsize_in_config

Both configuration values use millimeters, so the result is in detector pixels.
The pixel size and beam center are used to calculate geometry rather than
stored as separate datasets.

For a flat detector perpendicular to the incident beam, ``calc_from_coords``
uses beam-relative pixel coordinates ``cx``, ``cy``, distance ``D = detd``,
and radius ``R = ewald_rad``::

    L = sqrt(cx**2 + cy**2 + D**2)
    qx = R * cx / L
    qy = R * cy / L
    qz = R * (D / L - 1)

``cx``, ``cy``, and ``D`` must use the same units; use detector pixels for
compatibility with image assembly. With this convention, the direct beam has
``q = (0, 0, 0)`` and off-axis pixels have negative ``qz``.

Correction Factors
------------------

``corr`` is the multiplicative factor that maps model intensity to expected
pixel intensity. The standard geometry generator includes solid-angle and
polarization corrections. It does not store corrected photon counts.

By default, ``CDetector`` divides the loaded correction factors by their mean
over all pixels, giving mean correction one. This changes the in-memory array,
not the file. Use ``CDetector(fname, norm=False)`` when the original correction
normalization must be preserved.

Mask Values
-----------

.. list-table:: Stored pixel classifications
   :header-rows: 1
   :widths: 10 25 65

   * - Value
     - Classification
     - Reconstruction behavior
   * - ``0``
     - Good
     - Used for orientation determination and model updates.
   * - ``1``
     - Non-orientational
     - Excluded from orientation determination, but included in model updates.
   * - ``2``
     - Bad
     - Excluded from both orientation determination and model updates.

Value ``1`` is useful for valid photon data in regions with poor reciprocal-space
coverage, such as detector corners. Value ``2`` is appropriate for beamstop,
dead, or otherwise unreliable pixels. Keep their entries in the detector and
photon arrays unless both are reindexed together.

The stored ``mask`` is exposed as ``raw_mask`` in Python. Do not confuse it
with ``Detector.mask``, which is a derived Boolean display mask where ``True``
means usable. By default, image assembly includes both raw mask ``0`` and ``1``;
``keep_mask_1=False`` excludes the latter from display as well.

Optional Background
-------------------

``Detector.write()`` includes a ``background`` dataset if the object's
background array is populated. However, loading a detector does **not**
automatically load or enable this background. Use
``det.parse_background('background.h5')`` to load the root-level ``background``
dataset, with one value per pixel in detector order. The same file may hold
both geometry and background.

For reconstruction, configure ``background_file`` (or ``background_list`` for
multiple detectors) in the ``[emc]`` section to load the background explicitly.

Creating a Detector File
------------------------

For geometry specified in ``config.ini``, run::

    dragonfly.utils.make_detector -c config.ini

The output path comes from ``[make_detector] out_detector_file``; see
:doc:`config` and :doc:`simulator`.

For custom geometry, the following example creates an illustrative 4-by-4
flat detector. Replace its coordinates, distance, Ewald radius, and mask with
values appropriate for your experiment:

.. code-block:: python

   import numpy as np
   from dragonfly import Detector

   det = Detector()
   cx, cy = np.meshgrid(np.arange(4) - 1.5, np.arange(4) - 1.5,
                        indexing='ij')
   det.cx = cx.ravel()
   det.cy = cy.ravel()
   det.detd = 1000.0  # Detector-pixel units, not millimeters
   det.ewald_rad = 100.0  # Model voxel units
   det.calc_from_coords(pol='x')
   det.raw_mask = np.zeros(det.cx.size, dtype='u1')
   det.write('detector.h5')

   loaded = Detector('detector.h5')
   print(loaded.num_pix)  # 16

Historical ASCII Format
-----------------------

The original ``.dat`` format was a text file whose first line contained the
number of pixels, followed by one line per pixel with five whitespace-separated
columns::

    num_pix
    qx qy qz corr mask
    qx qy qz corr mask
    ...

Later ASCII files extended the first line to include distance and Ewald radius::

    num_pix detd ewald_rad

The coordinate, correction, and mask meanings match the HDF5 fields above.
The first pixel row corresponds to pixel index ``0``. This is the format
described in the old
`Dragonfly wiki <https://github.com/SerialXI/Dragonfly/wiki/Data-stream-simulator#make_detector>`_.

.. note::

   Current Python readers require HDF5. Renaming a text file to ``.h5`` does
   not convert it. To migrate an ASCII file, store its five columns in the
   corresponding HDF5 datasets and provide scalar ``detd`` and ``ewald_rad``.
   For a one-value historical header, these two values must be recovered from
   the experimental geometry or configuration, using the units described above.
