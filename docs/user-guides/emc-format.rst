EMC File Format
===============

Dragonfly stores diffraction patterns in EMC photon data files, used for both simulated and experimental data. This page describes the historical sparse binary format, its header extensions, and the alternative dense and HDF5 representations. For Python reading and writing examples, see :doc:`../api/emcfile`.

Sparse Photon Encoding
----------------------

High-resolution single-particle imaging patterns often contain few photons: most pixels receive zero or one photon. Rather than storing every pixel in each pattern, the sparse format stores only pixels with nonzero counts:

* For a pixel receiving exactly one photon, only its pixel index is stored.
* For a pixel receiving two or more photons, both its pixel index and its photon count are stored.
* Pixels absent from both lists have zero photons.

Pixel indices are zero-based, from ``0`` to ``num_pix - 1``, and refer to the pixel ordering in the detector geometry file, not assembled image coordinates.  The photon file does not contain detector geometry; its pixel count and ordering must match the detector file used for reconstruction.

See :doc:`detector-format` for the detector geometry datasets and mask values.

Binary Layout
------------------------

The sparse ``photons.emc`` format consists of six consecutive blocks.  All integers are signed 32-bit values (four bytes). The implementation uses native byte order, normally little-endian on current systems; the format does not specify a portable byte order.

The first block is a **1024-byte header**. Its first four bytes contain the number of patterns, ``num_data``, and the next four contain the number of pixels per pattern, ``num_pix``. In the historical format, the remaining 1016 bytes are zero-filled. Header extensions are described below.

Let ``N = num_data``, ``So = sum(ones)``, and ``Sm = sum(multi)``. Here an event means a nonzero pixel in a pattern: ``Sm`` counts multi-photon pixels, not the total number of photons in them.

.. list-table:: Sparse binary blocks in file order
   :header-rows: 1
   :widths: 10 20 20 50

   * - Block
     - Contents
     - Size (bytes)
     - Meaning
   * - 1
     - Header
     - 1024
     - Number of patterns and pixels, followed by reserved space.
   * - 2
     - ``ones``
     - ``4 * N``
     - Number of single-photon pixels in each pattern.
   * - 3
     - ``multi``
     - ``4 * N``
     - Number of multi-photon pixels in each pattern.
   * - 4
     - ``place_ones``
     - ``4 * So``
     - Pixel indices of all single-photon events.
   * - 5
     - ``place_multi``
     - ``4 * Sm``
     - Pixel indices of all multi-photon events.
   * - 6
     - ``count_multi``
     - ``4 * Sm``
     - Photon counts corresponding to the entries in ``place_multi``.

.. figure:: /_static/images/dataFormat.png
   :alt: Historical EMC sparse format with a 1024-byte header, two per-pattern count arrays, and three concatenated event arrays
   :width: 280px

   Original historical sparse layout from the
   `Dragonfly wiki <https://github.com/SerialXI/Dragonfly/wiki/Data-stream-simulator#make_data>`_.
   The header extensions below occupy previously reserved bytes; the five
   sparse data blocks are unchanged.

Each event array concatenates patterns in file order. For pattern ``d``, its single-photon entries begin at ``sum(ones[:d])`` and have length ``ones[d]``.  Its multi-photon entries begin at ``sum(multi[:d])`` and have length ``multi[d]`` in both ``place_multi`` and ``count_multi``. Thus, the per-pattern count arrays define frame boundaries without separators in the event arrays.

The sparse file size is ``1024 + 8*N + 4*So + 8*Sm`` bytes.

Current Header Extensions
-------------------------

Going forward, the header also carries explicit data type information and an EMC file signature in previously reserved space. The header remains **1024 bytes**, and sparse data still begin at byte offset **1024**.

.. list-table:: Extended binary header (zero-based byte offsets)
   :header-rows: 1
   :widths: 20 20 60

   * - Bytes
     - Type
     - Contents
   * - 0--3
     - 32-bit integer
     - ``num_data``: number of patterns.
   * - 4--7
     - 32-bit integer
     - ``num_pix``: number of pixels per pattern.
   * - 8--11
     - 32-bit integer
     - Data type: ``0`` sparse, ``1`` dense integer, ``2`` dense double.
   * - 12--19
     - Eight bytes
     - EMC signature: ``89 45 4D 43 0D 0A 1A 0A`` (hexadecimal), or
       ``b'\x89EMC\r\n\x1a\n'`` in Python.
   * - 20--1023
     - Reserved space
     - Normally zero-filled; may contain additional application-specific
       header values.

.. note::

   Historical sparse files have zero in the data type field, which already denotes sparse data, and no EMC signature. Current Dragonfly readers accept files without the signature, so historical sparse files remain readable without conversion.

Dense Binary Data
-----------------

For data type ``1`` or ``2``, the same 1024-byte header is followed by ``num_data * num_pix`` values, with all pixels of the first pattern followed by all pixels of the second, and so on:

* ``1``: signed 32-bit integers (four bytes per pixel).
* ``2``: 64-bit double-precision floating-point values (eight bytes per pixel).

Dense files do not contain the sparse count or event arrays. The binary ``EMCWriter`` writes sparse data (type ``0``).

HDF5 Storage
------------

Dragonfly also supports sparse photon data in HDF5 files, conventionally named ``.h5``. This representation has no 1024-byte EMC binary header. It contains:

* ``num_pix``: the number of pixels per pattern.
* ``place_ones``: one variable-length array of single-photon pixel indices per pattern.
* ``place_multi``: one variable-length array of multi-photon pixel indices per pattern.
* ``count_multi``: one variable-length array of corresponding photon counts per pattern.

The three event datasets have shape ``(num_data,)``, with each entry holding a variable-length array of 32-bit integers. Array lengths replace the binary format's separate ``ones`` and ``multi`` arrays.

``EMCWriter`` defaults to HDF5 output; pass ``hdf5=False`` to write the binary format. ``EMCReader`` detects HDF5 from the file contents rather than its extension. It can also read dense HDF5 datasets when a dataset name is supplied; see :doc:`../api/emcfile` for the API.

Browsing Versus Reconstruction Inputs
-------------------------------------

The reconstruction loader accepts the sparse and dense binary EMC formats described above, and the sparse EMC-HDF5 representation. Sparse and dense integer data are used with the photon-counting (Poisson) likelihood; dense double data use the Gaussian path and require a positive ``gaussian_sigma``.  See :doc:`config` for the noise-model settings.

``EMCReader`` additionally supports browsing an arbitrary dense HDF5 dataset through ``dset_list``. The reconstruction loader does not have an equivalent dataset-name option: it expects HDF5 photon data to have ``num_pix``, ``place_ones``, ``place_multi``, and ``count_multi``. Convert a dense HDF5 photon-count dataset with ``EMCWriter`` before reconstruction; a short example is provided in :doc:`experimental-quickstart`.
