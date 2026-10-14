EMC File I/O Module
===================

The EMC file I/O module provides classes for reading and writing EMC photon data files.

For the on-disk layout, header fields, and historical format diagram, see :doc:`../user-guides/emc-format`.

Usage Example
-------------

**Reading EMC files:**

.. code-block:: python

   from dragonfly import Detector, EMCReader

   det = Detector('detector.h5')
   emc = EMCReader('photons.emc', det)

   aframe = emc.get_frame(0) # Assembled masked array
   frame = emc.get_frame(0, raw=True) # Unassembled dense array
   # Sparse format (fastest to read)
   place_ones, place_multi, count_multi = emc.get_frame(0, sparse=True)

   powder = emc.get_powder() # Virtual powder sum (assembled)

**Writing EMC files:**

.. code-block:: python

   import numpy as np
   from dragonfly import EMCWriter

   # Write a random dataset with a dispersion of 1 ph/pixel on average
   num_pix = 10000
   with EMCWriter('output.emc', num_pix, hdf5=False) as emc:
       for i in range(100):
           emc.write_frame(np.random.poisson(1.0, num_pix))
           
**Writing subsets of data:**

.. code-block:: python
   
   from dragonfly import Detector, EMCReader, EMCWriter
   
   det = Detector('detector.h5')
   emc = EMCReader('original.emc', det)

   # Write every 10th frame to a new file
   with EMCWriter('every10.emc', det.num_pix, hdf5=False) as wemc:
       for i in range(0, emc.num_frames, 10):
           wemc.write_sparse_frame(*emc.get_frame(i, sparse=True))

   # Write dataset with 1/10th ph/pixel (randomly chosen) 
   with EMCWriter('lowsignal.emc', det.num_pix, hdf5=False) as wemc:
       for i in range(emc.num_frames):
           wemc.write_frame(emc.get_frame(i, raw=True), fraction=0.1)

.. automodule:: dragonfly.emcfile
   :members:
   :undoc-members:
   :show-inheritance:
