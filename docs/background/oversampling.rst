Oversampling in Single-Particle Imaging
=======================================

What is Oversampling?
---------------------

In coherent diffractive imaging, a speckle is the diffraction-space analogue of
a Shannon pixel. Its width is set mainly by the real-space size of the object.
If the detector samples each speckle with several detector pixels, the pattern
is oversampled.

Some oversampling is useful, but too much is wasteful. When a single speckle is
very wide in detector pixels, neighboring pixels carry highly redundant
information while the data volume and FFT cost continue to grow.

Why Dragonfly Warns About It
----------------------------

During simulation setup, Dragonfly estimates the reciprocal-space field of view
implied by the detector geometry and compares it with the size of the input
density map. If that combination would require an excessively padded FFT grid,
the code prints a warning that the oversampling ratio is high.

This warning is mainly about efficiency. A very large oversampling ratio means
the simulated diffraction pattern is being sampled much more finely than needed.
The extra pixels usually do not add useful information, but they do increase
memory use and runtime.

Oversampling in Simulations
---------------------------

The ``dragonfly.utils.make_intensities`` tool reads a real-space density map,
pads it, and Fourier transforms it to produce the 3D intensity used for EMC
simulations. The warning appears when the density size is too small relative to
the reciprocal-space sampling implied by the detector geometry. In practical
terms, the detector pixels are much finer than the speckle width.

Two common causes are:

* the detector is too far from the sample, which makes speckles spread across
  many detector pixels
* the detector pixel grid is unnecessarily fine for the chosen geometry

How To Reduce Excessive Oversampling
------------------------------------

If Dragonfly warns that the oversampling ratio is high, the usual remedies are:

#. Reduce the detector distance ``detd`` if that does not hide too much signal
   behind the central beamstop or detector hole.
#. Bin detector pixels by reducing ``detsize`` and increasing ``pixsize`` so the
   physical detector size stays consistent while each effective pixel covers a
   larger area.

Both changes make a speckle span fewer detector pixels, which reduces the amount
of zero-padding needed in the simulated intensity volume.

What EMC Does Not Do
--------------------

EMC reconstructs the reciprocal-space intensity volume. Recovering a real-space
density or structure still requires phase retrieval or another downstream
inversion method.
