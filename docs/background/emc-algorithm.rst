The EMC Algorithm
=================

Overview
--------

Dragonfly implements the Expand-Maximize-Compress (EMC) algorithm for
single-particle imaging. The input is a set of diffraction patterns recorded
from nominally identical particles in unknown orientations. The main output is
a 3D reciprocal-space intensity volume whose central sections are consistent
with those measured 2D patterns.

Algorithm Principle
-------------------

The reconstruction problem has two coupled unknowns:

#. the 3D intensity distribution
#. the orientation of every measured frame

EMC addresses this by alternating between estimating how well each frame matches
many candidate orientations and updating the 3D model from those weighted
matches. No orientation labels are required in advance.

The Three Steps
---------------

Each EMC iteration consists of:

**1. Expand**

Sample the current 3D model on the orientation grid used for the search.
For each frame, Dragonfly compares the observed photon counts with the model
prediction for every tested orientation.

**2. Maximize**

Convert those comparison scores into orientation probabilities and accumulate
their weighted contributions back into reciprocal space. This is the step that
updates the model toward the maximum-likelihood solution.

**3. Compress**

Merge the expanded orientation-specific information back onto the regular 3D
grid that stores the reconstructed intensity volume. The next iteration starts
again from this compressed representation.

What Dragonfly Reconstructs
---------------------------

Dragonfly reconstructs reciprocal-space intensity, not a real-space electron
density map. In a typical workflow:

#. EMC estimates a 3D intensity volume.
#. That intensity can then be used in later analysis or phase retrieval.
#. Recovering a real-space density requires additional steps beyond EMC itself.

Orientation Sampling
--------------------

Orientations are sampled on a quaternion grid. In Dragonfly this is controlled
primarily by ``num_div``. A larger ``num_div`` means finer angular sampling and
more orientations to test per iteration.

For 3D reconstructions, the number of sampled orientations is:

.. math::

    10\,(5\,\mathrm{num\_div}^3 + \mathrm{num\_div})

This count grows quickly, so increasing ``num_div`` improves angular resolution
at a significant computational cost.

Deterministic Annealing
-----------------------

Dragonfly also supports deterministic annealing through the ``beta`` and
``beta_schedule`` settings. Early iterations can use broader orientation
probabilities; later iterations gradually sharpen them. This often helps the
reconstruction avoid poor local optima, especially for cleaner data.

Convergence and Practical Use
-----------------------------

In practice, users monitor quantities such as log-likelihood, RMS model change,
and the visual stability of the reconstructed volume. EMC is iterative and its
behavior depends strongly on data quality, detector geometry, masking, angular
sampling, and annealing settings.

The algorithm is designed for large sets of sparse diffraction patterns, which
is why Dragonfly parallelizes the work with MPI and OpenMP.

References
----------

* Loh, D. M. & Elser, V. (2009). Phys. Rev. E 80, 026705.
* Ayyer et al. (2016). J. Appl. Cryst. 49, 1320-1335.
