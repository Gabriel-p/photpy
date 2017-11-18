
Crowded field photometry
========================

.. todo::
   Not finished.

There are three steps for obtaining your final (matched, combined, and
transformed) magnitudes for your observed field:

* `Find stars`_: identify the coordinates of the stars in the
  observed field.
* `PSF photometry on detected stars`_: obtain aperture photometry for
  each of the selected standard star frames.
* `Define and solve the transformation equations`_: define and fit the
  instrumental to standard transformation equations.


Find stars
-----------------------

The ``find_stars`` script uses one of two methods to detect stars in your
observed field: IRAFStarFinder or DAOStarFinder.

The input folder is ``input/fits_find``

The output folder is created as ``output/fits_find``

.. warning::

   This script assumes that the ``fitstats`` script was previously used and
   hence that the ``fitstats.dat`` file exists in the ``output/fits_find``
   folder.