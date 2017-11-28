
Crowded field photometry
========================

.. todo::
   Not finished.

There are three steps for obtaining your final (matched, combined, and
transformed) magnitudes for your observed field:

* `Find stars`_: identify the coordinates of the stars in the
  observed fields.
* `PSF photometry on detected stars`_: obtain PSF photometry for
  each of the observed fields.
* `Match and transform`_: match stars for all filters and exposures,
  and transform their magnitudes to the standard system.


Find stars
----------

The ``find_stars`` script uses one of two methods to detect stars in your
observed field: IRAFStarFinder or DAOStarFinder.

The input folder is ``input/fits_find``

The output folder is created as ``output/fits_find``

.. warning::

   This script assumes that the ``fitstats`` script was previously used and
   hence that the ``fitstats.dat`` file exists in the ``output/fits_find``
   folder.



PSF photometry on detected stars
--------------------------------

.. todo::
   Not finished.


Match and transform
-------------------

.. todo::
   Not finished.