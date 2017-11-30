
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

.. todo::
   Not finished.

.. warning::

   This script assumes that the ``fitstats`` script was previously used and
   hence that the ``fitstats.dat`` file exists in the ``output/fits_find``
   folder.


``find_stars`` parameters

* ``fits_find``        field
* ``filter_proc``      U
* ``find_method``      DAO
* ``thresh_find``      5.
* ``round_method``     AND
* ``round_max``        1.
* ``do_plots_F``       y


The ``find_stars`` script uses one of two methods to detect stars in your
observed field: IRAFStarFinder or DAOStarFinder.

The input folder is ``input/fits_find``

The output folder is created as ``output/fits_find``



PSF photometry on detected stars
--------------------------------

.. todo::
   Not finished.

``psf_phot`` parameters

* ``fits_psf``         field/stk_fcd0044.fits
* ``filter_psf``       V
* ``fitshape``         3.
* ``group_sep``        2.
* ``psf_thresh``       5.
* ``niters``           1
* ``do_plots_G``       y



Match and transform
-------------------

.. todo::
   Not finished.

``match_transf`` parameters

* ``match_folder``     field
* ``load_format``      allstar, list
* ``#load_format``      default
* ``ref_frame``        none
* ``maxrad``           5.
* ``method``           mean
* ``do_plots_H``       y