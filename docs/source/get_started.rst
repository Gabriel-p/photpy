
Getting started
===============

.. warning::
   Not finished.


Folder structure
----------------

The code expects standard and field observed frames to exist separately in
corresponding folders within the ``photpy/input/`` folder. The folder structure
must look like this:


.. code-block:: none

    input/
    │
    ├── standards/
    │   ├── stand_001.fits
    │   ├── stand_002.fits
    │   ├── stand_003.fits
    │   └── ...
    │
    ├── field/
    │   ├── field_100.fits
    │   ├── field_101.fits
    │   ├── field_102.fits
    │   └── ...


The names of the ``.fits`` files are not important, as long as their headers
contain all the required information. See the next section to learn of to check
if your headers are correct.


Fixing your headers
-------------------

The required information that needs to be present in the header of each .fits
file is:

* ``Gain``
* ``Read noise``
* ``Filter``
* ``Exposure time``

These keys should all be present and equal throughout all files that will be
processed. To display the header of a .fits file you can use the following code:

.. code-block:: python

    from astropy.io import fits

    # Load .fits file
    hdulist = fits.open(image_file)
    # Extract header data
    hdr = hdulist[0].header
    # Display header data
    for k, v in zip(*[hdr.keys(), hdr.values()]):
        print(k, v)


.. _secinput:

Input parameters
----------------


The required input information for all the scripts is listed in its ``.pars``
file. This file can be accessed within the ``/tasks`` folder, or simply filled
when the script is called. A description of each required parameter is presented
below.

.. code-block:: none

    ff_proc        Name of the .fits file to be processed, or a folder
                   containing more than one .fits files.
    do_plots       Flag to determine whether the output plot is produced.
                   Accepted inputs are y/n.
    dmax           Maximum flux value of a non-saturated star.
    thresh_level   Threshold detection level in units of the sky's STDDEV, used
                   by DAOStarFinder.
    fwhm_init      Initial estimate of the FWHM, used by DAOStarFinder.
    max_stars      Maximum number of bright unsaturated stars used to estimate
                   the average FWHM of stars in the frame.
    ellip_max      Maximum accepted ellipticity value.
    fwhm_min       Minimum accepted FWHM value.
    gain_key       Header key for the gain value.
    rdnoise_key    Header key for the noise value.
    filter_key     Header key for the filter's name.
    exp_key        Header key for the exposure time of the frame.



Extract data from your observed frames
--------------------------------------

The ``fitstats`` script is used to estimate the FWHM, sky mean, and sky standard
deviation for your observed set of standard and field frames.
Once executed, it will go through all the files defined as input 
(see :ref:`Input <secinput>` section) and automatically process them.

The steps followed by the script are:

1. Estimate the sky's mean and standard deviation values using the
   `sigma_clipped_stats`__ function.
2. Find candidate stars in the frame through the `DAOStarFinder`__ class.
   Only bright, unsaturated stars are selected.
3. Extract FWHM values for each of the stars selected in the above step,
   using IRAF's `psfmeasure`__ task. Those stars with large ellipticities or
   suspiciously small FWHMs are rejected.
4. Remove outliers with large FWHM values.
5. Obtain mean and standard deviation FWHM values for each frame processed.
6. Save date to files and plot.

The script generates the following output files (where ``xxxxx`` is the name of
the .fits file processed):

* ``xxxxx`` **.coo**: output data with x,y coordinates, `FWHM`, ellipticity,
  and relative magnitude values of the stars selected in the  .fits file.

.. parsed-literal::
    # x      y        FWHM   Ellip  Mag
    2635.46  847.5    5.076  0.02   3.23
    130.46   3820.8   4.788  0.04   1.91
    3848.14  2100.48  5.224  0.04   2.24
    3858.27  108.83   4.468  0.12   4.26
    ...

* ``xxxxx`` **.png**: output image showing the analysis performed on each
  .fits file processed.

.. image:: _figs/fitstats.png
   :width: 95%

* **fitstats.dat**: output file that contains the relevant data found after
  the analysis of either the single .fits file processed, or all the .fits files
  in the processed folder.

.. parsed-literal::
     # image           filter  exposure    Sky_mean  Sky_STDDEV  FWHM_(N_stars)  FWHM_(mean)  FWHM_(std) 
     stk_2153.fits          U      20.0        1.96        3.48              46         4.73        0.70 
     stk_2085.fits          U     250.0       19.36        5.50              14         5.33        0.11 
     stk_2151.fits          U      20.0        1.96        3.48              49         4.31        0.62 
     ....


Align your images
-----------------

xxxxxx



.. __: http://docs.astropy.org/en/stable/api/astropy.stats.sigma_clipped_stats.html
.. __: http://photutils.readthedocs.io/en/stable/api/photutils.DAOStarFinder.html
.. __: http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?psfmeasure