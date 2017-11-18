
Standard stars photometry
=========================

.. todo::
   Not finished.

This section describes how to obtain instrumental magnitudes for your standards,
and set up the transformation from instrumental to standard. The steps are

* `Identify standard stars`_: identify the coordinates of the stars in the
  observed standard field.
* `Aperture photometry on standard stars`_: obtain aperture photometry for
  each of the selected standard star frames.
* `Define and solve the transformation equations`_: define and fit the
  instrumental to standard transformation equations.


Identify standard stars
-----------------------

* ``ref_id_std``:       stk_2078_crop.fits
* ``landolt_fld``:      pg1323
* ``scale_min``:        0.1
* ``scale_max``:        10.0
* ``rot_min``:          0.0
* ``rot_max``:          5.0
* ``do_plots_C``:       y


The ``id_standard`` script allows the automatic identification of standard
stars in your observed frame, requiring only minimal information.

.. warning::

   This script assumes three things: that you are using `Landolt standards`_,
   that **all** your standard frames are aligned, and that there are **at
   least** three standard stars in your observed Landolt frame.

Landolt equatorial standards are the only set supported currently by the code.
The alignment of your observed standard frames is required since the coordinates
of the standard stars will be transformed to the observed system using a single
observed frame.
The requirement of three standard stars minimum is related to the algorithm
used by the script to perform the transformation to the observed coordinates
system.

The user selects a reference (aligned) observed standards frame, and informs the
script to which Landolt frame it corresponds. The script will then search for
the proper scaling, rotation, and ``x, y`` translation that allows it to
identify which stars in your observed frame are standards.

An image is produced marking the stars in the observed frame that where
identified as standards.

.. image:: _figs/std_id.png
   :width: 95%

An ``Lframe_obs.coo`` file is also generated (where ``Lframe`` is the name of
the Landolt frame processed, given in the ``landolt_fld`` keyword onf the input
parameters file) containing the Landolt stars information, as well as their
``x, y`` coordinates in the observed system:

.. parsed-literal::
         ID      x      y       V     BV      UB      VR      RI      VI      x_obs      y_obs 
  pg1323-86  211.0  158.3  13.481  -0.14  -0.681  -0.048  -0.078  -0.127   1872.768   1847.208 
 pg1323-86A  162.5  137.5  13.591  0.393  -0.019   0.252   0.252   0.506   1521.456   1695.408 
 pg1323-86B  158.1  128.0  13.406  0.761   0.265   0.426   0.407   0.833   1489.735   1626.427


Aperture photometry on standard stars
-------------------------------------

.. todo::
   Not finished.

* ``stnd_obs_fields``:  pg1323, stk_2082_crop.fits, stk_2127_crop.fits, stk_2129_crop.fits, stk_2131_crop.fits
* ``aperture``:         15
* ``annulus_in``:       20
* ``annulus_out``:      25
* ``do_plots_D``:       y


The ``aperphot_standard`` automatically performs aperture photometry on your
observed standard frames, for the identified stars in one.

.. warning::
  We assume that you will work with a **single** aperture radius value for the
  standards from all the nights for all the filters.

The selected aperture value should be large enough to contain as much light from
your observed standards as possible, but at the same time small enough to
minimize contamination from bad pixels and other stars.
The default convention is to use an aperture radius that is ``~4.5*FWHM`` of a
stellar image. This is, for a ``FWHM=3 px`` you'll use an aperture around
``14-15 px``.

Assuming you've already run the ``fitstats`` script, the median FWHM values for
your standard frames are stored in the ``fitstats.dat`` file. 

.. warning::
  If you have more than one exposure per filter for your standard frame, at this
  point you need to select only one. An ideal frame should allow the detection
  of all the standard stars in it. This means no over-exposed saturated stars,
  and no under-exposed undetectable stars.

This script reads as many standard .fits files as you input in the
``stnd_obs_fields`` keyword in the ``params_input.dat`` file. After that, the
corresponding the ``Lframe_obs.coo`` file for this Landolt frame is read (the
script ``id_standards`` generates this file). From this file we read the
``x, y`` coordinates for each standard star in the observed frame, along with
its calibrated (Landolt) photometry.

Circular aperture photometry is performed on each frame for each standard star.
Their calculated instrumental magnitudes are corrected for zero airmass.

Final zero airmass magnitudes are stored in the ``stnd_aperphot.dat`` file,
along with Landolt default colors and magnitudes for  each standard star in
each observed standard field:

.. parsed-literal::
  Filt   Stnd_field    ID                 file   exp_t       A       ZA_mag      Col      Mag  
     I       pg1323    86   stk_2082_crop.fits    20.0   1.071       15.358   -0.127   13.608  
     I       pg1323   86A   stk_2082_crop.fits    20.0   1.071       14.682    0.506   13.085  
     I       pg1323   86B   stk_2082_crop.fits    20.0   1.071       14.269    0.833   12.573  
     I       pg1323   86C   stk_2082_crop.fits    20.0   1.071       14.956    0.759   13.244  



Define and solve the transformation equations
---------------------------------------------

.. todo::
   Not finished.

* ``R^2_min``:          0.98
* ``RMSE_max``:         0.05
* ``extin_coeffs``:     U .49, B .27, V .12, I .02
* ``do_plots_E``:       y

.. warning::
  The extinction coefficients for your observed filters are assumed to be known.

.. warning::
  The V filter is assumed to be present among your observed filters.

This script will obtain instrumental magnitudes for all the observed
standard stars, for each selected frame. The transformation equations are used
to put these magnitudes on the standard system.

.. _Landolt standards: http://www.eso.org/sci/observing/tools/standards/Landolt.html