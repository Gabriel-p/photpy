
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

.. todo::
   Not finished.

**We assume that your standard frames are already properly aligned.**
This can be done using the ``align_crop`` script. This step is necessary, since
the coordinates of the standard stars will be transformed to the observed system
using a single observed frame.

.. note::
  If you have more than one exposure per filter for your standard frame, at this
  point **you need to select only one**. An ideal frame should allow the
  detection of all the standard stars in it. This means no over-exposed
  saturated stars, and no under-exposed undetectable stars.

The ``id_standard`` script allows the automatic identification of standard
stars in your observed frame, requiring minimal information.


Aperture photometry on standard stars
-------------------------------------

.. todo::
   Not finished.

**We assume that you will work with a single aperture radius value for the
standards from all the nights for all the filters.**

The aperture value should be large enough to contain as much light from
your observed standards as possible, but at the same time small enough to
minimize contamination from bad pixels and other stars.
The default convention is to use an aperture radius that is ``~4.5*FWHM`` of a
stellar image. This is, for a ``FWHM=3 px`` you'll use an aperture around
``14-15 px``.

Assuming you've already run the ``getdata`` script, these values are stored in
the ``fwhm_final.dat`` file for each of the observed filters. You can quickly
obtain an average estimate by running the ``data_avrg`` script on the parent
folder that contains the filter sub-folders.



Define and solve the transformation equations
---------------------------------------------

.. todo::
   Not finished.

**We assume you are using** `Landolt standards`_. These are the only standards
supported so far by the code.

The previous step generated instrumental magnitudes for all the observed
standard stars, for each selected frame. The transformation equations are used
to put these magnitudes on the standard system.




.. _Landolt standards: http://www.eso.org/sci/observing/tools/standards/Landolt.html