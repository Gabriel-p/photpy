
Standard stars photometry and reduction
=======================================

.. todo::
   Not finished.

Obtaining aperture photometry of your standards
-----------------------------------------------

Standard stars provide a good example of relatively uncrowded photometry, and in
this section we will describe hot to obtain instrumental for your standards.
The basic steps are

- `Picking an aperture size`_: decide what aperture size you whish to use for
  measuring your standards (**this should be the same for all the frames**). At
  the same time pick a sky annulus.
- `Setting things up`_: set up various parameters (**fwhm, threshold, sky
  annulus, aperture radius**) to have correct values.
- `Doing it: aperture photometry at last`_: for each frame

   1. Identify the standard star(s) either interactively using a cursor or by
      using an automatic star finding algorithm.

   2. Obtain aperture photometry for each of your standard stars.


Picking an aperture size
------------------------

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


Setting things up
-----------------

.. todo::
   Not finished.


Doing it: aperture photometry at last
-------------------------------------

.. todo::
   Not finished.

