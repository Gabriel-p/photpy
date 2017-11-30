
Introduction
============

.. todo::
   Not finished.

Requirements
------------

- matplotlib, numpy, astropy, photutils
- PyRAF, IRAF (temporary until photutils allows replacing the ``psfmeasure``
  tasks still used from IRAF)
- `imexam`_, `ginga`_ (optional)


Installing
----------

Download the `Anaconda`_ distribution, or its smaller `Miniconda`_ alternative
(recommended). Follow the steps given in those links to install.
Once the `conda`_ package manager is installed, add the `Astroconda`_ channel
to be able to easily install PyRAF + IRAF:


.. code-block:: bash

  $ conda config --add channels http://ssb.stsci.edu/astroconda

Create a dedicated environment with all the necessary dependences:

.. code-block:: bash

  $ conda create -n iraf27 python=2.7 pyraf iraf matplotlib numpy astropy photutils

Access the environment with:

.. code-block:: bash

  $ source activate iraf27

Generate a ``login.cl`` file and ``uparm/`` directory:

.. code-block:: bash

  $ mkiraf

and select ``xgterm`` as the terminal type when prompted.

The ``pyraf`` command should now correctly load an IRAF session.


Useful links
............

-  `PSF Photometry in Crowded Fields with Photutils`_
-  `Photometry overview`_
-  `CCD Gain Lab: The Theory`_
-  `Using Python for Astronomical Data Analysis in the Era of JWST`_


.. _Anaconda: https://www.continuum.io/downloads
.. _Miniconda: https://conda.io/miniconda.html
.. _conda: https://conda.io/docs/
.. _Astroconda: http://astroconda.readthedocs.io/en/latest/index.html


.. _imexam: http://imexam.readthedocs.io/en/latest/index.html
.. _ginga: http://ejeschke.github.io/ginga/
.. _IRAF + DS9: http://www.astronomy.ohio-state.edu/~khan/iraf/%20iraf_step_by_step_installation_64bit
.. _PSF Photometry in Crowded Fields with Photutils: https://github.com/astropy/photutils-datasets/blob/master/notebooks/ArtificialCrowdedFieldPSFPhotometry.ipynb
.. _Photometry overview: http://telvsn.fcaglp.unlp.edu.ar/normativas/charlas/%20seminario_baume.pdf
.. _`CCD Gain Lab: The Theory`: http://www.astro.umd.edu/~veilleux/ASTR310/%20fall06/ccd_theory.pdf
.. _Using Python for Astronomical Data Analysis in the Era of JWST: http://www.astrobetter.com/blog/2016/09/26/using-python-for-astronomical-data-analysis-in-the-era-of-jwst/