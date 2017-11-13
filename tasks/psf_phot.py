
import read_pars_file as rpf

import os
from os.path import join, isfile
import sys
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import ascii, fits
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from photutils.psf import IterativelySubtractedPSFPhotometry


def in_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    out_path = pars['mypath'].replace('tasks', 'output/field')

    fits_list = []
    for file in os.listdir(in_out_path):
        f = join(in_out_path, file)
        if isfile(f):
            if f.endswith('_crop.fits'):
                fits_list.append(f)

    if not fits_list:
        print("No '*_crop.fits' files found in 'output/standards' folder."
              " Exit.")
        sys.exit()

    return pars, out_path


def readData(pars, in_out_path, imname):
    """
    """
    # Load .fits file.
    hdulist = fits.open(join(in_out_path, imname))
    # Extract header and data.
    hdr, hdu_data = hdulist[0].header, hdulist[0].data
    filt, exp_time, airmass = hdr[pars['filter_key']],\
        hdr[pars['exposure_key']], hdr[pars['airmass_key']]
    print("  Filter {}, Exp time {}, Airmass {}".format(
        filt, exp_time, airmass))

    return hdu_data


def main(hdu_data, fwhm, xy_cents):
    """
    """
    daogroup = DAOGroup(2.0 * sigma_psf * gaussian_sigma_to_fwhm)
    mmm_bkg = MMMBackground()
    # fitter = LevMarLSQFitter()
    psf_model = IntegratedGaussianPRF(sigma=sigma_psf)

    photometry = IterativelySubtractedPSFPhotometry(
        group_maker=daogroup, bkg_estimator=mmm_bkg,
        psf_model=psf_model, fitter=LevMarLSQFitter(), niters=1,
        fitshape=(11, 11))

    result_tab = photometry(image=hdu_data, positions=xy_cents)
    residual_image = photometry.get_residual_image()

    plt.subplot(1, 2, 1)
    plt.imshow(hdu_data, cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower')
    plt.title('Simulated data')
    plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
    plt.subplot(1, 2, 2)
    plt.imshow(residual_image, cmap='viridis', aspect=1,
               interpolation='nearest', origin='lower')
    plt.title('Residual Image')
    plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
    plt.show()


if __name__ == '__main__':
    pars, in_out_path = in_params()
    hdu_data = readData()
    sigma_psf = -1.
    main()

# PSF Photometry
# https://photutils.readthedocs.io/en/stable/psf.html
# PSF Photometry in Crowded Fields with Photutils
# https://github.com/astropy/photutils-datasets/blob/master/notebooks/ArtificialCrowdedFieldPSFPhotometry.ipynb
# Uncertainties
# issue:
# https://github.com/pyDANDIA/pyDANDIA/issues/10
# example:
# https://gist.github.com/larrybradley/e23d905cd69ca44d032541f69c2b5433