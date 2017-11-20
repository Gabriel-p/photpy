
import os
import sys
from os.path import join, realpath, dirname
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
from astropy.stats import gaussian_sigma_to_fwhm

from photutils.utils import cutout_footprint
from photutils.background import MADStdBackgroundRMS
from photutils.psf import IntegratedGaussianPRF
from photutils.psf import DAOPhotPSFPhotometry
from photutils import DAOStarFinder
from photutils import CircularAperture

import imexam
# from imexam.imexamine import Imexamine
# from imexam.math_helper import gfwhm


os.environ['XPA_METHOD'] = "local"

# Load data.
mypath = realpath(join(os.getcwd(), dirname(__file__)))
image_file = mypath + '/input/standards/filt_B/stk_2080.fits'

# viewer = imexam.connect(viewer='ginga')
viewer = imexam.connect()
viewer.load_fits(image_file)
viewer.imexam()
# viewer.close()

# Load .fits file
hdulist = fits.open(image_file)
# Image data.
hdu_data = hdulist[0].data
# Extract header.
hdr = hdulist[0].header
hdulist.close()

import pdb; pdb.set_trace()  # breakpoint 674d8c4a //


# Crop image
crop = cutout_footprint(hdu_data, (2100, 1800), (500, 1100))
hdu_crop = crop[0]

bkgrms = MADStdBackgroundRMS()
std = bkgrms(hdu_crop)
thresh = 50. * std
sigma_psf = 5.
fwhm_sigma = sigma_psf * gaussian_sigma_to_fwhm
fitshape = int(3 * np.ceil(fwhm_sigma) // 2 * 2 + 1)

# DAOStarFinder
stfind = DAOStarFinder(threshold=thresh, fwhm=fwhm_sigma)
sources = stfind(hdu_crop)
print(sources)

positions = (sources['xcentroid'], sources['ycentroid'])
# Make a list of locations as a tuple of (x, y, ID)
starlist = list()
for point in zip(*positions):
    starlist.append((point[0], point[1]))

# get an object up with your data attached
plots = Imexamine()
plots.set_data(hdu_crop)

# make a function to calculate the fwhm
# g = lambda x: x * np.sqrt(8.0 * np.log(2.))

# fwhm = gfwhm(stddev)[0]  #to get just the x fwhm

# sys.stdout = open(os.devnull, "w")
# plots.line_fit_pars["center"][0] = False

results = []
for star in starlist:
    x, y = star
    gauss_x = plots.line_fit(x, y, genplot=False)
    gauss_y = plots.column_fit(x, y, genplot=False)
    # print(g(gauss.stddev), gfwhm(gauss.stddev)[0])
    results.append([x, y, gauss_x.stddev, gauss_x.mean, gauss_x.amplitude,
                    gfwhm(gauss_x.stddev)[0], gfwhm(gauss_y.stddev)[0]])
# sys.stdout = sys.__stdout__

for _ in results:
    print(_)

# # DAOPhotPSFPhotometry
# psf_model = IntegratedGaussianPRF(sigma=sigma_psf)
# photometry = DAOPhotPSFPhotometry(
#     crit_separation=10. * fwhm_sigma, threshold=thresh, fwhm=fwhm_sigma,
#     psf_model=psf_model, fitshape=fitshape, niters=1)
# result_tab = photometry(image=hdu_crop)
# print(result_tab)

median, std = np.median(hdu_crop), np.std(hdu_crop)

plt.subplot(1, 2, 1)
plt.title('DAOStarFinder')
plt.imshow(hdu_crop, cmap='viridis', aspect=1, interpolation='nearest',
           origin='lower', vmin=0., vmax=median + std)
positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r=4.)
apertures.plot(color='red', lw=1.5)
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)

# plt.subplot(1, 2, 2)
# plt.title('DAOPhotPSFPhotometry')
# plt.imshow(hdu_crop, cmap='viridis', aspect=1, interpolation='nearest',
#            origin='lower', vmin=0., vmax=median + std)
# positions = (result_tab['x_fit'], result_tab['y_fit'])
# apertures = CircularAperture(positions, r=4.)
# apertures.plot(color='red', lw=1.5)
# plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)

plt.show()
