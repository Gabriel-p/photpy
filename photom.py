
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from photutils import CircularAperture
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

# from matplotlib.colors import LogNorm

from photutils.detection import IRAFStarFinder
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.fitting import LinearLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from photutils.psf import DAOPhotPSFPhotometry

"""
http://docs.astropy.org/en/stable/io/fits/#
http://www.astropy.org/astropy-tutorials/FITS-images.html
"""

# Load test dataset.
# from photutils import datasets
# hdu = datasets.load_star_image()
# hdu_data = hdu.data
# Load real data.
image_file = 'stk_fcd0048.fits'  # 'stk_2061.fits'

# hdulist = fits.open(image_file)
# print(hdulist.info())
# prihdr = hdulist[0].header
# print(prihdr)
# print(prihdr.keys())
# scidata = hdulist[0].data
# hdulist.close()

# Extract data from FITS file.
hdu_data = fits.getdata(image_file)

#
# Simple PSF photometry example
bkgrms = MADStdBackgroundRMS()
std = bkgrms(hdu_data)
plt.subplot(1, 2, 1)
sigma_psf = 2.0
iraffind = IRAFStarFinder(
    threshold=3.5 * std, fwhm=sigma_psf * gaussian_sigma_to_fwhm,
    minsep_fwhm=0.01, roundhi=5.0, roundlo=-5.0,
    sharplo=0.0, sharphi=2.0)
daofind = DAOStarFinder(fwhm=3.0, threshold=5. * std)

daogroup = DAOGroup(2.0 * sigma_psf * gaussian_sigma_to_fwhm)
mmm_bkg = MMMBackground()
# fitter = LevMarLSQFitter()
psf_model = IntegratedGaussianPRF(sigma=sigma_psf)
daophot_photometry = DAOPhotPSFPhotometry(finder=iraffind,  # daofind,
                                          group_maker=daogroup,
                                          bkg_estimator=mmm_bkg,  # None,
                                          psf_model=psf_model,
                                          fitter=LevMarLSQFitter(),
                                          niters=1, fitshape=(11, 11))
print('Performing PSF photometry.')
result_tab, residual_image = daophot_photometry(image=hdu_data)

plt.subplot(1, 2, 1)
plt.imshow(hdu_data, cmap='viridis', aspect=1, interpolation='nearest',
           origin='lower', vmin=0., vmax=50)
plt.title('Observed data')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.subplot(1, 2, 2)
plt.imshow(residual_image, cmap='viridis', aspect=1, interpolation='nearest',
           origin='lower', vmin=0., vmax=50)
plt.title('Residual Image')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.show()

import pdb; pdb.set_trace()  # breakpoint 51fb6845 //


f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
# ax1.imshow(hdu_data, origin='lower', cmap='gray')

# # Coarse background subtraction.
# median, std = np.median(hdu_data), np.std(hdu_data)
# img_no_bckg = hdu_data - median
# print('1) Median +- std: ', np.median(hdu_data), np.std(hdu_data))
# ax2.imshow(img_no_bckg.data, origin='lower', cmap='gray', vmin=0,
#            vmax=median + std)

print(np.mean(hdu_data), np.median(hdu_data), np.std(hdu_data),
      np.max(hdu_data), np.min(hdu_data))

# Better background subtraction.
mean, median, std = sigma_clipped_stats(hdu_data, sigma=3.0, iters=2)
img_no_bckg = hdu_data - median
print('Median +- std: ', median, std)
ax3.imshow(img_no_bckg.data, origin='lower', cmap='Greys', vmin=0.,
           vmax=median + std)

# # Best background subtraction.
# from photutils import make_source_mask
# mask = make_source_mask(hdu_data, snr=2, npixels=5, dilate_size=11)
# mean, median, std = sigma_clipped_stats(hdu_data, sigma=3.0, mask=mask)
# img_no_bckg = hdu_data - median
# print('3) Median +- std: ', median, std)
# ax4.imshow(img_no_bckg.data, origin='lower', cmap='gray', vmin=0,
#            vmax=median + std)

# # Variable background estimation.
# from photutils.background import Background2D
# bkg = Background2D(hdu_data, (50, 50), filter_size=(3, 3),
#                    method='median')
# # Background image.
# ax3.imshow(bkg.background, origin='lower', cmap='Greys_r')
# # Background-subtracted image.
# norm = ImageNormalize(stretch=SqrtStretch())
# ax4.imshow(hdu_data - bkg.background, norm=norm, origin='lower',
#            cmap='Greys_r')

#
daofind = DAOStarFinder(fwhm=3.0, threshold=5. * std)
sources = daofind(img_no_bckg)
print(sources)
positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r=4.)
norm = ImageNormalize(stretch=SqrtStretch())
ax4.imshow(img_no_bckg.data, cmap='Greys', origin='lower',
           vmin=0., vmax=median + std)
apertures.plot(color='red', lw=1.5, alpha=0.5)

plt.tight_layout()
plt.show()

for x, y in zip(*[sources['xcentroid'], sources['ycentroid']]):
    print(x, y)
    import pdb; pdb.set_trace()  # breakpoint d845a590 //

#
#
##############################################################################
image_file = 'stk_2061.fits'

# hdulist = fits.open(image_file)

# print(hdulist.info())

# prihdr = hdulist[0].header
# # print(prihdr)
# print(prihdr.keys())

# scidata = hdulist[0].data

# hdulist.close()

image_data = fits.getdata(image_file)
print(type(image_data))
print(image_data.shape)

print('Min:', np.min(image_data))
print('Max:', np.max(image_data))
print('Mean:', np.mean(image_data))
print('Stdev:', np.std(image_data))

# # plt.imshow(image_data, cmap='gray', norm=LogNorm())
plt.imshow(image_data, cmap='gray', vmin=0, vmax=30)
plt.colorbar()

# img_arr = np.array(image_data.flat)
# print(len(img_arr))
# img_filt = img_arr[(img_arr >= 0.) & (img_arr <= 100.)]
# print(len(img_filt))
# histogram = plt.hist(img_filt, 50)

plt.show()
