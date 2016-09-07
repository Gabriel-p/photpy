
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from photutils import datasets
from astropy.stats import sigma_clipped_stats
from photutils import make_source_mask
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.background import Background2D
# from matplotlib.colors import LogNorm

"""
http://docs.astropy.org/en/stable/io/fits/#
http://www.astropy.org/astropy-tutorials/FITS-images.html
"""

# Load test dataset.
# hdu = datasets.load_star_image()
# hdu_data = hdu.data
# Load real data.
image_file = 'stk_2135.fits'  # 'stk_fcd0048.fits' # 'stk_2061.fits'

# hdulist = fits.open(image_file)

hdulist = fits.open(image_file)
print(hdulist.info())
prihdr = hdulist[0].header
print(prihdr)
print(prihdr.keys())
scidata = hdulist[0].data
import pdb; pdb.set_trace()  # breakpoint 7cc0842b //
hdulist.close()

# Extract data from FITS file.
hdu_data = fits.getdata(image_file)

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
ax1.imshow(hdu_data, origin='lower', cmap='gray')

# Coarse background subtraction.
# image = hdu_data[500:700, 500:700].astype(float)
median, std = np.median(hdu_data), np.std(hdu_data)
img_no_bckg = hdu_data - median
print('Median +- std: ', np.median(hdu_data), np.std(hdu_data))
ax2.imshow(img_no_bckg.data, origin='lower', cmap='gray', vmin=0,
           vmax=median + std)

# Better background subtraction.
mean, median, std = sigma_clipped_stats(hdu_data, sigma=3.0, iters=5)
img_no_bckg = hdu_data - median
print('Median +- std: ', median, std)
ax3.imshow(img_no_bckg.data, origin='lower', cmap='gray', vmin=0,
           vmax=median + std)

# Best background subtraction.
mask = make_source_mask(hdu_data, snr=2, npixels=5, dilate_size=11)
mean, median, std = sigma_clipped_stats(hdu_data, sigma=3.0, mask=mask)
img_no_bckg = hdu_data - median
print('Median +- std: ', median, std)
ax4.imshow(img_no_bckg.data, origin='lower', cmap='gray', vmin=0,
           vmax=median + std)

# # Variable background estimation.
# bkg = Background2D(hdu_data, (50, 50), filter_size=(3, 3),
#                    method='median')
# # Background image.
# ax3.imshow(bkg.background, origin='lower', cmap='Greys_r')
# # Background-subtracted image.
# norm = ImageNormalize(stretch=SqrtStretch())
# ax4.imshow(hdu_data - bkg.background, norm=norm, origin='lower',
#            cmap='Greys_r')

plt.tight_layout()
plt.show()
import pdb; pdb.set_trace()  # breakpoint 3e914842 //

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
