
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from photutils import datasets
# import photutils
# from matplotlib.colors import LogNorm

"""
http://docs.astropy.org/en/stable/io/fits/#
http://www.astropy.org/astropy-tutorials/FITS-images.html
"""

# Load test dataset.
# hdu = datasets.load_star_image()
# hdu_data = hdu.data
image_file = 'stk_2061.fits'
hdu_data = fits.getdata(image_file)

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
ax1.imshow(hdu_data, origin='lower', cmap='gray', vmin=0, vmax=30)

# Coarse background subtraction.
# image = hdu_data[500:700, 500:700].astype(float)
img_no_bckg = hdu_data - np.median(hdu_data)
ax2.imshow(img_no_bckg.data, origin='lower', cmap='gray', vmin=0, vmax=30)

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
