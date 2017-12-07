
import astropy.units as u
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt


als_f = '/home/gabriel/Github/photpy/input/field/stk_2059.als.1'
data_iraf = ascii.read(als_f, format='daophot')

# ascii.write(
#     data_iraf['XCENTER', 'YCENTER'],
#     'stk_2059.src', names=['x_0', 'y_0'], format='fixed_width',
#     delimiter=' ',
#     formats={'x_0': '%10.4f', 'y_0': '%10.4f'},
#     fill_values=[(ascii.masked, 'nan')], overwrite=True)


def calibrate_magnitudes(tab, itime=15., zmag=25.):
    tab['cal_mags'] = (zmag - 2.5 * np.log10(tab['flux_fit'] / itime)) * u.mag
    return tab


apert_corr = -0.210260869565  # 0.187956521739


def readData(phot_run):
    if phot_run == 'fix_coord':
        psf_f = '/home/gabriel/Github/photpy/output/field/filt_I/stk_2059_psf_1.out'
    else:
        psf_f = '/home/gabriel/Github/photpy/output/field/filt_I/stk_2059_psf.out'
    data_photut = ascii.read(psf_f)
    data_photut = calibrate_magnitudes(data_photut)
    data_photut.sort(['id'])

    print(len(data_iraf), len(data_photut))

    # from scipy.spatial import distance
    # c1 = np.array([data_iraf['XCENTER'], data_iraf['YCENTER']])
    # c2 = np.array([data_photut['x_fit'], data_photut['y_fit']])
    # d = distance.cdist(c1.T, c2.T)

    # min_idxs = d.argmin(axis=1)
    # data_photut = data_photut[min_idxs]

    # for i, (x, y) in enumerate(data_iraf['XCENTER', 'YCENTER']):
    #     xp, yp = data_photut['x_fit', 'y_fit'][i]
    #     d = np.sqrt((x-xp)**2 + (y-yp)**2)
    #     print(d, x, y, xp, yp)

    print("Sum of distances:")
    print(sum(
        np.sqrt((data_iraf['XCENTER'] - data_photut['x_fit'])**2 +
                (data_iraf['YCENTER'] - data_photut['y_fit'])**2)))

    med = np.nanmedian(data_iraf['MAG'] - data_photut['cal_mags'])
    print("Median: ", med)
    med_corr = med - apert_corr

    return data_photut, med_corr


data_photut_fix, med_corr_fix = readData('fix_coord')

plt.subplot(211)
plt.ylim(med_corr_fix - .5, med_corr_fix + .5)
plt.grid()
plt.title("Median (IRAF-photutils_fix)={:.4f}".format(med_corr_fix))
plt.scatter(
    data_iraf['MAG'],
    data_iraf['MAG'] - (data_photut_fix['cal_mags'] + apert_corr), s=5)
plt.axhline(y=med_corr_fix, color='r', linestyle='--', lw=.8)
# plt.xlabel("IRAF mags (zmag=25.)")
plt.ylabel("(IRAF - photutil_fix) mags")

data_photut, med_corr = readData('')

plt.subplot(212)
plt.ylim(med_corr - .5, med_corr + .5)
plt.grid()
plt.title("Median (IRAF-photutils)={:.4f}".format(med_corr))
plt.scatter(
    data_iraf['MAG'],
    data_iraf['MAG'] - (data_photut['cal_mags'] + apert_corr), s=5)
plt.axhline(y=med_corr, color='r', linestyle='--', lw=.8)
plt.xlabel("IRAF mags (zmag=25.)")
plt.ylabel("(IRAF - photutil) mags")

plt.show()
