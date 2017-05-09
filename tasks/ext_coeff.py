
import read_pars_file as rpf
from fit_standard import instrumMags
from hlpr import st_fwhm_select

import os
from os.path import join, isfile
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table


def in_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    in_out_path = pars['mypath'].replace('tasks', 'output/standards')

    fits_list = []
    for file in os.listdir(in_out_path):
        f = join(in_out_path, file)
        if isfile(f):
            if f.endswith('.fits'):
                fits_list.append(f)

    return pars, fits_list, in_out_path


def main():
    """
    Rough sketch of script to determine extinction coefficients.
    """
    pars, fits_list, in_out_path = in_params()
    f_key, exp_key, air_key = pars['filter_key'], pars['exposure_key'],\
        pars['airmass_key']

    # Perform a Daofind on a selected .fits frame.
    hdulist = fits.open(join(in_out_path, 'stk_2083_crop.fits'))
    hdu_data = hdulist[0].data
    psf_select = st_fwhm_select(
        float(pars['dmax']), int(pars['max_stars']),
        float(pars['thresh_level']), float(pars['fwhm_init']), 35.,
        hdu_data)[0]
    psf_select['ID'], psf_select['x_obs'], psf_select['y_obs'] =\
        psf_select['id'], psf_select['xcentroid'], psf_select['ycentroid']
    # This table holds the IDs and coordinates of all stars detected.
    landolt_fl = Table([
        psf_select['ID'], psf_select['x_obs'], psf_select['y_obs']])

    filters = {'U': [], 'B': [], 'V': [], 'R': [], 'I': []}
    for imname in fits_list:
        fname = imname.replace(pars['mypath'].replace('tasks', 'output'), '')

        # Load .fits file.
        hdulist = fits.open(imname)
        # Extract header and data.
        hdr, hdu_data = hdulist[0].header, hdulist[0].data
        filt, exp_time, airmass = hdr[f_key], hdr[exp_key], hdr[air_key]

        print("Aperture photometry on: {}".format(fname))
        photu = instrumMags(
            landolt_fl, hdu_data, exp_time, float(pars['aperture']),
            float(pars['annulus_in']), float(pars['annulus_out']))

        stars = []
        for i, inst_mag in enumerate(photu['cal_mags'].value):
            stars.append([airmass, inst_mag, landolt_fl['ID'][i]])

        # Group frames by filter.
        filters[filt].append(stars)

    for filt, fdata in filters.iteritems():
        if fdata:
            print("Filter {}".format(filt))
            # Transform list into into:
            # airm_mag = [st1, st2, st3, ...]
            # stX = [airmass, instrum_magnitudes]
            airm_mag = [[y for y in zip(*x)] for x in zip(*fdata)]

            plt.style.use('seaborn-darkgrid')
            slopes = []
            for st in airm_mag:
                # plt.scatter(st[0], st[1])
                m, c, r_value, p_value, std_err = linregress(st[0], st[1])
                slopes.append(m)
                print("  Star {}, K={:.3f}".format(st[2][0], m))

            # Obtain median of extinction coefficients for this filter.
            median_K = np.nanmedian(slopes)
            print("Median K: {:.3f}".format(median_K))
            # Plot histogram of coefficients.
            plt.xlim(0., 2 * median_K)
            slopes = np.array(slopes)
            plt.hist(slopes[~np.isnan(slopes)], bins=50)
            plt.show()


if __name__ == '__main__':
    main()

# https://arxiv.org/PS_cache/arxiv/pdf/0906/0906.3014v1.pdf

# v3 = +0.16, b3 = +0.25, i3 = +0.08, u3 = +0.45
