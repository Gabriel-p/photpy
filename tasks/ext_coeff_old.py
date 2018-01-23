
import read_pars_file as rpf
from aperphot_standards import instrumMags
from hlpr import st_fwhm_select

import os
from os.path import join, isfile
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt

from astropy.io import fits, ascii
from astropy.table import Table


def in_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    in_path = pars['mypath'].replace('tasks', 'input/')
    ou_path = pars['mypath'].replace('tasks', 'output/')

    # TODO hardcoded
    # Filed should be in : '..input/field/filt_X' folders
    fits_list = [[], [], [], []]
    for i, fold in enumerate(['filt_U', 'filt_B', 'filt_V', 'filt_I']):
        for file in os.listdir(join(in_path, 'field', fold)):
            f = join(in_path, 'field', fold, file)
            if isfile(f):
                if f.endswith('.fits'):
                    fits_list[i].append(f)

    # TODO hardcoded
    # .mch files with longest exposure as reference: U B V I
    field_mch = ['stk0085.fits', 'stk1108.mch', 'stk3103.fits', 'stk1115.fits']
    mch_files = []
    for f in field_mch:
        f = join(ou_path, f)
        mch_files.append(f)

    return pars, fits_list, mch_files


def read_mch():
    ascii.read()


def main():
    """
    Rough sketch of script to determine extinction coefficients.
    """
    pars, fits_list, mch_files = in_params()
    f_key, exp_key, air_key = pars['filter_key'], pars['exposure_key'],\
        pars['airmass_key']

    for i, filt in enumerate(['U', 'B', 'V', 'I']):
        read_mch()

        ffile = 'stk_2083_crop.fits'
        print("Obtain standard stars coordinates from: {}".format(ffile))
        # Perform a Daofind on a selected .fits frame.
        hdulist = fits.open(join(in_out_path, ffile))
        hdu_data = hdulist[0].data
        psf_select = st_fwhm_select(
            float(pars['dmax']), int(pars['max_stars']),
            float(pars['thresh_fit']), float(pars['fwhm_init']), 35.,
            hdu_data)[0]
        psf_select['ID'], psf_select['x_obs'], psf_select['y_obs'] =\
            psf_select['id'], psf_select['xcentroid'], psf_select['ycentroid']
        # This table holds the IDs and coordinates of all stars detected.
        positions = Table([
            psf_select['ID'], psf_select['x_obs'], psf_select['y_obs']])

        filters = {'U': [], 'B': [], 'V': [], 'R': [], 'I': []}
        for imname in fits_list:
            f_name = imname.split('/')[-1].split('.')[0]

            # Load .fits file.
            hdulist = fits.open(imname)
            # Extract header and data.
            hdr, hdu_data = hdulist[0].header, hdulist[0].data
            filt, exp_time, airmass = hdr[f_key], hdr[exp_key], hdr[air_key]

            print("Aperture photometry on: {}".format(f_name))
            photu = instrumMags(
                f_name, positions, hdu_data, exp_time, float(pars['aperture']),
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
                print("Filter {}, median K: {:.3f}".format(filt, median_K))
                # Plot coefficients.
                plt.title("Filter {}".format(filt))
                plt.ylim(0., 2 * median_K)
                slopes = np.array(slopes)
                slopes = slopes[~np.isnan(slopes)]
                plt.scatter(range(len(slopes)), slopes)
                plt.axhline(y=median_K, color='r')
                plt.show()


if __name__ == '__main__':
    main()

# https://arxiv.org/PS_cache/arxiv/pdf/0906/0906.3014v1.pdf
# v3 = +0.16, b3 = +0.25, i3 = +0.08, u3 = +0.45

# From an old 'noche.ans' file when BO14 was processed.
# U .49, B .27, V .12, I .02
