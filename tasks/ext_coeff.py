
import read_pars_file as rpf
from aperphot_standards import read_standard_coo, instrumMags

import os
from os.path import join, isfile
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from astropy.io import fits


def in_params(in_fold, ldlt_fold):
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    in_path = pars['mypath'].replace('tasks', 'input/' + in_fold)
    stnd_path = pars['mypath'].replace('tasks', ldlt_fold)

    fits_list = []
    for file in os.listdir(in_path):
        f = join(in_path, file)
        if isfile(f):
            if f.endswith('.fits'):
                fits_list.append(f)

    return pars, in_path, stnd_path, fits_list


def main(in_fold, ldlt_fold, stnd_fl):
    """
    Rough sketch of script to determine extinction coefficients.
    """
    pars, in_path, stnd_path, fits_list = in_params(in_fold, ldlt_fold)
    f_key, exp_key, air_key = pars['filter_key'], pars['exposure_key'],\
        pars['airmass_key']

    print("Obtain standard stars coordinates from: {}".format(stnd_fl))
    # Read data for this Landolt field.
    landolt_fl = read_standard_coo(stnd_path, stnd_fl)

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
            f_name, landolt_fl, hdu_data, exp_time, float(pars['aperture']),
            float(pars['annulus_in']), float(pars['annulus_out']))

        stars = []
        for i, inst_mag in enumerate(photu['cal_mags'].value):
            print(filt, landolt_fl['ID'][i], airmass, inst_mag)
            stars.append([airmass, inst_mag, landolt_fl['ID'][i]])

        # Group frames by filter.
        filters[filt].append(stars)

    for filt, fdata in filters.iteritems():
        if fdata:
            print("Filter {}".format(filt))

            # airm, mag = [], []
            # for f in fdata:
            #     airm += [_[0] for _ in f]
            #     mag += [_[1] for _ in f]

            # m = linregress(airm, mag)[0]
            # print(" K: {:.3f}".format(m))

            # if not np.isnan(m):
            #     plt.style.use('seaborn-darkgrid')
            #     plt.title("Filter {}".format(filt))
            #     plt.scatter(airm, mag)
            #     plt.show()

            # Transform list into into:
            # airm_mag = [st1, st2, st3, ...]
            # stX = [airmass, instrum_magnitudes]
            airm_mag = [[y for y in zip(*x)] for x in zip(*fdata)]

            plt.style.use('seaborn-darkgrid')
            slopes = []
            for st in airm_mag:
                m, c, r_value, p_value, std_err = linregress(st[0], st[1])
                slopes.append(m)
                print("  Star {}, K={:.3f}".format(st[2][0], m))
                plt.scatter(st[0], st[1])
                plt.show()

            # Obtain median of extinction coefficients for this filter.
            median_K = np.nanmedian(slopes)
            print("Filter {}, median K: {:.3f}".format(filt, median_K))
            # Plot coefficients.
            if not np.isnan(median_K):
                plt.title("Filter {}".format(filt))
                plt.ylim(0., 2 * median_K)
                slopes = np.array(slopes)
                slopes = slopes[~np.isnan(slopes)]
                plt.scatter(range(len(slopes)), slopes)
                plt.axhline(y=median_K, color='r')
                plt.show()


if __name__ == '__main__':
    """
    To run this script we need *at least* two observations for each filter
    with different airmass values.
    """

    # Name of the Landolt frame.
    stnd_fl = 'SA95'
    # Path of folder where .fits files are stored.
    in_path = 'standards_95'
    # Path to the _match.coo file for this frame.
    ldlt_fold = 'output'
    main(in_path, ldlt_fold, stnd_fl)

# https://arxiv.org/PS_cache/arxiv/pdf/0906/0906.3014v1.pdf
# v3 = +0.16, b3 = +0.25, i3 = +0.08, u3 = +0.45

# From an old 'noche.ans' file when BO14 was processed.
# U .49, B .27, V .12, I .02
