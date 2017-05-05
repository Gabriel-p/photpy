
import read_pars_file as rpf
from fit_standard import read_standard_coo, instrumMags

import os
from os.path import join, isfile
import sys
import numpy as np
from scipy.stats import linregress
import gc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText

import astropy.units as u
from astropy.io import ascii, fits
from astropy.table import Table, Column

from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry


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
    """
    pars, fits_list, in_out_path = in_params()
    f_key, exp_key, air_key = pars['filter_key'], pars['exposure_key'],\
        pars['airmass_key']
    landolt_fl = read_standard_coo(in_out_path, pars['landolt_fld'])

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
            cc = [zip(*_) for _ in fdata]
            airm_mag = [[[], [], []] for _ in fdata[0]]
            for a in cc:
                aa = zip(*a)
                for i, st in enumerate(aa):
                    airm_mag[i][0].append(st[0])
                    airm_mag[i][1].append(st[1])
                    airm_mag[i][2].append(st[2])

            yl = []
            for st in airm_mag:
                plt.scatter(st[0], st[1])
                m, c, r_value, p_value, std_err = linregress(st[0], st[1])
                print("Star {}, m={:.3f}".format(st[2][0], m))
                yl.append(max(st[1]))
                yl.append(min(st[1]))
            plt.ylim(max(yl) + .05, min(yl) - .05)
            plt.show()

    import pdb; pdb.set_trace()  # breakpoint 1aa181a1 //


if __name__ == '__main__':
    main()

# https://arxiv.org/PS_cache/arxiv/pdf/0906/0906.3014v1.pdf
