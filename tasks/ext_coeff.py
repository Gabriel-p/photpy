
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
    f_key, exp_key, air_key = pars['filter_key'][0], pars['exposure_key'][0],\
        pars['airmass_key'][0]
    landolt_fl = read_standard_coo(in_out_path, pars['landolt_fld'][0])

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
            landolt_fl, hdu_data, exp_time, float(pars['aperture'][0]),
            float(pars['annulus_in'][0]), float(pars['annulus_out'][0]))

        stars = []
        for inst_mag in photu['cal_mags'].value:
            stars.append([inst_mag, airmass])

        # Group frames by filter.
        filters[filt].append(stars)

    import pdb; pdb.set_trace()  # breakpoint 1aa181a1 //



if __name__ == '__main__':
    main()
