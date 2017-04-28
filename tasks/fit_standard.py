
import os
from os.path import join, realpath, dirname, isfile
import sys
import numpy as np
from scipy.stats import linregress
import gc
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
from matplotlib.ticker import NullFormatter

import astropy.units as u
from astropy.io import ascii, fits
from astropy.table import Table

from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry


def read_params():
    """
    Read parameter values.
    """
    mypath = realpath(join(os.getcwd(), dirname(__file__)))
    pars_f = join(mypath.replace('tasks', ''), 'params_input.dat')
    if not isfile(pars_f):
        print("Parameters file is missing. Exit.")
        sys.exit()

    pars = {}
    with open(pars_f, 'r') as f:
        for line in f:
            if not line.startswith('#') and line != '\n':
                key, value = line.replace('\n', '').split()
                pars[key] = value

    in_out_path = mypath.replace('tasks', 'output/standards')

    fits_list = []
    if os.path.isdir(in_out_path):
        for file in os.listdir(in_out_path):
            f = join(in_out_path, file)
            if isfile(f):
                if f.endswith('_crop.fits'):
                    fits_list.append(f)

    if not fits_list:
        print("No '*_crop.fits' files found in 'output/standards' folder."
              " Exit.")
        sys.exit()

    return mypath, pars, fits_list, in_out_path


def read_standard_coo(in_out_path, landolt_fld):
    """
    Read .coo file with standard stars coordinates in the system of the
    observed frame.
    """
    f = join(in_out_path, landolt_fld + '_obs.coo')
    landolt_fl = ascii.read(f)

    return landolt_fl


def landoltZeroFilters(landolt_fl, filt, airmass, extin_coeffs):
    """
    Separate into individual filters, and correct for airmass (i.e.,
    instrumental. magnitude at zero airmass)
    """
    fe = [_.split(',') for _ in extin_coeffs.split(';')]
    f_idx = zip(*fe)[0].index(filt)
    obs_f, ext = zip(*fe)[0][f_idx], map(float, zip(*fe)[1])[f_idx]

    if obs_f == 'U':
        inst_mag = landolt_fl['UB'] + landolt_fl['BV'] + landolt_fl['V']
    elif obs_f == 'B':
        inst_mag = landolt_fl['BV'] + landolt_fl['V']
    elif obs_f == 'V':
        inst_mag = landolt_fl['V']
    elif obs_f == 'R':
        inst_mag = landolt_fl['V'] - landolt_fl['VR']
    elif obs_f == 'I':
        inst_mag = landolt_fl['V'] - landolt_fl['VI']

    zero_airmass_f = inst_mag - ext * airmass

    return zero_airmass_f


def calibrate_magnitudes(tab, itime=1., zmag=25.):
    tab['cal_mags'] = (zmag - 2.5 * np.log10(tab['flux_fit'] / itime)) * u.mag
    return tab


def photutils_phot(landolt_fl, hdu_data, exp_time, pars):
    """
    Perform aperture photometry for all 'landolt_fl' standard stars observed in
    the  'hdu_data' file.
    """
    aper_rad, annulus_in, annulus_out = float(pars['aperture']),\
        float(pars['annulus_in']), float(pars['annulus_out'])

    # Coordinates from observed frame.
    positions = zip(*[landolt_fl['x_obs'], landolt_fl['y_obs']])
    apertures = CircularAperture(positions, r=aper_rad)
    annulus_apertures = CircularAnnulus(
        positions, r_in=annulus_in, r_out=annulus_out)

    apers = [apertures, annulus_apertures]
    phot_table = aperture_photometry(hdu_data, apers)
    bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area()
    bkg_sum = bkg_mean * apertures.area()
    phot_table['flux_fit'] = phot_table['aperture_sum_0'] - bkg_sum
    phot_table = calibrate_magnitudes(phot_table, itime=exp_time)

    return phot_table


def main():
    """
    """
    mypath, pars, fits_list, in_out_path = read_params()

    landolt_fl = read_standard_coo(in_out_path, pars['landolt_fld'])

    # For each _crop.fits observed (aligned and cropped) standard file.
    for imname in fits_list:
        print("Fits file: {}".format(imname))

        # Load .fits file.
        hdulist = fits.open(imname)
        # Extract header and data.
        hdr, hdu_data = hdulist[0].header, hdulist[0].data
        filt, exp_time, airmass = hdr[pars['filter_key']],\
            hdr[pars['exposure_key']], hdr[pars['airmass_key']]

        zero_airmass_f = landoltZeroFilters(
            landolt_fl, filt, airmass, pars['extin_coeffs'])

        photu_data = photutils_phot(landolt_fl, hdu_data, exp_time, pars)
        print(photu_data['cal_mags'])

        m, c, r_value, p_value, std_err = linregress(
            zero_airmass_f, photu_data['cal_mags'])
        print(r_value**2, std_err)
        # regression line
        line = m * zero_airmass_f + c
        print(line)
        plt.plot(
            zero_airmass_f, line, 'r-', zero_airmass_f,
            photu_data['cal_mags'], 'o')
        # plt.scatter(zero_airmass_f, photu_data['cal_mags'])
        plt.show()

    print("\nFinished.")


if __name__ == '__main__':
    main()
