
import read_pars_file as rpf

import os
from os.path import join, isfile
import sys
import numpy as np

import astropy.units as u
from astropy.io import ascii, fits
from astropy.table import Table, vstack

from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry


def in_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    in_path = join(pars['mypath'].replace('tasks', 'input'))
    out_path = in_path.replace('input', 'output')

    landolt_fld, match_fldr = [], []
    for line in pars['stnd_obs_fields']:
        landolt_fld.append(line[0])
        match_fldr.append(line[1])
    pars['landolt_fld'], pars['match_fldr'] = landolt_fld, match_fldr

    # Generate list of fits files for each input folder.
    fits_list = []
    for folder in pars['match_fldr']:
        folder = folder[1:] if folder.startswith('/') else folder
        f_path = join(pars['mypath'].replace('tasks', 'input'), folder)

        list_temp = []
        if os.path.isdir(f_path):
            for file in os.listdir(f_path):
                f = join(f_path, file)
                if isfile(f):
                    if f.endswith('.fits'):
                        list_temp.append(f)
        if not list_temp:
            print("{}\n No .fits files found in match folder. Exit.".format(
                f_path))
            sys.exit()

        # Store list for this folder.
        print("Files found in '{}' folder:".format(folder))
        for fit in list_temp:
            print(" * {}".format(fit.replace(f_path, '')[1:]))
        fits_list.append(list_temp)

    return pars, fits_list, in_path, out_path


def read_standard_coo(in_path, landolt_fld):
    """
    Read _match.coo file created by 'match', with calibrated photometric
    data on Landolt standard stars, and their coordinates in the system of
    the observed frame.
    """
    f = join(in_path, landolt_fld + '.mch')
    # Change 'nan' values for '-999.9' for sources not detected.
    landolt_fl = ascii.read(f, fill_values=[('nan', '0', 'x_obs', 'y_obs')])
    landolt_fl['x_obs'].fill_value = -999.9
    landolt_fl['y_obs'].fill_value = -999.9
    landolt_fl = landolt_fl.filled()

    return landolt_fl


def standardMagCol(landolt_fl, filt):
    """
    Separate into individual filters. Return single filter, and color term
    used in the fitting process for that filter.
    """
    if filt == 'U':
        standard_mag = landolt_fl['UB'] + landolt_fl['BV'] + landolt_fl['V']
        standard_col = landolt_fl['UB']
    elif filt == 'B':
        standard_mag = landolt_fl['BV'] + landolt_fl['V']
        standard_col = landolt_fl['BV']
    elif filt == 'V':
        standard_mag = landolt_fl['V']
        standard_col = landolt_fl['BV']
    elif filt == 'R':
        standard_mag = landolt_fl['V'] - landolt_fl['VR']
        standard_col = landolt_fl['VR']
    elif filt == 'I':
        standard_mag = landolt_fl['V'] - landolt_fl['VI']
        standard_col = landolt_fl['VI']

    return list(standard_mag), list(standard_col)


def calibrate_magnitudes(tab, itime=1., zmag=25.):
    tab['cal_mags'] = (zmag - 2.5 * np.log10(tab['flux_fit'] / itime)) * u.mag
    return tab


def instrumMags(f_name, tfilt, hdu_data, exp_time, aper_rad, annulus_in,
                annulus_out):
    """
    Perform aperture photometry for all 'tfilt' standard stars observed in
    the 'hdu_data' file.
    """
    # Coordinates from observed frame.
    positions = zip(*[tfilt['x_obs'], tfilt['y_obs']])
    apertures = CircularAperture(positions, r=aper_rad)
    annulus_apertures = CircularAnnulus(
        positions, r_in=annulus_in, r_out=annulus_out)

    apers = [apertures, annulus_apertures]
    # TODO obtain errors for aperture photometry.

    phot_table = aperture_photometry(hdu_data, apers)
    bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area()
    bkg_sum = bkg_mean * apertures.area()
    phot_table['flux_fit'] = phot_table['aperture_sum_0'] - bkg_sum
    phot_table = calibrate_magnitudes(phot_table, itime=exp_time)

    # Change weird fluxes for sources outside the frame (-999.9, -999.9)
    # to 'nan'.
    for i, (x, y) in enumerate(phot_table['xcenter', 'ycenter']):
        if x.value < -999. and y.value < -999.:
            phot_table['cal_mags'][i] = np.nan

    return phot_table


def zeroAirmass(phot_table, extin_coeffs, filt, airmass):
    """
    Correct for airmass, i.e. instrumental magnitude at zero airmass.
    """
    # Identify correct index for this filter's extinction coefficient.
    f_idx = extin_coeffs.index(filt) + 1
    # Extinction coefficient.
    ext = float(extin_coeffs[f_idx])
    # Obtain zero airmass instrumental magnitude for this filter.
    phot_table['instZA'] = phot_table['cal_mags'] - (ext * airmass) * u.mag

    return phot_table


def writeAperPhot(out_path, filters):
    """
    """
    tables = []
    for v in filters.values():
        tables.append(Table(zip(*v)))
    aper_phot = Table(
        vstack(tables),
        names=('Filt', 'Stnd_field', 'ID', 'file', 'exp_t', 'A', 'ZA_mag',
               'Col', 'Mag'))

    ascii.write(
        aper_phot, out_path + '/stnd_aperphot.dat',
        format='fixed_width', delimiter=' ', formats={'ZA_mag': '%10.4f'},
        fill_values=[(ascii.masked, 'nan')], overwrite=True)


def main():
    """
    Performs aperture photometry on selected standard fields.

    Requires the '_match.coo' file generated by the 'match' script.

    Returns zero airmass corrected instrumental magnitudes for each filter.
    """
    pars, fits_list, in_path, out_path = in_params()

    filters = {'U': [], 'B': [], 'V': [], 'R': [], 'I': []}

    for proc_gr, stnd_fl in enumerate(pars['landolt_fld']):
        print("\nPerform aperture photometry and fit transformation\n"
              "equations for the standard field: {}\n".format(stnd_fl))

        # Read data for this Landolt field.
        landolt_fl = read_standard_coo(out_path, stnd_fl)

        # For each observed .fits standard file, from '.mch' file.
        for fr in fits_list[proc_gr]:
            f_name = fr.split('/')[-1].split('.')[0]
            tfilt = landolt_fl[landolt_fl['frame'] == f_name]

            print("Aperture photometry on: {}".format(f_name))

            # Load .fits file.
            hdulist = fits.open(join(in_path, fr))
            # Extract header and data.
            hdr, hdu_data = hdulist[0].header, hdulist[0].data
            filt, exp_time, airmass = hdr[pars['filter_key']],\
                hdr[pars['exposure_key']], hdr[pars['airmass_key']]
            print("  Filter {}, Exp time {}, Airmass {}".format(
                filt, exp_time, airmass))

            # Obtain instrumental magnitudes for the standard stars in the
            # defined Landolt field, in this observed frame.
            photu = instrumMags(
                f_name, tfilt, hdu_data, exp_time,
                float(pars['aperture']), float(pars['annulus_in']),
                float(pars['annulus_out']))

            print("  Correct instrumental magnitudes for zero airmass.")
            photu = zeroAirmass(photu, pars['extin_coeffs'][0], filt, airmass)

            # Extract data for this filter.
            stand_mag, stand_col = standardMagCol(tfilt, filt)
            # Group frames by filter.
            # Filt  Stnd_field ID  file  exp_t  A  ZA_mag  Col  Mag
            for i, ID in enumerate(tfilt['ID']):
                filters[filt].append(
                    [filt, stnd_fl, ID, f_name, str(exp_time), airmass,
                     photu['instZA'][i].value, stand_col[i], stand_mag[i]])

    # Remove not observed filters from dictionary.
    filters = {k: v for k, v in filters.iteritems() if v}
    if 'V' not in filters.keys():
        print("  WARNING: Filter V is missing.")
        if 'B' not in filters.keys():
            print("  WARNING: Filter B is missing.")

    print("\nWrite final aperture photometry to output file.")
    writeAperPhot(out_path, filters)


if __name__ == '__main__':
    main()
