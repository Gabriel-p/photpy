
import read_pars_file as rpf
from hlpr import rmNaNrows

import os
from os.path import join, isdir
import sys
import numpy as np
import warnings
from astropy.io import ascii
from astropy.table import Table, hstack
import matplotlib.pyplot as plt


def main():

    pars = in_params()

    # Read cross-matched photometry.
    group_phot = readPhot(pars['files_list'])

    # Combine magnitudes for each filter, for each exposure time.
    avrg_phot, id_coords = avrgMags(
        group_phot, pars['airmasses'], pars['extin_coeffs'][0], pars['method'])

    # Transform combined magnitudes to the standard system.
    stand_phot = standardCalib(avrg_phot)

    # Create final output file and make final plot.
    writeToFile(pars['in_out_path'], stand_phot, id_coords)
    if pars['do_plots_I'] == 'y':
        make_plots(pars['in_out_path'])


def in_params():
    """
    Read and prepare input parameter values.

    Returns
    -------
    pars : dictionary
        All parameters.
    """
    pars = rpf.main()

    pars['in_out_path'] = pars['mypath'].replace(
        'tasks', 'output/' + pars['transf_folder'])

    f_list = []
    if isdir(pars['in_out_path']):
        for file in os.listdir(pars['in_out_path']):
            if file.startswith('filter_') and file.endswith('.mag'):
                f_list.append(join(pars['in_out_path'], file))
    else:
        print("{}\n is not a folder. Exit.".format(pars['in_out_path']))
        sys.exit()

    # TODO this is a convoluted way to read this info and it is HARDCODED
    files_data = ascii.read(join(pars['in_out_path'], 'list'))
    airmass = {'U': {}, 'B': {}, 'V': {}, 'I': {}, 'R': {}}
    for fname, filt, expT, K, x0, y0 in files_data:
        airmass[filt].update({round(expT, 2): K})
    pars['airmasses'] = airmass

    pars['files_list'] = f_list
    return pars


def readPhot(f_list):
    """
    # TODO
    """

    group_phot = {'U': {}, 'B': {}, 'V': {}, 'I': {}, 'R': {}}
    for f in f_list:
        filt = f.split('_')[-1].replace('.mag', '')
        # if filt == 'U':  # TODO remove
        group_phot[filt] = ascii.read(f)

    group_phot = {k: v for k, v in group_phot.items() if v}

    return group_phot


def avrgMags(group_phot, airmasses, extin_coeffs, method):
    """
    Combine magnitudes for all cross-matched stars in all frames. Available
    methods are: #TODO

    Parameters
    ----------
    group_phot : dict
        Dictionary of astropy Tables(), one per observed filter with columns
        ordered as 'x_XXX  y_XXX  mag_XXX  emag_XXX', where 'XXX' represents
        the exposure time.
    method : string
        Selected method to obtain the final magnitudes and sigmas for each
        observed filter.

    Returns
    -------
    avrg_phot : dictionary
        Combined magnitudes for each filter, using the selected method.
    id_coords : dictionary
        IDs and averaged coordinates for each filter.

    """
    # Process each observed filter.
    avrg_phot, xcen_all, ycen_all = {}, [], []
    for filt, fvals in group_phot.iteritems():
        print("Processing filter: {}".format(filt))
        # Extract coordinates.
        for col in fvals.itercols():
            if col.name.startswith('x_'):
                xcen_all.append(col)
            elif col.name.startswith('y_'):
                ycen_all.append(col)

        # Extract magnitudes and exposure times.
        mags, emags, exp_t = [], [], []
        for col in fvals.itercols():
            if col.name.startswith('mag_'):
                expTime = col.name.split('_')[1]
                exp_t.append(expTime)
                # Correct to zero airmass magnitudes.
                # K = airmasses[filt][round(float(expTime), 2)]
                # print(" To zero airmass, expT={}, K={}".format(expTime, K))
                # zAmag = magtoZA(extin_coeffs, filt, K, col)
                # mags.append(zAmag)  # TODO uncmmt
                mags.append(col)
            elif col.name.startswith('emag_'):
                emags.append(col)
        mags = np.asarray(mags)

        # We expect some RuntimeWarnings in this block, so suppress them.
        # with warnings.catch_warnings():
        #     warnings.simplefilter("ignore", category=RuntimeWarning)

        if method == 'mean':
            # Index of the frame with the longest exposure: master frame.
            m_fr_id = [float(_) for _ in exp_t].index(
                max([float(_) for _ in exp_t]))
            print("Reference frame selected: expT={}".format(exp_t[m_fr_id]))
            # Median magnitude difference between master frame and all others.
            dmag = np.nanmedian(mags - mags[m_fr_id], axis=1)

            # Convert magnitudes and their sigmas to flux.
            flux = mag2flux(mags - dmag[:, None])
            flux = mag2flux(mags)
            eflux = emag2eflux(emags, flux)

            # Average flux, and convert back to magnitudes.
            flux_mean = np.nanmean(zip(*flux), axis=1)
            mag_mean = flux2mag(flux_mean)

            # Idem for errors (sigmas)
            # Obtain variances.
            flux_vari = np.array(zip(*np.array(eflux) ** 2))
            # Count non-nan values
            non_nans = (~np.isnan(flux_vari)).sum(1)
            # Replace 0 count with np.nan
            non_nans = np.where(non_nans == 0, np.nan, non_nans)
            # Sigma for the flux mean, obtained as:
            # e_f = sqrt(sum(var ** 2)) / N =
            #     = sqrt(mean(var) / N)
            eflux_mean = np.sqrt(
                (1. / non_nans) * np.nanmean(flux_vari, axis=1))
            # Back to magnitudes.
            emag_mean = eflux2emag(eflux_mean, flux_mean)

        elif method == 'stetson':
            print("Method: Stetson")
            # Index of the frame with the longest exposure: master frame.
            m_fr_id = [float(_) for _ in exp_t].index(
                max([float(_) for _ in exp_t]))
            print("Reference frame selected: expT={}".format(exp_t[m_fr_id]))

            # Identify master frame as the one with the longest exposure.
            flux_master = mag2flux(mags[m_fr_id])
            N_frames = len(mags)

            # TODO check order
            # Median magnitude difference between master frame and all others.
            dmag = np.nanmedian(mags - mags[m_fr_id], axis=1)
            # Flux in master's system.
            flux_m = mag2flux(mags - dmag[:, None])
            # Variance of the flux in master's system.
            eflux = emag2eflux(emags, flux_m)
            var_flux_m = np.array(zip(*np.array(eflux) ** 2))

            mag_mean = []
            for i, st_f_m in enumerate(flux_master):
                print('1', i, st_f_m)
                B = np.inf
                while abs(B) >= 0.00005 * st_f_m:
                    sum_m, sumwt = 0., 0.
                    for j in range(N_frames):
                        st_f, st_f_var = flux_m[j][i], var_flux_m[j][i]
                        print('2', B, st_f, st_f_var)
                        if ~np.isnan(st_f):
                            m_diff = st_f - st_f_m
                            WT = (1. / st_f_var) *\
                                (4. / (4. + (m_diff**2 / st_f_var)))
                            sum_m += m_diff * WT
                            sumwt += WT
                            print('3', j, m_diff, sum_m, sumwt)

                    B = sum_m / sumwt
                    st_f_m = st_f_m + B
                    print('4', B)
                print('5', st_f_m)
                mag_mean.append(st_f_m)
                import pdb; pdb.set_trace()  # breakpoint 1d8e8ee4 //


        # Add photometric data for this filter.
        avrg_phot.update({filt: [mag_mean, emag_mean]})

    # with warnings.catch_warnings():
    #     warnings.simplefilter("ignore", category=RuntimeWarning)
    # Median for (x, y) coordinates.
    xmedian = np.nanmedian(zip(*xcen_all), axis=1)
    ymedian = np.nanmedian(zip(*ycen_all), axis=1)
    # Store IDs and averaged coordinate data for this filter.
    id_coords = [group_phot['V']['ID'], xmedian, ymedian]

    return avrg_phot, id_coords


def magtoZA(extCoeffs, filt, K, mag):
    """
    Convert instrumental magnitudes to zero airmass magnitudes.
    This needs to be done here because when the frames are cross-matched,
    their magnitudes are averaged.
    """
    ec_idx = int(extCoeffs.index(filt) + 1)
    # Correct to zero airmass magnitudes. This needs to be done here because
    zAmag = mag - K * float(extCoeffs[ec_idx])

    return zAmag


def mag2flux(mag, zmag=15.):
    """
    Convert magnitudes to flux.

    Parameters
    ----------
    mag : array
        Magnitude values.
    zmag : float
        Arbitrary (fixed) constant.

    Returns
    -------
    array
        Fluxes.

    """
    return 10 ** ((zmag - mag) / 2.5)


def emag2eflux(emags, flux):
    """
    Convert magnitude sigmas into flux sigmas.

    Parameters
    ----------
    emags : array
        Magnitude sigma values.
    flux : array
        Flux values.

    Returns
    -------
    array
        Flux sigmas.

    Notes
    -----
    The float used is:
        1.0857362 = 2.5 * log10(e)

    """
    return emags * flux * 1.0857362


def flux2mag(flux, zmag=15.):
    """
    Convert fluxes into flux magnitudes.

    Parameters
    ----------
    flux : array
        Flux values.
    zmag : float
        Arbitrary (fixed) constant.

    Returns
    -------
    array
        Magnitudes.

    """
    return -2.5 * np.log10(flux) + zmag


def eflux2emag(eflux, flux_mean):
    """
    Convert flux sigmas into magnitude sigmas.

    Parameters
    ----------
    eflux : array
        Flux sigmas.
    flux_mean : array
        Flux mean values.

    Returns
    -------
    array
        Magnitude sigmas.

    Notes
    -----
    The float used is:
        1.0857362 = 2.5 * log10(e)

    """
    return eflux * 1.0857362 / flux_mean


def inst2cal(avrg_phot, ab_coeffs, col):
    """
    Transform an instrumental (zero airmass) color and its standard deviation
    to the calibrated system, using coefficients previously obtained.

    Parameters
    ----------
    avrg_phot : dictionary
        See avrgMags()
    ab_coeffs : dictionary
        See standardCalib()
    col : string
        Identifies the color to process.

    Returns
    -------
    cal_col : array
        Calibrated color.
    cal_sig : array
        Calibrated standard deviation for the processed color.

    """
    f1, f2 = col

    # Instrumental (zero airmass) magnitudes.
    f1_inst, f2_inst = avrg_phot[f1][0], avrg_phot[f2][0]
    # Instrumental magnitudes color.
    f12_inst = f1_inst - f2_inst
    # To calibrated standard system.
    cal_col = ab_coeffs['a'][col] * f12_inst + ab_coeffs['b'][col]

    # Instrumental magnitude sigmas.
    ef1_inst, ef2_inst = avrg_phot[f1][1], avrg_phot[f2][1]
    # Instrumental color sigma.
    ef12_inst = ef1_inst + ef2_inst
    # Calibrated color sigma.
    cal_sig = np.sqrt(
        (ab_coeffs['sa'][col] * cal_col) ** 2 +
        (ab_coeffs['a'][col] * ef12_inst) ** 2 +
        ab_coeffs['sb'][col] ** 2)

    return cal_col, cal_sig


def inst2calMag(avrg_phot, ab_coeffs, col_cal, col_sig, mag):
    """
    Transform an instrumental (zero airmass) magnitude and its standard
    deviation to the calibrated system, using coefficients previously obtained.

    Parameters
    ----------
    avrg_phot : dictionary
        See avrgMags()
    ab_coeffs : dictionary
        See standardCalib()
    col_cal : array
        Calibrated color (color term in transformation equation)
    col_sig : array
        Calibrated standard deviation for the color.
    mag : string
        Identifies the magnitude to process.

    Returns
    -------
    cal_mag : array
        Calibrated magnitude.
    cal_sig : array
        Calibrated standard deviation for the processed magnitude.

    """
    cal_mag = avrg_phot[mag][0] + ab_coeffs['a'][mag] * col_cal +\
        ab_coeffs['b'][mag]
    # Use the  calibrated standard deviation for the selected color.
    cal_sig = np.sqrt(
        avrg_phot[mag][1] ** 2 + (ab_coeffs['a'][mag] * col_sig) ** 2 +
        (ab_coeffs['sa'][mag] * col_cal) ** 2 + ab_coeffs['sb'][mag] ** 2)

    return cal_mag, cal_sig


def standardCalib(avrg_phot, ab_coeffs={}):
    """
    Transform instrumental (zero airmass) magnitudes into the standard V
    magnitude and standard colors.

    !!! --> Assumes that the V and B filters are present. <-- !!!

    Parameters
    ----------
    avrg_phot : dictionary
        See avrgMags()
    ab_coeffs : dictionary
        TODO

    Returns
    -------
    stand_phot : dictionary
        Calibrated standard V magnitude and colors, and their standard
        deviations.

    """
    # ab_coeffs = {
    #     'a': {'V': -0.066, 'BV': 1.22, 'UB': 0.98, 'VI': 0.9},
    #     'sa': {'V': 0.03, 'BV': 0.01, 'UB': 0.01, 'VI': 0.01},
    #     'b': {'V': -1.574, 'BV': -0.109, 'UB': -1.993, 'VI': 0.081},
    #     'sb': {'V': 0.03, 'BV': 0.01, 'UB': 0.01, 'VI': 0.01}
    # }

    # DELETE
    # # Values for the .als files used for testing.
    # u1, u3 = 3.964132, -0.135271
    # b1, b3 = 2.086563, 0.01751888
    # v1, v3 = 1.445104, -0.02613502
    # i1, i3 = 1.823921, 0.09020784
    # aV, bV = -v3, -v1
    # aBV, bBV = 1. / (1. - v3 + b3), (v1 - b1) / (1. - v3 + b3)
    # aUB, bUB = 1. / (1. - b3 + u3), (b1 - u1) / (1. - b3 + u3)
    # aVI, bVI = 1. / (1. - i3 + v3), (i1 - v1) / (1. - i3 + v3)
    # ab_coeffs = {
    #     'a': {'V': aV, 'BV': aBV, 'UB': aUB, 'VI': aVI},
    #     'sa': {'V': 0.03, 'BV': 0.01, 'UB': 0.01, 'VI': 0.01},
    #     'b': {'V': bV, 'BV': bBV, 'UB': bUB, 'VI': bVI},
    #     'sb': {'V': 0.03, 'BV': 0.01, 'UB': 0.01, 'VI': 0.01}
    # }
    # DELETE

    ab_coeffs = {
        'a': {'V': -0.065, 'BV': 1.218, 'UB': .962, 'VI': .901},
        'sa': {'V': 0.01, 'BV': 0.01, 'UB': 0.01, 'VI': 0.01},
        'b': {'V': -1.645, 'BV': .006, 'UB': -1.939, 'VI': .095},
        'sb': {'V': 0.01, 'BV': 0.01, 'UB': 0.01, 'VI': 0.01}
    }

    # Obtain the calibrated BV color and its standard deviation.
    BV_cal, BV_sig = inst2cal(avrg_phot, ab_coeffs, 'BV')
    # Add to dictionary.
    stand_phot = {'BV': BV_cal, 'eBV': BV_sig}

    # For the V magnitude use the obtained calibrated BV color.
    V_cal, V_sig = inst2calMag(avrg_phot, ab_coeffs, BV_cal, BV_sig, 'V')
    stand_phot.update({'V': V_cal, 'eV': V_sig})

    # Transform the rest of the filters, if they are available.
    if 'U' in avrg_phot.keys():
        UB_cal, UB_sig = inst2cal(avrg_phot, ab_coeffs, 'UB')
        stand_phot.update({'UB': UB_cal, 'eUB': UB_sig})
    if 'I' in avrg_phot.keys():
        VI_cal, VI_sig = inst2cal(avrg_phot, ab_coeffs, 'VI')
        stand_phot.update({'VI': VI_cal, 'eVI': VI_sig})
    if 'R' in avrg_phot.keys():
        VR_cal, VR_sig = inst2cal(avrg_phot, ab_coeffs, 'VR')
        stand_phot.update({'VR': VR_cal, 'eVR': VR_sig})

    return stand_phot


def writeToFile(in_out_path, stand_phot, id_coords):
    """
    TODO
    """
    print("\nWriting 'final_phot.dat' file.")
    t1 = Table(id_coords, names=('ID', 'xmedian', 'ymedian'))
    t2 = Table(stand_phot)
    t3 = hstack([t1, t2])
    tab = rmNaNrows(t3, 3)
    ascii.write(
        tab, in_out_path + '/final_phot.dat',
        format='fixed_width', delimiter=' ',
        formats={_: '%10.4f' for _ in tab.keys()[1:]},
        fill_values=[(ascii.masked, 'nan')], overwrite=True)


def make_plots(in_out_path):
    """
    TODO
    """


if __name__ == '__main__':
    main()
