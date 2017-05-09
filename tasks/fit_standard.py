
import read_pars_file as rpf

import os
from os.path import join, isfile
import sys
import numpy as np
from scipy.stats import linregress
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
            if f.endswith('_crop.fits'):
                if f.split('/')[-1] in pars['stnd_obs_fields'][1:]:
                    fits_list.append(f)

    if not fits_list:
        print("No '*_crop.fits' files found in 'output/standards' folder."
              " Exit.")
        sys.exit()

    return pars, fits_list, in_out_path


def read_standard_coo(in_out_path, landolt_fld):
    """
    Read .coo file with standard stars coordinates in the system of the
    observed frame.
    """
    f = join(in_out_path, landolt_fld + '_obs.coo')
    landolt_fl = ascii.read(f)

    return landolt_fl


def standardMagnitude(landolt_fl, filt):
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

    return np.asarray(standard_mag), np.asarray(standard_col)


def calibrate_magnitudes(tab, itime=1., zmag=25.):
    tab['cal_mags'] = (zmag - 2.5 * np.log10(tab['flux_fit'] / itime)) * u.mag
    return tab


def instrumMags(landolt_fl, hdu_data, exp_time, aper_rad, annulus_in,
                annulus_out):
    """
    Perform aperture photometry for all 'landolt_fl' standard stars observed in
    the 'hdu_data' file.
    """
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


def zeroAirmass(phot_table, extin_coeffs, filt, airmass):
    """
    Correct for airmass, i.e. instrumental. magnitude at zero airmass.
    """
    f_idx = extin_coeffs.index(filt) + 1
    ext = float(extin_coeffs[f_idx])
    phot_table['instZA'] = phot_table['cal_mags'] - (ext * airmass) * u.mag

    return phot_table


def rmse(targets, predictions):
    return np.sqrt(((predictions - targets) ** 2).mean())


def redchisqg(ydata, ymod, deg=2, sd=.05):
    """
    Returns the reduced chi-square error statistic for an arbitrary model,
    chisq/nu, where nu is the number of degrees of freedom. If individual
    standard deviations (array sd) are supplied, then the chi-square error
    statistic is computed as the sum of squared errors divided by the standard
    deviations. See http://en.wikipedia.org/wiki/Goodness_of_fit for reference.
    ydata : data
    ymod : model evaluated at the same x points as ydata
    """
    print((ydata - ymod) / sd) ** 2
    chisq = np.sum(((ydata - ymod) / sd) ** 2)
    print("  chisq: ", chisq)

    # Number of degrees of freedom assuming 2 free parameters
    nu = len(ydata) - 1. - deg
    print("  nu:", nu)

    return chisq / nu


def distPoint2Line(m, c, x, y):
    """
    Distance from (x, y) point, to line with equation:

    y = m*x + c

    http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
    """
    d = abs(m * np.array(x) + c - np.array(y)) / np.sqrt(m ** 2 + 1.)
    d_sort = sorted(zip(d, x, y))
    d, x, y = zip(*d_sort)[:3]

    return x, y


def regressRjctOutliers(x, y, chi_min=.95, RMSE_max=.05):
    """
    Perform a linear regression fit, rejecting outliers until the conditions
    of abs(1 - red_chi)<chi_min_delta and RMSE<RMSE_max are met.
    """
    m, c, r_value, p_value, std_err = linregress(x, y)
    predictions = m * np.array(x) + c
    red_chisq = redchisqg(y, predictions)
    print("Red Chi^2: {:.3f}".format(red_chisq))
    chi, RMSE = r_value**2, rmse(y, predictions)
    print("N, m, c, R^2, RMSE: {}, {:.3f}, {:.3f}, {:.3f}, {:.3f}".format(
        len(x), m, c, chi, RMSE))

    x_accpt, y_accpt, x_rjct, y_rjct = x[:], y[:], [], []
    while chi < chi_min or RMSE > RMSE_max:
        x_dsort, y_dsort = distPoint2Line(m, c, x_accpt, y_accpt)
        x_accpt, y_accpt = x_dsort[:-1], y_dsort[:-1]
        x_rjct.append(x_dsort[-1])
        y_rjct.append(y_dsort[-1])

        m, c, r_value, p_value, std_err = linregress(x_accpt, y_accpt)
        predictions = m * np.array(x_accpt) + c
        red_chisq = redchisqg(y_accpt, predictions)
        print("Red Chi^2: {:.3f}".format(red_chisq))
        # STD = np.std(y - predictions)
        chi, RMSE = r_value**2, rmse(y_accpt, predictions)
        print("N, m, c, R^2, RMSE: {}, {:.3f}, {:.3f}, {:.3f}, {:.3f}".format(
            len(x_accpt), m, c, chi, RMSE))

    if len(x_accpt) == 2:
        print("  WARNING: only two stars were used in the fit.")

    xy_accpt, xy_rjct = [x_accpt, y_accpt], [x_rjct, y_rjct]

    return xy_accpt, xy_rjct, m, c, chi, RMSE


def fitTransfEquations(filters):
    """
    Solve transformation equation to the Landolt system.
    """
    # Remove not observed filters from dictionary.
    filters = {k: v for k, v in filters.iteritems() if v}

    # Assume that the 'V' filter was used.
    stand_V, stand_BV, instZA_V = filters['V']
    x_fit, y_fit = stand_BV, (stand_V - instZA_V)
    # Find slope and/or intersection of linear regression.
    V_accpt, V_rjct, Vm, Vc, Vchi, VRMSE = regressRjctOutliers(x_fit, y_fit)
    print("V: {:.2f}, {:.2f}".format(Vm, Vc))
    filt_col_data = [["V", V_accpt, V_rjct, Vm, Vc, Vchi, VRMSE]]

    if 'B' in filters.keys():
        stand_B, stand_BV, instZA_B = filters['B']
        x_fit, y_fit = (instZA_B - instZA_V), stand_BV
        BV_accpt, BV_rjct, BVm, BVc, BVchi, BVRMSE = regressRjctOutliers(
            x_fit, y_fit)
        print("BV: {:.2f}, {:.2f}".format(BVm, BVc))
        filt_col_data.append(
            ["BV", BV_accpt, BV_rjct, BVm, BVc, BVchi, BVRMSE])

        if 'U' in filters.keys():
            stand_U, stand_UB, instZA_U = filters['U']
            x_fit, y_fit = (instZA_U - instZA_B), stand_UB
            UB_accpt, UB_rjct, UBm, UBc, UBchi, UBRMSE = regressRjctOutliers(
                x_fit, y_fit)
            print("UB: {:.2f}, {:.2f}".format(UBm, UBc))
            filt_col_data.append(
                ["UB", UB_accpt, UB_rjct, UBm, UBc, UBchi, UBRMSE])

    if 'I' in filters.keys():
        stand_I, stand_VI, instZA_I = filters['I']
        x_fit, y_fit = (instZA_V - instZA_I), stand_VI
        VI_accpt, VI_rjct, VIm, VIc, VIchi, VIRMSE = regressRjctOutliers(
            x_fit, y_fit)
        print("VI: {:.2f}, {:.2f}".format(VIm, VIc))
        filt_col_data.append(
            ["VI", VI_accpt, VI_rjct, VIm, VIc, VIchi, VIRMSE])

    if 'R' in filters.keys():
        stand_R, stand_VR, instZA_R = filters['R']
        x_fit, y_fit = (instZA_V - instZA_R), stand_VR
        VR_accpt, VR_rjct, VRm, VRc, VRchi, VRRMSE = regressRjctOutliers(
            x_fit, y_fit)
        print("VR: {:.2f}, {:.2f}".format(VRm, VRc))
        filt_col_data.append(
            ["VR", VR_accpt, VR_rjct, VRm, VRc, VRchi, VRRMSE])

    return filt_col_data


def writeTransfCoeffs(mypath, filt_col_data):
    """
    """
    f_path = join(
        mypath.replace('tasks', 'output/standards'), 'fit_coeffs.dat')
    with open(f_path, 'w') as f:
        hdr = """#\n# Instrumental zero airmass magnitude.\n# F_I^(0A) =""" +\
            """ F_I  - K * X  ; X: airmass, K: extinction coefficient\n""" +\
            """#\n# Standard color.\n# (F_1 - F_2)_L = c_1 """ +\
            """(F_1 - F_2)_I^(0A) + c_2\n#\n# Standard magnitude.\n""" +\
            """# V_L = V_I^(0A) + c_1 * (B - V)_L + c_2\n#\n# ID  N_a""" +\
            """   N_r       c_1      c_2  red_Chi     RMSE\n"""
        f.write(hdr)

    with open(f_path, 'a') as f:
        for fc in filt_col_data:
            ln = "{:2}      {}     {}".format(
                fc[0], len(fc[1][0]), len(fc[2][0]))
            ln += "    {:6.3f}   {:6.3f}   {:6.3f}   {:6.3f}\n".format(
                *map(float, fc[3:]))
            f.write(ln)


def make_plot(mypath, filt_col_data):
    """
    """
    print("\nPlotting.")
    # print(plt.style.available)
    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(10, 25))
    gs = gridspec.GridSpec(5, 2)

    i = 0
    for data in filt_col_data:
        if data[0] == "V":
            xlbl, ylbl = r'$(B-V)_L$', r'$(V_L-V^{0A}_{I})$'
        elif data[0] == "BV":
            xlbl, ylbl = r'$(B-V)^{0A}_{I}$', r'$(B-V)_L$'
        elif data[0] == "UB":
            xlbl, ylbl = r'$(U-B)^{0A}_{I}$', r'$(U-B)_L$'
        elif data[0] == "VI":
            xlbl, ylbl = r'$(V-I)^{0A}_{I}$', r'$(V-I)_L$'
        elif data[0] == "VR":
            xlbl, ylbl = r'$(V-R)^{0A}_{I}$', r'$(V-R)_L$'

        X_accpt, X_rjct, m, c, chi, RMSE = data[1:]
        x_accpt, y_accpt = X_accpt
        x_rjct, y_rjct = X_rjct

        x_r = np.arange(min(x_accpt), max(x_accpt) * 1.1, .01)
        predictions = m * x_r + c

        ax = fig.add_subplot(gs[0 + i * 2])
        sign = '+' if c >= 0. else '-'
        ax.set_title("{} = {:.3f} {} {} {:.3f}".format(
            ylbl, m, xlbl, sign, abs(c)), fontsize=9)
        plt.xlabel(xlbl)
        plt.ylabel(ylbl)
        ax.plot(x_r, predictions, 'b-', zorder=1)
        plt.scatter(x_accpt, y_accpt, c='g', zorder=4)
        plt.scatter(x_rjct, y_rjct, c='r', marker='x', zorder=4)

        # t1 = "m, c: {:.3f}, {:.3f}\n".format(m, c)
        t1 = r"$R^{{2}}={:.3f}$".format(chi) + "\n"
        t2 = r"$RMSE={:.3f}$".format(RMSE)
        loc = 2 if i != 0 else 3
        txt = AnchoredText(t1 + t2, loc=loc, prop=dict(size=9))
        txt.patch.set(boxstyle='square,pad=0.', alpha=0.75)
        ax.add_artist(txt)

        ax = fig.add_subplot(gs[1 + i * 2])
        plt.xlabel(ylbl)
        plt.ylabel("residuals (Landolt - pred)")
        resid_accpt = y_accpt - (m * np.asarray(x_accpt) + c)
        plt.scatter(y_accpt, resid_accpt, c='g', zorder=4)
        ax.axhline(0., linestyle=':', color='b', zorder=1)
        t1 = r"$\sigma={:.3f}$".format(np.std(resid_accpt))
        txt = AnchoredText(t1, loc=2, prop=dict(size=9))
        txt.patch.set(boxstyle='square,pad=0.', alpha=0.75)
        ax.add_artist(txt)

        i += 1

    fig.tight_layout()
    fn = join(mypath.replace('tasks', 'output/standards'), 'fitstand.png')
    plt.savefig(fn, dpi=150, bbox_inches='tight')
    # Close to release memory.
    plt.clf()
    plt.close()


def main():
    """
    """
    pars, fits_list, in_out_path = in_params()
    landolt_fl = read_standard_coo(in_out_path, pars['stnd_obs_fields'][0])

    filters = {'U': [], 'B': [], 'V': [], 'R': [], 'I': []}
    # For each _crop.fits observed (aligned and cropped) standard file.
    for imname in fits_list:
        fname = imname.replace(pars['mypath'].replace('tasks', 'output'), '')

        # Load .fits file.
        hdulist = fits.open(imname)
        # Extract header and data.
        hdr, hdu_data = hdulist[0].header, hdulist[0].data
        filt, exp_time, airmass = hdr[pars['filter_key']],\
            hdr[pars['exposure_key']], hdr[pars['airmass_key']]

        stand_mag, stand_col = standardMagnitude(landolt_fl, filt)
        print("Aperture photometry on: {}".format(fname))
        photu = instrumMags(
            landolt_fl, hdu_data, exp_time, float(pars['aperture']),
            float(pars['annulus_in']), float(pars['annulus_out']))
        photu = zeroAirmass(photu, pars['extin_coeffs'], filt, airmass)

        # Group frames by filter.
        filters[filt] = [stand_mag, stand_col, photu['instZA'].value]

    filt_col_data = fitTransfEquations(filters)

    writeTransfCoeffs(pars['mypath'], filt_col_data)

    if pars['do_plots_D'] == 'y':
        make_plot(pars['mypath'], filt_col_data)

    print("\nFinished.")


if __name__ == '__main__':
    main()
