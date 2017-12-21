
import read_pars_file as rpf

import os
from os.path import join, isfile
import sys
import numpy as np
from scipy.optimize import leastsq, least_squares, curve_fit
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText

from astropy.io import ascii


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
                fits_list.append(f)

    if not fits_list:
        print("No '*_crop.fits' files found in 'output/standards' folder."
              " Exit.")
        sys.exit()

    return pars, in_out_path


def readData(in_out_path):
    """
    """
    f = join(in_out_path, "stnd_aperphot.dat")
    stndData = ascii.read(f)

    for f_id in ['V', 'B']:
        if f_id not in stndData['Filt']:
            print("Filter {} is missing.".format(f_id))
            sys.exit()

    # Group by filters
    f_grouped = stndData.group_by('Filt')

    return f_grouped


def rmse(targets, predictions):
    return np.sqrt(((predictions - targets) ** 2).mean())


def redchisqg(ydata, ymod, deg=2, sd=.01):
    """
    Returns the reduced chi-square error statistic for an arbitrary model,
    chisq/nu, where nu is the number of degrees of freedom. If individual
    standard deviations (array sd) are supplied, then the chi-square error
    statistic is computed as the sum of squared errors divided by the standard
    deviations. See http://en.wikipedia.org/wiki/Goodness_of_fit for reference.
    ydata : data
    ymod : model evaluated at the same x points as ydata
    """
    chisq = np.sum(((ydata - ymod) / sd) ** 2)

    # Number of degrees of freedom assuming 2 free parameters
    nu = len(ydata) - deg

    return chisq / nu


def f(x, a, b):
    return a * x + b


def residual(p, x, y):
    return y - f(x, *p)


def distPoint2Line(m, c, x, y, z):
    """
    Distance from (x, y) point, to line with equation:

    y = m*x + c

    http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
    """
    d = abs(m * np.array(x) + c - np.array(y)) / np.sqrt(m ** 2 + 1.)
    d_sort = sorted(zip(d, x, y, z))
    x, y, z = zip(*d_sort)[1:]

    return x, y, z


def regressRjctOutliers(x, y, z, R2_min=.97, RMSE_max=.03):
    """
    Perform a linear regression fit, rejecting outliers until the conditions
    of abs(1 - red_chi)<chi_min_delta and RMSE<RMSE_max are met.
    """

    # TODO: read chi_min and RMSE_max
    # TODO: finish reduced chi function
    # TODO: errors for fit coefficients?

    from scipy.stats import linregress
    m, c, r_value, p_value, std_err = linregress(x, y)

    predictions = m * np.array(x) + c
    red_chisq = redchisqg(y, predictions)
    R2, RMSE = r_value**2, rmse(y, predictions)

    x_accpt, y_accpt, z_accpt, x_rjct, y_rjct, z_rjct =\
        x[:], y[:], z[:], [], [], []
    while R2 < R2_min or RMSE > RMSE_max:
        x_dsort, y_dsort, z_dsort = distPoint2Line(
            m, c, x_accpt, y_accpt, z_accpt)
        x_accpt, y_accpt, z_accpt = x_dsort[:-1], y_dsort[:-1], z_dsort[:-1]
        x_rjct.append(x_dsort[-1])
        y_rjct.append(y_dsort[-1])
        z_rjct.append(z_dsort[-1])
        print("   Rjct: {}, ({:.3f}, {:.3f})".format(
            z_dsort[-1], x_dsort[-1], y_dsort[-1]))

        m, c, r_value, p_value, std_err = linregress(x_accpt, y_accpt)
        predictions = m * np.array(x_accpt) + c
        red_chisq = redchisqg(y_accpt, predictions)
        # STD = np.std(y - predictions)
        R2, RMSE = r_value**2, rmse(y_accpt, predictions)

    # TODO optional ways of fitting
    coeff0 = [1., 0.]
    coeff, cov, infodict, mesg, ier = leastsq(
        residual, coeff0, args=(np.array(x_accpt), np.array(y_accpt)),
        full_output=True)
    print("leastsq: m={}, c={}".format(*coeff))

    sol = least_squares(
        residual, coeff0, args=(np.array(x_accpt), np.array(y_accpt)))
    print("least_squares: m={}, c={}".format(*sol.x))

    popt, pcov = curve_fit(f, x_accpt, y_accpt, coeff0)
    print("curve_fit: m={}, c={}".format(*popt))
    # TODO optional ways of fitting

    print("  N, m, c: {}, {:.3f}, {:.3f}".format(len(x_accpt), m, c))
    print("  Red_Chi^2: {:.3f}".format(red_chisq))
    print("  R^2: {:.3f}".format(R2))
    print("  RMSE: {:.3f}".format(RMSE))

    if len(x_accpt) == 2:
        print("  WARNING: only two stars were used in the fit.")

    xy_accpt, xy_rjct = [x_accpt, y_accpt], [x_rjct, y_rjct]

    return xy_accpt, xy_rjct, m, c, R2, RMSE


def extractData(filters, f_id):
    """
    """
    mask = filters.groups.keys['Filt'] == f_id
    f_dat = filters.groups[mask]
    stand_mag, stand_col, instZA_mag = f_dat['Mag'], f_dat['Col'],\
        f_dat['ZA_mag']
    # Flatten extra data used to identify rejected stars
    extra_data = [
        _[0] + ' ' + _[1] + '-' + _[2] + ' (' + str(_[3]) + ')' for _ in
        zip(*[f_dat['Filt'], f_dat['Stnd_field'], f_dat['ID'], f_dat['file']])]

    return stand_mag, stand_col, instZA_mag, extra_data


def fitTransfEquations(filters):
    """
    Solve transformation equation to the Landolt system.
    """
    # Assume that the 'V' filter was used.
    stand_V, stand_BV, instZA_V, z = extractData(filters, 'V')
    x_fit, y_fit = stand_BV, (stand_V - instZA_V)
    # Find slope and/or intersection of linear regression.
    print("Filter V")
    V_accpt, V_rjct, Vm, Vc, Vchi, VRMSE = regressRjctOutliers(
        x_fit.flatten(), y_fit.flatten(), z)
    filt_col_data = [["V", V_accpt, V_rjct, Vm, Vc, Vchi, VRMSE]]

    if 'B' in filters['Filt']:
        _, stand_BV, instZA_B, z = extractData(filters, 'B')
        x_fit, y_fit = (instZA_B - instZA_V), stand_BV
        print("Filter B")
        BV_accpt, BV_rjct, BVm, BVc, BVchi, BVRMSE = regressRjctOutliers(
            x_fit.flatten(), y_fit.flatten(), z)
        filt_col_data.append(
            ["BV", BV_accpt, BV_rjct, BVm, BVc, BVchi, BVRMSE])

        if 'U' in filters['Filt']:
            _, stand_UB, instZA_U, z = extractData(filters, 'U')
            x_fit, y_fit = (instZA_U - instZA_B), stand_UB
            print("Filter U")
            UB_accpt, UB_rjct, UBm, UBc, UBchi, UBRMSE = regressRjctOutliers(
                x_fit.flatten(), y_fit.flatten(), z)
            filt_col_data.append(
                ["UB", UB_accpt, UB_rjct, UBm, UBc, UBchi, UBRMSE])

    if 'I' in filters['Filt']:
        _, stand_VI, instZA_I, z = extractData(filters, 'I')
        x_fit, y_fit = (instZA_V - instZA_I), stand_VI
        print("Filter I")
        VI_accpt, VI_rjct, VIm, VIc, VIchi, VIRMSE = regressRjctOutliers(
            x_fit.flatten(), y_fit.flatten(), z)
        filt_col_data.append(
            ["VI", VI_accpt, VI_rjct, VIm, VIc, VIchi, VIRMSE])

    if 'R' in filters['Filt']:
        _, stand_VR, instZA_R, z = extractData(filters, 'R')
        x_fit, y_fit = (instZA_V - instZA_R), stand_VR
        print("Filter R")
        VR_accpt, VR_rjct, VRm, VRc, VRchi, VRRMSE = regressRjctOutliers(
            x_fit.flatten(), y_fit.flatten(), z)
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
            """   N_r       c_1      c_2  R^2     RMSE\n"""
        f.write(hdr)
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
        if data[0] == "V":
            xlbl, x_scatter = ylbl, y_accpt
        else:
            x_scatter = x_accpt
        plt.xlabel(xlbl)
        plt.ylabel("residuals (Landolt - pred)")
        resid_accpt = y_accpt - (m * np.asarray(x_accpt) + c)
        plt.scatter(x_scatter, resid_accpt, c='g', zorder=4)
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

    !!! --> Assumes that the V and B filters are present. <-- !!!

    """
    pars, in_out_path = in_params()

    filters = readData(in_out_path)

    print("Obtain transformation coefficients.")
    filt_col_data = fitTransfEquations(filters)

    print("Write 'fit_coeffs.dat' output file.")
    writeTransfCoeffs(pars['mypath'], filt_col_data)

    if pars['do_plots_E'] == 'y':
        make_plot(pars['mypath'], filt_col_data)

    print("\nFinished.")


if __name__ == '__main__':
    main()
