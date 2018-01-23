
import read_pars_file as rpf

from os.path import join
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText
import astropy.units as u
from astropy.io import ascii


def in_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    in_out_path = pars['mypath'].replace('tasks', 'output/')

    return pars, in_out_path


def readData(in_out_path):
    """
    """
    f = join(in_out_path, "stnd_aperphot.dat")
    stndData = ascii.read(f)

    # Group by filters
    f_grouped = stndData.group_by('Filt')

    return f_grouped


def zeroAirmass(phot_table, extin_coeffs):
    """
    Correct for airmass, i.e. instrumental magnitude at zero airmass.
    """
    # Identify correct index for this filter's extinction coefficient.
    ext = []
    for filt in phot_table['Filt']:
        f_idx = extin_coeffs.index(filt) + 1
        # Extinction coefficient.
        ext.append(float(extin_coeffs[f_idx]))

    # Obtain zero airmass instrumental magnitude for this filter.
    phot_table['ZA_mag'] = (
        phot_table['mag'] - (ext * phot_table['A'])) * u.mag

    return phot_table


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


def linearFunc(p, x):
    """Linear model."""
    a, b = p
    return a * x + b


def linear_func(x, a, b):
    """Linear model."""
    return a * x + b


def residual(p, x, y):
    return y - linear_func(x, *p)


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
    of R2>R2_min and RMSE<RMSE_max are met.
    """
    # Filter out nan values carried from ZA mag.
    x, y, z = x[~np.isnan(y)], y[~np.isnan(y)],\
        np.array(z)[~np.isnan(y)].tolist()
    x, y, z = x[~np.isnan(x)], y[~np.isnan(x)],\
        np.array(z)[~np.isnan(x)].tolist()

    # TODO: read chi_min and RMSE_max
    # TODO: finish reduced chi function
    # TODO: errors for fit coefficients?

    from scipy.stats import linregress
    m, c, r_value, p_value, std_err = linregress(x, y)
    predictions = m * np.array(x) + c
    # red_chisq = redchisqg(y, predictions)
    R2, RMSE = r_value**2, rmse(y, predictions)

    x_accpt, y_accpt, z_accpt = x[:], y[:], z[:]
    x_rjct, y_rjct, z_rjct = [], [], []
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
        # red_chisq = redchisqg(y_accpt, predictions)
        # STD = np.std(y - predictions)
        R2, RMSE = r_value**2, rmse(y_accpt, predictions)
    print("linregress   : m={:.3f}, c={:.3f}".format(m, c))
    print("  R^2: {:.3f}".format(R2))
    print("  RMSE: {:.3f}".format(RMSE))
    # print("  Red_Chi^2: {:.3f}".format(red_chisq))

    # TODO optional ways of fitting

    from scipy.optimize import leastsq
    coeff0 = [1., 0.]
    coeff, cov, infodict, mesg, ier = leastsq(
        residual, coeff0, args=(np.array(x_accpt), np.array(y_accpt)),
        full_output=True)
    print("leastsq      : m={:.3f}, c={:.3f}".format(*coeff))

    from scipy.optimize import least_squares
    sol = least_squares(
        residual, coeff0, args=(np.array(x_accpt), np.array(y_accpt)))
    print("least_squares: m={:.3f}, c={:.3f}".format(*sol.x))

    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(linear_func, x_accpt, y_accpt, coeff0)
    print("curve_fit    : m={:.3f}, c={:.3f}".format(*popt))

    from scipy.odr import Model, RealData, ODR
    # Create a model for fitting.
    linear_model = Model(linearFunc)
    # Create a RealData object using our initiated data from above.
    data = RealData(x_accpt, y_accpt)  # , sx=xerr, sy=yerr)
    # Set up ODR with the model and data.
    odr = ODR(data, linear_model, beta0=[0., 1.])
    # Run the regression.
    out = odr.run()
    # Estimated parameter values
    beta = out.beta
    print("ODR          : m={:.3f}, c={:.3f}".format(*beta))
    # Standard errors of the estimated parameters
    std = out.sd_beta
    print(" Standard errors: {}, {}".format(*std))
    # Covariance matrix of the estimated parameters
    cov = out.cov_beta
    stddev = np.sqrt(np.diag(cov))
    print(" Squared diagonal covariance: {}".format(stddev))

    import statsmodels.formula.api as smf
    from astropy.table import Table
    df = Table([x_accpt, y_accpt], names=('x', 'y'))
    mod = smf.ols(formula='y ~ x', data=df)
    res = mod.fit()
    print("OLS          : m={:.3f}, c={:.3f}".format(
        res.params['x'], res.params['Intercept']))
    print(" Standard errors: {}, {}".format(
        res.bse['x'], res.bse['Intercept']))

    # TODO optional ways of fitting

    if len(x_accpt) == 2:
        print("  WARNING: only two stars were used in the fit.")

    xy_accpt, xy_rjct = [x_accpt, y_accpt], [x_rjct, y_rjct]

    return xy_accpt, xy_rjct, m, c, R2, RMSE


def extractData(filters, f_id):
    """
    """
    mask = filters.groups.keys['Filt'] == f_id
    f_dat = filters.groups[mask]
    stand_mag, stand_col, instZA_mag = f_dat['Mag_L'], f_dat['Col_L'],\
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
    filt_col_data = []
    for filt in ['U', 'B', 'V', 'R', 'I']:
        if filt in filters['Filt']:
            print("\nFilter {}".format(filt))
            # Assume that the 'V' filter was used.
            stand_mag, stand_col, instZA_mag, z = extractData(filters, filt)
            x_fit, y_fit = stand_col, (stand_mag - instZA_mag)
            # Find slope and/or intersection of linear regression.
            st_accpt, st_rjct, m, c, R2, RMSE = regressRjctOutliers(
                x_fit.flatten(), y_fit.flatten(), z)
            filt_col_data.append([filt, st_accpt, st_rjct, c, m, R2, RMSE])
        else:
            print("\nFilter missing {}".format(filt))

    return filt_col_data


def writeTransfCoeffs(mypath, filt_col_data):
    """
    """
    f_path = join(
        mypath.replace('tasks', 'output/standards'), 'fit_coeffs.dat')
    with open(f_path, 'w') as f:
        hdr = """#\n# Instrumental zero airmass magnitude.\n# m_I^(0A) =""" +\
            """ m_I  - K * X  ; X: airmass, K: extinction coefficient\n""" +\
            """#\n# Standard magnitude.\n""" +\
            """# m_L = m_I^(0A) + c_1  + c_2 * col_L\n#\n# ID   N_a""" +\
            """    N_r       c_1      c_2      R^2     RMSE\n"""
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
        elif data[0] == "B":
            xlbl, ylbl = r'$(B-V)_L$', r'$(B_L-B^{0A}_{I})$'
        elif data[0] == "U":
            xlbl, ylbl = r'$(U-B)_L$', r'$(U_L-U^{0A}_{I})$'
        elif data[0] == "I":
            xlbl, ylbl = r'$(V-I)_L$', r'$(I_L-I^{0A}_{I})$'
        elif data[0] == "R":
            xlbl, ylbl = r'$(V-R)_L$', r'$(R_L-R^{0A}_{I})$'

        X_accpt, X_rjct, c, m, chi, RMSE = data[1:]
        x_accpt, y_accpt = X_accpt
        x_rjct, y_rjct = X_rjct

        x_r = np.arange(min(x_accpt), max(x_accpt) * 1.1, .01)
        predictions = m * x_r + c

        ax = fig.add_subplot(gs[0 + i * 2])
        sign = '+' if m >= 0. else '-'
        ax.set_title("{} = {:.3f} {} {:.3f} {}".format(
            ylbl, c, sign, abs(m), xlbl), fontsize=9)
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
        xlbl, x_scatter = ylbl, y_accpt
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

    print("Correct instrumental magnitudes for zero airmass.")
    filters = zeroAirmass(filters, pars['extin_coeffs'][0])

    print("Obtain transformation coefficients.")
    filt_col_data = fitTransfEquations(filters)

    print("Write 'fit_coeffs.dat' output file.")
    writeTransfCoeffs(pars['mypath'], filt_col_data)

    if pars['do_plots_E'] == 'y':
        make_plot(pars['mypath'], filt_col_data)

    print("\nFinished.")


if __name__ == '__main__':
    main()
