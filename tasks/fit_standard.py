
import read_pars_file as rpf

import os
from os.path import join
import numpy as np
from scipy.stats import linregress
from sklearn.linear_model import RANSACRegressor, TheilSenRegressor
from sklearn.metrics import r2_score, mean_squared_error
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table


def in_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    # This output path is assumed (hardcoded)
    out_path = pars['mypath'].replace('tasks', 'output/')
    apert_file = join(out_path, pars['apert_file_in'])
    file_out = join(out_path, pars['fit_file_out'])
    # Final image path+name
    file_out_ext = file_out.split('/')[-1].split('.')[-1]
    img_out = file_out.replace(file_out_ext, 'png')

    return pars, apert_file, file_out, img_out


def readData(apert_file):
    """
    """
    stndData = ascii.read(apert_file)
    # Group by filters
    f_grouped = stndData.group_by('Filt')

    return f_grouped


def zeroAirmass(phot_table, extin_coeffs):
    """
    Correct for airmass, i.e. instrumental magnitude at zero airmass.
    """
    # Identify correct index for this filter's extinction coefficient.
    ext, K_filts = [], {}
    for filt in phot_table['Filt']:
        f_idx = extin_coeffs.index(filt) + 1
        # Extinction coefficient.
        ext.append(float(extin_coeffs[f_idx]))
        K_filts[filt] = float(extin_coeffs[f_idx])

    # Obtain zero airmass instrumental magnitude for this filter.
    phot_table['ZA_mag'] = (
        phot_table['mag'] - (ext * phot_table['A'])) * u.mag

    return phot_table, K_filts


def rmse(targets, predictions):
    return np.sqrt(((predictions - targets) ** 2).mean())


def redchisq(ydata, ymod, sd=None, deg=2):
    """
    Returns the squared root of the reduced chi-square error statistic
    chisq/nu, where nu is the number of degrees of freedom.

    See http://en.wikipedia.org/wiki/Goodness_of_fit for reference.
    Source: http://astropython.blogspot.com.ar/2012/02/
    computing-chi-squared-and-reduced-chi.html

    ydata : data
    ymod : model evaluated at the same x points as ydata
    """
    if sd is None:
        chisq = np.nansum(((ydata - ymod)) ** 2)
    else:
        chisq = np.nansum(((ydata - ymod) / sd) ** 2)

    # Number of degrees of freedom assuming 'deg' free parameters
    non_nans = (~np.isnan(ydata + ymod)).sum()
    nu = non_nans - deg

    return np.sqrt(chisq / nu)


def linearFunc(p, x):
    """Linear model."""
    a, b = p
    return a * x + b


def linear_func(x, a, b):
    """Linear model."""
    return a * x + b


def residual(p, x, y):
    return y - linear_func(x, *p)


def distPoint2Line(m, c, x, y):
    """
    Distance from (x, y) point, to line with equation:

    y = m*x + c

    http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html

    The ravel() is necessary for the 'theilsen' mode.

    Returns the index of the element with largest distance.
    """
    d = abs(m * np.array(x).ravel() +
            c - np.array(y).ravel()) / np.sqrt(m ** 2 + 1.)
    return d.argmax()


def regressRjctOutliers(x, y, z, mode, R2_min, RMSE_max):
    """
    Perform a linear regression fit, rejecting outliers until the conditions
    of R2>R2_min and RMSE<RMSE_max are met.

    # TODO errors for fit coefficients
    # TODO finish reduced chi value

    """
    # Filter out nan values carried from ZA mag.
    x, y, z = x[~np.isnan(y)], y[~np.isnan(y)],\
        np.array(z)[~np.isnan(y)].tolist()
    x, y, z = x[~np.isnan(x)], y[~np.isnan(x)],\
        np.array(z)[~np.isnan(x)].tolist()

    if mode == 'linregress':
        m, c, r_value, p_value, std_err = linregress(x, y)
        predictions = m * np.array(x) + c
        red_chisq = redchisq(y, predictions)
        R2, RMSE = r_value**2, rmse(y, predictions)

        x_accpt, y_accpt, z_accpt = x.tolist()[:], y.tolist()[:], z[:]
        x_rjct, y_rjct, z_rjct = [], [], []
        while R2 <= R2_min or RMSE >= RMSE_max:
            idx_rm = distPoint2Line(m, c, x_accpt, y_accpt)
            x_rjct.append(x_accpt[idx_rm])
            y_rjct.append(y_accpt[idx_rm])
            z_rjct.append(z_accpt[idx_rm])
            # print("   Rjct: {}, ({:.3f}, {:.3f})".format(
            #     z_accpt[idx_rm], x_accpt[idx_rm], y_accpt[idx_rm]))
            del x_accpt[idx_rm]
            del y_accpt[idx_rm]
            del z_accpt[idx_rm]

            m, c, r_value, p_value, std_err = linregress(x_accpt, y_accpt)
            predictions = m * np.array(x_accpt) + c
            red_chisq = redchisq(y_accpt, predictions)
            R2, RMSE = r_value**2, rmse(y_accpt, predictions)

        # TODO other ways of fitting

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
        linModel = Model(linearFunc)
        # Create a RealData object using our initiated data from above.
        data = RealData(x_accpt, y_accpt)  # , sx=xerr, sy=yerr)
        # Set up ODR with the model and data.
        odr = ODR(data, linModel, beta0=[0., 1.])
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

        # TODO other ways of fitting

    elif mode == 'RANSAC':
        x, y = np.array(x).reshape(-1, 1), np.array(y).reshape(-1, 1)
        # Robustly fit linear model with RANSAC algorithm
        ransac = RANSACRegressor()
        R2, RMSE, _iter, iter_max = R2_min * .5, RMSE_max * 2., 0, 1000
        res, R2_old = [], R2
        while (R2 <= R2_min or RMSE >= RMSE_max) and _iter < iter_max:
            ransac.fit(x, y)
            # Estimated coefficients
            m, c = float(ransac.estimator_.coef_),\
                float(ransac.estimator_.intercept_)
            inlier_mask = ransac.inlier_mask_
            outlier_mask = np.logical_not(inlier_mask)
            x_accpt, y_accpt = x[inlier_mask], y[inlier_mask]
            x_rjct, y_rjct = x[outlier_mask], y[outlier_mask]
            y_predict = m * x_accpt + c
            red_chisq = redchisq(y_accpt, y_predict)
            R2 = r2_score(y_accpt, y_predict)
            RMSE = np.sqrt(mean_squared_error(y_accpt, y_predict))
            if R2 > R2_old:
                res = [m, c, x_accpt, y_accpt, x_rjct, y_rjct, R2, RMSE]
                R2_old = R2
            _iter += 1
        if _iter == iter_max:
            print(" WARNING: R2 and/or RMSE precision was not achieved.")
            m, c, x_accpt, y_accpt, x_rjct, y_rjct, R2, RMSE = res

    elif mode == 'theilsen':
        theilsen = TheilSenRegressor()
        x, y = np.array(x).reshape(-1, 1), y.ravel()
        x_accpt, y_accpt = x.tolist()[:], y.tolist()[:]
        x_rjct, y_rjct = [], []
        R2, RMSE, _iter, iter_max = R2_min * .5, RMSE_max * 2., 0, 1000
        res, R2_old = [], R2
        while (R2 <= R2_min or RMSE >= RMSE_max) and _iter < iter_max:
            theilsen.fit(x_accpt, y_accpt)
            # Estimated coefficients
            m, c = float(theilsen.coef_), float(theilsen.intercept_)
            # Remove point furthest away from fit line
            idx_rm = distPoint2Line(m, c, x_accpt, y_accpt)
            x_rjct.append(x_accpt[idx_rm][0])
            y_rjct.append(y_accpt[idx_rm])
            del x_accpt[idx_rm]
            del y_accpt[idx_rm]
            # Stats
            y_predict = m * np.array(x_accpt) + c
            red_chisq = redchisq(y_accpt, y_predict)
            R2 = r2_score(y_accpt, y_predict)
            RMSE = np.sqrt(mean_squared_error(y_accpt, y_predict))

        # Flatten list before storing.
        x_accpt = [_[0] for _ in x_accpt]

    print("m={:.3f}, c={:.3f}".format(m, c))
    # print("  red_chi: {:.3f}".format(red_chisq))
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
    stand_mag, stand_col, instZA_mag = f_dat['Mag_L'], f_dat['Col_L'],\
        f_dat['ZA_mag']
    # Flatten extra data used to identify rejected stars
    extra_data = [
        _[0] + ' ' + _[1] + '-' + _[2] + ' (' + str(_[3]) + ')' for _ in
        zip(*[f_dat['Filt'], f_dat['Stnd_field'], f_dat['ID'], f_dat['file']])]

    return stand_mag, stand_col, instZA_mag, extra_data


def fitTransfEquations(filters, mode, R2_min, RMSE_max):
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
                x_fit.flatten(), y_fit.flatten(), z, mode, R2_min, RMSE_max)
            filt_col_data.append([filt, st_accpt, st_rjct, c, m, R2, RMSE])
        else:
            print("\nFilter missing {}".format(filt))

    return filt_col_data


def writeTransfCoeffs(file_out, filt_col_data, K_filts):
    """
    """
    with open(file_out, 'w') as f:
        hdr = """#\n# Instrumental zero airmass magnitude.\n# m_I^(0A) =""" +\
            """ m_I  - K * X  ; X: airmass, K: extinction coefficient\n""" +\
            """#\n# Standard magnitude.\n""" +\
            """# m_L = m_I  - K * X + c_1  + c_2 * col_L\n#\n"""
        f.write(hdr)

    with open(file_out, mode='a') as f:
        t_data = []
        for data in filt_col_data:
            filt, accpt, rjct, c, m, r2, RMSE = data
            K = K_filts[filt]
            t_data.append([
                filt, len(accpt[0]), len(rjct[0]), K, c, m, r2, RMSE])

        tt = Table(
            zip(*t_data), names=(
                'ID', 'N_a', 'N_r', 'K', 'c_1', 'c_2', 'R^2', 'RMSE'))
        f.seek(0, os.SEEK_END)
        ascii.write(
            tt, f, format='fixed_width', delimiter='',
            formats={
                'K': '%8.3f', 'c_1': '%8.3f', 'c_2': '%8.3f', 'R^2': '%8.3f',
                'RMSE': '%8.3f'})


def make_plot(img_out, filt_col_data, K_filts):
    """
    """
    print("\nPlotting.")
    # print(plt.style.available)
    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(10, 25))
    gs = gridspec.GridSpec(5, 2)

    i = 0
    for data in filt_col_data:

        filt, X_accpt, X_rjct, c, m, r2, RMSE = data
        x_accpt, y_accpt = X_accpt
        x_rjct, y_rjct = X_rjct

        if filt == "V":
            xlbl, ylbl, K = r'$(B-V)_L$', r'$(V_L-V^{0A}_{I})$', K_filts[filt]
        elif filt == "B":
            xlbl, ylbl, K = r'$(B-V)_L$', r'$(B_L-B^{0A}_{I})$', K_filts[filt]
        elif filt == "U":
            xlbl, ylbl, K = r'$(U-B)_L$', r'$(U_L-U^{0A}_{I})$', K_filts[filt]
        elif filt == "I":
            xlbl, ylbl, K = r'$(V-I)_L$', r'$(I_L-I^{0A}_{I})$', K_filts[filt]
        elif filt == "R":
            xlbl, ylbl, K = r'$(V-R)_L$', r'$(R_L-R^{0A}_{I})$', K_filts[filt]

        x_r = np.arange(min(x_accpt), max(x_accpt) * 1.1, .01)
        predictions = m * x_r + c

        ax = fig.add_subplot(gs[0 + i * 2])
        sign = '+' if m >= 0. else '-'
        ax.set_title("{} = {:.3f} {} {:.3f} {} ; K={:.3f}".format(
            ylbl, c, sign, abs(m), xlbl, K), fontsize=9)
        plt.xlabel(xlbl)
        plt.ylabel(ylbl)
        ax.plot(x_r, predictions, 'b-', zorder=1)
        plt.scatter(x_accpt, y_accpt, c='g', zorder=4)
        plt.scatter(x_rjct, y_rjct, c='r', marker='x', zorder=4)

        # t1 = "m, c: {:.3f}, {:.3f}\n".format(m, c)
        t1 = r"$R^{{2}}={:.3f}$".format(r2) + "\n"
        t2 = r"$RMSE={:.3f}$".format(RMSE) + "\n"
        t3 = r"$N_a={:.0f},\;N_r={:.0f}$".format(len(x_accpt), len(x_rjct))
        txt = AnchoredText(t1 + t2 + t3, loc=2, prop=dict(size=9))
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
    plt.savefig(img_out, dpi=150, bbox_inches='tight')
    # Close to release memory.
    plt.clf()
    plt.close()


def main():
    """

    !!! --> Assumes that the V and B filters are present. <-- !!!

    """
    pars, apert_file, file_out, img_out = in_params()

    filters = readData(apert_file)

    print("Correct instrumental magnitudes for zero airmass.")
    filters, K_filts = zeroAirmass(filters, pars['extin_coeffs'][0])

    mode = pars['fit_mode']
    print("Obtain transformation coefficients ({}).".format(mode))
    R2_min, RMSE_max = map(float, [pars['R^2_min'], pars['RMSE_max']])
    filt_col_data = fitTransfEquations(filters, mode, R2_min, RMSE_max)

    print("\nWrite output file.")
    writeTransfCoeffs(file_out, filt_col_data, K_filts)

    if pars['do_plots_E'] == 'y':
        make_plot(img_out, filt_col_data, K_filts)

    print("\nFinished.")


if __name__ == '__main__':
    main()
