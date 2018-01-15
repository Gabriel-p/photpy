
"""
FB post with lots of info:
https://www.facebook.com/groups/astro.r/permalink/1497636676998999/

One more:
https://www.facebook.com/groups/123898011017097/permalink/1783251948415020/

Bayesian fitting package:
https://github.com/dokester/BayesicFitting
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import Model, RealData, ODR
from scipy import random

# import random


# random.seed(9001)
# np.random.seed(117)


def linear_func(p, x):
    """Linear model."""
    a, b = p
    return a * x + b


def getData(c=1.):
    # Initiate some data, giving some randomness using random.random().
    N = 5
    x = np.random.uniform(-.5, 3., N)  # color
    a, b = np.random.uniform(-.5, 2., 2)
    print("Real parameters (a, b): ({:.3f}, {:.3f})".format(a, b))
    y = linear_func((a, b), x) + np.random.uniform(-.1, .1, N)
    # np.array([i**2 + 7random.random() for i in x])
    xerr = c * np.random.uniform(.05, .1, N)
    yerr = c * np.random.uniform(.05, .1, N)

    return x, y, xerr, yerr, a, b


def fitModel(x, y, xerr, yerr, f_type, f_name):
    # Create a model for fitting.
    linear_model = Model(linear_func)
    # Create a RealData object using our initiated data from above.
    data = RealData(x, y, sx=xerr, sy=yerr)
    # Set up ODR with the model and data.
    odr = ODR(data, linear_model, beta0=[0., 1.])

    # Fit type
    print("\n* Type: {}".format(f_name))
    odr.set_job(fit_type=f_type)
    # Run the regression.
    out = odr.run()

    print("Reason: {}".format(out.stopreason[0]))
    # Estimated parameter values
    beta = out.beta
    print("Parameters (a, b): {}, {}".format(*beta))
    # Standard errors of the estimated parameters
    std = out.sd_beta
    print("sd_beta: {}, {}".format(*std))
    print("sd_beta * sqrt(N): {}, {}".format(*std * np.sqrt(len(x))))
    # v This is equivalent to out.sd_beta
    # print("sqrt(diag cov * res_var): {}, {}".format(
    #     *np.sqrt(np.diag(out.cov_beta * out.res_var))))

    # print("sd_beta / sqrt(res_var): {}, {}".format(
    #     *std / np.sqrt(out.res_var)))
    # Covariance matrix of the estimated parameters
    cov = out.cov_beta
    stddev = np.sqrt(np.diag(cov))
    print("sqrt(diag_cov): {}, {}".format(*stddev))

    return out


def bces(y1, y1err, y2, y2err, cerr):
    """
    BCES fitting

    Does the entire regression calculation for 4 slopes:
    OLS(Y|X), OLS(X|Y), bisector, orthogonal.
    Fitting form: Y=AX+B.

    Usage:

    >>> a,b,aerr,berr,covab=bces(x,xerr,y,yerr,cov)

    Output:

    - a,b : best-fit parameters a,b of the linear regression
    - aerr,berr : the standard deviations in a,b
    - covab : the covariance between a and b (e.g. for plotting confidence
      bands)

    Arguments:

    - x,y : data
    - xerr,yerr: measurement errors affecting x and y
    - cov : covariance between the measurement errors
    (all are arrays)

    v1 Mar 2012: ported from bces_regress.f. Added covariance output.
    Rodrigo Nemmen, http://goo.gl/8S1Oo
    """
    # Arrays holding the code main results for each method:
    # Elements: 0-Y|X, 1-X|Y, 2-bisector, 3-orthogonal
    a, b, avar, bvar, covarxiz =\
        np.zeros(4), np.zeros(4), np.zeros(4), np.zeros(4), np.zeros(4)
    # Lists holding the xi and zeta arrays for each method above
    xi, zeta = [], []

    # Calculate sigma's for datapoints using length of conf. intervals
    sig11var = np.mean(y1err**2)
    sig22var = np.mean(y2err**2)
    sig12var = np.mean(cerr)

    # Covariance of Y1 (X) and Y2 (Y)
    covar_y1y2 = np.mean((y1 - y1.mean()) * (y2 - y2.mean()))

    # Compute the regression slopes
    a[0] = (covar_y1y2 - sig12var) / (y1.var() - sig11var)  # Y|X
    a[1] = (y2.var() - sig22var) / (covar_y1y2 - sig12var)  # X|Y
    a[2] = (a[0] * a[1] - 1.0 + np.sqrt((1.0 + a[0] ** 2) *
            (1. + a[1] ** 2))) / (a[0] + a[1])  # bisector

    if covar_y1y2 < 0:
        sign = -1.
    else:
        sign = 1.
    a[3] = .5 * ((a[1] - (1. / a[0])) + sign * np.sqrt(4. +
                 (a[1] - (1. / a[0])) ** 2))  # orthogonal

    # Compute intercepts
    for i in range(4):
        b[i] = y2.mean() - a[i] * y1.mean()

    # Set up variables to calculate standard deviations of slope/intercept
    # Y|X
    xi.append(((y1 - y1.mean()) * (y2 - a[0] * y1 - b[0]) +
               a[0] * y1err ** 2) / (y1.var() - sig11var))
    # X|Y
    xi.append(
        ((y2 - y2.mean()) * (y2 - a[1] * y1 - b[1]) - y2err ** 2) / covar_y1y2)
    # bisector
    xi.append(
        xi[0] * (1. + a[1] ** 2) * a[2] /
        ((a[0] + a[1]) * np.sqrt((1. + a[0] ** 2) * (1. + a[1] ** 2))) +
        xi[1] * (1. + a[0] ** 2) * a[2] /
        ((a[0] + a[1]) * np.sqrt((1. + a[0] ** 2) * (1. + a[1] ** 2))))
    # orthogonal
    xi.append(
        xi[0] * a[3] / (a[0] ** 2 * np.sqrt(4. + (a[1] - 1. / a[0]) ** 2)) +
        xi[1] * a[3] / np.sqrt(4. + (a[1] - 1. / a[0]) ** 2))

    for i in range(4):
        zeta.append(y2 - a[i] * y1 - y1.mean() * xi[i])

    for i in range(4):
        # Calculate variance for all a and b
        avar[i] = xi[i].var() / xi[i].size
        bvar[i] = zeta[i].var() / zeta[i].size

        # Sample covariance obtained from xi and zeta (paragraph after equation
        # 15 in AB96)
        covarxiz[i] = np.mean(
            (xi[i] - xi[i].mean()) * (zeta[i] - zeta[i].mean()))

    # Covariance between a and b (equation after eq. 15 in AB96)
    covar_ab = covarxiz / y1.size

    return a, b, np.sqrt(avar), np.sqrt(bvar), covar_ab


def bootstrap(v):
    """
    Constructs Monte Carlo simulated data set using the
    Bootstrap algorithm.

    Usage:

    >> > bootstrap(x)

    where x is either an array or a list of arrays. If it is a
    list, the code returns the corresponding list of bootstrapped
    arrays assuming that the same position in these arrays map the
    same "physical" object.

    Rodrigo Nemmen, http://goo.gl/8S1Oo

    """
    if type(v) == list:
        vboot = []    # list of boostrapped arrays
        n = v[0].size
        iran = random.randint(0, n, n)    # Array of random indexes
        for x in v:
            vboot.append(x[iran])
    else:    # if v is an array, not a list of arrays
        n = v.size
        iran = random.randint(0, n, n)    # Array of random indexes
        vboot = v[iran]

    return vboot


def bcesboot(y1, y1err, y2, y2err, cerr, nsim=10000):
    """
    Does the BCES with bootstrapping.

    Usage:

    >> > a, b, aerr, berr, covab = bcesboot(x, xerr, y, yerr, cov, nsim)

    : param x, y: data
    : param xerr, yerr: measurement errors affecting x and y
    : param cov: covariance between the measurement errors(all are arrays)
    : param nsim: number of Monte Carlo simulations(bootstraps)

    : returns: a, b - - best-fit parameters a, b of the linear regression
    : returns: aerr, berr - - the standard deviations in a, b
    : returns: covab - - the covariance between a and b(e.g. for plotting
      confidence bands)

    .. note:: this method is definitely not nearly as fast as bces_regress.f.
      Needs to be optimized. Maybe adapt the fortran routine using f2python?

    v1 Mar 2012: ported from bces_regress.f. Added covariance output.
    Rodrigo Nemmen, http: // goo.gl/8S1Oo

    """

    # print("Bootstrapping progress:")
    print("* Type: BCES (N={})".format(nsim))

    """
    My convention for storing the results of the bces code below as
    matrixes for processing later are as follow:

          simulation\method  y|x x|y bisector orthogonal
              sim0           ...
    Am =      sim1           ...
              sim2           ...
              sim3           ...
    """
    for i in range(nsim):
        [y1sim, y1errsim, y2sim, y2errsim, cerrsim] = bootstrap(
            [y1, y1err, y2, y2err, cerr])

        asim, bsim, errasim, errbsim, covabsim = bces(
            y1sim, y1errsim, y2sim, y2errsim, cerrsim)

        if i == 0:
            # Initialize the matrixes
            am, bm = asim.copy(), bsim.copy()
        else:
            am = np.vstack((am, asim))
            bm = np.vstack((bm, bsim))

    # Bootstrapping results
    a = np.array([am[:, 0].mean(), am[:, 1].mean(),
                  am[:, 2].mean(), am[:, 3].mean()])
    b = np.array([bm[:, 0].mean(), bm[:, 1].mean(),
                  bm[:, 2].mean(), bm[:, 3].mean()])

    # Error from unbiased sample variances
    erra, errb, covab = np.zeros(4), np.zeros(4), np.zeros(4)
    for i in range(4):
        erra[i] = np.sqrt(
            1. / (nsim - 1) * (np.sum(am[:, i] ** 2) - nsim *
                               (am[:, i].mean()) ** 2))
        errb[i] = np.sqrt(
            1. / (nsim - 1) * (np.sum(bm[:, i] ** 2) - nsim *
                               (bm[:, i].mean()) ** 2))
        covab[i] = 1. / (nsim - 1) *\
            (np.sum(am[:, i] * bm[:, i]) - nsim * am[:, i].mean() *
             bm[:, i].mean())

    return a, b, erra, errb, covab


for c in np.arange(1., 10., 1.):
    print("\n c (errors length) = {}".format(c))

    x, y, xerr, yerr, a, b = getData(c)
    # Plot
    y_fit = linear_func((a, b), x)
    plt.errorbar(
        x, y, xerr=xerr, yerr=yerr, linestyle='None', marker='x')
    plt.plot(x, y_fit, label="Actual fit")

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.odr.ODR.set_job.html#scipy.odr.ODR.set_job
    # f_type
    # 0 -> explicit ODR
    # 1 -> implicit ODR
    # 2 -> ordinary least-squares

    # for f_type, f_name in enumerate(
    #         ["Explicit ODR", "Implicit ODR", "Ordinary least-squares"]):
    #     if f_type != 1:
    #         out = fitModel(x, y, xerr, yerr, f_type, f_name)
    #         # Plot
    #         y_fit = linear_func(out.beta, x)
    #         plt.errorbar(
    #             x, y, xerr=xerr, yerr=yerr, linestyle='None', marker='x')
    #         plt.plot(x, y_fit, label=f_name)

    # BCES
    N_b = 10000
    cov = np.zeros(len(x))
    a, b, erra, errb, covab = bcesboot(x, xerr, y, yerr, cov, N_b)
    # Selects the desired BCES method for plotting
    i = 0
    print("Parameters (a, b): {}, {}".format(a[i], b[i]))
    print("Standard deviations: {}, {}".format(erra[i], errb[i]))
    # print(covab[i])
    y_fit = linear_func([a[i], b[i]], x)
    plt.plot(x, y_fit, label="BCES (N_b={})".format(N_b))

    plt.legend()
    plt.show()
