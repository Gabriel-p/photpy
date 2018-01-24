
import read_pars_file as rpf

from astropy.io import ascii
from astropy.table import Column
from os.path import join
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def in_params(apert_fl):
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()
    apert_fl = join(pars['mypath'].replace('tasks', 'output/'), apert_fl)
    return apert_fl


def readPhot(apert_fl):
    """
    """
    phot = ascii.read(apert_fl)

    filts = {'U': [], 'B': [], 'V': [], 'R': [], 'I': []}
    for f in filts.keys():
        # Select stars for this filter.
        mask = phot['Filt'] == f
        f_phot = phot[mask]
        # Remove nan magnitudes
        mask = ~np.isnan(f_phot['mag'])
        f_phot = f_phot[mask]
        # Combine 'Stnd_field' and 'ID' for unique ids
        c = Column(
            [_ + f_phot['ID'][i] for i, _ in enumerate(f_phot['Stnd_field'])])
        f_phot.add_column(c, name='Unq_ID')
        # Use set() to extract unique ids
        unq_stars = list(set(f_phot['Unq_ID']))

        for st in unq_stars:
            # Isolate this star in this filter.
            mask = f_phot['Unq_ID'] == st
            filts[f].append([f_phot[mask]['A'], f_phot[mask]['mag'], st])

    return filts


def reject_outliers(data, m=2.):
    """
    https://stackoverflow.com/a/16562028/1391441
    """
    # distance from data median
    d = np.abs(data - np.median(data))
    # median of distances to data median
    mdev = np.median(d)
    s = d / mdev if mdev else 0.
    mask = s < m
    return np.array(mask)


def main(apert_fl):
    """
    Rough sketch of script to determine extinction coefficients.
    """
    apert_fl = in_params(apert_fl)

    # Read aperture photometry for each matched star.
    filters = readPhot(apert_fl)

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 2)
    f_sbp = {'U': 0, 'B': 1, 'V': 2, 'I': 3}
    i = 0

    for filt, fdata in filters.iteritems():
        if fdata:
            n_stars = sum([1. for _ in fdata if len(_) > 1.])
            print("Filter {}, N_stars={:.0f}".format(filt, n_stars))
            K_median_outl = []
            outliers_lims = np.arange(1., 5., .5)
            for m_outliers in outliers_lims:
                plt.style.use('seaborn-darkgrid')
                slopes = []
                for st in fdata:

                    # If this star was detected more than once.
                    if len(st[1]) > 1:
                        mask = reject_outliers(st[1], m=m_outliers)
                        # If there's more than one star after outlier rejection
                        if sum(mask) > 1:
                            X, Y = st[0][mask], st[1][mask]
                        else:
                            X, Y = st[0], st[1]

                        m, c, r_value, p_value, m_err = linregress(X, Y)
                        slopes.append(m)
                        # print("  Star {}, K={:.3f}".format(st[2], m))

                    # # Plot
                    # plt.title("{}: K={:.3f}, m_err={:.3f}".format(
                    #     st[2], m, m_err))
                    # plt.xlabel('airmass')
                    # plt.ylabel('instrumental mag')
                    # plt.scatter(st[0][mask], st[1][mask])
                    # plt.scatter(st[0][~mask], st[1][~mask], c='r')
                    # fit_fn = np.poly1d([m, c])
                    # X = [min(st[0]), max(st[0])]
                    # plt.plot(X, fit_fn(X), '--k')
                    # plt.gca().invert_yaxis()
                    # plt.show()

                # Obtain median of extinction coefficients for this filter.
                median_K = np.nanmedian(slopes)
                print("  m_out: {:.1f} median K: {:.3f}".format(
                    m_outliers, median_K))
                # # Plot coefficients.
                # if not np.isnan(median_K):
                #     plt.title("Filter {}, med={:.3f}".format(filt, median_K))
                #     plt.ylim(0., 2 * median_K)
                #     slopes = np.array(slopes)
                #     slopes = slopes[~np.isnan(slopes)]
                #     plt.scatter(range(len(slopes)), slopes)
                #     plt.axhline(y=median_K, color='r')
                #     plt.show()

                K_median_outl.append(median_K)

            # Obtain median of extinction coefficients for this filter.
            K_median = np.nanmedian(K_median_outl)
            print(" K_{}: {:.3f}".format(filt, K_median))

            # Plot coefficients.
            fig.add_subplot(gs[f_sbp[filt]])
            i += 1
            plt.title(r"$K_{}={:.3f}\;(N_{{stars}}={})$".format(
                filt, K_median, int(n_stars)), y=.93)
            plt.xlabel("Outlier detection threshold")
            plt.ylabel("median(K_outl)")
            plt.ylim(0., 2 * K_median)
            K_median_outl = np.array(K_median_outl)
            K_median_outl = K_median_outl[~np.isnan(K_median_outl)]
            plt.scatter(outliers_lims, K_median_outl)
            plt.axhline(y=K_median, color='r')
            # plt.show()

    fig.tight_layout()
    plt.savefig('ext_coeffs.png', dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    """
    To run this script we need *at least* two observations for each filter
    with different airmass values.
    """

    # Name of the .apert file containing the matched stars and their aperture
    # photometry.
    apert_fl = 'k_coeffs.apert'
    main(apert_fl)

# https://arxiv.org/PS_cache/arxiv/pdf/0906/0906.3014v1.pdf
# v3 = +0.16, b3 = +0.25, i3 = +0.08, u3 = +0.45

# From an old 'noche.ans' file when BO14 was processed.
# U .49, B .27, V .12, I .02
