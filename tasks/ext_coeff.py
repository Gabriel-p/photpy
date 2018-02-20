
import read_pars_file as rpf

from astropy.io import ascii
from astropy.table import Column
from os.path import join
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def in_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()
    out_path = pars['mypath'].replace('tasks', 'output/')

    # Name of the .apert file containing the matched stars and their aperture
    # photometry.
    apert_fl = join(out_path, pars['apert_file_ext'])

    return pars, apert_fl, out_path


def readPhot(apert_fl):
    """
    Read input .apert file with aperture photometry.
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
            # Isolate this star in this filter and store.
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


def makePlot(filters, out_path):
    """
    """

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 2)
    plt.style.use('seaborn-darkgrid')

    f_sbp = {'U': 0, 'B': 1, 'V': 2, 'I': 3}
    for filt, data in filters.iteritems():
        n_stars, outliers_lims, K_median_outl, K_median, filt_slopes = data
        fig.add_subplot(gs[f_sbp[filt]])

        plt.title(r"$K_{}={:.3f}\;(N_{{stars}}={})$".format(
            filt, K_median, int(n_stars)), y=.93)

        if len(outliers_lims) > 1:
            fig.suptitle('Outlier rejection')
            # Plot coefficients.
            plt.xlabel("Outlier detection threshold")
            plt.ylabel(r"$K_{median}^{outl}$")
            plt.ylim(0., 2 * K_median)
            K_median_outl = np.array(K_median_outl)
            K_median_outl = K_median_outl[~np.isnan(K_median_outl)]
            plt.scatter(outliers_lims, K_median_outl)
            plt.axhline(y=K_median, color='r')

            f = join(out_path, 'ext_coeffs.png')
        else:
            fig.suptitle('No outlier rejection')
            plt.xlabel("Unique stars")
            plt.ylabel(r"$K_{median}$")
            std = np.std(filt_slopes[filt])
            plt.ylim(K_median - std, K_median + std)
            plt.scatter(range(len(filt_slopes[filt])), filt_slopes[filt])
            plt.axhline(y=K_median, color='r')

            f = join(out_path, 'ext_coeffs_no_outlier.png')

    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    plt.savefig(f, dpi=150, bbox_inches='tight')


def main():
    """
    Determine extinction coefficients from an existing aperture photometry
    file.
    """
    pars, apert_fl, out_path = in_params()

    # Read aperture photometry for each matched star.
    phot_data = readPhot(apert_fl)

    filters = {'U': [], 'B': [], 'V': [], 'R': [], 'I': []}
    for filt, fdata in phot_data.iteritems():
        if fdata:
            n_stars = sum([1. for _ in fdata if len(_) > 1.])
            print("Filter {}, N unique stars={:.0f}".format(filt, n_stars))

            if pars['outlier_reject'] == 'yes':
                outl_vals = map(float, pars['outliers_lims'][0])
                outliers_lims = np.arange(*outl_vals)
            else:
                outliers_lims = [np.inf]

            filt_slopes = {'U': [], 'B': [], 'V': [], 'R': [], 'I': []}
            K_median_outl = []
            for i, m_outliers in enumerate(outliers_lims):
                slopes = []
                for st in fdata:

                    # If this star was detected more than once.
                    if len(st[1]) > 1:

                        X, Y = st[0], st[1]
                        mask = reject_outliers(st[1], m=m_outliers)
                        # If there's more than one star after outlier
                        # rejection
                        if sum(mask) > 1:
                            X, Y = st[0][mask], st[1][mask]

                        m, c, r_value, p_value, m_err = linregress(X, Y)
                        slopes.append(m)
                        # print("  Star {}, K={:.3f}".format(st[2], m))
                    else:
                        if i == 0:
                            print("  Single detection for star: {}".format(
                                st[2]))

                    print("  star {}: K={:.3f}, m_err={:.3f}".format(
                        st[2], m, m_err))

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
                print(" m_out: {:.1f} median K: {:.3f}".format(
                    m_outliers, median_K))
                filt_slopes[filt] = slopes

                K_median_outl.append(median_K)

            # Obtain median of extinction coefficients for this filter.
            K_median = np.nanmedian(K_median_outl)
            print(" K_{}: {:.3f}".format(filt, K_median))
            # Store data for plotting.
            filters[filt] = [
                n_stars, outliers_lims, K_median_outl, K_median, filt_slopes]

    # Remove not observed filters from dictionary.
    filters = {k: v for k, v in filters.iteritems() if v}
    makePlot(filters, out_path)


if __name__ == '__main__':
    main()
