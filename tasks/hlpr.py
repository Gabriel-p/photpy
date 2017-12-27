
import numpy as np
from astropy.stats import biweight_location, biweight_midvariance
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from photutils import DAOStarFinder


def bckg_data(hdr, hdu_data, gain_key, rdnoise_key, method):
    """
    Estimate sky's mean, median, and standard deviation.
    Two methods available:

    * Bivariate weight (fast & approximated)
    * Sigma clipping (slow & accurate)

    Usually BW is close enough to the SC values to justify using it.
    """
    gain = hdr[gain_key]
    rdnoise = hdr[rdnoise_key]
    print("Gain, Rdnoise: {}, {}".format(gain, rdnoise))
    print("\nEstimating background mean, median, and STDDEV.")

    if method == 'BW':
        sky_mean = biweight_location(hdu_data)
        sky_median = np.median(hdu_data)
        sky_std = biweight_midvariance(hdu_data)
    elif method == 'SC':
        sky_mean, sky_median, sky_std = sigma_clipped_stats(
            hdu_data, sigma=3.0, iters=5)

    print("Mean, median, and STDDEV: {:.2f}, {:.2f}, {:.2f}".format(
        sky_mean, sky_median, sky_std))
    print("STDDEV estimated from GAIN, RDNOISE, and median: {:.2f}".format(
        np.sqrt(sky_median * gain + rdnoise ** 2) / gain))

    return sky_mean, sky_median, sky_std


def st_fwhm_select(
        dmax, max_stars, thresh_level, fwhm_init, std, hdu_data):
    """
    Find stars in image with DAOStarFinder, filter out saturated stars (with
    peak values greater than the maximum accepted in 'dmax'), order by the
    largest (most negative) magnitude, and select a maximum of 'max_stars'.
    """
    print("\nFinding stars.")
    thresh = thresh_level * std
    print('Threshold, initial FWHM: {:.1f}, {:.1f}'.format(thresh, fwhm_init))
    stfind = DAOStarFinder(threshold=thresh, fwhm=fwhm_init)
    sources = stfind(hdu_data)
    print("Sources found: {}".format(len(sources)))
    mask = sources['peak'] < dmax
    sour_no_satur = sources[mask]
    print("Non-saturated sources found: {}".format(len(sour_no_satur)))
    sour_no_satur.sort('mag')
    psf_select = sour_no_satur[:int(max_stars)]

    return psf_select, len(sources), len(sour_no_satur)


def psf_filter(fwhm_min, ellip_max, psf_data, out_f=2.):
    """
    Reject stars with a large ellipticity (> ellip_max), a very low
    FWHM (<fwhm_min), or large FWHMs (outliers).
    """
    fwhm_min_rjct = psf_data[psf_data['FWHM'] <= fwhm_min]
    print("Stars with FWHM<{:.1f} rejected: {}".format(
        fwhm_min, len(fwhm_min_rjct)))
    fwhm_min_accpt = psf_data[psf_data['FWHM'] > fwhm_min]
    fwhm_estim = fwhm_min_accpt[fwhm_min_accpt['Ellip'] <= ellip_max]
    ellip_rjct = fwhm_min_accpt[fwhm_min_accpt['Ellip'] > ellip_max]
    print("Stars with ellip>{:.1f} rejected: {}".format(
        ellip_max, len(ellip_rjct)))

    # Remove duplicates, if any.
    fwhm_estim_no_dups = list(set(np.array(fwhm_estim).tolist()))
    print("Duplicated stars rejected: {}".format(
        len(fwhm_estim) - len(fwhm_estim_no_dups)))

    fwhm_accptd, fwhm_outl = [], []
    if fwhm_estim_no_dups:
        fwhm_median = np.median(zip(*fwhm_estim_no_dups)[2])
        for st in fwhm_estim_no_dups:
            if st[2] <= out_f * fwhm_median:
                fwhm_accptd.append(st)
            else:
                fwhm_outl.append(st)
        print("Outliers with FWHM>{:.2f} rejected: {}".format(
            out_f * fwhm_median, len(fwhm_outl)))

    if fwhm_accptd:
        print("\nFinal number of accepted stars: {}".format(len(fwhm_accptd)))
        fwhm_mean, fwhm_std = np.mean(zip(*fwhm_accptd)[2]),\
            np.std(zip(*fwhm_accptd)[2])
        print("Mean FWHM +/- std: {:.2f}, {:.2f}".format(fwhm_mean, fwhm_std))
    else:
        print("  WARNING: all stars were rejected.\n"
              "  Could not obtain a mean FWHM.")
        fwhm_mean, fwhm_std = 0., 0.

    return fwhm_min_rjct, ellip_rjct, fwhm_accptd, fwhm_mean, fwhm_std,\
        fwhm_outl


def rmNaNrows(tab, not_cols):
    """
    Remove from 'tab' all those rows that contain only NaN values.

    Parameters
    ----------
    tab : class astropy.table
        All cross-matched stars for a given filter.
    not_cols : int
        Leave out this many columns when searching for all NaN values in a row.
        This is used because otherwise the 'ID' column (for example) would
        count towards that row containing one non-NaN value.

    Returns
    -------
    tab : class astropy.table
        Same table minus rows with all NaN values.

    """
    # Convert to pandas dataframe.
    tab_df = tab.to_pandas()
    # Find rows with *all* nan values (~ means 'not').
    nan_idx = ~tab_df[tab.keys()[not_cols:]].isnull().all(1)
    # Filter out all nan rows and transform back to Table.
    tab = Table.from_pandas(tab_df[nan_idx])

    return tab
