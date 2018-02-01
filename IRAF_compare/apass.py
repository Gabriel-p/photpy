
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table


def photRead(final_phot, col_IDs):
    """
    Select a file with photometry to read and compare with APASS.
    """
    # Final calibrated photometry
    print("Read final photometry")
    phot = ascii.read(final_phot, fill_values=('INDEF', np.nan))

    # plt.style.use('seaborn-darkgrid')

    # plt.subplot(131)
    # plt.xlabel("(B-V)")
    # plt.ylabel("V")
    # plt.scatter(phot['col6'], phot['col4'], s=5)
    # plt.gca().invert_yaxis()

    # plt.subplot(132)
    # plt.xlabel("(V-I)")
    # plt.ylabel("V")
    # plt.scatter(phot['col10'], phot['col4'], s=5)
    # plt.gca().invert_yaxis()

    # plt.subplot(133)
    # plt.xlabel("(B-V)")
    # plt.ylabel("U-B")
    # plt.scatter(phot['col6'], phot['col8'], s=5)
    # plt.gca().invert_yaxis()

    # plt.show()
    id_x, id_y, id_v, id_bv = col_IDs
    x_p, y_p, v_p, bv_p = phot[id_x], phot[id_y], phot[id_v], phot[id_bv]
    b_p = bv_p + v_p

    return x_p, y_p, v_p, bv_p, b_p


def astrometryFeed(f_id, x_p, y_p, v_p, regs_filt):
    """
    Create file with the proper format to feed astrometry.net
    """
    t = Table([x_p, y_p, v_p], names=('x', 'y', 'V'))
    t.sort('V')

    # Define central region limits.
    xmin, xmax, ymin, ymax = regs_filt
    mask = [xmin < t['x'], t['x'] < xmax, ymin < t['y'], t['y'] < ymax]
    total_mask = reduce(np.logical_and, mask)
    xm, ym = t['x'][total_mask], t['y'][total_mask]

    ascii.write(
        [xm, ym], 'output/' + f_id + "_astrometry.dat",
        delimiter=' ', format='fixed_width_no_header', overwrite=True)


def apassAstroRead(astro_cross, apass_reg):
    """
    Read APASS data.
    """
    apass = ascii.read(apass_reg, fill_values=('NA', np.nan))

    # astrometry.net cross-matched data
    hdul = fits.open(astro_cross)
    cr_m_data = hdul[1].data

    return apass, cr_m_data


def px2Eq(x_p, y_p, cr_m_data):
    """
    Transform pixels to (ra, dec) using the correlated astrometry.net file.
    """
    m_ra, h_ra = curve_fit(f, cr_m_data['field_x'], cr_m_data['field_ra'])[0]
    print("RA transf: ra = {:.5f} * x + {:.5f}".format(m_ra, h_ra))
    m_de, h_de = curve_fit(f, cr_m_data['field_y'], cr_m_data['field_dec'])[0]
    print("DEC transf: dec = {:.5f} * y + {:.5f}".format(m_de, h_de))
    x_p, y_p = m_ra * x_p + h_ra, m_de * y_p + h_de
    # plt.subplot(121)
    # plt.scatter(cr_m_data['field_x'], cr_m_data['field_ra'])
    # plt.plot(
    #     [min(cr_m_data['field_x']), max(cr_m_data['field_x'])],
    #     [f(min(cr_m_data['field_x']), m_ra, h_ra),
    #      f(max(cr_m_data['field_x']), m_ra, h_ra)], c='r')
    # plt.subplot(122)
    # plt.scatter(cr_m_data['field_y'], cr_m_data['field_dec'])
    # plt.plot(
    #     [min(cr_m_data['field_y']), max(cr_m_data['field_y'])],
    #     [f(min(cr_m_data['field_y']), m_de, h_de),
    #      f(max(cr_m_data['field_y']), m_de, h_de)], c='r')
    # plt.show()

    return x_p, y_p


def centerFilter(x_p, y_p, v_p, b_p, bv_p, apass, mag_max, mag_min):
    """
    Center APASS frame, filter data according to V range.
    """
    # Center frame for APASS data with proper range.
    xmin, xmax, ymin, ymax = x_p.min(), x_p.max(), y_p.min(), y_p.max()
    ra_c, de_c = .5 * (xmin + xmax), .5 * (ymin + ymax)
    ra_l, de_l = .5 * (xmax - xmin), .5 * (ymax - ymin)

    # Filter APASS frame to match the observed frame.
    mask = [apass['radeg'] < ra_c + ra_l, ra_c - ra_l < apass['radeg'],
            apass['decdeg'] < de_c + de_l, apass['decdeg'] > de_c - de_l,
            mag_min < apass['Johnson_V'], apass['Johnson_V'] < mag_max]
    total_mask = reduce(np.logical_and, mask)
    ra_apass = apass['radeg'][total_mask]
    deg_apass = apass['decdeg'][total_mask]
    V_apass = apass['Johnson_V'][total_mask]
    B_apass = apass['Johnson_B'][total_mask]
    BV_apass = B_apass - V_apass
    x, y = ra_apass, deg_apass
    print("Max APASS V: {:.1f}".format(max(V_apass)))

    # Filter observed data to the fixed magnitude range.
    mask = [mag_min < v_p, v_p < mag_max]  # , ev_p < .03
    mask = reduce(np.logical_and, mask)
    print("\nMag limit for myphot: {}".format(mag_max))
    x_i, y_i, v_i = x_p[mask], y_p[mask], v_p[mask]
    b_i, bv_i = b_p[mask], bv_p[mask]

    print("APASS stars: {}, myphot stars: {}".format(len(x), len(x_i)))

    return x, y, x_i, y_i, V_apass, B_apass, BV_apass, v_i, b_i, bv_i


def closestStar(x_fr1, y_fr1, x_fr2, y_fr2):
    """
    For every star in fr1, find the closest star in fr2.

    Parameters
    ----------
    x_fr1 : list
       x coordinates for stars in the reference frame.
    y_fr1 : list
       y coordinates for stars in the reference frame.
    x_fr2 : list
       x coordinates for stars in the processed frame.
    y_fr2 : list
       y coordinates for stars in the processed frame.

    Returns
    -------
    min_dist_idx : numpy array
        Index to the processed star closest to the reference star, for each
        reference star:
        * fr2[min_dist_idx[i]]: closest star in fr2 to the ith star in fr1.
        Also the index of the minimum distance in dist[i], i.e.: distance to
        the closest processed star to the ith reference star:
        * dist[i][min_dist_idx[i]]: distance between these two stars.
    min_dists : list
        Minimum distance for each star in the reference frame to a star in the
        processed frame.

    Notes
    -----
    len(fr1) = len(dist) = len(min_dist_idx)

    """
    fr1 = np.array(zip(*[x_fr1, y_fr1]))
    fr2 = np.array(zip(*[x_fr2, y_fr2]))
    min_dists, min_dist_idx = cKDTree(fr2).query(fr1, 1)

    return min_dist_idx, min_dists


def matchStars(
        x_apass, y_apass, x_iraf, y_iraf, V_apass, B_apass, BV_apass, v_iraf,
        b_iraf, bv_iraf):
    """
    """
    min_dist_idx, min_dists = closestStar(x_apass, y_apass, x_iraf, y_iraf)

    print("Match tolerance: {} arcsec".format(N_tol))
    rad = .00028 * N_tol
    x_a, y_a, x_i, y_i = [], [], [], []
    V_a_f, B_a_f, BV_a_f, V_i_f, B_i_f, BV_i_f = [], [], [], [], [], []
    for st1_i, st2_i in enumerate(min_dist_idx):
        d = min_dists[st1_i]
        if d < rad:  # and abs(V_apass[st1_i] - v_i[st2_i]) < .5:
            # print("St1, St2: d={:.5f}".format(d))
            # print(" ({:.2f}, {:.2f}) ; ({:.2f}, {:.2f})".format(
            #     x[st1_i], y[st1_i], x_i[st2_i], y_i[st2_i]))
            # print(" V_1={:.2f}, V_2={:.2f}".format(
            #     V_apass[st1_i], v_i[st2_i]))
            x_a.append(x_apass[st1_i])
            y_a.append(y_apass[st1_i])
            V_a_f.append(V_apass[st1_i])
            B_a_f.append(B_apass[st1_i])
            BV_a_f.append(BV_apass[st1_i])

            x_i.append(x_iraf[st2_i])
            y_i.append(y_iraf[st2_i])
            V_i_f.append(v_iraf[st2_i])
            B_i_f.append(b_iraf[st2_i])
            BV_i_f.append(bv_iraf[st2_i])

    print("Matched stars: {}".format(len(x_apass)))

    V_a_f, B_a_f, BV_a_f, V_i_f, B_i_f, BV_i_f =\
        np.array(V_a_f), np.array(B_a_f), np.array(BV_a_f),\
        np.array(V_i_f), np.array(B_i_f), np.array(BV_i_f)
    return x_a, y_a, x_i, y_i, V_a_f, B_a_f, BV_a_f, V_i_f, B_i_f, BV_i_f


def f(x, m, h):
    """'straight line"""
    return m * x + h


def star_size(mag, N=None, min_m=None):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    # Scale factor.
    if N is None:
        N = len(mag)
    if min_m is None:
        min_m = np.nanmin(mag)
        # print("min mag used: {}".format(min_m))
    factor = 500. * (1 - 1 / (1 + 150 / N ** 0.85))
    return 0.1 + factor * 10 ** ((np.array(mag) - min_m) / -2.5)


def makePlot(
    f_id, V_min, V_max, N_tol, x, y, V_apass, x_i, y_i, v_i,
    x1, y1, x2, y2, V_a_f, B_a_f, BV_a_f, V_i_f, B_i_f, BV_i_f):
    """
    """
    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(12, 18)

    plt.subplot(gs[0:6, 0:6])
    plt.title(r"APASS ($N$={}, $rad_{{match}}$={} arcsec)".format(
        len(x), N_tol))
    plt.gca().invert_xaxis()
    plt.xlabel("ra")
    plt.ylabel("dec")
    plt.scatter(x, y, s=star_size(V_apass), c='r')
    plt.scatter(x1, y1, s=star_size(V_a_f))

    plt.subplot(gs[0:6, 6:12])
    plt.title(r"IRAF (N={}, $N_{{matched}}$={})".format(
        len(x_i), len(x1)))
    plt.gca().invert_xaxis()
    plt.xlabel("ra")
    plt.ylabel("dec")
    plt.scatter(x_i, y_i, s=star_size(v_i), c='r')
    plt.scatter(x2, y2, s=star_size(V_i_f))

    plt.subplot(gs[0:3, 12:18])
    plt.title("{:.1f} < V < {:.1f}".format(V_min, V_max))
    plt.xlabel(r"$V_{{APASS}}$")
    plt.ylabel(r"$V_{{IRAF}}$")
    plt.scatter(V_a_f, V_i_f, s=4)
    plt.plot([min(V_a_f), max(V_a_f)], [min(V_a_f), max(V_a_f)], c='r')
    plt.xlim(min(V_a_f), max(V_a_f))
    plt.ylim(min(V_a_f), max(V_a_f))

    plt.subplot(gs[3:6, 12:18])
    plt.xlabel(r"$B_{{APASS}}$")
    plt.ylabel(r"$B_{{IRAF}}$")
    plt.scatter(B_a_f, B_i_f, s=4)
    plt.plot([min(B_a_f), max(B_a_f)], [min(B_a_f), max(B_a_f)], c='r')
    plt.xlim(min(B_a_f), max(B_a_f))
    plt.ylim(min(B_a_f), max(B_a_f))

    Vmed, Vmean = np.nanmedian(V_a_f - V_i_f), np.nanmean(V_a_f - V_i_f)
    print("median(V_APASS-V_IRAF): {:.4f}".format(Vmed))
    print("mean(V_APASS-V_IRAF): {:.4f}".format(Vmean))
    plt.subplot(gs[6:12, 0:6])
    plt.ylim(-.5, .5)
    plt.title("V median diff: {:.4f}".format(Vmed))
    plt.xlabel(r"$V_{{APASS}}$")
    plt.ylabel(r"$V_{{APASS}}-V_{{IRAF}}$")
    plt.scatter(V_a_f, V_a_f - V_i_f, s=4)
    plt.axhline(y=Vmed, c='r')
    plt.axhline(y=Vmean, ls='--', c='g')

    Bmed, Bmean = np.nanmedian(B_a_f - B_i_f), np.nanmean(B_a_f - B_i_f)
    print("median(B_APASS-B_IRAF): {:.4f}".format(Bmed))
    print("mean(B_APASS-B_IRAF): {:.4f}".format(Bmean))
    plt.subplot(gs[6:12, 6:12])
    plt.ylim(-.5, .5)
    plt.title("B median diff: {:.4f}".format(Bmed))
    plt.xlabel(r"$B_{{APASS}}$")
    plt.ylabel(r"$B_{{APASS}}-B_{{IRAF}}$")
    plt.scatter(B_a_f, B_a_f - B_i_f, s=4)
    plt.axhline(y=Bmed, c='r')
    plt.axhline(y=Bmean, ls='--', c='g')

    plt.subplot(gs[6:12, 12:18])
    BVmed = np.nanmedian(BV_a_f - BV_i_f)
    plt.title("BV median diff: {:.4f}".format(BVmed))
    plt.scatter(BV_a_f, V_a_f, s=7, label="APASS")
    plt.scatter(BV_i_f, V_i_f, s=7, label="IRAF")
    plt.xlim(min(BV_a_f) - .3, max(BV_a_f) + .3)
    plt.ylim(min(V_a_f) - .5, max(V_a_f) + .25)
    plt.xlabel(r"$(B-V)$")
    plt.ylabel(r"$V$")
    plt.gca().invert_yaxis()
    plt.legend()

    fig.tight_layout()
    plt.savefig(
        'output/apass_' + f_id + '.png', dpi=300, bbox_inches='tight')


def main(
    f_id, astro_cross, apass_reg, final_phot, col_IDs, V_min, V_max, N_tol,
    astrom_gen, regs_filt):
    """
    """
    # Read cluster photometry.
    x_p, y_p, v_p, bv_p, b_p = photRead(final_phot, col_IDs)

    if astrom_gen is True:
        # astrometry.net file
        astrometryFeed(f_id, x_p, y_p, v_p, regs_filt)
        return

    # Read APASS data and astrometry.net correlated coordinates.
    apass, cr_m_data = apassAstroRead(astro_cross, apass_reg)

    # Pixels to RA,DEC
    x_p, y_p = px2Eq(x_p, y_p, cr_m_data)

    # Center APASS and filter data
    x_apass, y_apass, x_iraf, y_iraf, V_apass, B_apass, BV_apass,\
        v_iraf, b_iraf, bv_iraf =\
        centerFilter(x_p, y_p, v_p, b_p, bv_p, apass, V_max, V_min)

    # Find stars within match tolerance.
    x_a, y_a, x_i, y_i, V_a_f, B_a_f, BV_a_f, V_i_f, B_i_f, BV_i_f =\
        matchStars(
            x_apass, y_apass, x_iraf, y_iraf, V_apass, B_apass, BV_apass,
            v_iraf, b_iraf, bv_iraf)

    makePlot(
        f_id, V_min, V_max, N_tol, x_apass, y_apass, V_apass,
        x_iraf, y_iraf, v_iraf, x_a, y_a, x_i, y_i, V_a_f, B_a_f, BV_a_f,
        V_i_f, B_i_f, BV_i_f)


if __name__ == '__main__':

    # Identify cluster to process.
    clusters = ['rup44']

    for f_id in clusters:
        # Path to astrometry-net cross-matched file.
        astro_cross = 'input/' + f_id + "_corr.fits"

        # Path to APASS region
        apass_reg = 'input/' + f_id + "_apass.csv"

        # Path to final photometry file.
        final_phot = 'input/' + f_id + "_final.dat"
        # Column IDs in 'final_phot' file: x, y, V, BV
        col_IDs = ['col2', 'col3', 'col4', 'col6']

        # Magnitude limit for the photometry.
        V_min, V_max = 7., 16.
        # Tolerance in arcsec for the cross-match
        N_tol = 15

        # Create file to feed astrometry.net?
        astrom_gen = True
        # xmin, xmax, ymin, ymax
        regs_filt = [1400., 2800., 1400., 2800.]

        main(
            f_id, astro_cross, apass_reg, final_phot, col_IDs,
            V_min, V_max, N_tol, astrom_gen, regs_filt)
