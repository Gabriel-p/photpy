
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import ascii
from astropy.io import fits


def photRead():
    """
    Select a file with photometry to read and compare with APASS.
    """

    # # .als file
    # N_als = '42'
    # print("{} .als.1 V file".format(N_als))
    # phot = ascii.read('2_PSF/stk_fcd00' + N_als + '.als.1', format='daophot')
    # x_p, y_p, v_p, = phot['XCENTER'], phot['YCENTER'], phot['MAG']
    # b_p = []

    # NEW .als file
    N_als = '42'
    print("NEW {} .als.1 V file".format(N_als))
    phot = ascii.read(
        '2_PSF/test/' + N_als + '/stk_fcd00' + N_als + '.als_NEW2.1',
        format='daophot')
    x_p, y_p, v_p = phot['XCENTER'], phot['YCENTER'], phot['MAG']
    # To zero airmass
    K, XV = .136, 1.078
    v_p = v_p - K * XV
    b_p = []

    # # .mag DAOM V file
    # print("DAOM V .mag file")
    # phot = ascii.read('3_DAOM/vfilter.mag', fill_values=('INDEF', np.nan))
    # x_p, y_p, v_p = phot['col2'], phot['col3'], phot['col4']
    # b_p = []

    # # .mag DAOM V file
    # print("DAOM V .mag file, no I filter")
    # phot = ascii.read(
    #     '3_DAOM/no_I_filt/vfilter.mag', fill_values=('INDEF', np.nan))
    # x_p, y_p, v_p = phot['col2'], phot['col3'], phot['col4']
    # b_p = []

    # # .mag DAOM .obs file
    # print("DAOM .obs file")
    # phot = ascii.read('3_DAOM/daom.obs', fill_values=('INDEF', np.nan))
    # x_p, y_p, v_p, b_p = phot['col2'], phot['col3'], phot['col4'], phot['col6']
    # v1, v2, v3, v4, b1, b2, b3, b4 = 1.606477, 0.136, 0.03963335, 0.03036864,\
    #     1.686928, 0.232, -0.1767316, 0.06384443
    # XV, XB = 1.075, 1.072
    # BV = ((b_p - v_p) - b1 + v1 - b2*XB + v2*XV)/(1.+b3-v3+b4*XB-v4*XV)
    # v_p = v_p - v1 - v2*XV - v3*BV - v4*BV*XV

    # # .mag DAOM B file
    # print("DAOM B .mag file")
    # phot = ascii.read('3_DAOM/bfilter.mag', fill_values=('INDEF', np.nan))
    # x_p, y_p, v_p = phot['col2'], phot['col3'], phot['col4']
    # b_p = []

    # # Final calibrated photometry
    # print("Final photometry")
    # phot = ascii.read(
    #     'phot_compare/BH73_IRAF_4_2_coeffs.dat', fill_values=('INDEF', np.nan))
    # x_p, y_p, v_p, bv_p = phot['x'], phot['y'], phot['V'], phot['BV']
    # b_p = bv_p + v_p

    # # Final calibrated photometry
    # print("Final photometry using BO14 coeffs")
    # phot = ascii.read(
    #     '4_INVERTFIT/bo14_ans/BH73_BO14_coeffs.dat',
    #     fill_values=('INDEF', np.nan))
    # x_p, y_p, v_p, bv_p = phot['x'], phot['y'], phot['V'], phot['BV']
    # b_p = bv_p + v_p

    # # Final calibrated photometry
    # print("Final photometry using all standard frames")
    # phot = ascii.read(
    #     '4_INVERTFIT/all_standards/bh73_final.dat',
    #     fill_values=('INDEF', np.nan))
    # x_p, y_p, v_p, bv_p = phot['x'], phot['y'], phot['V'], phot['BV']
    # b_p = bv_p + v_p

    return x_p, y_p, v_p, b_p


def transfABCD(x, y, ABCD):
    """
    Apply transformation equations to obtain new (xt, yt) coordinates.
    """
    A, B, C, D = ABCD
    xt = A + x * C + y * D
    yt = B + y * C + x * D

    return xt, yt


def solveABCD(x1, y1, x2, y2):
    """
    Obtain new A, B, C, D parameters, solving the linear equations:

    1*A + 0*B + x2*C + y2*D = x1
    0*A + 1*B + y2*C + x2*D = y1
    """
    N = len(x1)
    l1 = np.array([np.ones(N), np.zeros(N), x2, y2])
    l2 = np.array([np.zeros(N), np.ones(N), y2, x2])
    M1 = np.vstack([l1.T, l2.T])
    M2 = np.concatenate([x1, y1])
    A, B, C, D = np.linalg.lstsq(M1, M2)[0]
    print("New A,B,C,D parameters: "
          "({:.3f}, {:.3f}, {:.3f}, {:.3f})".format(A, B, C, D))

    return A, B, C, D


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
    mag_sel, N_tol, x, y, V_apass, x_i, y_i, v_i, x1, y1, mag1, x2, y2,
    mag2):
    """
    """
    mag1, mag2 = np.array(mag1), np.array(mag2)
    med = np.nanmedian(mag1 - mag2)
    print("median({}_APASS-V_DAOM): {:.4f}".format(mag_sel, med))
    print("mean({}_APASS-V_DAOM): {:.4f}".format(
        mag_sel, np.nanmean(mag1 - mag2)))

    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 2)

    # plt.subplot(221)
    fig.add_subplot(gs[0])
    plt.title(r"APASS {} ($N$={}, $rad_{{match}}$={} arcsec)".format(
        mag_sel, len(x), N_tol))
    plt.gca().invert_xaxis()
    plt.xlabel("ra")
    plt.ylabel("dec")
    plt.scatter(x, y, s=star_size(V_apass), c='r')
    plt.scatter(x1, y1, s=star_size(mag1))

    # plt.subplot(222)
    fig.add_subplot(gs[1])
    plt.title(r"DAOM {} ($N_{{matched}}$={})".format(mag_sel, len(x1)))
    plt.gca().invert_xaxis()
    plt.xlabel("ra")
    plt.ylabel("dec")
    plt.scatter(x_i, y_i, s=star_size(v_i), c='r')
    plt.scatter(x2, y2, s=star_size(mag2))

    # plt.subplot(223)
    fig.add_subplot(gs[2])
    plt.ylim(med - 2. * med, med + 2. * med)
    # plt.ylim(-.5, .5)
    plt.title("Median diff: {:.4f}".format(med))
    plt.xlabel(r"${}_{{APASS}}$".format(mag_sel))
    plt.ylabel(r"${}_{{APASS}}-V_{{DAOM}}$".format(mag_sel))
    plt.scatter(mag1, mag1 - mag2, s=4)
    plt.axhline(y=med, c='r')

    # plt.subplot(224)
    fig.add_subplot(gs[3])
    plt.xlabel(r"${}_{{APASS}}$".format(mag_sel))
    plt.ylabel(r"${}_{{DAOM}}$".format(mag_sel))
    plt.scatter(mag1, mag2, s=4)
    plt.plot([0., 100.], [0., 100.], c='r')
    plt.xlim(9., 19.)
    plt.ylim(9., 19.)

    # plt.subplot(121)
    # plt.gca().invert_xaxis()
    # plt.scatter(x, y, s=star_size(V_apass))
    # plt.subplot(122)
    # plt.gca().invert_xaxis()
    # plt.scatter(x_i, y_i, s=star_size(np.array(v_i)))
    # plt.show()

    fig.tight_layout()
    plt.savefig(
        'output/apass_compare_' + mag_sel + '.png',
        dpi=300, bbox_inches='tight')


def main():
    """
    """
    input_apass_f = 'input/APASS/bh73_apass.csv'
    # Read APASS data.
    apass = ascii.read(input_apass_f, fill_values=('NA', np.nan))

    # Manually select values to filter the APASS frame to match the observed
    # frame.
    # Center
    # RA: 142.9833 , DEC: -50.2167
    mask = [apass['radeg'] < 143.35, 142.59 < apass['radeg'],
            apass['decdeg'] < -49.98, apass['decdeg'] > -50.47]
    total_mask = reduce(np.logical_and, mask)
    ra_apass = apass['radeg'][total_mask]
    deg_apass = apass['decdeg'][total_mask]
    V_apass = apass['Johnson_V'][total_mask]
    B_apass = apass['Johnson_B'][total_mask]
    # BV_apass = apass['BV'][total_mask]

    x, y = ra_apass, deg_apass

    # ra_cent = (max(ra_apass) + min(ra_apass)) / 2.
    # dec_cent = (max(deg_apass) + min(deg_apass)) / 2.
    # x = (np.array(ra_apass) - ra_cent) * np.cos(np.deg2rad(dec_cent))
    # y = (np.array(deg_apass) - dec_cent)

    # # Flip
    # x = -x + 1.
    # x = ((x - min(x)) / (max(x) - min(x))) * 4096
    # y = ((y - min(y)) / (max(y) - min(y))) * 4096

    x_p, y_p, v_p, b_p = photRead()
    # Save photometry to feed astrometry.net
    # from astropy.table import Table
    # t = Table([x_p, y_p, v_p])
    # t.sort('V')
    # x_p, y_p, v_p = t['x'], t['y'], t['V']
    # mask = [x_p < 2500., 1500. < x_p, y_p < 2500., y_p > 1500.]
    # total_mask = reduce(np.logical_and, mask)
    # xm, ym = x_p[total_mask], y_p[total_mask]
    # ascii.write([xm, ym], "bh73_astrometry.dat")

    # Transform pixels to (ra, dec)
    # astrometry.net cross-matched data
    hdul = fits.open('input/APASS/corr.fits')
    data = hdul[1].data
    m_ra, h_ra = curve_fit(f, data['field_x'], data['field_ra'])[0]
    print("RA transf: ra = {:.5f} * x + {:.5f}".format(m_ra, h_ra))
    m_de, h_de = curve_fit(f, data['field_y'], data['field_dec'])[0]
    print("DEC transf: dec = {:.5f} * y + {:.5f}".format(m_de, h_de))
    x_p, y_p = m_ra * x_p + h_ra, m_de * y_p + h_de

    # Select magnitude to analyze
    mag_sel = 'V'
    # Magnitude limit for the photometry.
    mag_l = 20.

    mask = v_p < mag_l
    print("\nMag limit for myphot: {}".format(mag_l))
    if mag_sel == 'V':
        x_i, y_i, v_i = x_p[mask], y_p[mask], v_p[mask]
    else:
        x_i, y_i, v_i, b_i = x_p[mask], y_p[mask], v_p[mask], b_p[mask]

    # # Initial translations
    # ABCD = (0., 0., 1., 0.)
    # x_i, y_i = transfABCD(x_i, y_i, ABCD)

    print("APASS stars: {}, myphot stars: {}".format(len(x), len(x_i)))

    min_dist_idx, min_dists = closestStar(x, y, x_i, y_i)

    N_tol = 7
    print("Match tolerance: {} arcsec".format(N_tol))
    rad = .00028 * N_tol
    x1, y1, x2, y2, mag1, mag2 = [], [], [], [], [], []
    for st1_i, st2_i in enumerate(min_dist_idx):
        d = min_dists[st1_i]
        if d < rad:  # and abs(V_apass[st1_i] - v_i[st2_i]) < .5:
            # print("St1, St2: d={:.5f}".format(d))
            # print(" ({:.2f}, {:.2f}) ; ({:.2f}, {:.2f})".format(
            #     x[st1_i], y[st1_i], x_i[st2_i], y_i[st2_i]))
            # print(" V_1={:.2f}, V_2={:.2f}".format(
            #     V_apass[st1_i], v_i[st2_i]))
            x1.append(x[st1_i])
            y1.append(y[st1_i])
            if mag_sel == 'V':
                mag1.append(V_apass[st1_i])
            else:
                mag1.append(B_apass[st1_i])

            x2.append(x_i[st2_i])
            y2.append(y_i[st2_i])
            if mag_sel == 'V':
                mag2.append(v_i[st2_i])
            else:
                mag2.append(b_i[st2_i])

    # solveABCD(x1, y1, x2, y2)

    # x1, y1, x2, y2 = [np.array(_) for _ in [x1, y1, x2, y2]]
    # print(np.median(x1 - x2), np.median(y1 - y2))

    print("Matched stars: {}".format(len(x1)))

    makePlot(
        mag_sel, N_tol, x, y, V_apass, x_i, y_i, v_i, x1, y1, mag1, x2, y2,
        mag2)


if __name__ == '__main__':
    main()
