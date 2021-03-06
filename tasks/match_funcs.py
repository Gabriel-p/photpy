
from hlpr import bckg_data, st_fwhm_select
from aperphot_standards import calibrate_magnitudes

import numpy as np
import itertools
import math
import functools
from scipy.spatial.distance import pdist, cdist
from scipy.spatial import cKDTree

from photutils.utils import cutout_footprint
from photutils import detect_threshold, detect_sources
from photutils import source_properties
from astropy.convolution import Gaussian2DKernel


def getTriangles(set_X, X_combs):
    """
    # TODO
    Inefficient way of obtaining the lengths of each triangle's side.
    Normalized so that the minimum length is 1 and sorted so that smaller
    sides are presented first.
    """
    triang, tr_not_scaled = [], []
    for p0, p1, p2 in X_combs:
        d1 = np.sqrt((set_X[p0][0] - set_X[p1][0]) ** 2 +
                     (set_X[p0][1] - set_X[p1][1]) ** 2)
        d2 = np.sqrt((set_X[p0][0] - set_X[p2][0]) ** 2 +
                     (set_X[p0][1] - set_X[p2][1]) ** 2)
        d3 = np.sqrt((set_X[p1][0] - set_X[p2][0]) ** 2 +
                     (set_X[p1][1] - set_X[p2][1]) ** 2)
        d_min = min(d1, d2, d3)
        if d_min != 0.:
            d_unsort = [d1 / d_min, d2 / d_min, d3 / d_min]
        else:
            print("WARNING: duplicate star found.")
            d_unsort = [d1, d2, d3]
        # Sorting is important so that the minimum length side is stored first.
        triang.append(sorted(d_unsort))
        # These are used to obtain the scale between frames.
        tr_not_scaled.append(sorted([d1, d2, d3]))

    return np.asarray(triang), tr_not_scaled


def indices_merged_arr_generic_using_cp(arr):
    """
    Based on cartesian_product
    http://stackoverflow.com/a/11146645/190597 (senderle)
    """
    shape = arr.shape
    arrays = [np.arange(s, dtype='int') for s in shape]
    broadcastable = np.ix_(*arrays)
    broadcasted = np.broadcast_arrays(*broadcastable)
    rows, cols = functools.reduce(np.multiply, broadcasted[0].shape),\
        len(broadcasted) + 1
    out = np.empty(rows * cols, dtype=arr.dtype)
    start, end = rows, 2 * rows
    for a in broadcasted:
        out[start:end] = a.reshape(-1)
        start, end = end, end + rows
    out[0:rows] = arr.flatten()
    return out.reshape(cols, rows).T


def lengthDiff(a, b):
    """
    Sum of the absolute valued difference of the lengths of each normalized
    triangle (a, b).

    Source: https://stackoverflow.com/a/48285624/1391441
    which in turn is based on:
    https://stackoverflow.com/a/46135435/3293881 (@unutbu)
    """
    val = np.abs(a[:, 1, None] - b[:, 1]) + np.abs(a[:, 2, None] - b[:, 2])

    # Add triangle indexes.
    diff_idxs = indices_merged_arr_generic_using_cp(val)

    return diff_idxs


def findRotAngle(A_pts, B_pts):
    '''
    Find the appropriate rotation between two equivalent triangles, when one
    of them is also translated and scaled.

    Uses the function to obtain the angle between two vectors defined in:
    http://stackoverflow.com/a/13226141/1391441

    The angle returned is in degrees, but we do not know if it is the correct
    rotation angle or the (360 - angle) value.
    '''

    # Distances between each point,in the order given by the dict below.
    dA, dB = pdist(A_pts), pdist(B_pts)
    idx = {"0": (0, 1), "1": (0, 2), "2": (1, 2)}

    # Indexes the points that define the shortest and longest segments for
    # each triangle.
    dA_short, dA_long = np.argmin(dA), np.argmax(dA)
    dB_short, dB_long = np.argmin(dB), np.argmax(dB)
    idxA_shr, idxA_lng = idx[str(dA_short)], idx[str(dA_long)]
    idxB_shr, idxB_lng = idx[str(dB_short)], idx[str(dB_long)]

    # Identify the point in each triangle that belongs to both the shortest
    # and longest segment. This represents the same corner in both triangles,
    # i.e.: the same (assuming the match is right) star.
    A_corner, B_corner = [], []
    for idx_sh in idxA_shr:
        for _ in [A_pts[idxA_lng[0]], A_pts[idxA_lng[1]]]:
            if list(A_pts[idx_sh]) == list(_):
                A_corner = list(_)
    for idx_sh in idxB_shr:
        for _ in [B_pts[idxB_lng[0]], B_pts[idxB_lng[1]]]:
            if B_pts[idx_sh] == _:
                B_corner = list(_)
    # print("A_corner", A_corner)
    # print("B_corner", B_corner)

    # Move the shortest vector (segment) to the (0., 0.) origin, by
    # subtracting the corner point found above.
    if list(A_pts[idxA_shr[0]]) != A_corner:
        vect_A = A_pts[idxA_shr[0]] - np.array(A_corner)
    else:
        vect_A = A_pts[idxA_shr[1]] - np.array(A_corner)

    if list(B_pts[idxB_shr[0]]) != B_corner:
        vect_B = B_pts[idxB_shr[0]] - np.array(B_corner)
    else:
        vect_B = B_pts[idxB_shr[1]] - np.array(B_corner)

    # Use the dot product definition to find the angle between these two
    # translated vectors.
    x1, y1 = vect_A
    x2, y2 = vect_B
    inner_product = x1 * x2 + y1 * y2
    len1 = math.hypot(x1, y1)
    len2 = math.hypot(x2, y2)
    rot_ang = np.rad2deg(math.acos(inner_product / (len1 * len2)))

    return rot_ang, np.asarray(A_corner), np.asarray(B_corner)


def rotatePoints(center, xy, angle):
    """
    Rotates points in 'xy' around 'center'. Angle is in degrees.
    Rotation is counter-clockwise

    http://stackoverflow.com/a/20024348/1391441
    """
    angle = np.radians(angle)
    xy_rot = xy[0] - center[0], xy[1] - center[1]
    xy_rot = (xy_rot[0] * np.cos(angle) - xy_rot[1] * np.sin(angle),
              xy_rot[0] * np.sin(angle) + xy_rot[1] * np.cos(angle))
    xy_rot = xy_rot[0] + center[0], xy_rot[1] + center[1]

    return np.asarray(xy_rot)


def standard2observed(
        xy_ref, scale, rot_angle, ref_cent, obs_cent, obs_tr_match):
    """
    Transform standard stars coordinates to the observed frame coordinates.
    """

    # Move reference stars so that the selected star 'ref_cent' is the new
    # origin.
    ref_trans = xy_ref - ref_cent
    # Apply scaling.
    ref_scaled = ref_trans * scale
    # Move scaled reference stars to the observed star 'obs_cent' as new
    # center.
    ref_obs = (ref_scaled + obs_cent).T

    # Rotate standard stars. Check both possible rotation angles.
    xy_rot1 = rotatePoints(obs_cent, ref_obs, rot_angle)
    xy_rot2 = rotatePoints(obs_cent, ref_obs, -rot_angle)

    # Find the rotation angle that produces the minimum summed distance of
    # the rotated stars to the three manually selected stars.
    d1 = np.min(cdist(obs_tr_match, zip(*xy_rot1)), axis=1).sum()
    d2 = np.min(cdist(obs_tr_match, zip(*xy_rot2)), axis=1).sum()

    # Select rotated standard stars.
    xy_rot = xy_rot1 if d1 < d2 else xy_rot2
    # Pass proper rotation angle.
    rot_angle = rot_angle if d1 < d2 else -rot_angle

    return xy_rot.T, rot_angle


def scaleTransRot(
    A_pts, B_pts, A_combs, B_combs, A_triang, A_tr_not_scaled, B_triang,
    B_tr_not_scaled, mtoler, scale_range, rot_range, trans_range, hdu_data):
    """
    For each normalized triangle in A, compare with each normalized triangle
    in B. Find the differences between their sides, sum their absolute values,
    and select the two triangles with the smallest sum of absolute differences.
    """

    # Store indexes of each triangle vs triangle comparison, along with the
    # sums of the absolute valued differences between the normalized lengths
    # of each triangle vs triangle.
    tr_sum = lengthDiff(A_triang, B_triang)
    # Put smallest sums first.
    tr_sum = tr_sum[tr_sum[:, 0].argsort()]

    mdist_old = 1.e6
    match_flag = None
    for i, tr_sij in enumerate(tr_sum):
        # Index of the triangles in A and B.
        A_idx, B_idx = int(tr_sij[1]), int(tr_sij[2])

        # Scale between observed/reference triangle.
        scale = np.mean(
            np.array(B_tr_not_scaled[B_idx]) / A_tr_not_scaled[A_idx])

        # Accepted range of scaling.
        if scale_range[0] <= scale <= scale_range[1]:

            # Indexes of points of the best match triangles in each set.
            A_idx_pts, B_idx_pts = A_combs[A_idx], B_combs[B_idx]
            # print('Triangle std {} matches triangle obs {}'.format(
            #     A_idx_pts, B_idx_pts))

            # Matched points in A and B.
            A_tr_match = [A_pts[_] for _ in A_idx_pts]
            B_tr_match = [B_pts[_] for _ in B_idx_pts]

            # Rotation angle between triangles.
            rot_angle, ref_cent, obs_cent = findRotAngle(
                A_tr_match, B_tr_match)

            tr_x, tr_y = ref_cent - obs_cent
            if rot_range[0] <= rot_angle <= rot_range[1] and\
                    trans_range[0][0] <= tr_x <= trans_range[0][1] and\
                    trans_range[1][0] <= tr_y <= trans_range[1][1]:

                # Apply translation, scaling, and rotation to reference
                # coordinates.
                xy_transf, rot_angle = standard2observed(
                    A_pts, scale, rot_angle, ref_cent, obs_cent, B_tr_match)

                # Mean distance between detected sources in the observed
                # frame, and the closest transformed standard stars.
                mdist = np.mean(cdist(np.array(B_pts), xy_transf).min(axis=1))
                if mdist < mdist_old:
                    print(" Mean match dist: {:.3f} ({:.0f}% sols, "
                          "{:.3f})".format(
                              mdist, 100. * i / len(tr_sum), tr_sij[0]))
                    mdist_old = mdist
                    # Update values to pass.
                    A_tr_match_f, B_tr_match_f, scale_f, rot_angle_f,\
                        xy_transf_f = A_tr_match, B_tr_match, scale,\
                        rot_angle, xy_transf
                    match_flag = False

                # Accept mtoler error for each matched star.
                if mdist < mtoler * len(xy_transf):
                    # print("Mean dist: {:.3f}".format(mdist))
                    # # print(cdist(np.array(B_pts), xy_transf).min(axis=1))
                    # print("sum: {:.3f}, sc: {:.2f}, rot: {:.2f}".format(
                    #     tr_sij[0], scale, rot_angle))
                    # print("Ref triang match: {}".format(A_tr_match))
                    # print("Obs triang match: {}".format(B_tr_match))
                    # print("A star, B star: {}, {}".format(
                    #     ref_cent, obs_cent))

                    # fig, ax = plt.subplots(1)
                    # interval = ZScaleInterval()
                    # zmin, zmax = interval.get_limits(hdu_data)
                    # ax.imshow(
                    #     hdu_data, cmap='Greys', aspect=1,
                    #     interpolation='nearest',
                    #     origin='lower', vmin=zmin, vmax=zmax)
                    # # Transformed standard stars
                    # ax.scatter(*xy_transf.T, marker='s', edgecolor='g',
                    #            facecolor='', lw=.7, s=60., zorder=4)
                    # # Observed stars marked by user
                    # ax.scatter(*zip(*B_pts), marker='P', edgecolor='r',
                    #            facecolor='', lw=.7, s=60., zorder=4)
                    # plt.show()
                    match_flag = True
                    break

    if match_flag is None:
        print("  WARNING: no match found within range limits.")
        A_tr_match_f, B_tr_match_f, scale_f, rot_angle_f, tr_x, tr_y,\
            xy_transf_f = [], [], np.nan, np.nan, np.nan, np.nan, []
    elif match_flag is False:
        print("  WARNING: match found is outside 'match_toler' value.")
    elif match_flag is True:
        pass

    return A_tr_match_f, B_tr_match_f, scale_f, rot_angle_f, [tr_x, tr_y],\
        xy_transf_f


def triangleMatch(A_pts, B_pts, mtoler, scale_range, rot_range, trans_range, hdu_data):
    """
    Given two sets of 2D points, generate all possible combinations of three
    points for each (triangles), normalize the lengths, and identify the best
    match triangle in each set.
    """
    # Find all possible triangles.
    print("Generating all triangle combinations.")
    A_combs = list(itertools.combinations(range(len(A_pts)), 3))
    B_combs = list(itertools.combinations(range(len(B_pts)), 3))

    # Obtain normalized triangles.
    print("Normalizing triangles.")
    A_triang, A_tr_not_scaled = getTriangles(A_pts, A_combs)
    B_triang, B_tr_not_scaled = getTriangles(B_pts, B_combs)

    # Index of the (A, B) triangles with the smallest difference.
    print("Find best matching triangles ({} vs {}).".format(
        len(A_triang), len(B_triang)))
    A_tr_match, B_tr_match, scale, rot_angle, xy_shift, xy_transf =\
        scaleTransRot(
            A_pts, B_pts, A_combs, B_combs, A_triang, A_tr_not_scaled,
            B_triang, B_tr_not_scaled, mtoler, scale_range, rot_range,
            trans_range, hdu_data)

    return A_tr_match, B_tr_match, scale, rot_angle, xy_shift, xy_transf


def xyTrans(max_shift, xy_ref, mags_ref, xy_dtct, mags_dtct, mtoler):
    """
    Average minimal (Euclidean) distance from points in 'xy_dtct' to points in
    'xy_ref', until the stopping condition.
    """
    x_dtct, y_dtct = np.array(xy_dtct)[:, 0], np.array(xy_dtct)[:, 1]
    xmin, xmax = max_shift[0]
    ymin, ymax = max_shift[1]
    x0, y0, x0_old, y0_old = 0., 0., np.inf, np.inf

    # The 'runs' value is tied to the 'tol' value so that after this
    # many iterations, the limits will have been reduced by tol^runs%.
    # 'li' is the resolution of the shift grid.
    runs, tol, li = 25, .7, 10
    for _ in range(runs):
        # Inner 'for' block values.
        x0_in, y0_in, d_in = x0, y0, np.inf
        for sx in np.linspace(xmin, xmax, li):
            # Shifts in y.
            for sy in np.linspace(ymin, ymax, li):
                # Apply possible x,y translation
                xy_shifted = np.array([x_dtct + sx + x0, y_dtct + sy + y0]).T

                # Minimum distance for each shifted star to the
                # closest star in xy_ref. The indexes are in the sense:
                # xy_ref[min_dist_idx[IDX]] ~ xy_shifted[IDX] (closest star)
                min_dists, min_dist_idx = cKDTree(xy_ref).query(xy_shifted, 1)

                # Matched star's magnitudes absolute difference
                d_mag = np.abs(mags_dtct - mags_ref[min_dist_idx])
                # Pixel distances weighted by mag differences.
                d = np.median(min_dists * d_mag)

                # Update inner 'for' block vars with better shifted values.
                if d < d_in:
                    x0_in, y0_in, d_in = sx + x0, sy + y0, d

        # import matplotlib.pyplot as plt
        # print("d_in={:.2f}, d_px={:.2f}".format(d_in, d_px))
        # plt.scatter(
        #     np.array(xy_ref).T[0][:100] + x0_in,
        #     np.array(xy_ref).T[1][:100] + y0_in, c='r', s=5)
        # plt.scatter(*zip(*xy_dtct[:100]), c='g', s=5)
        # plt.grid()
        # plt.show()

        # Decrease shift limits by tol%
        xmin, xmax = xmin * tol, xmax * tol
        ymin, ymax = ymin * tol, ymax * tol

        # Update final values
        x0, y0 = x0_in, y0_in
        print("  Best x,y shift found: ({:.2f}, {:.2f})".format(x0, y0))
        if abs(x0 - x0_old) < mtoler and abs(y0 - y0_old) < mtoler:
            print("  Tolerance achieved.")
            break
        else:
            x0_old, y0_old = x0, y0

    if _ == runs:
        print("  WARNING: match did not converge to the requested tolerance.")
    if not (max_shift[0][0] < x0 < max_shift[0][1]):
        print("  WARNING: x shift found is outside the limits.")
    if not (max_shift[1][0] < y0 < max_shift[1][1]):
        print("  WARNING: y shift found is outside the limits.")

    return [x0, y0]


def reCenter(hdu_data, positions, side=30):
    """
    Find better center coordinates for selected stars.
    """
    xy_cent = []
    for x0, y0 in positions:
        # Check that position falls inside of the observed frame.
        if 0. < x0 < hdu_data.shape[1] and 0. < y0 < hdu_data.shape[0]:
            crop = cutout_footprint(hdu_data, (x0, y0), side)[0]

            threshold = detect_threshold(crop, snr=10)
            sigma = 3.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))   # FWHM = 3
            kernel = Gaussian2DKernel(sigma)
            kernel.normalize()
            segm = detect_sources(crop, threshold, npixels=5,
                                  filter_kernel=kernel)
            try:
                tbl = source_properties(crop, segm).to_table()
                # Index of the brightest detected source
                idx_b = tbl['source_sum'].argmax()
                # Index of the closest source.
                sh, d_old, idx_c = side / 2., 1.e6, 0
                for i, xy1 in enumerate(tbl['xcentroid', 'ycentroid']):
                    x, y = xy1[0].value, xy1[1].value
                    d = np.sqrt((sh - x) ** 2 + (sh - y) ** 2)
                    tbl['source_sum']
                    if d < d_old:
                        idx_c = i
                        d_old = 1. * d
                # Choose the brightest star except when its distance to the
                # center of the region is N or more that of the closest star.
                N = 4.
                if np.sqrt(
                        (sh - tbl['xcentroid'][idx_b].value) ** 2 +
                        (sh - tbl['ycentroid'][idx_b].value) ** 2) > N * d_old:
                    idx = idx_c
                else:
                    idx = idx_b
                xy1 = tbl['xcentroid'][idx].value, tbl['ycentroid'][idx].value
                x, y = (x0 + xy1[0] - side * .5, y0 + xy1[1] - side * .5)
            except ValueError:
                # "SourceCatalog contains no sources" may happen when no
                # stars are found within the cropped region.
                x, y = x0, y0
        else:
            x, y = x0, y0

        xy_cent.append((x, y))

    return np.asarray(xy_cent)


def autoSrcDetect(pars, hdulist):
    """
    """
    # Background estimation.
    hdu_data = hdulist[0].data
    hdr = hdulist[0].header
    sky_mean, sky_median, sky_std = bckg_data(
        hdr, hdu_data, pars['gain_key'], pars['rdnoise_key'],
        pars['sky_method'])

    # Stars selection.
    psf_select = st_fwhm_select(
        float(pars['dmax']), 1000000, float(pars['thresh_fit']),
        float(pars['fwhm_init']), sky_std, hdu_data)[0]

    psf_select.rename_column('flux', 'flux_fit')
    psf_select = calibrate_magnitudes(psf_select, hdr[pars['exposure_key']])

    # Filter by min/max x,y limits.
    xmi, xma = pars['min_x-max_x'][0]
    ymi, yma = pars['min_y-max_y'][0]
    xmi = float(xmi) if xmi != 'min' else min(psf_select['xcentroid'])
    xma = float(xma) if xma != 'max' else max(psf_select['xcentroid'])
    ymi = float(ymi) if ymi != 'min' else min(psf_select['ycentroid'])
    yma = float(yma) if yma != 'max' else max(psf_select['ycentroid'])

    print(" Selection (x,y) limits: "
          "({:.0f}, {:.0f}) ; ({:.0f}, {:.0f})".format(
              xmi, xma, ymi, yma))
    stars_filter = []
    for st in psf_select:
        if xmi <= st['xcentroid'] <= xma and\
                ymi <= st['ycentroid'] <= yma:
            stars_filter.append(
                [st['xcentroid'], st['ycentroid'], st['cal_mags']])
    # Filter by max number of stars.
    stars_filter = stars_filter[:int(pars['max_stars_match'])]

    id_pick = [str(_) for _ in range(int(pars['max_stars_match']))]
    zip_stars_filter = np.array(zip(*stars_filter))
    xy_cent = list(zip(*zip_stars_filter[:2]))
    xy_mags = list(zip_stars_filter[2])

    return id_pick, xy_cent, xy_mags


def allInFrame(id_ref, hdu_data, xy_transf):
    """
    Assign nan to reference stars located outside the limits
    of the observed frame.
    """
    xy_all, id_all = [], []
    for i, st in enumerate(xy_transf):
        if 0. < st[0] < hdu_data.shape[1] and\
                0. < st[1] < hdu_data.shape[0]:
            xy_all.append(st)
            id_all.append(id_ref[i])
        else:
            xy_all.append([np.nan, np.nan])
            id_all.append(id_ref[i])
    xy_all = np.array(xy_all)

    return id_all, xy_all
