
import read_pars_file as rpf

import landolt_fields
import os
from os.path import join, isfile, exists
import sys
import itertools
import math
import functools
import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.distance import cdist

from hlpr import bckg_data, st_fwhm_select

from astropy.table import Table, hstack
from astropy.io import ascii, fits
from photutils import centroid_2dg
from photutils.utils import cutout_footprint

from astropy.visualization import ZScaleInterval
from photutils import CircularAperture
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def in_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    in_path = join(pars['mypath'].replace('tasks', 'input'))
    out_path = in_path.replace('input', 'output')

    landolt_fld, ref_frame, match_fldr = [], [], []
    for line in pars['ref_fld_fldr']:
        # Identify if this is a Landolt field or an observed frame.
        if not line[0].endswith(".fits"):
            # A Landolt field was set. We are matching standard frames.
            landolt_fld.append(line[0])
            ref_frame.append('--')
            match_fldr.append(line[1])
        else:
            landolt_fld.append('--')
            ref_frame.append(line[0])
            match_fldr.append(line[1])

    pars['landolt_fld'], pars['ref_frame'], pars['match_fldr'] =\
        landolt_fld, ref_frame, match_fldr

    # Generate full path to reference image, and check existence.
    for i, ref_im in enumerate(pars['ref_frame']):
        if ref_im != '--':
            ref_frame = join(in_path, pars['match_fldr'][i], ref_im)
            if not os.path.isfile(ref_frame):
                print("{}\n Reference frame is not present. Exit.".format(
                    ref_frame))
                sys.exit()
            else:
                # Store full path.
                pars['ref_frame'][i] = ref_frame

    # Generate list of fits files for each input folder.
    fits_list = []
    for folder in pars['match_fldr']:
        folder = folder[1:] if folder.startswith('/') else folder
        in_path = join(pars['mypath'].replace('tasks', 'input'), folder)

        list_temp = []
        if os.path.isdir(in_path):
            for file in os.listdir(in_path):
                f = join(in_path, file)
                if isfile(f):
                    if f.endswith('.fits'):
                        list_temp.append(f)
        if not list_temp:
            print("{}\n No .fits files found in match folder. Exit.".format(
                in_path))
            sys.exit()

        # Store list for this folder.
        print("Files found in '{}' folder:".format(folder))
        for fit in list_temp:
            print(" * {}".format(fit.replace(in_path, '')[1:]))
        fits_list.append(list_temp)

    return pars, fits_list, out_path


def zoom_fun(event, ax):
    """
    Source: https://gist.github.com/tacaswell/3144287
    """
    base_scale = 2.
    # get the current x and y limits
    cur_xlim = ax.get_xlim()
    cur_ylim = ax.get_ylim()
    # set the range
    cur_xrange = (cur_xlim[1] - cur_xlim[0]) * .5
    cur_yrange = (cur_ylim[1] - cur_ylim[0]) * .5
    # get event location
    xdata, ydata = event.xdata, event.ydata
    if event.button == 'up':
        # deal with zoom in
        scale_factor = 1 / base_scale
    elif event.button == 'down':
        # deal with zoom out
        scale_factor = base_scale
    else:
        # deal with something that should never happen
        scale_factor = 1
        print event.button
    # set new limits
    ax.set_xlim([xdata - cur_xrange * scale_factor,
                 xdata + cur_xrange * scale_factor])
    ax.set_ylim([ydata - cur_yrange * scale_factor,
                 ydata + cur_yrange * scale_factor])
    # force re-draw
    ax.figure.canvas.draw()


def reCenter(hdu_data, positions, side=30):
    """
    Find better center coordinates for selected stars.
    """

    # from astropy.visualization import ZScaleInterval
    # from photutils import CircularAperture
    # interval = ZScaleInterval()
    # zmin, zmax = interval.get_limits(hdu_data)
    # plt.imshow(hdu_data, cmap='viridis', aspect=1, interpolation='nearest',
    #            origin='lower', vmin=zmin, vmax=zmax)
    # fwhm_mean = 3.
    # apertures = CircularAperture(positions, r=2. * fwhm_mean)
    # apertures.plot(color='red', lw=1.)
    # plt.show()

    xy_cent = []
    for x, y in positions:
        # Check that new position falls inside of the observed frame.
        if 0. < x < hdu_data.shape[1] and 0. < y < hdu_data.shape[0]:
            crop = cutout_footprint(hdu_data, (x, y), side)[0]
            xy = centroid_2dg(crop)
            xy_cent.append((x + xy[0] - side * .5, y + xy[1] - side * .5))
            # median, std = np.median(crop), np.std(crop)
            # plt.imshow(
            #     crop, cmap='viridis', aspect=1, interpolation='nearest',
            #     origin='lower', vmin=0., vmax=median + std)
            # # PLot new center.
            # plt.scatter(*xy, c='g')
            # plt.show()

    return np.asarray(xy_cent)


def coo_ref_frame(fr, id_selec, xy_selec, pars, proc_grps):
    """
    Find sources on 'fr' fits file. This uses parameter values from other
    tasks:

    General: gain_key, rdnoise_key, dmax
    fitstats: sky_method, fwhm_init, thresh_fit

    """
    answ = 'n'
    if xy_selec:
        answ = raw_input("Use existing coordinates? (y/n): ")

    # Load .fits file.
    hdulist = fits.open(fr)
    hdu_data = hdulist[0].data
    hdr = hdulist[0].header

    if answ != 'y':

        # # Background estimation.
        # sky_mean, sky_median, sky_std = bckg_data(
        #     hdr, hdu_data, pars['gain_key'], pars['rdnoise_key'],
        #     pars['sky_method'])

        # # Stars selection.
        # psf_select, _, _ = st_fwhm_select(
        #     float(pars['dmax']), 1000000, float(pars['thresh_fit']),
        #     float(pars['fwhm_init']), sky_std, hdu_data)

        # # Filter by min/max x,y limits.
        # xmin, xmax = pars['min_x-max_x'][0]
        # ymin, ymax = pars['min_y-max_y'][0]
        # xmin = float(xmin) if xmin != 'min' else min(psf_select['xcentroid'])
        # xmax = float(xmax) if xmax != 'max' else max(psf_select['xcentroid'])
        # ymin = float(ymin) if ymin != 'min' else min(psf_select['ycentroid'])
        # ymax = float(ymax) if ymax != 'max' else max(psf_select['ycentroid'])
        # print(" Select stars within the (x,y) limits: "
        #       "({:.0f}, {:.0f}) ; ({:.0f}, {:.0f})".format(
        #           xmin, xmax, ymin, ymax))
        # stars_filter = []
        # for st in psf_select:
        #     if xmin <= st['xcentroid'] <= xmax and\
        #             ymin <= st['ycentroid'] <= ymax:
        #         stars_filter.append(
        #             [st['xcentroid'], st['ycentroid'], st['flux']])
        # # Filter by max number of stars.
        # stars_filter = stars_filter[:int(pars['max_stars_match'])]

        # zip_stars_filter = np.array(zip(*stars_filter))
        # xy_sources = list(zip(*zip_stars_filter[:2]))
        # print(" Selected detected sources: {}".format(len(xy_sources)))
        # print(xy_sources)

        # Ref stars selection
        landolt_field_img = join(
            pars['mypath'], 'landolt', pars['landolt_fld'][proc_grps] + '.gif')
        # Manual selection of standard stars in the observed frame.
        fig, (ax1, ax2) = plt.subplots(1, 2)
        interval = ZScaleInterval()
        zmin, zmax = interval.get_limits(hdu_data)
        ax1.imshow(plt.imread(landolt_field_img))
        ax2.imshow(
            hdu_data, cmap='Greys', aspect=1, interpolation='nearest',
            origin='lower', vmin=zmin, vmax=zmax)

        def onclick(event, ax):
            if event.button == 1:
                apertures = CircularAperture(
                    (event.xdata, event.ydata), r=10.)
                apertures.plot(color='green', lw=1.)
                ax.figure.canvas.draw()
                ref_id = raw_input(" ID of selected star: ")
                print(" {} added to list: ({:.2f}, {:.2f})".format(
                    ref_id, event.xdata, event.ydata))
                id_selec.append(ref_id)
                xy_selec.append((event.xdata, event.ydata))
                if len(id_selec) == 3:
                    plt.close()
            elif event.button == 2:
                print("scroll click")
            elif event.button == 3:
                print("right click")
            else:
                pass

        # Mouse click / scroll zoom events.
        fig.canvas.mpl_connect(
            'scroll_event', lambda event: zoom_fun(event, ax2))
        fig.canvas.mpl_connect(
            'button_press_event', lambda event: onclick(event, ax2))
        print("\nSelect three reference stars covering "
              "as much frame as possible.")
        plt.show()

        id_selec = ['190', '263', '101']
        xy_selec = [[1999.0899459817479, 1932.7875835241143], [824.1312661197503, 3343.7059611550444], [2328.7031396812936, 65.03048532537136]]

        # Re-center coordinates.
        d, runs = 100., 0
        while d > .5 and runs < 10:
            xy_cent = reCenter(hdu_data, xy_selec)
            d = np.mean(abs(xy_cent - np.array(xy_selec)))
            xy_selec = xy_cent
            runs += 1

        xy_selec = xy_selec.tolist()
        print(xy_selec)

    print(" Stars selected for match: {}".format(len(xy_selec)))

    return hdu_data, id_selec, xy_selec


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


def findRotAngle(A_pts, B_rot):
    '''
    Find the appropriate rotation between two equivalent triangles, when one
    of them is also translated and scaled.

    Uses the function to obtain the angle between two vectors defined in:
    http://stackoverflow.com/a/13226141/1391441

    The angle returned is in degrees, but we do not know if it is the correct
    rotation angle or the (360 - angle) value.
    '''

    # Distances between each point,in the order given by the dict below.
    dA, dB = pdist(A_pts), pdist(B_rot)
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
        for _ in [B_rot[idxB_lng[0]], B_rot[idxB_lng[1]]]:
            if B_rot[idx_sh] == _:
                B_corner = list(_)
    # print("A_corner", A_corner)
    # print("B_corner", B_corner)

    # Move the shortest vector (segment) to the (0., 0.) origin, by
    # subtracting the corner point found above.
    if list(A_pts[idxA_shr[0]]) != A_corner:
        vect_A = A_pts[idxA_shr[0]] - np.array(A_corner)
    else:
        vect_A = A_pts[idxA_shr[1]] - np.array(A_corner)

    if list(B_rot[idxB_shr[0]]) != B_corner:
        vect_B = B_rot[idxB_shr[0]] - np.array(B_corner)
    else:
        vect_B = B_rot[idxB_shr[1]] - np.array(B_corner)

    # Use the dot product definition to find the angle between these two
    # translated vectors.
    x1, y1 = vect_A
    x2, y2 = vect_B
    inner_product = x1 * x2 + y1 * y2
    len1 = math.hypot(x1, y1)
    len2 = math.hypot(x2, y2)
    rot_ang = np.rad2deg(math.acos(inner_product / (len1 * len2)))

    return rot_ang, np.asarray(A_corner), np.asarray(B_corner)


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


def standard2observed(xy_ref, scale, rot_angle, ref_cent, obs_cent):
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
    # xy_rot2 = rotatePoints(obs_cent, ref_obs, 360. - rot_angle)

    # # Find the rotation angle that produces the minimum summed distance of
    # # the rotated stars to the standard stars.
    # d1 = np.min(cdist(xy_ref, zip(*xy_rot1)), axis=1).sum()
    # d2 = np.min(cdist(xy_ref, zip(*xy_rot2)), axis=1).sum()
    # # Select rotated standard stars.
    # xy_rot = xy_rot1 if d1 < d2 else xy_rot2

    return xy_rot1.T


def scaleTransRot(
    A_pts, B_pts, A_combs, B_combs, A_triang, A_tr_not_scaled, B_triang,
    B_tr_not_scaled, mtoler, scale_range, rot_range, hdu_data):
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
    A_tr_match_f, B_tr_match_f, scale_f, rot_angle_f, xy_transf_f =\
        [], [], 0., 0., []
    match_flag = False
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

            if rot_range[0] <= rot_angle <= rot_range[1]:

                # Apply translation, scaling, and rotation to reference
                # coordinates.
                xy_transf = standard2observed(
                    A_pts, scale, rot_angle, ref_cent, obs_cent)

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

                # Accept half a pixel error for each matched star.
                if mdist < mtoler * len(xy_transf):
                    print("Mean dist: {:.3f}".format(mdist))
                    # print(cdist(np.array(B_pts), xy_transf).min(axis=1))
                    print("sum: {:.3f}, sc: {:.2f}, rot: {:.2f}".format(
                        tr_sij[0], scale, rot_angle))
                    print("Ref triang match: {}".format(A_tr_match))
                    print("Obs triang match: {}".format(B_tr_match))
                    print("A star, B star: {}, {}".format(
                        ref_cent, obs_cent))

                    fig, ax = plt.subplots(1)
                    interval = ZScaleInterval()
                    zmin, zmax = interval.get_limits(hdu_data)
                    ax.imshow(
                        hdu_data, cmap='Greys', aspect=1,
                        interpolation='nearest',
                        origin='lower', vmin=zmin, vmax=zmax)
                    # Transformed standard stars
                    ax.scatter(*xy_transf.T, marker='s', edgecolor='g',
                               facecolor='', lw=.7, s=60., zorder=4)
                    # Observed stars marked by user
                    ax.scatter(*zip(*B_pts), marker='P', edgecolor='r',
                               facecolor='', lw=.7, s=60., zorder=4)
                    plt.show()
                    match_flag = True
                    break

    if match_flag is False:
        print("  WARNING: no match found within minimum mean match distance.")

    return A_tr_match_f, B_tr_match_f, scale_f, rot_angle_f, xy_transf_f


def triangleMatch(A_pts, B_pts, mtoler, scale_range, rot_range, hdu_data):
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
    A_tr_match, B_tr_match, scale, rot_angle, xy_transf =\
        scaleTransRot(
            A_pts, B_pts, A_combs, B_combs, A_triang, A_tr_not_scaled,
            B_triang, B_tr_not_scaled, mtoler, scale_range, rot_range, hdu_data)

    return A_tr_match, B_tr_match, scale, rot_angle, xy_transf


def outFileHeader(out_data_file):
    """
    Create output stats file with a header.
    """
    with open(out_data_file, mode='w') as f:
        f.write(
            "#  frame      x_obs      y_obs   ref    ID          x          "
            "y      V      BV      UB     VR     RI     VI     e_V    e_BV    "
            "e_UB    e_VR    e_RI    e_VI\n")


def make_out_file(
    f_name, id_inframe, xy_inframe, img_id, landolt_t, out_data_file):
    """
    Write coordinates of reference stars in the observed frame system.
    """
    landolt_in = Table(dtype=landolt_t.dtype)
    for r_id in id_inframe:
        i = landolt_t['ID'].tolist().index(r_id)
        landolt_in.add_row(landolt_t[i])

    landolt_f = Table({'reference': [img_id for _ in landolt_in['ID']]})
    ref_transf = [[f_name, xy[0], xy[1]] for xy in xy_inframe]
    obs_tbl = Table(zip(*ref_transf), names=('frame', 'x_obs', 'y_obs'))
    tt = hstack([obs_tbl, landolt_f, landolt_in])

    with open(out_data_file, mode='a') as f:
        # Some platforms don't automatically seek to end when files opened
        # in append mode
        f.seek(0, os.SEEK_END)
        ascii.write(
            tt, f, format='fixed_width_no_header', delimiter='',
            formats={'x_obs': '%9.3f', 'y_obs': '%9.3f', 'x': '%9.3f',
                     'y': '%9.3f'})


def posFinder(xy_inframe, max_x, min_x, max_y, min_y):
    """
    Finds the offset position to place the IDs of the standard stars in the
    observed frame, such that they don't overlap with other stars or text.
    """
    # Ranges in x,y
    x_rang, y_rang = max_x - min_x, max_y - min_y
    # This value sets the minimum spacing allowed.
    txt_sep = max(x_rang, y_rang) * .025
    # Initiate list with the coordinates of the points themselves.
    used_positions, xy_offset = xy_inframe[:], []
    # For every standard star positioned in the observed frame.
    for i, xy in enumerate(xy_inframe):
        point_done = False
        # Factor that defines the offset separation.
        for j in np.arange(1., 10, .5):
            rang_perc = .025 * j
            # For each possible combination of the offsets.
            for s in [[1., 1.], [1., -1.], [-1., 1.], [-1., -1.]]:
                offset_x = s[0] * rang_perc * x_rang
                offset_y = s[1] * rang_perc * y_rang
                # Offset coordinates.
                xi_off, yi_off = xy[0] + offset_x, xy[1] + offset_y
                # Check that the coordinates are inside the figure.
                if min_x < xi_off < max_x and min_y < yi_off < max_y:
                    # Distance between the position of the text and all
                    # positions already used.
                    d = cdist(
                        np.array([[xi_off, yi_off]]), np.array(used_positions))
                    min_d = np.min(d, axis=1)
                    # Store if it is far away enough.
                    if min_d > txt_sep:
                        used_positions.append([xi_off, yi_off])
                        xy_offset.append([xi_off, yi_off])
                        point_done = True
                        break
            if point_done:
                break
        if not point_done:
            # If no suitable position was found for this point, assign an
            # offset close to its coordinates.
            offset_x, offset_y = .025 * x_rang, .025 * y_rang
            xi_off, yi_off = xy[0] + offset_x, xy[1] + offset_y
            used_positions.append([xi_off, yi_off])
            xy_offset.append([xi_off, yi_off])

    return xy_offset


def make_plot(
    f_name, hdu_data, out_plot_file, std_tr_match, obs_tr_match,
    xy_inframe, id_inframe, landolt_field_img):
    """
    Make plots.
    """
    print("Plotting.")
    fig = plt.figure(figsize=(20, 20))
    gs = gridspec.GridSpec(10, 10)

    ax1 = plt.subplot(gs[0:4, 0:4])
    ax1.set_title("Standard/Reference frame")
    land_img = ax1.imshow(plt.imread(landolt_field_img))
    # Extract maximum y axis value
    max_y = float(len(land_img.get_array()))
    x, y = zip(*std_tr_match)
    y_inv = max_y - np.array(y)
    ax1.scatter(x, y_inv, marker='s', edgecolor='r', facecolor='', lw=1.,
                s=60)

    ax2 = plt.subplot(gs[0:4, 4:8])
    ax2.set_aspect('auto')
    xmin, ymin = xy_inframe.T.min(axis=1)
    xmax, ymax = xy_inframe.T.max(axis=1)
    ax2.set_xlim(max(0., xmin - xmax * .2),
                 min(hdu_data.shape[1], xmax + xmax * .2))
    ax2.set_ylim(max(0., ymin - ymax * .2),
                 min(hdu_data.shape[0], ymax + ymax * .2))
    ax2.set_title("Observed frame ({})".format(f_name))
    ax2.grid(lw=1., ls='--', color='grey', zorder=1)
    interval = ZScaleInterval()
    zmin, zmax = interval.get_limits(hdu_data)
    ax2.imshow(hdu_data, cmap='Greys', aspect=1, interpolation='nearest',
               origin='lower', vmin=zmin, vmax=zmax)
    # standard stars in coordinates of the observed frame.
    ax2.scatter(*zip(*xy_inframe), marker='s', edgecolor='g',
                facecolor='', lw=.7, s=60., zorder=4)

    # Define offsets.
    xy_offset = posFinder(xy_inframe.tolist(), xmax, xmin, ymax, ymin)
    for i, txt in enumerate(id_inframe):
        ax2.annotate(
            txt, xy=(xy_inframe[i][0], xy_inframe[i][1]),
            xytext=(xy_offset[i][0], xy_offset[i][1]),
            fontsize=10, arrowprops=dict(
                arrowstyle="-", color='b',
                connectionstyle="angle3,angleA=90,angleB=0"))
    ax2.scatter(*zip(*obs_tr_match), marker='s', edgecolor='r',
                facecolor='', lw=.7, s=60., zorder=6)

    fig.tight_layout()
    plt.savefig(out_plot_file, dpi=150, bbox_inches='tight')
    plt.clf()
    plt.close()


def main():
    """
    This algorithm expects at least three stars detected in the observed field.
    """
    pars, fits_list, out_path = in_params()

    # Process each defined group.
    for proc_grps, ref_frame in enumerate(pars['ref_frame']):

        id_selec, xy_selec = [], []

        if pars['landolt_fld'][proc_grps] != '--':
            print("\nLandolt field: {}".format(pars['landolt_fld'][proc_grps]))
            # Selected standard field.
            landolt_t = landolt_fields.main(pars['landolt_fld'][proc_grps])
            # Sort putting brightest stars at the top.
            landolt_t.sort('V')
            # Extract (x,y) coords and IDs.
            xy_ref = zip(*[landolt_t['x'], landolt_t['y']])
            id_ref = landolt_t['ID']
            # TODO
            # # Use only selected number of brightest stars.
            # xy_ref = xy_ref[:int(pars['max_stars_match'])]
            # id_ref = id_ref[:int(pars['max_stars_match'])]
            print(" Stars selected for match: {}".format(len(xy_ref)))
            # Path to Landolt image.
            landolt_field_img = join(
                pars['mypath'], 'landolt',
                pars['landolt_fld'][proc_grps] + '.gif')
            # ID for final image
            img_id = pars['landolt_fld'][proc_grps]

        else:
            # TODO finish
            print("\nReference frame: {}\n".format(ref_frame.split('/')[-1]))
            _, xy_sources, _ = coo_ref_frame(
                ref_frame, xy_selec, pars, proc_grps)
            # ID for final image/file.
            img_id = ref_frame.split('/')[-1].split('.')[0]
            landolt_t, id_ref = [], []

        # Generate output subdir if it doesn't exist.
        out_folder = join(out_path, pars['match_fldr'][proc_grps])
        if not exists(out_folder):
            os.makedirs(out_folder)

        # Write output file header
        out_data_file = join(out_path, img_id + "_match.coo")
        outFileHeader(out_data_file)

        # Process each fits file in list.
        for fr in fits_list[proc_grps]:
            f_name = fr.split('/')[-1].split('.')[0]
            # Name of final image.
            out_img = join(out_folder, img_id + '_' + f_name + "_obs.png")

            # Skip reference frame.
            if f_name != img_id:

                print("\nMatching frame: {}".format(f_name))

                # Coordinates from observed frame.
                hdu_data, id_selec, xy_selec = coo_ref_frame(
                    fr, id_selec, xy_selec, pars, proc_grps)

                # # Scale and rotation ranges, and match tolerance (in pixels).
                # scale_range = (
                #     float(pars['scale_min']), float(pars['scale_max']))
                # rot_range = (float(pars['rot_min']), float(pars['rot_max']))
                # mtoler = float(pars['match_toler'])

                # print("\nFinding scale, translation, and rotation.")
                # std_tr_match, obs_tr_match, scale, rot_angle, xy_transf =\
                #     triangleMatch(
                #         xy_ref, xy_selec, mtoler, scale_range,
                #         rot_range, hdu_data)

                xy_ref_sel = []
                for r_id in id_selec:
                    i = id_ref.tolist().index(r_id)
                    xy_ref_sel.append(xy_ref[i])

                # Matched stars in ref and obs
                std_tr_match, obs_tr_match = xy_ref_sel, xy_selec

                _, A_tr_not_scaled = getTriangles(xy_ref_sel, [(0, 1, 2)])
                _, B_tr_not_scaled = getTriangles(xy_selec, [(0, 1, 2)])

                # Obtain scale.
                scale = np.mean(np.array(B_tr_not_scaled[0]) / A_tr_not_scaled[0])
                print(np.array(B_tr_not_scaled[0]) / A_tr_not_scaled[0])
                # Rotation angle between triangles.
                rot_angle, ref_cent, obs_cent = findRotAngle(
                    xy_ref_sel, xy_selec)

                print("Scale: {:.2f}, Rot: {:.2f}".format(scale, rot_angle))

                # Apply translation, scaling, and rotation to reference
                # coordinates.
                xy_transf = standard2observed(
                    xy_ref, scale, rot_angle, ref_cent, obs_cent)

                # Remove reference stars located outside the observed frame.
                # x axis
                xy_inframe, id_inframe = [], []
                for i, st in enumerate(xy_transf):
                    if 0. < st[0] < hdu_data.shape[1] and\
                            0. < st[1] < hdu_data.shape[0]:
                        xy_inframe.append(st)
                        id_inframe.append(id_ref[i])
                xy_inframe = np.array(xy_inframe)

                print("Re-center final coordinates.")
                d, runs = 100., 0
                while d > .5 and runs < 10:
                    xy_cent = reCenter(hdu_data, xy_inframe, side=20)
                    d = np.mean(abs(xy_cent - np.array(xy_inframe)))
                    xy_inframe = xy_cent
                    runs += 1

                if pars['do_plots_C'] == 'y':
                    make_plot(
                        f_name, hdu_data, out_img, std_tr_match,
                        obs_tr_match, xy_inframe, id_inframe,
                        landolt_field_img)

            else:
                xy_inframe = []

            print("Write final .coo file.")
            make_out_file(
                f_name, id_inframe, xy_inframe, img_id, landolt_t,
                out_data_file)

    print("\nFinished.")


if __name__ == '__main__':
    main()
