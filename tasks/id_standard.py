
import read_pars_file as rpf

import landolt_fields
import os
from os.path import join
import sys
import itertools
import math

import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.distance import cdist

from hlpr import bckg_data, st_fwhm_select

from astropy.table import Table, hstack
from astropy.io import ascii, fits
from photutils import centroid_2dg
from photutils.utils import cutout_footprint

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def in_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    in_path = join(pars['mypath'].replace('tasks', 'output/standards'))
    ref_id_std = join(in_path, pars['ref_id_std'])
    if not os.path.isfile(ref_id_std):
        print("{}\n Reference frame is not present. Exit.".format(
            ref_id_std))
        sys.exit()
    else:
        pars['ref_id_std'] = ref_id_std

    # Create path to output folder
    out_path = in_path.replace('input', 'output')

    return pars, out_path


def relative_mag(flux):
    """
    """
    return -2.5 * np.log10(flux / max(flux))


def coo_ref_frame(pars):
    """
    """
    # Load .fits file.
    hdulist = fits.open(pars['ref_id_std'])
    hdr = hdulist[0].header
    hdu_data = hdulist[0].data

    # Background estimation.
    sky_mean, sky_std = bckg_data(
        hdr, hdu_data, pars['gain_key'], pars['rdnoise_key'],
        pars['sky_method'])

    # Stars selection.
    psf_select, all_sources, n_not_satur = st_fwhm_select(
        float(pars['dmax']), int(pars['max_stars']),
        float(pars['thresh_level']), float(pars['fwhm_init']),
        sky_std, hdu_data)

    xy_obs = zip(*[psf_select['xcentroid'], psf_select['ycentroid']])
    obs_mag = relative_mag(psf_select['flux'])

    return xy_obs, obs_mag, hdu_data


def getTriangles(set_X, X_combs):
    """
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

    return triang, tr_not_scaled


def findRotAngle(A_pts, B_rot):
    '''
    Find the appropriate rotation between two equivalent triangles, when one
    of them is also translated and scaled.

    Uses the function to obtain the angle between two vectors defined in:
    http://stackoverflow.com/a/13226141/1391441

    The angle returned is in degrees, but we do not know it is the correct
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
    # and longest segment. This represents the same corner in both triangles.
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

    return rot_ang


def scalerotTriangles(
        A_pts, B_pts, A_combs, B_combs, A_triang, A_tr_not_scaled, B_triang,
        B_tr_not_scaled, scale_range, rot_range):
    """
    For each normalized triangle in A, compare with each normalized triangle
    in B. Find the differences between their sides, sum their absolute values,
    and select the two triangles with the smallest sum of absolute differences.
    """
    tr_sum = []
    for i, A_tr in enumerate(A_triang):
        for j, B_tr in enumerate(B_triang):
            # Absolute value of lengths differences.
            tr_diff = abs(np.array(A_tr) - np.array(B_tr))
            # Sum the differences
            tr_sum.append([sum(tr_diff), i, j])

    # Put smallest sums first.
    tr_sum.sort()

    for tr_sij in tr_sum:
        # Index of the triangles in A and B with the smallest sum of absolute
        # length differences.
        A_idx, B_idx = tr_sij[1], tr_sij[2]
        # print("Smallest difference: {}".format(min(tr_sum)))

        # Scale between observed stars and standard stars triangle.
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
            # print("Std: {}".format(A_tr_match))
            # print("Obs: {}".format(B_tr_match))

            # Rotation angle between best match triangles.
            rot_angle = findRotAngle(A_tr_match, B_tr_match)

            if rot_range[0] <= rot_angle <= rot_range[1]:
                print("Scale (obs/std): {:.2f}".format(scale))
                print("Rotation: {:.2f} deg".format(rot_angle))
                break
        #     else:
        #         print("  Rotation angle out of range")
        # else:
        #     print("  Scale out of range")

    return A_tr_match, B_tr_match, scale, rot_angle


def triangleMatch(A_pts, B_pts, scale_range, rot_range):
    """
    Given two sets of 2D points, generate all possible combinations of three
    points for each (triangles), normalize the lengths, and identify the best
    match triangle in each set.
    """
    print("Generating all triangle combinations.")
    # Find all possible triangles.
    A_combs = list(itertools.combinations(range(len(A_pts)), 3))
    B_combs = list(itertools.combinations(range(len(B_pts)), 3))

    # Obtain normalized triangles.
    A_triang, A_tr_not_scaled = getTriangles(A_pts, A_combs)
    B_triang, B_tr_not_scaled = getTriangles(B_pts, B_combs)

    # Index of the (A, B) triangles with the smallest difference.
    A_tr_match, B_tr_match, scale, rot_angle = scalerotTriangles(
        A_pts, B_pts, A_combs, B_combs, A_triang, A_tr_not_scaled, B_triang,
        B_tr_not_scaled, scale_range, rot_range)

    return A_tr_match, B_tr_match, scale, rot_angle


def centTriangle(tr_match):
    """
    Obtain centroid of triangle.
    """
    x_sum, y_sum = np.sum(tr_match, axis=0)
    return np.array([x_sum / 3., y_sum / 3.])


def scalePoints(xy_center, xy, scale):
    """
    Scaled xy points.

    http://codereview.stackexchange.com/q/159183/35351
    """
    delta_x, delta_y = xy_center[0] - xy[0], xy_center[1] - xy[1]
    x_scale = xy_center[0] - scale * delta_x
    y_scale = xy_center[1] - scale * delta_y
    xy_scale = [x_scale, y_scale]

    return xy_scale


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

    return xy_rot


def standard2observed(xy_std, std_tr_match, obs_tr_match, scale, rot_angle):
    """
    Transform standard stars coordinates to the observed frame coordinates.
    """

    # Centroids for the best match triangles.
    std_tr_cent, obs_tr_cent = centTriangle(std_tr_match),\
        centTriangle(obs_tr_match)
    print("Translation (x, y): ({:.2f}, {:.2f})".format(
        *(obs_tr_cent - std_tr_cent)))

    # Move all standard stars to observed triangle's center.
    std_trans = xy_std - (std_tr_cent - obs_tr_cent)
    # Scale standard stars to observed stars scale.
    xy_scale = scalePoints(obs_tr_cent, zip(*std_trans), scale)

    # Rotate standard stars. Check both possible rotation angles.
    xy_rot1 = rotatePoints(obs_tr_cent, xy_scale, rot_angle)
    xy_rot2 = rotatePoints(obs_tr_cent, xy_scale, 360. - rot_angle)

    # Find the rotation angle that produces the minimum summed distance of
    # the rotated stars to the standard stars.
    d1 = np.min(cdist(xy_std, zip(*xy_rot1)), axis=1).sum()
    d2 = np.min(cdist(xy_std, zip(*xy_rot2)), axis=1).sum()
    # Select rotated standard stars.
    xy_rot = xy_rot1 if d1 < d2 else xy_rot2

    return xy_rot


def star_size(mag, max_mag):
    '''
    Convert "relative magnitudes" from psfmeasure into intensities.

    http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?psfmeasure
    '''
    # Scale factor.
    factor = 10. * (1 - 1 / (1 + 150 / len(mag) ** 0.85))
    return 0.1 + factor * 10 ** ((np.array(mag) / -2.5)) * max_mag


def posFinder(xy_rot, max_x, min_x, max_y, min_y):
    """
    Finds the offset position to place the IDs of the standard stars in the
    observed frame, such that they don't overlap with other stars or text.
    """
    # Ranges in x,y
    x_rang, y_rang = max_x - min_x, max_y - min_y
    # This value sets the minimum spacing allowed.
    txt_sep = max(x_rang, y_rang) * .025
    # Initiate list with the coordinates of the points themselves.
    used_positions, xy_offset = xy_rot[:], []
    # For every standard star positioned in the observed frame.
    for i, xy in enumerate(xy_rot):
        point_done = False
        # Factor that defined the offset separation.
        for j in np.arange(1., 10, .5):
            rang_perc = .025 * j
            # For each possible combinatin of the offsets.
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
        coo_file, out_plot_file, xy_obs, obs_mag, std_tr_match, obs_tr_match,
        xy_rot, id_std, landolt_field_img):
    """
    Make plots.
    """
    print("Plotting.")
    fig = plt.figure(figsize=(20, 20))
    gs = gridspec.GridSpec(10, 10)

    ax1 = plt.subplot(gs[0:4, 0:4])
    # ax1.set_aspect('auto')
    ax1.set_title("Standard frame")
    land_img = ax1.imshow(plt.imread(landolt_field_img))
    # Extract maximum y axis value
    max_y = float(len(land_img.get_array()))
    x, y = zip(*std_tr_match)
    y_inv = max_y - np.array(y)
    ax1.scatter(x, y_inv, marker='s', edgecolor='r', facecolor='', lw=1.,
                s=60)

    ax2 = plt.subplot(gs[0:4, 4:8])
    ax2.set_aspect('auto')
    ax2.set_title("Observed frame ({})".format(coo_file))
    ax2.grid(lw=1., ls='--', color='grey', zorder=1)
    st_sizes_arr = star_size(obs_mag, max(obs_mag))
    x_obs, y_obs = zip(*xy_obs)
    ax2.scatter(x_obs, y_obs, c='k', s=st_sizes_arr, zorder=3)
    # Define offsets.
    max_x, min_x, max_y, min_y = max(x_obs), min(x_obs), max(y_obs), min(y_obs)
    x, y = xy_rot
    ax2.scatter(x, y, marker='s', edgecolor='g',
                facecolor='', lw=.7, s=st_sizes_arr + 50., zorder=4)
    xy_offset = posFinder(zip(*xy_rot), max_x, min_x, max_y, min_y)
    for i, txt in enumerate(id_std):
        ax2.annotate(
            txt, xy=(x[i], y[i]), xytext=(xy_offset[i][0], xy_offset[i][1]),
            fontsize=10, arrowprops=dict(
                arrowstyle="-", connectionstyle="angle3,angleA=90,angleB=0"))
    ax2.scatter(*zip(*obs_tr_match), marker='s', edgecolor='r',
                facecolor='', lw=.7, s=st_sizes_arr + 50., zorder=6)

    fig.tight_layout()
    plt.savefig(out_plot_file, dpi=150, bbox_inches='tight')
    plt.clf()
    plt.close()


def reCenter(hdu_data, positions, side=30):
    """
    Find better center coordinates for standards in the observed frame.
    """
    xy_cent = [[], []]
    for x, y in zip(*positions):
        crop = cutout_footprint(hdu_data, (x, y), side)[0]
        xy = centroid_2dg(crop)
        xy_cent[0].append(x + xy[0] - side * .5)
        xy_cent[1].append(y + xy[1] - side * .5)
        # median, std = np.median(crop), np.std(crop)
        # plt.imshow(crop, cmap='viridis', aspect=1, interpolation='nearest',
        #            origin='lower', vmin=0., vmax=median + std)
        # plt.show()

    return xy_cent


def make_out_file(landolt_fld, out_data_file, landolt_t, xy_rot):
    """
    Write coordinates of standard stars in the observed frame system.
    """
    landolt_f = Table({'Landolt': [landolt_fld for _ in landolt_t['ID']]})
    xy_obs = Table(xy_rot, names=('x_obs', 'y_obs'))
    tt = hstack([landolt_f, landolt_t, xy_obs])

    ascii.write(tt, out_data_file, format='fixed_width', delimiter='',
                formats={'x_obs': '%9.3f', 'y_obs': '%9.3f'}, overwrite=True)


def main():
    """
    This algorithm expects at least three standard stars observed in the
    observed field.
    """
    pars, out_path = in_params()
    f_name = pars['ref_id_std'].replace(pars['mypath'].replace(
        'tasks', 'output/standards/'), '')
    print("Reference frame: {}".format(f_name))

    # Selected standard field.
    landolt_t = landolt_fields.main(pars['landolt_fld'])
    xy_std = zip(*[landolt_t['x'], landolt_t['y']])
    id_std = landolt_t['ID']

    # Coordinates from observed frame.
    xy_obs, obs_mag, hdu_data = coo_ref_frame(pars)

    # Best match triangles, scale, and rotation angle between them.
    scale_range = (float(pars['scale_min']), float(pars['scale_max']))
    rot_range = (float(pars['rot_min']), float(pars['rot_max']))
    print("\nFinding scale, translation, and rotation.")
    std_tr_match, obs_tr_match, scale, rot_angle = triangleMatch(
        xy_std, xy_obs, scale_range, rot_range)

    # Move, scale, and rotate standard stars to observed frame.
    xy_rot = standard2observed(
        xy_std, std_tr_match, obs_tr_match, scale, rot_angle)

    print("Find better center coordinates.")
    xy_cent = reCenter(hdu_data, xy_rot)
    out_data_file = join(out_path, pars['landolt_fld'] + "_obs.coo")
    make_out_file(pars['landolt_fld'], out_data_file, landolt_t, xy_cent)

    if pars['do_plots_C'] == 'y':
        out_plot_file = join(out_path, pars['landolt_fld'] + "_obs.png")
        landolt_field_img = join(
            pars['mypath'], 'landolt', pars['landolt_fld'] + '.gif')
        make_plot(
            f_name, out_plot_file, xy_obs, obs_mag, std_tr_match,
            obs_tr_match, xy_rot, id_std, landolt_field_img)

    print("\nFinished.")


if __name__ == '__main__':
    main()
