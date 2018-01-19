
import numpy as np
import itertools
import math
import functools
from scipy.spatial.distance import pdist
from scipy.spatial.distance import cdist


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

    # TODO is this necessary?
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
            B_triang, B_tr_not_scaled, mtoler, scale_range, rot_range,
            hdu_data)

    return A_tr_match, B_tr_match, scale, rot_angle, xy_transf
