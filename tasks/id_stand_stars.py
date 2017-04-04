
import itertools
import math
from scipy.spatial.distance import pdist
from scipy.spatial.distance import cdist
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt


def getTriangles(set_X, X_combs):
    """
    Inefficient way of obtaining the lengths of each triangle's side.
    Normalized so that the minimum length is 1 and sorted so that smaler
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
        print("Scale (obs/std): {:.2f}".format(scale))

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
            print("Rotation: {:.2f} deg".format(rot_angle))

            if rot_range[0] <= rot_angle <= rot_range[1]:
                break
            else:
                print("  Rotation angle out of range")
        else:
            print("  Scale out of range")

    return A_tr_match, B_tr_match, scale, rot_angle


def triangleMatch(A_pts, B_pts, scale_range, rot_range):
    """
    Given two sets of 2D points, generate all possible combinations of three
    points for each (triangles), normalize the lenghts, and identify the best
    match triangle in each set.
    """
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
    print("Translation (x, y): {}".format(obs_tr_cent - std_tr_cent))

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


def make_plot(xy_std, xy_obs, std_tr_match, obs_tr_match, xy_rot):
    """
    """
    ax1 = plt.subplot(131)
    ax1.set_title("Standard frame")
    ax1.scatter(*zip(*xy_std), c='k')
    ax1.scatter(*zip(*std_tr_match), marker='s', edgecolor='g',
                facecolor='', lw=2., s=70)

    ax2 = plt.subplot(132)
    ax2.set_title("Observed frame")
    ax2.scatter(*zip(*xy_obs), c='r')
    ax2.scatter(*zip(*obs_tr_match), marker='s', edgecolor='g',
                facecolor='', lw=2., s=70)

    ax3 = plt.subplot(133)
    ax3.set_title("Standard stars in observed frame")
    # ax3.scatter(*obs_tr_cent, marker='x', c='b')
    ax3.scatter(*zip(*obs_tr_match), edgecolor='r', s=70, lw=1., facecolor='')
    ax3.scatter(*xy_rot, marker='s', edgecolor='k', s=150, lw=1., facecolor='')
    plt.show()


def main():
    """
    This algorithm expects at least three standard stars observed in the
    observed field.
    """

    # Load Landolt standard file. Notice that y axis is inverted.
    # plt.imshow(plt.imread("../tasks/landolt/pg1323-086.gif"))
    # plt.show()

    # Dictionary of stars for this standard field.
    pg1323 = {
        '86': (211., 158.3), '86A': (162.5, 137.5), '86B': (158.1, 128.),
        '86C': (160.1, 171.2), '86D': (89.6, 133.7)}
    # Coordinates of stars for this standard field.
    xy_std = [_ for _ in pg1323.values()]

    # Coordinates from observed frame.
    fname = '../output/standards/filt_I/stk_2083.coo'
    t = Table.read(fname, format='ascii')
    xy_obs = zip(*[t['x'], t['y']])

    scale_range = (0.1, 10.)
    rot_range = (-5., 5.)

    # Best match triangles, scale, and rotation angle between them.
    std_tr_match, obs_tr_match, scale, rot_angle = triangleMatch(
        xy_std, xy_obs, scale_range, rot_range)

    # Move, scale, and rotate standard stars to observed frame.
    xy_rot = standard2observed(
        xy_std, std_tr_match, obs_tr_match, scale, rot_angle)

    make_plot(xy_std, xy_obs, std_tr_match, obs_tr_match, xy_rot)


if __name__ == '__main__':
    main()
