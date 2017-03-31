
import itertools
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt


def getTriangles(set_X, X_combs):
    """
    Inefficient way of obtaining the lengths of each triangle's side.
    Normalized so that the minimum length is 1.
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
        d_unsort = [d1 / d_min, d2 / d_min, d3 / d_min]
        triang.append(sorted(d_unsort))
        # These are used to obtain the scale between frames.
        tr_not_scaled.append(sorted([d1, d2, d3]))

    return triang, tr_not_scaled


def sumTriangles(A_triang, B_triang):
    """
    For each normalized triangle in A, compare with each normalized triangle
    in B. find the differences between their sides, sum their absolute values,
    and select the two triangles with the smallest sum of absolute differences.
    """
    tr_sum, tr_idx = [], []
    for i, A_tr in enumerate(A_triang):
        for j, B_tr in enumerate(B_triang):
            # Absolute value of lengths differences.
            tr_diff = abs(np.array(A_tr) - np.array(B_tr))
            # Sum the differences
            tr_sum.append(sum(tr_diff))
            tr_idx.append([i, j])

    # Index of the triangles in A and B with the smallest sum of absolute
    # length differences.
    tr_idx_min = tr_idx[tr_sum.index(min(tr_sum))]
    A_idx, B_idx = tr_idx_min[0], tr_idx_min[1]
    # print("Smallest difference: {}".format(min(tr_sum)))

    return A_idx, B_idx


def centTriangle(tr_match):
    """
    Obtain centroid of triangle.
    """
    x_sum, y_sum = np.sum(tr_match, axis=0)
    return (x_sum / 3., y_sum / 3.)


def scalePoints(xy_center, xy, scale):
    """
    Scaled xy points.

    http://codereview.stackexchange.com/q/159183/35351
    """
    delta_x, delta_y = xy_center[0] - xy[0], xy_center[1] - xy[1]
    x_scale = xy_center[0] - scale * delta_x
    y_scale = xy_center[1] - scale * delta_y

    return x_scale, y_scale


def rotatePoints(center, x, y, angle):
    """
    Rotates points in 'xy' around 'center'. Angle is in degrees.
    Rotation is counter-clockwise

    http://stackoverflow.com/a/20024348/1391441
    """
    angle = np.radians(angle)
    xy_rot = x - center[0], y - center[1]
    xy_rot = (xy_rot[0] * np.cos(angle) - xy_rot[1] * np.sin(angle),
              xy_rot[0] * np.sin(angle) + xy_rot[1] * np.cos(angle))
    xy_rot = xy_rot[0] + center[0], xy_rot[1] + center[1]

    return xy_rot


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
    fname = '../output/standards/filt_U/stk_2085.coo'
    t = Table.read(fname, format='ascii')
    xy_obs = zip(*[t['x'], t['y']])

    # Find all possible triangles.
    std_combs = list(itertools.combinations(range(len(xy_std)), 3))
    obs_combs = list(itertools.combinations(range(len(xy_obs)), 3))

    # Obtain normalized triangles.
    std_triang, std_tr_not_scaled = getTriangles(xy_std, std_combs)
    obs_triang, obs_tr_not_scaled = getTriangles(xy_obs, obs_combs)

    # Index of the triangles with the smallest difference.
    std_idx, obs_idx = sumTriangles(std_triang, obs_triang)

    # Indexes of points in A and B of the best match triangles.
    std_idx_pts, obs_idx_pts = std_combs[std_idx], obs_combs[obs_idx]
    print 'triangle A %s matches triangle B %s' % (std_idx_pts, obs_idx_pts)

    # Matched points in A and B.
    std_tr_match, obs_tr_match = [xy_std[_] for _ in std_idx_pts],\
        [xy_obs[_] for _ in obs_idx_pts]
    print "Std:", std_tr_match
    print "Obs:", obs_tr_match
    obs_tr_match_copy = obs_tr_match[:]

    # Centroids for the best match triangles.
    std_tr_cent, obs_tr_cent = centTriangle(std_tr_match),\
        centTriangle(obs_tr_match)
    # Move observed stars triangle to standard triangle's center.
    obs_tr_match -= np.array(obs_tr_cent) - np.array(std_tr_cent)

    # Scale between standard stars and observed stars triangle.
    scale = np.mean(
        np.array(std_tr_not_scaled[std_idx]) / obs_tr_not_scaled[obs_idx])
    print("Scale (std / obs): {}".format(scale))

    x_scale, y_scale = scalePoints(std_tr_cent, zip(*obs_tr_match), scale)

    ax1 = plt.subplot(131)
    ax1.set_title("Standard frame")
    ax1.scatter(*zip(*xy_std), c='k')
    ax1.scatter(*zip(*std_tr_match), marker='s', edgecolor='g',
                facecolor='', lw=2., s=70)

    ax2 = plt.subplot(132)
    ax2.set_title("Observed frame")
    ax2.scatter(*zip(*xy_obs), c='r')
    ax2.scatter(*zip(*obs_tr_match_copy), marker='s', edgecolor='g',
                facecolor='', lw=2., s=70)

    ax3 = plt.subplot(133)
    ax3.set_title("Matched stars in standard frame")
    ax3.scatter(*std_tr_cent, marker='x', c='b')
    ax3.scatter(*zip(*std_tr_match), edgecolor='k', s=70, lw=1., facecolor='')
    ax3.scatter(x_scale, y_scale, edgecolor='r', s=70, lw=1., facecolor='')
    plt.show()


if __name__ == '__main__':
    main()
