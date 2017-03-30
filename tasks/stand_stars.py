
from astropy.table import Table
import numpy as np
from scipy.spatial import distance
import math
from collections import OrderedDict
import matplotlib.pyplot as plt


def scalePoints(xy_center, delta_x, delta_y, scale):
    """
    Scaled xy points.

    http://codereview.stackexchange.com/q/159183/35351
    """
    x_scale = xy_center[0] - scale * delta_x
    y_scale = xy_center[1] - scale * delta_y

    return x_scale, y_scale


def rotatePoints(center, x, y, angle):
    """
    Rotates points in 'xy' around 'center'. Angle is in degrees.
    Rotation is counter-clockwise

    http://stackoverflow.com/a/20024348/1391441
    """
    angle = math.radians(angle)
    xy_rot = x - center[0], y - center[1]
    xy_rot = (xy_rot[0] * math.cos(angle) - xy_rot[1] * math.sin(angle),
              xy_rot[0] * math.sin(angle) + xy_rot[1] * math.cos(angle))
    xy_rot = xy_rot[0] + center[0], xy_rot[1] + center[1]

    return xy_rot


def distancePoints(xy_st, x_transl, y_transl):
    """
    Find the sum of the minimum distance of stars in the Landolt frame to
    stars in the observed frame (scaled, rotates, and moved)
    """
    d = distance.cdist(xy_st, zip(*[x_transl, y_transl]), 'euclidean')
    # Sum of all minimal distances.
    d_sum = np.sum(np.min(d, axis=1))

    return d_sum


def match_frames(
        xy_st, xy_center, delta_x, delta_y, tol, sc_min, sc_max, sc_step,
        ang_min, ang_max, ang_step, xmin, xmax, xstep, ymin, ymax, ystep):
    """
    Process all possible solutions in the defined ranges.
    """
    N = len(xy_st)
    # Ranges
    sc_range = np.arange(sc_min, sc_max, sc_step)
    ang_range = np.arange(ang_min, ang_max, ang_step)
    x_range = np.arange(xmin, xmax, xstep)
    y_range = np.arange(ymin, ymax, ystep)
    print("Total solutions: {:.2e}".format(
          np.prod([len(_) for _ in [sc_range, ang_range, x_range, y_range]])))

    d_sum_all, xy_transl_all = [], []
    # Zoom scale.
    for scale in sc_range:
        # Scaled points.
        x_scale, y_scale = scalePoints(xy_center, delta_x, delta_y, scale)
        for ang in ang_range:
            # Rotated points.
            xy_rot = rotatePoints(xy_center, x_scale, y_scale, ang)
            for x_tr in x_range:
                x_transl = xy_rot[0] + x_tr
                for y_tr in y_range:
                    y_transl = xy_rot[1] + y_tr

                    # Find minimum distance sum.
                    d_sum = distancePoints(xy_st, x_transl, y_transl)

                    d_sum_all.append([d_sum, scale, ang, x_tr, y_tr])
                    xy_transl_all.append([x_transl, y_transl])

                    # Condition to break out if given tolerance is achieved.
                    if d_sum <= tol * N:
                        print("Stars matched:", scale, ang, x_tr, y_tr)
                        return d_sum_all, xy_transl_all

        # Print best solution found so far.
        i_min = zip(*d_sum_all)[0].index(min(zip(*d_sum_all)[0]))
        print(d_sum_all[i_min])

    print("WARNING: no suitable match could be found.")
    return d_sum_all, xy_transl_all


def main():
    """
    """
    # Load Landolt standard file. Notice that y axis is inverted.
    # plt.imshow(plt.imread("../tasks/landolt/pg1323-086.gif"))
    # plt.show()

    # Coordinates of stars for this standard field.
    pg1323 = OrderedDict((
        ['86', (211., 158.3)], ['86A', (162.5, 137.5)], ['86B', (158.1, 128.)],
        ['86C', (160.1, 171.2)], ['86D', (89.6, 133.7)]))
    xy_st = np.array([_ for _ in pg1323.values()])
    print("Landolt frame x,y range:", max(xy_st[0]) - min(xy_st[0]),
          max(xy_st[1]) - min(xy_st[1]))

    # Coordinates from observed frame.
    fname = '../output/standards/filt_U/stk_2085.coo'
    t = Table.read(fname, format='ascii')
    x, y = t['x'], t['y']

    # Center of observed xy points, defined as the center of the minimal
    # rectangle that contains all points.
    xy_center = [(min(x) + max(x)) * .5, (min(y) + max(y)) * .5]
    # Difference between the center coordinates and the xy points.
    delta_x, delta_y = xy_center[0] - x, xy_center[1] - y

    # Move Landolt frame to the center of the obserced image. This makes the
    # reasonable guess that the observed standards field was observed centering
    # the standard stars in it.
    landolt_f_cent = [(min(xy_st[0]) + max(xy_st[0])) * .5,
                      (min(xy_st[1]) + max(xy_st[1])) * .5]
    landolt_delta = np.array(xy_center) - np.array(landolt_f_cent)
    xy_st = xy_st + landolt_delta

    # Tolerance in pixels for match.
    tol = 1.
    # Boundaries.
    sc_min, sc_max, sc_step = .01, 1., .01
    ang_min, ang_max, ang_step = 0., 10., 1.
    xmin, xmax, xstep = -50., 50., 1.
    ymin, ymax, ystep = -50., 50., 1.

    # Find proper zomm + rotation + translation.
    d_sum_all, xy_transl_all = match_frames(
        xy_st, xy_center, delta_x, delta_y, tol, sc_min, sc_max, sc_step,
        ang_min, ang_max, ang_step, xmin, xmax, xstep, ymin, ymax, ystep)

    i_min = zip(*d_sum_all)[0].index(min(zip(*d_sum_all)[0]))
    print(d_sum_all[i_min])

    ax = plt.subplot(111)
    ax.scatter(*xy_center, marker='x', c='g')
    ax.scatter(*zip(*xy_st), c='k', s=10)
    ax.scatter(xy_transl_all[i_min][0], xy_transl_all[i_min][1], c='r', s=10)
    plt.show()


if __name__ == '__main__':
    main()
