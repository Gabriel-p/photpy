
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def rotatePoints(center, xy, angle):
    """
    Rotates points in 'xy' around 'center'. Angle is in degrees.
    Rotation is counter-clockwise

    http://stackoverflow.com/a/20024348/1391441
    """
    angle = math.radians(angle)
    xy_rot = xy[0] - center[0], xy[1] - center[1]
    xy_rot = (xy_rot[0] * math.cos(angle) - xy_rot[1] * math.sin(angle),
              xy_rot[0] * math.sin(angle) + xy_rot[1] * math.cos(angle))
    xy_rot = xy_rot[0] + center[0], xy_rot[1] + center[1]

    return xy_rot


N = 5
xy = [np.random.uniform(0., 1000., 2) for _ in range(N)]
x, y = zip(*xy)

# Center of xy points, defined as the center of the minimal rectangle that
# contains all points.
xy_center = ((min(x) + max(x)) * .5, (min(y) + max(y)) * .5)

# Defference between the center coordinates and the xy points.
delta_x, delta_y = xy_center[0] - x, xy_center[1] - y

# Zoom scale.
scale = 1.5

# Scaled xy points.
x_scale = xy_center[0] - scale * delta_x
y_scale = xy_center[1] - scale * delta_y

# Rotated scaled points.
angle = 5.
xy_rot = rotatePoints(xy_center, zip([x_scale, y_scale]), angle)

ax = plt.subplot(111)
# Center.
ax.scatter(*xy_center, marker='x', c='g')
# Square: bottom left corner, width, height
ax.add_patch(
    patches.Rectangle(
        (min(x), min(y)), (max(x) - min(x)), (max(y) - min(y)), fill=False))
# Original xy points.
ax.scatter(x, y, c='b')
# Zoomed points.
ax.scatter(x_scale, y_scale, c='r')
# Rotated
ax.scatter(*xy_rot, c='g')

plt.show()
