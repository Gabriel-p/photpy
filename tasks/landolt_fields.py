
from astropy.table import Table


def main(field):
    """
    Dictionary of Landolt standard stars.

    http://www.eso.org/sci/observing/tools/standards/Landolt.html
    """

    # Load Landolt standard file. Notice that y axis is inverted.
    # import matplotlib.pyplot as plt
    # plt.imshow(plt.imread("landolt/pg1323.gif"))
    # plt.show()

    #  x  y  V  B-V    U-B    V-R    R-I    V-I
    data_rows = [
        ('86', 211., 158.3, 13.481, -0.14, -0.681, -0.048, -0.078, -0.127),
        ('86A', 162.5, 137.5, 13.591, 0.393, -0.019, 0.252, 0.252, 0.506),
        ('86B', 158.1, 128., 13.406, 0.761, 0.265, 0.426, 0.407, 0.833),
        ('86C', 160.1, 171.2, 14.003, 0.707, 0.245, 0.395, 0.363, 0.759),
        ('86D', 89.6, 133.7, 12.080, 0.587, 0.005, 0.346, 0.335, 0.684)]
    pg1323 = Table(
        rows=data_rows,
        names=('ID', 'x', 'y', 'V', 'BV', 'UB', 'VR', 'RI', 'VI'))

    all_fields = {'pg1323': pg1323}

    return all_fields[field]


if __name__ == '__main__':
    main('pg1323')
