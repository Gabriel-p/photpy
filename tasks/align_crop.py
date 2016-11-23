
import numpy as np
from scipy.spatial import distance


def main():

    # Generate main figure, and shifted sub-figures with fewer points.
    pts, trs, figs = gen_data()

    # Obtain optimal shifts for each sub-figure.
    shifts = []
    for i, f in enumerate(figs):
        print('\nTrue x,y shift: {}'.format(str(trs[i])))
        d = avrg_dist(pts, f)
        print('Found x,y shift: {}'.format((d[0], d[1])))
        # I have to add instead of substract
        print('Shift dist: {:.2f}'.format(
            np.sqrt((trs[i][0] + d[0]) ** 2 + (trs[i][1] + d[1]) ** 2)))
        shifts.append([d[0], d[1]])

    make_plot(pts, figs, shifts)


if __name__ == "__main__":
    main()
