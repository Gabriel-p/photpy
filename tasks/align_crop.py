
import os
import sys
from os.path import join, realpath, dirname
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from itertools import cycle

from getdata import st_fwhm_select
import psfmeasure

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval
from photutils import CircularAperture
from photutils.utils import cutout_footprint


def create_pars_file(pars_f, pars_list=None):
    """
    """
    # Default values.
    if pars_list is None:
        pars_list = [
            'None', 'None', 'True', '5.', '3.', '60000.', '0.15', '1.5', '100',
            '0.', '0.', '-1.', '0.05']
    with open(pars_f, 'w') as f:
        f.write(
            "# Default parameters for the align_crop.py script\n#\n"
            "ff_proc {}\nref_im {}\ndo_plots {}\nthresh_level {}\n"
            "fwhm_init {}\ndmax {}\nellip_max {}\nfwhm_min {}\n"
            "max_shift_stars {}\nx_init_shift {}\ny_init_shift {}\n"
            "max_shift {}\ntol {}\n".format(*pars_list))
    return


def read_params():
    """
    """
    pars = {}
    mypath = realpath(join(os.getcwd(), dirname(__file__)))
    pars_f = join(mypath, '.align_crop.pars')
    if not os.path.isfile(pars_f):
        print("Parameters file missing. Create it.")
        create_pars_file(pars_f)

    with open(pars_f, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                key, value = line.replace('\n', '').split(' ')
                pars[key] = value

    return mypath, pars_f, pars


def get_params(mypath, pars_f, pars):
    """
    """
    pars_list = []

    answ = raw_input("\nDir or file to process ({}): ".format(
        pars['ff_proc']))
    fname = str(answ) if answ is not '' else pars['ff_proc']
    fname = fname[1:] if fname.startswith('/') else fname
    r_path = join(mypath.replace('/tasks', ''), fname)
    fits_list = []
    if os.path.isdir(r_path):
        for subdir, dirs, files in os.walk(r_path):
            for file in files:
                if file.endswith('.fits'):
                    fits_list.append(os.path.join(subdir, file))
    elif os.path.isfile(r_path):
        print("  Single .fits file found in folder.")
        sys.exit()
    else:
        print("{}\nis neither a folder nor a file. Exit.".format(str(r_path)))
        sys.exit()
    pars['ff_proc'] = fname
    pars_list.append(pars['ff_proc'])

    answ = raw_input("Reference image ({}): ".format(pars['ref_im']))
    if answ in ['None', 'none']:
        pars['ref_im'] = 'none'
    elif answ is not '':
        pars['ref_im'] = str(answ)
        try:
            os.path.isfile(join(r_path, pars['ref_im']))
        except:
            print("Reference image does not exist. Exit.")
    pars_list.append(pars['ref_im'])

    answ = raw_input("Create plots? (y/n) ({}): ".format(
        pars['do_plots']))
    pars['do_plots'] = False if answ in ('n', 'N', 'no', 'NO') else True
    pars_list.append(pars['do_plots'])

    answ = raw_input("Threshold level above STDDEV ({}): ".format(
        pars['thresh_level']))
    pars['thresh_level'] = float(answ) if answ is not '' else\
        float(pars['thresh_level'])
    pars_list.append(pars['thresh_level'])

    answ = raw_input("Initial FWHM ({}): ".format(pars['fwhm_init']))
    pars['fwhm_init'] = float(answ) if answ is not '' else\
        float(pars['fwhm_init'])
    pars_list.append(pars['fwhm_init'])

    answ = raw_input("Maximum flux ({}): ".format(pars['dmax']))
    pars['dmax'] = float(answ) if answ is not '' else float(pars['dmax'])
    pars_list.append(pars['dmax'])

    answ = raw_input("Maximum ellipticity ({}): ".format(pars['ellip_max']))
    pars['ellip_max'] = float(answ) if answ is not '' else\
        float(pars['ellip_max'])
    pars_list.append(pars['ellip_max'])

    answ = raw_input("Minimum FWHM ({}): ".format(pars['fwhm_min']))
    pars['fwhm_min'] = float(answ) if answ is not '' else\
        float(pars['fwhm_min'])
    pars_list.append(pars['fwhm_min'])

    answ = raw_input("Max number of stars used to obtain "
                     "shift ({}): ".format(pars['max_shift_stars']))
    pars['max_shift_stars'] = int(answ) if answ is not '' else\
        int(pars['max_shift_stars'])
    pars_list.append(pars['max_shift_stars'])

    answ = raw_input("Initial estimated shift in x ({}): ".format(
        pars['x_init_shift']))
    pars['x_init_shift'] = float(answ) if answ is not '' else\
        float(pars['x_init_shift'])
    pars_list.append(pars['x_init_shift'])

    answ = raw_input("Initial estimated shift in y ({}): ".format(
        pars['y_init_shift']))
    pars['y_init_shift'] = float(answ) if answ is not '' else\
        float(pars['y_init_shift'])
    pars_list.append(pars['y_init_shift'])

    answ = raw_input("Maximum shift allowed ({}): ".format(
        pars['max_shift']))
    pars['max_shift'] = float(answ) if answ is not '' else\
        float(pars['max_shift'])
    pars_list.append(pars['max_shift'])

    answ = raw_input("Matching tolerance ({}): ".format(pars['tol']))
    pars['tol'] = float(answ) if answ is not '' else float(pars['tol'])
    pars_list.append(pars['tol'])

    # Write values to file.
    create_pars_file(pars_f, pars_list)

    return r_path, fits_list, pars


def avrg_dist(init_shift, max_shift, tol, ref, f):
    """
    Average minimal (Euclidean) distance from points in 'f' to points in
    'ref', until the stopping condition.

    Source: http://codereview.stackexchange.com/a/134918/35351
    """
    # Shift.
    x0, y0 = init_shift[0], init_shift[1]
    if max_shift < 1.:
        # Use the full length of both sides.
        max_shift = max(max(ref[0]), max(ref[1]))

    while max_shift >= 1.:
        dists = []
        li = 10
        # Shift in x
        for sx in np.linspace(-1. * max_shift, max_shift, li):
            # Shift in y
            for sy in np.linspace(-1. * max_shift, max_shift, li):
                # Apply possible x,y translation
                sf = [f[0] + sx + x0, f[1] + sy + y0]
                # Average minimal distance.
                d = np.median(
                    np.min(distance.cdist(zip(*sf), zip(*ref)), axis=1))
                # Store x,y shift values, and the average minimal distance.
                dists.append([sx + x0, sy + y0, d])

        # Index of the x,y shifts that resulted in the average minimal
        # distance.
        min_idx = np.array(zip(*dists)[2]).argmin()
        x0, y0 = dists[min_idx][0], dists[min_idx][1]
        # Decrease max shift by X%
        max_shift = (1. - tol) * max_shift
        # Break condition: minimum accuracy reached.
        if max_shift < 1.:
            if dists[min_idx][2] > 5.:
                print("  WARNING: match is probably wrong.")

    return dists[min_idx]


def overlap_reg(hdu, shifts):
    """
    """
    # Height and width (h=y, w=x)
    h, w = np.shape(hdu)
    # Obtain edges of overlap.
    lef = max(zip(*shifts)[0])
    bot = max(zip(*shifts)[1])
    rig = min(zip(*shifts)[0]) + w
    top = min(zip(*shifts)[1]) + h
    # Width and height of overlap.
    xbox, ybox = rig - lef, top - bot
    # Center of overlap.
    xcen, ycen = (rig + lef) / 2., (top + bot) / 2.
    print("\nOverlapping area\nCenter: {:.2f}, {:.2f}".format(xcen, ycen))
    print("Box: {:.2f}, {:.2f}".format(xbox, ybox))

    return lef, bot, rig, top, xcen, ycen, xbox, ybox


def make_sub_plot(ax, hdu, fname, f_list, shifts, overlap, n_ref):
    # Overlap info.
    lef, bot, rig, top, xcen, ycen, xbox, ybox = overlap
    # Re-center
    xcs, ycs = xcen - shifts[0], ycen - shifts[1]
    # Crop image: (xc, yc), (y length, x length)
    hdu_crop = cutout_footprint(hdu, (xcs, ycs), (ybox, xbox))[0]
    # hdu_crop.writeto('crop_' + str(i) + '_fig.fits')
    # Zscale
    interval = ZScaleInterval()
    zmin, zmax = interval.get_limits(hdu_crop)
    plt.imshow(hdu_crop, cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower', vmin=zmin, vmax=zmax)
    ax.set_title('{}{} ({} stars)'.format(n_ref, fname, len(f_list[0])),
                 fontsize=8)
    positions = (f_list[0] - lef + shifts[0], f_list[1] - bot + shifts[1])
    apertures = CircularAperture(positions, r=20.)
    apertures.plot(color='r', lw=0.5)
    # ax.tick_params(axis='both', which='major', labelsize=8)


def make_plots(mypath, hdu, ref_i, fnames, f_list, shifts, overlap):
    """
    Make plots.
    """
    print("\nPlotting.")
    fig = plt.figure(figsize=(20, 20))
    p = int(np.sqrt(len(hdu))) + 1
    gs = gridspec.GridSpec(p, p)
    fn = join(mypath.replace('tasks', ''), '/'.join(fnames[0].split('/')[:-1]))

    ax = fig.add_subplot(gs[0])
    # Height and width (h=y, w=x)
    h, w = np.shape(hdu[ref_i])
    lef, bot, rig, top, xcen, ycen, xbox, ybox = overlap
    plt.scatter(w / 2., h / 2.)
    cols = ['r', 'b', 'g', 'm', 'c']
    col_cyc = cycle(cols)
    for s_xy in shifts:
        ax.add_patch(
            Rectangle((s_xy[0], s_xy[1]), h, w, fill=None, alpha=1,
                      color=next(col_cyc)))
    ax.add_patch(Rectangle((0., 0.), h, w, fill=None, alpha=1, lw=2.,
                           color='k'))
    # Overlap area.
    ax.add_patch(
        Rectangle((lef, bot), rig - lef, top - bot, fill=1, alpha=.25,
                  color=next(col_cyc)))
    plt.scatter(xcen, ycen, c='r', marker='x', zorder=5)
    ax.set_title("Overlap region", fontsize=8)

    # Selected stars.
    ax = fig.add_subplot(gs[1])
    for i, fl in enumerate(f_list):
        positions = (fl[0] - lef + shifts[i][0], fl[1] - bot + shifts[i][1])
        plt.scatter(
            positions[0], positions[1], label=fnames[i].split('/')[-1], lw=0.5,
            s=20., facecolors='none', edgecolors=next(col_cyc))
    plt.legend(loc='upper right', fontsize=8, scatterpoints=1, framealpha=0.85)
    ax.set_title("Stars used for alignment", fontsize=8)
    ax.set_xlim(0., w)
    ax.set_ylim(0., h)

    # Aligned and cropped frames.
    for i, hdu_data in enumerate(hdu):
        ax = fig.add_subplot(gs[i + 2])
        n_ref = '' if i != ref_i else 'Reference - '
        make_sub_plot(
            ax, hdu_data, fnames[i], f_list[i], shifts[i], overlap, n_ref)

    fig.tight_layout()
    plt.savefig(fn + '/align_crop.png', dpi=150, bbox_inches='tight')
    plt.clf()
    plt.close()


def main():
    """
    """
    mypath, pars_f, pars = read_params()
    r_path, fits_list, pars = get_params(mypath, pars_f, pars)

    f_list, l_f, hdu = [], [], []
    # For each .fits image in the root folder.
    for i, imname in enumerate(fits_list):
        print("\nFile: {}".format(
            imname.replace(mypath.replace('tasks', ''), "")))

        # Extract data.
        hdulist = fits.open(imname)
        hdu_data = hdulist[0].data

        # Background estimation.
        sky_mean, sky_median, sky_std = sigma_clipped_stats(
            hdu_data, sigma=3.0, iters=2)

        # Stars selection.
        psf_select, all_sources, n_not_satur = st_fwhm_select(
            pars['dmax'], pars['max_shift_stars'], pars['thresh_level'],
            pars['fwhm_init'], sky_std, hdu_data)

        # FWHM selection.
        fwhm_estim, psfmeasure_estim, fwhm_min_rjct, ellip_rjct =\
            psfmeasure.main(
                pars['dmax'], pars['ellip_max'], pars['fwhm_min'],
                psf_select, imname, hdu_data)
        print("Suitable sources found: {}".format(len(fwhm_estim)))

        if fwhm_estim:
            l_f.append(all_sources)
            hdu.append(hdu_data)
            f_list.append(np.array([zip(*fwhm_estim)[0], zip(*fwhm_estim)[1]]))
        else:
            print("\nWARNING: no stars left after rejecting\nby min FWHM"
                  "and max ellipticity.")

    # Identify reference frame.
    if pars['ref_im'] not in ['none', 'None', None]:
        ref_i = [_.split('/')[-1] for _ in fits_list].index(pars['ref_im'])
    else:
        ref_i = l_f.index(max(l_f))
    ref = f_list[ref_i]
    print('\nReference image: {}'.format(fits_list[ref_i].split('/')[-1]))

    # Find x,y shifts (to reference frame) for alignment.
    shifts, fnames = [], []
    for i, f in enumerate(f_list):
        file = fits_list[i].replace(mypath.replace('tasks', ''), "")
        fnames.append(file)
        if i != ref_i:
            print("\nFile: {}".format(file))
            # Obtain shifts
            d = avrg_dist((pars['x_init_shift'], pars['y_init_shift']),
                          pars['max_shift'], pars['tol'], ref, f)
            print("Reg shifted by: {:.2f}, {:.2f}".format(d[0], d[1]))
            print("Median average distance: {:.2f}".format(d[2]))
            shifts.append([d[0], d[1]])
        else:
            # Shift for reference frame.
            shifts.append([0., 0.])

    # Obtain overlapping region.
    overlap = overlap_reg(hdu_data, shifts)

    if pars['do_plots']:
        make_plots(mypath, hdu, ref_i, fnames, f_list, shifts, overlap)


if __name__ == "__main__":
    main()
