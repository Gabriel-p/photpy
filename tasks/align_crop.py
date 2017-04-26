
import os
import sys
from os.path import join, realpath, dirname
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from itertools import cycle
import datetime

from hlpr import st_fwhm_select

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval
from photutils import CircularAperture
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D


def read_params():
    """
    Read parameter values.
    """
    mypath = realpath(join(os.getcwd(), dirname(__file__)))
    pars_f = join(mypath.replace('tasks', ''), 'params_input.dat')
    if not os.path.isfile(pars_f):
        print("Parameters file is missing. Exit.")
        sys.exit()

    pars = {}
    with open(pars_f, 'r') as f:
        for line in f:
            if not line.startswith('#') and line != '\n':
                key, value = line.replace('\n', '').split()
                pars[key] = value

    fname = pars['ff_crop']
    fname = fname[1:] if fname.startswith('/') else fname
    in_path = join(mypath.replace('tasks', 'input'), fname)

    if pars['ref_im'] not in ['none', 'None', None]:
        ref_im_f = join(in_path, pars['ref_im'])
        if not os.path.isfile(ref_im_f):
            print("{}\nreference image is not a file. Exit.".format(ref_im_f))
            sys.exit()
        else:
            pars['ref_im'] = ref_im_f

    fits_list = []
    if os.path.isdir(in_path):
        for subdir, dirs, files in os.walk(in_path):
            for file in files:
                if file.endswith('.fits'):
                    fits_list.append(os.path.join(subdir, file))
    else:
        print("{}\nis not a folder. Exit.".format(in_path))
        sys.exit()

    return mypath, fits_list, pars


def get_coords_data(imname, pars, hdu_data):
    """
    """
    print("Obtaining stars coordinates for aligning.")
    coords_flag = False
    if pars['read_coords'] == 'y':
        try:
            fn = imname.replace('input', 'output').replace('fits', 'coo')
            fwhm_accptd = []
            with open(fn, 'r') as f:
                for l in f:
                    if not l.startswith('#'):
                        fwhm_accptd.append(map(float, l.split()))
            all_sources = len(fwhm_accptd)
        except IOError:
            coords_flag = True
    else:
        coords_flag = True

    if coords_flag:
        # Background estimation.
        sky_mean, sky_median, sky_std = sigma_clipped_stats(
            hdu_data, sigma=3.0, iters=2)

        # Stars selection.
        fwhm_accptd, all_sources, n_not_satur = st_fwhm_select(
            float(pars['dmax']), int(pars['max_stars']),
            float(pars['thresh_level']), float(pars['fwhm_init']),
            sky_std, hdu_data)

    return fwhm_accptd, all_sources


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


def crop_frames(mypath, fnames, pars, hdu, hdr, shifts, overlap):
    """
    """
    # Overlap info.
    lef, bot, rig, top, xcen, ycen, xbox, ybox = overlap
    # Crop image: (xc, yc), (y length, x length)
    hdu_crop = []
    for i, frame in enumerate(hdu):
        # Re-center
        xcs, ycs = xcen - shifts[i][0], ycen - shifts[i][1]
        # Crop image
        frame_crop = Cutout2D(frame, (xcs, ycs), (ybox, xbox), wcs=WCS(hdr[i]))
        # Save for plotting.
        hdu_crop.append(frame_crop.data)
        # Save cropped fits.
        if pars['crop_save'] is 'y':
            # Cropped WCS
            wcs_cropped = frame_crop.wcs
            # Update WCS in header
            hdr[i].update(wcs_cropped.to_header())
            # Add comment to header
            hdr[i]['COMMENT'] = "= Cropped fits file ({}).".format(
                datetime.date.today())
            # Write cropped frame to new fits file.
            crop_name = join(
                mypath.replace('tasks', ''),
                fnames[i].replace('input', 'output').replace(
                    '.fits', '_crop.fits'))
            try:
                os.remove(crop_name)
            except OSError:
                pass
            fits.writeto(crop_name, frame_crop.data, hdr[i])

    return hdu_crop


def make_sub_plot(ax, frame, fname, f_list, shifts, overlap, n_ref):
    lef, bot, rig, top, xcen, ycen, xbox, ybox = overlap
    # Zscale
    interval = ZScaleInterval()
    zmin, zmax = interval.get_limits(frame)
    plt.imshow(frame, cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower', vmin=zmin, vmax=zmax)
    ax.set_title('{}{} ({} stars)'.format(
        n_ref, fname.replace('input/', ''), len(f_list[0])), fontsize=8)
    positions = (f_list[0] - lef + shifts[0], f_list[1] - bot + shifts[1])
    apertures = CircularAperture(positions, r=20.)
    apertures.plot(color='r', lw=0.5)


def make_plots(
        mypath, out_folder, hdu, ref_i, fnames, f_list, shifts, overlap,
        hdu_crop):
    """
    Make plots.
    """
    print("\nPlotting.")
    fig = plt.figure(figsize=(20, 20))
    p = int(np.sqrt(len(hdu))) + 1
    gs = gridspec.GridSpec(p, p)

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
    maxl = 10
    for i, fl in enumerate(f_list):
        positions = (fl[0] - lef + shifts[i][0], fl[1] - bot + shifts[i][1])
        lbl = fnames[i].split('/')[-1]
        plt.scatter(
            positions[0], positions[1], label=(i // maxl) * "_" + lbl,
            lw=0.5, s=20., facecolors='none', edgecolors=next(col_cyc))
    plt.legend(loc='upper right', fontsize=8, scatterpoints=1, framealpha=0.85)
    ax.set_title("Stars used for alignment", fontsize=8)
    ax.set_xlim(0., w)
    ax.set_ylim(0., h)

    # Aligned and cropped frames.
    ax = fig.add_subplot(gs[2])
    make_sub_plot(
        ax, hdu_crop[ref_i], fnames[ref_i], f_list[ref_i], shifts[ref_i],
        overlap, 'Reference - ')
    for _ in [hdu_crop, fnames, f_list, shifts]:
        del _[ref_i]
    for i, frame in enumerate(hdu_crop):
        ax = fig.add_subplot(gs[i + 3])
        make_sub_plot(
            ax, frame, fnames[i], f_list[i], shifts[i], overlap, '')

    fig.tight_layout()
    fn = join(mypath.replace('tasks', 'output'), out_folder, 'align_crop.png')
    plt.savefig(fn, dpi=150, bbox_inches='tight')
    plt.clf()
    plt.close()


def main():
    """
    """
    mypath, fits_list, pars = read_params()

    f_list, l_f, hdu, hdr = [], [], [], []
    # For each .fits image in the root folder.
    for i, imname in enumerate(fits_list):
        print("\nFile: {}".format(
            imname.replace(mypath.replace('tasks', 'input'), "")))

        # Extract frame data.
        hdulist = fits.open(imname)
        hdu_data = hdulist[0].data
        # Header
        header = hdulist[0].header
        hdulist.close()

        # Obtain stars coordinates.
        coords_data, n_sources = get_coords_data(imname, pars, hdu_data)

        if coords_data:
            l_f.append(n_sources)
            hdu.append(hdu_data)
            hdr.append(header)
            f_list.append(
                np.array([zip(*coords_data)[0], zip(*coords_data)[1]]))
        else:
            print("\nWARNING: no stars left after rejecting\nby min FWHM"
                  "and max ellipticity.")

    # Identify reference frame.
    if pars['ref_im'] not in ['none', 'None', None]:
        ref_i = fits_list.index(pars['ref_im'])
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
            print("\nFile: {}".format(file.replace('input', '')))
            # Obtain shifts
            d = avrg_dist(
                (float(pars['x_init_shift']), float(pars['y_init_shift'])),
                float(pars['max_shift']), float(pars['tolerance']), ref, f)
            print("Reg shifted by: {:.2f}, {:.2f}".format(d[0], d[1]))
            print("Median average distance: {:.2f}".format(d[2]))
            shifts.append([d[0], d[1]])
        else:
            # Shift for reference frame.
            shifts.append([0., 0.])

    # Obtain overlapping region.
    overlap = overlap_reg(hdu_data, shifts)

    # Crop frames.
    hdu_crop = crop_frames(mypath, fnames, pars, hdu, hdr, shifts, overlap)

    if pars['do_plots_B'] == 'y':
        make_plots(
            mypath, pars['ff_crop'], hdu, ref_i, fnames, f_list, shifts,
            overlap, hdu_crop)


if __name__ == "__main__":
    main()
