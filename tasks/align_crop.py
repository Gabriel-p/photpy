
import os
import sys
from os.path import join, realpath, dirname, isfile
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import Rectangle
from itertools import cycle
import datetime

from astropy.io import fits
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
    if not isfile(pars_f):
        print("Parameters file is missing. Exit.")
        sys.exit()

    pars = {}
    with open(pars_f, 'r') as f:
        for line in f:
            if not line.startswith('#') and line != '\n':
                key, value = line.replace('\n', '').split()
                pars[key] = value

    folder = pars['in_folder']
    folder = folder[1:] if folder.startswith('/') else folder
    in_path = join(mypath.replace('tasks', 'input'), folder)

    fits_list = []
    if os.path.isdir(in_path):
        for file in os.listdir(in_path):
            f = join(in_path, file)
            if isfile(f):
                if f.endswith('.fits'):
                    fits_list.append(f)
    else:
        print("{}\nis not a folder. Exit.".format(in_path))
        sys.exit()

    if pars['ref_align'] not in ['none', 'None', None]:
        ref_im_f = join(in_path, pars['ref_align'])
        if not os.path.isfile(ref_im_f):
            print("{}\nreference image is not present. Exit.".format(ref_im_f))
            sys.exit()
        else:
            pars['ref_align'] = ref_im_f

    # Create path to output folder
    out_path = in_path.replace('input', 'output')

    return mypath, fits_list, pars, out_path


def get_coords_data(imname):
    """
    Read stars coordinates from .coo files.
    """
    fn = imname.replace('input', 'output').replace('fits', 'coo')
    coords_data = []
    try:
        with open(fn, 'r') as f:
            for l in f:
                if not l.startswith('#'):
                    coords_data.append(map(float, l.split())[:2])
    except IOError:
        print("{}\nFile not found.".format(fn))
        sys.exit()

    if not coords_data:
        print("\nERROR: no stars in .coo file:\n{}".format(fn))
        sys.exit()

    return zip(*coords_data)


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

    d = dists[min_idx]
    print("Reg shifted by: {:.2f}, {:.2f}".format(d[0], d[1]))
    print("Median average distance: {:.2f}".format(d[2]))

    return d


def overlap_reg(h, w, shifts):
    """
    Obtain edges of overlap.
    """
    lef = max(zip(*shifts)[0])
    bot = max(zip(*shifts)[1])
    rig = min(zip(*shifts)[0]) + w
    top = min(zip(*shifts)[1]) + h
    # Width and height of overlap.
    xbox, ybox = rig - lef, top - bot
    # Center of overlap.
    xcen, ycen = (rig + lef) / 2., (top + bot) / 2.
    print("\nOverlapping area")
    print("Center: {:.2f}, {:.2f}".format(xcen, ycen))
    print("Box: {:.2f}, {:.2f}".format(xbox, ybox))

    return lef, bot, rig, top, xcen, ycen, xbox, ybox


def crop_frames(fits_list, pars, hdu, hdr, shifts, overlap):
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
            crop_name = join(fits_list[i].replace('input', 'output').replace(
                             '.fits', '_crop.fits'))
            try:
                os.remove(crop_name)
            except OSError:
                pass
            fits.writeto(crop_name, frame_crop.data, hdr[i])

    return hdu_crop


def make_sub_plot(ax, frame, fname, xy, shifts, overlap, n_ref):
    lef, bot, rig, top, xcen, ycen, xbox, ybox = overlap
    # Zscale
    interval = ZScaleInterval()
    zmin, zmax = interval.get_limits(frame)
    plt.imshow(frame, cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower', vmin=zmin, vmax=zmax)
    ax.set_title('{}{} ({} stars)'.format(
        n_ref, fname.split('/')[-1], len(xy[0])), fontsize=8)
    positions = (xy[0] - lef + shifts[0], xy[1] - bot + shifts[1])
    apertures = CircularAperture(positions, r=20.)
    apertures.plot(color='r', lw=0.5)
    txt = AnchoredText("Median dist: {:.2f}".format(shifts[2]), loc=1,
                       prop=dict(size=8))
    txt.patch.set(boxstyle='square,pad=0.', alpha=0.85)
    ax.add_artist(txt)


def make_plots(out_path, hdu, ref_i, fits_list, xy_coo, shifts, overlap,
               hdu_crop):
    """
    Make plots.
    """
    print("\nPlotting.")
    fig = plt.figure(figsize=(20, 20))
    p = int(round(np.sqrt(len(hdu)))) + 1
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
    for i, xy in enumerate(xy_coo):
        positions = (xy[0] - lef + shifts[i][0], xy[1] - bot + shifts[i][1])
        lbl = fits_list[i].split('/')[-1]
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
        ax, hdu_crop[ref_i], fits_list[ref_i], xy_coo[ref_i], shifts[ref_i],
        overlap, 'Reference: ')
    for _ in [hdu_crop, fits_list, xy_coo, shifts]:
        del _[ref_i]
    for i, frame in enumerate(hdu_crop):
        ax = fig.add_subplot(gs[i + 3])
        make_sub_plot(
            ax, frame, fits_list[i], xy_coo[i], shifts[i], overlap, '')

    fig.tight_layout()
    fn = join(out_path, 'align_crop.png')
    plt.savefig(fn, dpi=150, bbox_inches='tight')
    plt.clf()
    plt.close()


def main():
    """
    """
    mypath, fits_list, pars, out_path = read_params()

    print("Read coordinates from .coo files.")
    xy_coo, hdu, hdr, hw = [], [], [], [[], []]
    # For each .fits image in the root folder.
    for i, imname in enumerate(fits_list):
        # Extract frame data.
        hdulist = fits.open(imname)
        hdu_data = hdulist[0].data
        hdu.append(hdu_data)
        # Height and width (h=y, w=x)
        h, w = np.shape(hdu_data)
        hw[0].append(h)
        hw[1].append(w)
        hdr.append(hdulist[0].header)
        hdulist.close()

        # Obtain stars coordinates.
        xy_coo.append(get_coords_data(imname))

    # Check that all observed frames are of the same size.
    if hw[0][1:] == hw[0][:-1] and hw[1][1:] == hw[1][:-1]:
        h, w = hw[0][0], hw[1][0]
    else:
        print("ERROR: frames are not all of the same size.")
        sys.exit()

    # Identify reference frame.
    if pars['ref_align'] not in ['none', 'None', None]:
        ref_i = fits_list.index(pars['ref_align'])
    else:
        # Else, select the one with the most stars.
        N_coo = [len(_[0]) for _ in xy_coo]
        ref_i = N_coo.index(max(N_coo))
    ref = xy_coo[ref_i]
    print('Reference image: {}'.format(fits_list[ref_i].split('/')[-1]))

    # Find x,y shifts (to reference frame) for alignment.
    shifts = []
    for i, xy in enumerate(xy_coo):
        if i != ref_i:
            print("\nFile: {}".format(fits_list[i].split('/')[-1]))
            # Obtain shifts
            d = avrg_dist(
                (float(pars['x_init_shift']), float(pars['y_init_shift'])),
                float(pars['max_shift']), float(pars['tolerance']), ref, xy)
            shifts.append(d)
        else:
            # Shift for reference frame.
            shifts.append([0., 0., 0.])

    # Obtain overlapping region.
    overlap = overlap_reg(h, w, shifts)

    # Crop frames.
    hdu_crop = crop_frames(fits_list, pars, hdu, hdr, shifts, overlap)

    if pars['do_plots_B'] == 'y':
        make_plots(
            out_path, hdu, ref_i, fits_list, xy_coo, shifts, overlap, hdu_crop)


if __name__ == "__main__":
    main()
