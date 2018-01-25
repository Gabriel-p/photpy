
import read_pars_file as rpf
from match_funcs import triangleMatch, getTriangles, findRotAngle,\
    standard2observed, xyTrans
from plot_events import zoom

import landolt_fields
import os
from os.path import join, isfile, exists
from pathlib2 import Path
import sys
import numpy as np
from scipy.spatial.distance import cdist

from hlpr import bckg_data, st_fwhm_select

from astropy.table import Table, hstack
from astropy.io import ascii, fits
from photutils.utils import cutout_footprint
from photutils import detect_threshold, detect_sources
from photutils import source_properties
from astropy.convolution import Gaussian2DKernel

from astropy.visualization import ZScaleInterval
from photutils import CircularAperture
import matplotlib.pyplot as plt
import time
import matplotlib.gridspec as gridspec


def in_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    in_path = join(pars['mypath'].replace('tasks', 'input'))
    out_path = in_path.replace('input', 'output')

    landolt_fld, ref_frame, match_fldr = [], [], []
    for line in pars['ref_fld_fldr']:
        # Identify if this is a Landolt field or an observed frame.
        if not line[0].endswith(".fits"):
            # A Landolt field was set. We are matching standard frames.
            landolt_fld.append(line[0])
            ref_frame.append('--')
            match_fldr.append(line[1])
        else:
            landolt_fld.append('--')
            ref_frame.append(line[0])
            match_fldr.append(line[1])

    pars['landolt_fld'], pars['ref_frame'], pars['match_fldr'] =\
        landolt_fld, ref_frame, match_fldr

    # Generate full path to reference image, and check existence.
    for i, ref_im in enumerate(pars['ref_frame']):
        if ref_im != '--':
            ref_frame = join(in_path, pars['match_fldr'][i], ref_im)
            if not os.path.isfile(ref_frame):
                print("{}\n Reference frame is not present. Exit.".format(
                    ref_frame))
                sys.exit()
            else:
                # Store full path.
                pars['ref_frame'][i] = ref_frame

    # Generate list of fits files for each input folder.
    fits_list = []
    for folder in pars['match_fldr']:
        folder = folder.strip('/')
        in_path = join(pars['mypath'].replace('tasks', 'input'), folder)

        list_temp = []
        if os.path.isdir(in_path):
            for file in os.listdir(in_path):
                f = join(in_path, file)
                if isfile(f):
                    if f.endswith('.fits'):
                        list_temp.append(f)
        if not list_temp:
            print("{}\n No .fits files found in match folder. Exit.".format(
                in_path))
            sys.exit()

        # Store list for this folder.
        print("Files found in '{}' folder:".format(folder))
        for fit in list_temp:
            print(" * {}".format(fit.replace(in_path, '')[1:]))
        fits_list.append(list_temp)

    return pars, fits_list, out_path


def reCenter(hdu_data, positions, side=30):
    """
    Find better center coordinates for selected stars.
    """
    xy_cent = []
    for x0, y0 in positions:
        # Check that position falls inside of the observed frame.
        if 0. < x0 < hdu_data.shape[1] and 0. < y0 < hdu_data.shape[0]:
            crop = cutout_footprint(hdu_data, (x0, y0), side)[0]

            threshold = detect_threshold(crop, snr=10)
            sigma = 3.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))   # FWHM = 3
            kernel = Gaussian2DKernel(sigma)
            kernel.normalize()
            segm = detect_sources(crop, threshold, npixels=5,
                                  filter_kernel=kernel)
            try:
                tbl = source_properties(crop, segm).to_table()
                # Index of the brightest detected source
                idx_b = tbl['source_sum'].argmax()
                # Index of the closest source.
                sh, d_old, idx_c = side / 2., 1.e6, 0
                for i, xy1 in enumerate(tbl['xcentroid', 'ycentroid']):
                    x, y = xy1[0].value, xy1[1].value
                    d = np.sqrt((sh - x) ** 2 + (sh - y) ** 2)
                    tbl['source_sum']
                    if d < d_old:
                        idx_c = i
                        d_old = 1. * d
                # Choose the brightest star except when its distance to the
                # center of the region is N or more that of the closest star.
                N = 4.
                if np.sqrt(
                        (sh - tbl['xcentroid'][idx_b].value) ** 2 +
                        (sh - tbl['ycentroid'][idx_b].value) ** 2) > N * d_old:
                    idx = idx_c
                else:
                    idx = idx_b
                xy1 = tbl['xcentroid'][idx].value, tbl['ycentroid'][idx].value
                x, y = (x0 + xy1[0] - side * .5, y0 + xy1[1] - side * .5)
            except ValueError:
                # "SourceCatalog contains no sources" may happen when no
                # stars are found within the cropped region.
                x, y = x0, y0
        else:
            x, y = x0, y0

        xy_cent.append((x, y))

    return np.asarray(xy_cent)


def coo_ref_frame(fr, id_selec, xy_selec, pars, proc_grps):
    """
    Mark sources on 'fr' fits file. This uses parameter values from other
    tasks:

    General: gain_key, rdnoise_key, dmax
    fitstats: sky_method, fwhm_init, thresh_fit

    """
    # Load .fits file.
    hdulist = fits.open(fr)
    hdu_data = hdulist[0].data

    answ, idx_choice = 'n', 0
    if id_selec and pars['match_mode'] != 'xyshift':

        for i, xy in enumerate(xy_selec):
            # Re-center coordinates.
            xy = reCenter(hdu_data, xy, side=40).tolist()

            # Manual selection of reference stars in the observed frame.
            landolt_field_img = join(
                pars['mypath'], 'landolt',
                pars['landolt_fld'][proc_grps] + '.gif')
            # Manual selection of reference stars in the observed frame.
            fig, (ax1, ax2) = plt.subplots(1, 2)
            ax1.imshow(plt.imread(landolt_field_img))
            # fig, ax = plt.subplots()
            interval = ZScaleInterval()
            zmin, zmax = interval.get_limits(hdu_data)
            ax2.imshow(
                hdu_data, cmap='Greys', aspect=1, interpolation='nearest',
                origin='lower', vmin=zmin, vmax=zmax)
            apertures = CircularAperture(xy, r=10.)
            apertures.plot(color='red', lw=1.)
            plt.title("Coordinates ({}) over {}".format(
                i, fr.split('/')[-1].split('.')[0]))
            xmax, ymax = hdu_data.shape[1], hdu_data.shape[0]
            for j, txt in enumerate(id_selec[i]):
                xtxt, ytxt = xy[j][0] + .01 * xmax,\
                    xy[j][1] + .01 * ymax
                ax2.annotate(
                    txt, xy=(xy[j][0], xy[j][1]),
                    xytext=(xtxt, ytxt),
                    fontsize=15, color='r', arrowprops=dict(
                        arrowstyle="-", color='b',
                        connectionstyle="angle3,angleA=90,angleB=0"))
            fig.canvas.mpl_connect(
                'scroll_event', lambda event: zoom(event, [ax1, ax2]))
            plt.show()
            plt.close('all')

            while True:
                answ = raw_input(
                    "Use coordinates ({})? (y/n): ".format(i)).strip()
                if answ == 'y':
                    idx_choice = i
                    # Re-write wit re-centered coordinates.
                    xy_selec[idx_choice] = xy
                    break
                elif answ == 'n':
                    break
                else:
                    print("Wrong answer. Try again.")

    # Select new coordinates
    if answ != 'y':

        if pars['match_mode'] in ['auto', 'xyshift']:

            # Background estimation.
            hdr = hdulist[0].header
            sky_mean, sky_median, sky_std = bckg_data(
                hdr, hdu_data, pars['gain_key'], pars['rdnoise_key'],
                pars['sky_method'])

            # Stars selection.
            psf_select = st_fwhm_select(
                float(pars['dmax']), 1000000, float(pars['thresh_fit']),
                float(pars['fwhm_init']), sky_std, hdu_data)[0]

            # Filter by min/max x,y limits.
            xmi, xma = pars['min_x-max_x'][0]
            ymi, yma = pars['min_y-max_y'][0]
            xmi = float(xmi) if xmi != 'min' else min(psf_select['xcentroid'])
            xma = float(xma) if xma != 'max' else max(psf_select['xcentroid'])
            ymi = float(ymi) if ymi != 'min' else min(psf_select['ycentroid'])
            yma = float(yma) if yma != 'max' else max(psf_select['ycentroid'])
            print(" Selection (x,y) limits: "
                  "({:.0f}, {:.0f}) ; ({:.0f}, {:.0f})".format(
                      xmi, xma, ymi, yma))
            stars_filter = []
            for st in psf_select:
                if xmi <= st['xcentroid'] <= xma and\
                        ymi <= st['ycentroid'] <= yma:
                    stars_filter.append(
                        [st['xcentroid'], st['ycentroid'], st['flux']])
            # Filter by max number of stars.
            stars_filter = stars_filter[:int(pars['max_stars_match'])]

            zip_stars_filter = np.array(zip(*stars_filter))
            xy_cent = list(zip(*zip_stars_filter[:2]))
            id_pick = [str(_) for _ in range(int(pars['max_stars_match']))]

            # Select latest coordinates.
            idx_choice = -1
            id_selec.append(id_pick)
            xy_selec.append(xy_cent)

        else:
            id_pick, xy_pick = [], []
            # Ref stars selection
            landolt_field_img = join(
                pars['mypath'], 'landolt',
                pars['landolt_fld'][proc_grps] + '.gif')
            # Manual selection of reference stars in the observed frame.
            fig, (ax1, ax2) = plt.subplots(1, 2)
            ax1.imshow(plt.imread(landolt_field_img))
            interval = ZScaleInterval()
            zmin, zmax = interval.get_limits(hdu_data)
            ax2.imshow(
                hdu_data, cmap='Greys', aspect=1, interpolation='nearest',
                origin='lower', vmin=zmin, vmax=zmax)

            def onclick(event, ax):
                ax.time_onclick = time.time()

            def onrelease(event, ax):
                # Only clicks inside this axis.
                if event.inaxes == ax:
                    if event.button == 1 and\
                            ((time.time() - ax.time_onclick) < .1):
                        apertures = CircularAperture(
                            (event.xdata, event.ydata), r=10.)
                        apertures.plot(color='green', lw=1.)
                        ax.figure.canvas.draw()
                        ref_id = raw_input(" ID of selected star: ").strip()
                        print(" {} added to list: ({:.2f}, {:.2f})".format(
                            ref_id, event.xdata, event.ydata))
                        id_pick.append(ref_id)
                        xy_pick.append((event.xdata, event.ydata))
                        if len(id_pick) == 3:
                            # Exit when 3 stars have been identified
                            plt.close()
                    elif event.button == 2:
                        print("scroll click")
                    elif event.button == 3:
                        print("right click")
                    else:
                        pass

            # Mouse click / scroll zoom events.
            fig.canvas.mpl_connect(
                'scroll_event', lambda event: zoom(event, [ax1, ax2]))
            # fig.canvas.mpl_connect(
            #     'button_press_event', lambda event: onclick(event, ax2))
            fig.canvas.mpl_connect(
                'button_press_event', lambda event: onclick(event, ax2))
            fig.canvas.mpl_connect(
                'button_release_event', lambda event: onrelease(event, ax2))

            print(" (Select three reference stars covering\n"
                  "  as much frame as possible)")
            plt.show()
            plt.close('all')

            # Re-center coordinates.
            xy_cent = reCenter(hdu_data, xy_pick, side=20).tolist()

            # Select latest coordinates.
            idx_choice = -1
            id_selec.append(id_pick)
            xy_selec.append(xy_cent)

            print(id_pick)
            print(xy_cent)

    print(" Stars selected for match: {}".format(len(xy_selec[idx_choice])))

    return hdu_data, id_selec, xy_selec, idx_choice


def outFileHeader(out_data_file, ref_id):
    """
    Create output stats file with a header.
    """
    if not Path(out_data_file).exists():
        if ref_id != '--':
            with open(out_data_file, mode='w') as f:
                f.write(
                    "# frame       x_obs      y_obs       ref        ID"
                    "          x          y      V     BV      UB     VR"
                    "     RI     VI     e_V    e_BV    e_UB    e_VR    e_RI"
                    "    e_VI\n")
        else:
            with open(out_data_file, mode='w') as f:
                f.write(
                    "# frame          A         B         C         D         "
                    "E         F\n")


def make_out_file(
    mode, f_name, id_all, xy_all, img_id, ref_data, xy_shift, out_data_file):
    """
    Write coordinates of reference stars in the observed frame system.
    """
    if mode == 'manual':
        landolt_in = Table(dtype=ref_data.dtype)
        for r_id in id_all:
            i = ref_data['ID'].tolist().index(r_id)
            landolt_in.add_row(ref_data[i])

        landolt_f = Table({'ref': [img_id for _ in landolt_in['ID']]})
        ref_transf = [[f_name, xy[0], xy[1]] for xy in xy_all]
        obs_tbl = Table(zip(*ref_transf), names=('frame', 'x_obs', 'y_obs'))
        tt = hstack([obs_tbl, landolt_f, landolt_in])

        with open(out_data_file, mode='a') as f:
            # Some platforms don't automatically seek to end when files opened
            # in append mode
            f.seek(0, os.SEEK_END)
            ascii.write(
                tt, f, format='fixed_width_no_header', delimiter='',
                formats={'x_obs': '%9.3f', 'y_obs': '%9.3f', 'ref': '%8s',
                         'ID': '%8s', 'x': '%9.3f', 'y': '%9.3f'})

    elif mode == 'xyshift':
        with open(out_data_file, mode='a') as f:
            tt = Table(zip([f_name, xy_shift[0], xy_shift[1], 1., 0., 0., 1.]),
                       names=('frame', 'A', 'B', 'C', 'D', 'E', 'F'))
            f.seek(0, os.SEEK_END)
            ascii.write(
                tt, f, format='fixed_width_no_header', delimiter='',
                formats={
                    'A': '%8.2f', 'B': '%8.2f', 'C': '%8.2f',
                    'D': '%8.2f', 'E': '%8.2f', 'F': '%8.2f'})


def posFinder(xy_inframe, max_x, min_x, max_y, min_y):
    """
    Finds the offset position to place the IDs of the standard stars in the
    observed frame, such that they don't overlap with other stars or text.
    """
    # Ranges in x,y
    x_rang, y_rang = max_x - min_x, max_y - min_y
    # This value sets the minimum spacing allowed.
    txt_sep = max(x_rang, y_rang) * .025
    # Initiate list with the coordinates of the points themselves.
    used_positions, xy_offset = xy_inframe[:], []
    # For every standard star positioned in the observed frame.
    for i, xy in enumerate(xy_inframe):
        point_done = False
        # Factor that defines the offset separation.
        for j in np.arange(1., 10, .5):
            rang_perc = .025 * j
            # For each possible combination of the offsets.
            for s in [[1., 1.], [1., -1.], [-1., 1.], [-1., -1.]]:
                offset_x = s[0] * rang_perc * x_rang
                offset_y = s[1] * rang_perc * y_rang
                # Offset coordinates.
                xi_off, yi_off = xy[0] + offset_x, xy[1] + offset_y
                # Check that the coordinates are inside the figure.
                if min_x < xi_off < max_x and min_y < yi_off < max_y:
                    # Distance between the position of the text and all
                    # positions already used.
                    d = cdist(
                        np.array([[xi_off, yi_off]]), np.array(used_positions))
                    min_d = np.min(d, axis=1)
                    # Store if it is far away enough.
                    if min_d > txt_sep:
                        used_positions.append([xi_off, yi_off])
                        xy_offset.append([xi_off, yi_off])
                        point_done = True
                        break
            if point_done:
                break
        if not point_done:
            # If no suitable position was found for this point, assign an
            # offset close to its coordinates.
            offset_x, offset_y = .025 * x_rang, .025 * y_rang
            xi_off, yi_off = xy[0] + offset_x, xy[1] + offset_y
            used_positions.append([xi_off, yi_off])
            xy_offset.append([xi_off, yi_off])

    return xy_offset


def make_plot(
    f_name, hdu_data, out_plot_file, std_tr_match, obs_tr_match,
    xy_all, id_all, landolt_field_img):
    """
    Make plots.
    """
    id_inframe, xy_inframe = [], []
    for i, st in enumerate(xy_all):
        if not np.isnan(st[0]) and not np.isnan(st[1]):
            xy_inframe.append(st)
            id_inframe.append(id_all[i])
    xy_inframe = np.array(xy_inframe)

    print("Plotting.")
    fig = plt.figure(figsize=(20, 20))
    gs = gridspec.GridSpec(10, 10)

    ax1 = plt.subplot(gs[0:4, 0:4])
    ax1.set_title("Standard/Reference frame")
    land_img = ax1.imshow(plt.imread(landolt_field_img))
    # Extract maximum y axis value
    max_y = float(len(land_img.get_array()))
    x, y = zip(*std_tr_match)
    y_inv = max_y - np.array(y)
    ax1.scatter(x, y_inv, marker='s', edgecolor='r', facecolor='', lw=1.,
                s=60)

    ax2 = plt.subplot(gs[0:4, 4:8])
    ax2.set_aspect('auto')
    xmin, ymin = np.nanmin(xy_inframe.T, axis=1)
    xmax, ymax = np.nanmax(xy_inframe.T, axis=1)
    xmin_e, xmax_e = max(0., xmin - xmax * .2),\
        min(hdu_data.shape[1], xmax + xmax * .2)
    ymin_e, ymax_e = max(0., ymin - ymax * .2),\
        min(hdu_data.shape[0], ymax + ymax * .2)
    ax2.set_xlim(xmin_e, xmax_e)
    ax2.set_ylim(ymin_e, ymax_e)
    ax2.set_title("Observed frame ({})".format(f_name))
    ax2.grid(lw=1., ls='--', color='grey', zorder=1)
    interval = ZScaleInterval()
    zmin, zmax = interval.get_limits(hdu_data)
    ax2.imshow(hdu_data, cmap='Greys', aspect=1, interpolation='nearest',
               origin='lower', vmin=zmin, vmax=zmax)
    # standard stars in coordinates of the observed frame.
    ax2.scatter(*zip(*xy_inframe), marker='s', edgecolor='g',
                facecolor='', lw=.7, s=60., zorder=4)

    # Define offsets.
    xy_offset = posFinder(xy_inframe.tolist(), xmax, xmin, ymax, ymin)
    for i, txt in enumerate(id_inframe):
        ax2.annotate(
            txt, xy=(xy_inframe[i][0], xy_inframe[i][1]),
            xytext=(xy_offset[i][0], xy_offset[i][1]),
            fontsize=10, arrowprops=dict(
                arrowstyle="-", color='b',
                connectionstyle="angle3,angleA=90,angleB=0"))
    ax2.scatter(*zip(*obs_tr_match), marker='s', edgecolor='r',
                facecolor='', lw=.7, s=60., zorder=6)

    fig.tight_layout()
    plt.savefig(out_plot_file, dpi=150, bbox_inches='tight')
    plt.clf()
    plt.close('all')


def main():
    """
    This algorithm expects at least three stars detected in the observed field.
    """
    pars, fits_list, out_path = in_params()

    # Process each defined group.
    for proc_grps, ref_frame in enumerate(pars['ref_frame']):

        id_selec, xy_selec = [], []

        if pars['landolt_fld'][proc_grps] != '--':
            # Selected standard field.
            ref_data = landolt_fields.main(pars['landolt_fld'][proc_grps])
            # Sort putting brightest stars at the top.
            ref_data.sort('V')
            # Extract (x,y) coords and IDs.
            xy_ref = zip(*[ref_data['x'], ref_data['y']])
            id_ref = ref_data['ID']

            # Path to Landolt image.
            ref_field_img = join(
                pars['mypath'], 'landolt',
                pars['landolt_fld'][proc_grps] + '.gif')
            # ID for final image
            img_id = pars['landolt_fld'][proc_grps]

            print("\nLandolt field: {} ({} stars)".format(
                pars['landolt_fld'][proc_grps], len(id_ref)))
            print("IDs as: {}, {}, {}...".format(*id_ref[:3]))

        else:
            # TODO finish
            print("\nReference frame: {}".format(ref_frame.split('/')[-1]))
            _, id_selec, xy_selec, idx_choice = coo_ref_frame(
                ref_frame, id_selec, xy_selec, pars, proc_grps)
            id_ref, xy_ref = id_selec[idx_choice], xy_selec[idx_choice]
            # ID for final image/file.
            img_id = ref_frame.split('/')[-1].split('.')[0]
            ref_field_img = ref_frame
            # Dummy lists
            ref_data = []

        # Generate output subdir if it doesn't exist.
        out_folder = join(out_path, pars['match_fldr'][proc_grps])
        if not exists(out_folder):
            os.makedirs(out_folder)

        # Write output file header
        out_data_file = join(out_folder, img_id + ".mch")
        outFileHeader(out_data_file, pars['landolt_fld'][proc_grps])

        # Process each fits file in list.
        for fr in fits_list[proc_grps]:
            f_name = fr.split('/')[-1].split('.')[0]
            # Name of final image.
            out_img = join(out_folder, img_id + '_' + f_name + "_obs.png")

            # Skip reference frame.
            if f_name != img_id:

                print("\nMatching frame: {}".format(f_name))

                # Coordinates from observed frame.
                hdu_data, id_selec, xy_selec, idx_choice = coo_ref_frame(
                    fr, id_selec, xy_selec, pars, proc_grps)

                # Selected coordinates.
                idx_choice, xy_choice = id_selec[idx_choice],\
                    xy_selec[idx_choice]

                if pars['match_mode'] == 'xyshift':
                    xmi, xma = map(float, pars['xtr_min-xtr_max'][0])
                    ymi, yma = map(float, pars['ytr_min-ytr_max'][0])
                    max_shift = [[xmi, xma], [ymi, yma]]
                    mtoler = float(pars['match_toler'])
                    xy_shift = xyTrans(max_shift, xy_ref, xy_choice, mtoler)

                    scale, rot_angle = np.nan, np.nan

                elif pars['match_mode'] == 'auto':
                    # Scale and rotation ranges, and match tolerance.
                    scale_range = (
                        float(pars['scale_min']), float(pars['scale_max']))
                    rot_range = (
                        float(pars['rot_min']), float(pars['rot_max']))
                    trans_range = (
                        map(float, pars['xtr_min-xtr_max'][0]),
                        map(float, pars['ytr_min-ytr_max'][0]))
                    mtoler = float(pars['match_toler'])

                    print("\nFinding scale, translation, and rotation.")
                    std_tr_match, obs_tr_match, scale, rot_angle, xy_shift,\
                        xy_transf = triangleMatch(
                            xy_ref, xy_choice, mtoler, scale_range,
                            rot_range, trans_range, hdu_data)

                elif pars['match_mode'] == 'manual':

                    # Match selected reference stars to standard stars by
                    # the IDs given by the user.
                    xy_ref_sel = []
                    for r_id in idx_choice:
                        i = id_ref.tolist().index(r_id)
                        xy_ref_sel.append(xy_ref[i])

                    # Matched stars in ref and obs
                    std_tr_match, obs_tr_match = xy_ref_sel, xy_choice

                    if len(idx_choice) < 3:
                        print("  WARNING: Less than 3 stars selected,\n"
                              "  no match can be performed.")
                        xy_transf = []
                        for st_id in id_ref:
                            if st_id in idx_choice:
                                i = idx_choice.index(st_id)
                                xy_transf.append(
                                    [xy_choice[i][0], xy_choice[i][1]])
                            else:
                                xy_transf.append([np.nan, np.nan])

                        scale, rot_angle, xy_shift = np.nan, np.nan,\
                            [np.nan, np.nan]

                    else:
                        _, A_tr_not_scaled = getTriangles(
                            xy_ref_sel, [(0, 1, 2)])
                        _, B_tr_not_scaled = getTriangles(
                            xy_choice, [(0, 1, 2)])

                        # Obtain scale.
                        scale = np.mean(
                            np.array(B_tr_not_scaled[0]) / A_tr_not_scaled[0])
                        # Rotation angle between triangles.
                        rot_angle, ref_cent, obs_cent = findRotAngle(
                            xy_ref_sel, xy_choice)

                        # Apply translation, scaling, and rotation to reference
                        # coordinates.
                        xy_transf, rot_angle = standard2observed(
                            xy_ref, scale, rot_angle, ref_cent, obs_cent,
                            obs_tr_match)

                        xy_shift = ref_cent - obs_cent

                print("Scale: {:.2f}, Rot: {:.2f}, "
                      "Trans: ({:.2f}, {:.2f})".format(
                          scale, rot_angle, xy_shift[0], xy_shift[1]))

                xy_all, id_all = [], []
                if pars['match_mode'] != 'xyshift':

                    print("Re-center final coordinates.")
                    xy_transf = reCenter(hdu_data, xy_transf, side=40)

                    # Assign nan to reference stars located outside the limits
                    # of the observed frame.
                    for i, st in enumerate(xy_transf):
                        if 0. < st[0] < hdu_data.shape[1] and\
                                0. < st[1] < hdu_data.shape[0]:
                            xy_all.append(st)
                            id_all.append(id_ref[i])
                        else:
                            xy_all.append([np.nan, np.nan])
                            id_all.append(id_ref[i])
                    xy_all = np.array(xy_all)

                    if pars['do_plots_C'] == 'y':
                        make_plot(
                            f_name, hdu_data, out_img, std_tr_match,
                            obs_tr_match, xy_all, id_all,
                            ref_field_img)

            else:
                # TODO reference image values.
                xy_shift = [0., 0.]
                id_all, xy_all = [], []

            # Write final match file
            make_out_file(
                pars['match_mode'], f_name, id_all, xy_all, img_id,
                ref_data, xy_shift, out_data_file)

    print("\nFinished.")


if __name__ == '__main__':
    main()
