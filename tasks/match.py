
import read_pars_file as rpf
from match_funcs import triangleMatch, getTriangles, findRotAngle,\
    standard2observed, xyTrans, reCenter, autoSrcDetect
from plot_events import refDisplay, refStrSlct, drawCirclesIDs

import landolt_fields
import os
from os.path import join, isfile, exists
from pathlib2 import Path
import sys
import numpy as np
from scipy.spatial.distance import cdist

from astropy.table import Table, hstack
from astropy.io import ascii, fits
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
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


def srcSelect(
    f_name, hdu_data, ref_field_img, id_ref, xy_ref, id_selec, xy_selec):
    """
    Manually mark sources on 'fr' fits file.
    """
    answ = 'n'
    if id_selec:

        for i, xy in enumerate(xy_selec):
            # Re-center coordinates.
            xy = reCenter(hdu_data, xy, side=40).tolist()
            # Display frames.
            refDisplay(
                ref_field_img, f_name, hdu_data, id_ref, xy_ref, id_selec[i],
                xy)

            while True:
                answ = raw_input(
                    "Use coordinates ({})? (y/n): ".format(i)).strip()
                if answ == 'y':
                    idx_c = i
                    # Re-write with re-centered coordinates.
                    xy_selec[idx_c] = xy
                    break
                elif answ == 'n':
                    break
                else:
                    print("Wrong answer. Try again.")

    # Select new coordinates
    if answ != 'y':
        while True:
            # Manually pick three reference stars.
            id_pick, xy_pick = refStrSlct(
                ref_field_img, f_name, hdu_data, id_ref, xy_ref)
            if len(id_pick) < 1:
                print("ERROR: at least one star must be selected.")
            elif '' in id_pick:
                print("ERROR: all stars must be identified. Try again.")
            else:
                break

        # Re-center coordinates.
        xy_cent = reCenter(hdu_data, xy_pick, side=20).tolist()

        # Select latest coordinates.
        idx_c = -1
        id_selec.append(id_pick)
        xy_selec.append(xy_cent)

        print(id_pick)
        print(xy_cent)

    print(" Stars selected for match: {}".format(len(xy_selec[idx_c])))

    return id_selec, xy_selec, idx_c


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


def outFileStandard(f_name, ref_data, id_all, xy_all, img_id, out_data_file):
    """
    Write coordinates of reference stars in the observed frame system.
    """
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


def outFileFrame(f_name, scale, rot_angle, xy_shift, out_data_file):
    """
    # TODO transform scale and rotation to C, D E F values.
    """
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
    id_ref, xy_ref, xy_all, id_all, ref_field_img):
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
    interval = ZScaleInterval()

    ax1 = plt.subplot(gs[0:4, 0:4])
    ax1.set_title("Standard/Reference frame")

    if ref_field_img.endswith('.fits'):
        # Load .fits file.
        hdulist = fits.open(ref_field_img)
        hdu_data_ref = hdulist[0].data
        zmin, zmax = interval.get_limits(hdu_data_ref)
        ax1.imshow(
            hdu_data_ref, cmap='Greys', aspect=1, interpolation='nearest',
            origin='lower', vmin=zmin, vmax=zmax)
        drawCirclesIDs(ax1, hdu_data_ref, id_ref, xy_ref)
    else:
        land_img = ax1.imshow(plt.imread(ref_field_img))
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
        # Identify if this is a Landolt or an observed frame.
        ref_id = pars['landolt_fld'][proc_grps]

        id_selec, xy_selec = [], []

        if ref_id != '--':
            # Selected standard field.
            ref_data = landolt_fields.main(ref_id)
            # Sort putting brightest stars at the top.
            ref_data.sort('V')
            # Extract (x,y) coords and IDs.
            xy_ref = zip(*[ref_data['x'], ref_data['y']])
            id_ref = ref_data['ID'].tolist()

            # Path to Landolt image.
            ref_field_img = join(pars['mypath'], 'landolt', ref_id + '.gif')
            # ID for final image
            img_id = ref_id

            print("\nLandolt field: {} ({} stars)".format(
                pars['landolt_fld'][proc_grps], len(id_ref)))
            print("IDs as: {}, {}, {}...".format(*id_ref[:3]))

        else:
            print("\nReference frame: {}".format(ref_frame.split('/')[-1]))
            hdulist = fits.open(ref_frame)
            id_ref, xy_ref, mags_ref = autoSrcDetect(pars, hdulist)

            ref_field_img = ref_frame
            # ID for final image/file.
            img_id = ref_frame.split('/')[-1].split('.')[0]
            # Dummy lists
            ref_data = []

        # Generate output subdir if it doesn't exist.
        out_folder = join(out_path, pars['match_fldr'][proc_grps])
        if not exists(out_folder):
            os.makedirs(out_folder)

        # Write output file header
        out_data_file = join(out_folder, img_id + ".mch")
        outFileHeader(out_data_file, ref_id)

        # Process each fits file in list.
        for fr in fits_list[proc_grps]:
            f_name = fr.split('/')[-1].split('.')[0]
            # Name of final image.
            out_img = join(out_folder, img_id + '_' + f_name + "_obs.png")

            # Skip reference frame.
            if f_name != img_id:

                # Load .fits file.
                hdulist = fits.open(fr)
                hdu_data = hdulist[0].data

                print("\nMatching frame: {}".format(f_name))

                if pars['match_mode'] == 'auto':
                    # Scale and rotation ranges, and match tolerance.
                    scale_range = (
                        float(pars['scale_min']), float(pars['scale_max']))
                    rot_range = (
                        float(pars['rot_min']), float(pars['rot_max']))
                    trans_range = (
                        map(float, pars['xtr_min-xtr_max'][0]),
                        map(float, pars['ytr_min-ytr_max'][0]))
                    mtoler = float(pars['match_toler'])

                    _, xy_choice, _ = autoSrcDetect(pars, hdulist)

                    print("\nFinding scale, translation, and rotation.")
                    std_tr_match, obs_tr_match, scale, rot_angle, xy_shift,\
                        xy_transf = triangleMatch(
                            xy_ref, xy_choice, mtoler, scale_range,
                            rot_range, trans_range, hdu_data)

                elif pars['match_mode'] == 'xyshift':

                    id_dtct, xy_dtct, mags_dtct = autoSrcDetect(
                        pars, hdulist)

                    xmi, xma = map(float, pars['xtr_min-xtr_max'][0])
                    ymi, yma = map(float, pars['ytr_min-ytr_max'][0])
                    max_shift = [[xmi, xma], [ymi, yma]]
                    mtoler = float(pars['match_toler'])
                    xy_shift = xyTrans(
                        max_shift, xy_ref, mags_ref, xy_dtct, mags_dtct,
                        mtoler)

                    scale, rot_angle = np.nan, np.nan

                elif pars['match_mode'] == 'manual':

                    # Coordinates from observed frame.
                    id_selec, xy_selec, idx_c =\
                        srcSelect(
                            f_name, hdu_data, ref_field_img, id_ref, xy_ref,
                            id_selec, xy_selec)

                    # Selected IDs and coordinates.
                    id_choice, xy_choice = id_selec[idx_c], xy_selec[idx_c]

                    # Match selected reference stars to standard stars by
                    # the IDs given by the user.
                    xy_ref_sel = []
                    for r_id in id_choice:
                        i = id_ref.index(r_id)
                        xy_ref_sel.append(xy_ref[i])

                    # Matched stars in ref and obs
                    std_tr_match, obs_tr_match = xy_ref_sel, xy_choice

                    if len(id_choice) < 3:
                        print("  WARNING: Less than 3 stars selected,"
                              " no match can be performed.")
                        xy_transf = []
                        for st_id in id_ref:
                            if st_id in id_choice:
                                i = id_choice.index(st_id)
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
                # TODO plot xyshift outcome
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

                    if id_all:
                        if pars['do_plots_C'] == 'y':
                            make_plot(
                                f_name, hdu_data, out_img, std_tr_match,
                                obs_tr_match, id_ref, xy_ref, xy_all, id_all,
                                ref_field_img)
                    else:
                        print("  ERROR: no stars could be matched.")

            else:
                # TODO reference image values.
                xy_shift, scale, rot_angle = [0., 0.], 1., 0.

            # Write final match file
            if ref_id != '--':
                if id_all:
                    outFileStandard(
                        f_name, ref_data, id_all, xy_all, img_id,
                        out_data_file)
            else:
                outFileFrame(f_name, scale, rot_angle, xy_shift, out_data_file)

    print("\nFinished.")


if __name__ == '__main__':
    main()
