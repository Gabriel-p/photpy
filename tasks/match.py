
import read_pars_file as rpf
from match_funcs import autoSrcDetect, allInFrame
from plot_events import drawCirclesIDs

import match_auto
import match_manual
import match_xyshift

import landolt_fields
import os
from os.path import join, isfile, exists
from pathlib2 import Path
import numpy as np
from scipy.spatial.distance import cdist

from astropy.table import Table, hstack
from astropy.io import ascii, fits
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from photutils import CircularAperture


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
            if pars['match_mode'] == 'xyshift':
                raise ValueError(
                    "\n'xyshift' mode can not be used with Landolt fields.")

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
                raise ValueError(
                    "{}\n Reference frame is not present. Exit.".format(
                        ref_frame))
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
            raise ValueError(
                "{}\n No .fits files found in match folder. Exit.".format(
                    in_path))

        # Store list for this folder.
        print("Files found in '{}' folder:".format(folder))
        for fit in list_temp:
            print(" * {}".format(fit.replace(in_path, '')[1:]))
        fits_list.append(list_temp)

    # Put reference frames to the front  of the lists for each group.
    if pars['match_mode'] == 'xyshift':
        for i, fit_group in enumerate(fits_list):
            ref_f_gr = pars['ref_frame'][i]
            # Index of the reference frame for this group.
            j = fit_group.index(ref_f_gr)
            # Bring ref frame to the front.
            fit_group.insert(0, fit_group.pop(j))
            # Update this group.
            fits_list[i] = fit_group

    return pars, fits_list, out_path


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
                'frame': '%-15s',
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


def auto_manualPlot(
    ref_field_img, id_ref, xy_ref, f_name, out_plot_file, hdu_data,
    std_tr_match, obs_tr_match, xy_all, id_all):
    """
    Make 'auto' and 'manual' mode plots.
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


def xyshiftPlot(
    rf_name, hdu_data_ref, out_plot_file, f_name, hdu_data, xy_shift, xy_dtct):
    """
    """
    print("Plotting.")

    # Plot 100 brightest stars (they are already sorted by brightest)
    x_dtct, y_dtct = np.array(xy_dtct[:100])[:, 0],\
        np.array(xy_dtct[:100])[:, 1]
    x_shft, y_shft = x_dtct + xy_shift[0], y_dtct + xy_shift[1]

    fig = plt.figure(figsize=(20, 20))
    gs = gridspec.GridSpec(10, 10)
    interval = ZScaleInterval()

    ax1 = plt.subplot(gs[0:5, 0:5])
    ax1.set_title("Reference frame {}".format(rf_name))
    zmin, zmax = interval.get_limits(hdu_data_ref)
    ax1.imshow(
        hdu_data_ref, cmap='Greys', aspect=1, interpolation='nearest',
        origin='lower', vmin=zmin, vmax=zmax)
    apertures = CircularAperture(np.array([x_shft, y_shft]).T, r=10.)
    apertures.plot(ax=ax1, color='red', lw=1.)

    ax2 = plt.subplot(gs[0:5, 5:10])
    ax2.set_title("Matched frame {}, shift: ({:.2f}, {:.2f})".format(
        f_name, *xy_shift))
    zmin, zmax = interval.get_limits(hdu_data)
    ax2.imshow(hdu_data, cmap='Greys', aspect=1, interpolation='nearest',
               origin='lower', vmin=zmin, vmax=zmax)
    apertures = CircularAperture(np.array([x_dtct, y_dtct]).T, r=10.)
    apertures.plot(ax=ax2, color='green', lw=1.)

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
            rf_name = ref_frame.split('/')[-1].split('.')[0]
            print("\nReference frame: {}".format(rf_name))
            hdulist = fits.open(ref_frame)
            hdu_data_ref = hdulist[0].data
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
                    std_tr_match, obs_tr_match, scale, rot_angle, xy_shift,\
                        xy_transf = match_auto.main(pars, xy_ref, hdulist)

                    id_all, xy_all = allInFrame(id_ref, hdu_data, xy_transf)

                elif pars['match_mode'] == 'manual':
                    std_tr_match, obs_tr_match, scale, rot_angle, xy_shift,\
                        xy_transf = match_manual.main(
                            f_name, hdu_data, ref_field_img, id_ref, xy_ref,
                            id_selec, xy_selec)

                    id_all, xy_all = allInFrame(id_ref, hdu_data, xy_transf)

                elif pars['match_mode'] == 'xyshift':
                    scale, rot_angle, xy_shift, xy_dtct =\
                        match_xyshift.main(pars, hdulist, xy_ref, mags_ref)

                print("Scale: {:.2f}, Rot: {:.2f}, "
                      "Trans: ({:.2f}, {:.2f})".format(
                          scale, rot_angle, xy_shift[0], xy_shift[1]))

                if pars['do_plots_C'] == 'y':
                    if pars['match_mode'] in ['auto', 'manual']:
                        if id_all:
                            auto_manualPlot(
                                ref_field_img, id_ref, xy_ref, f_name, out_img,
                                hdu_data, std_tr_match, obs_tr_match, xy_all,
                                id_all)
                        else:
                            print("  ERROR: no stars could be matched.")
                    else:
                        xyshiftPlot(
                            rf_name, hdu_data_ref, out_img, f_name,
                            hdu_data, xy_shift, xy_dtct)

            else:
                # Reference image values.
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
