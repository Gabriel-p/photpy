
import read_pars_file as rpf

import os
from os.path import exists, join, isfile
import sys
from operator import itemgetter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
from matplotlib.ticker import NullFormatter

from hlpr import bckg_data, st_fwhm_select, psf_filter

from astropy.io import ascii, fits
from astropy.table import Table
from astropy.visualization import ZScaleInterval

from photutils import CircularAperture
from photutils.utils import cutout_footprint

import psfmeasure


def read_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    folder = pars['in_folder_A']
    folder = folder[1:] if folder.startswith('/') else folder
    in_path = join(pars['mypath'].replace('tasks', 'input'), folder)

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

    # Create path to output folder
    out_path = in_path.replace('input', 'output')

    return pars, fits_list, out_path


def headerCheck(hdr, gain_key, rdnoise_key, filter_key, exposure_key,
                airmass_key):
    """
    Check that all necessary header information is present.
    """
    header_flag = True
    passed_keys = [
        gain_key, rdnoise_key, filter_key, exposure_key, airmass_key]
    i = 0
    try:
        hdr[gain_key]
        i += 1
        hdr[rdnoise_key]
        i += 1
        hdr[filter_key]
        i += 1
        hdr[exposure_key]
        i += 1
        hdr[airmass_key]
    except KeyError:
        _key = ['Gain', 'Rdnoise', 'Filter', 'Exposure', 'Airmass']
        print("{} key '{}' is not present in header.".format(
            _key[i], passed_keys[i]))
        header_flag = False

    return header_flag


def isolate_stars(hdu_data, fwhm_accptd, crop_side):
    """
    Create crop centered at each star used to obtain the FWHM.
    """
    stars = []
    for st in fwhm_accptd:
        # Crop image
        crop = cutout_footprint(hdu_data, (st[0], st[1]), crop_side)
        # Skip stars placed on borders.
        if crop[0].shape == (crop_side, crop_side):
            stars.append(crop[0])
        else:
            print("Skip star on border (xc, yc)=({}, {})".format(st[0], st[1]))

    if stars:
        # Obtain average combined PSF.
        # Meshgrid of coordinates (0,1,...,N) times (0,1,...,N)
        stars = np.array(stars, dtype='float32')
        y, x = np.mgrid[:len(stars[0, :, 0]), :len(stars[0, 0, :])]
        # duplicating the grids
        xcoord, ycoord = np.array([x] * len(stars)), np.array([y] * len(stars))
        # compute histogram with coordinates as x,y
        psf_avrg = np.histogram2d(
            xcoord.ravel(), ycoord.ravel(),
            bins=[len(stars[0, 0, :]), len(stars[0, :, 0])],
            weights=stars.ravel())[0]
    else:
        print("  WARNING: no suitable stars found to obtain average PSF.")
        psf_avrg = []

    return stars, psf_avrg


def make_plots(
    hdu_data, crop_side, max_stars, all_sources, n_not_satur,
        fwhm_min_rjct, ellip_rjct, fwhm_accptd, fwhm_outl, fwhm_mean,
        fwhm_std, stars, psf_avrg, out_path, imname, filter_val):
    """
    Make plots.
    """
    print("\nPlotting.")
    fig = plt.figure(figsize=(20, 15))
    gs = gridspec.GridSpec(10, 12)

    ax = plt.subplot(gs[0:5, 0:5])
    interval = ZScaleInterval()
    zmin, zmax = interval.get_limits(hdu_data)
    plt.imshow(hdu_data, cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower', vmin=zmin, vmax=zmax)
    if fwhm_accptd:
        positions = (zip(*fwhm_accptd)[0], zip(*fwhm_accptd)[1])
        apertures = CircularAperture(positions, r=2. * fwhm_mean)
        apertures.plot(color='red', lw=1.)
    ax.set_title('{} sources detected ({} not saturated)'.format(
        all_sources, n_not_satur), fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=8)
    # plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)

    gs0 = gridspec.GridSpec(10, 12)
    gs0.update(bottom=0.25, left=0.15, right=0.88)
    ax = plt.subplot(gs0[0:2, 5:10])
    ax.set_xlim(-0.2, 2.2)
    x = np.arange(3)
    y = [len(fwhm_min_rjct), len(ellip_rjct), len(fwhm_outl)]
    y_perc = np.array(y) / float(max_stars)
    up = max(y_perc) * .03
    ax.set_ylim(0, max(y_perc) + 5. * up)
    ax.bar(x, y_perc, align='center', width=0.2, color=('g', 'b', 'k'),
           zorder=4)
    for xi, yi, l in zip(*[x, y_perc, list(map(str, y))]):
        ax.text(xi - len(l) * .02, yi + up, l, fontsize=9,
                bbox=dict(facecolor='w', edgecolor='w', alpha=.5))
    ax.set_xticks(x)
    ax.set_xticklabels(['FWHM min rjct', 'Ellip max rjct', 'Outliers rjct'],
                       fontsize=9)
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.set_title("As % of max FWHM stars ({})".format(max_stars), fontsize=9)

    gs3 = gridspec.GridSpec(10, 12)
    gs3.update(hspace=0.35, wspace=0.1)
    if fwhm_accptd:
        ax = plt.subplot(gs3[2:5, 8:10])
        ax.hist(
            zip(*fwhm_accptd)[2], histtype='step', lw=2.,
            range=[fwhm_mean - 2. * fwhm_std, fwhm_mean + 2. * fwhm_std],
            bins=max(1, int(len(fwhm_accptd) * 0.3)))
        ax.axvline(fwhm_mean, ls='-', lw=3.5, c='g')
        ax.axvline(fwhm_mean + fwhm_std, ls='--', lw=1.5, c='r')
        ax.axvline(fwhm_mean - fwhm_std, ls='--', lw=1.5, c='r')
        ax.set_title('FWHM distribution ({} stars)'.format(
            len(fwhm_accptd)), fontsize=9)
        ax.tick_params(axis='x', which='major', labelsize=8)
        ax.tick_params(axis='y', which='major', labelleft='off')
        ax.set_xlim(fwhm_mean - 3. * fwhm_std, fwhm_mean + 3. * fwhm_std)
    elif fwhm_min_rjct:
        ax = plt.subplot(gs3[2:5, 8:10])
        ax.hist(zip(*fwhm_min_rjct)[2], histtype='step', lw=2.)
        ax.set_title('FWHM < min distribution ({} stars)'.format(
            len(fwhm_min_rjct)), fontsize=11)
        ax.tick_params(axis='x', which='major', labelsize=8)
        ax.tick_params(axis='y', which='major', labelleft='off')
    elif ellip_rjct:
        ax = plt.subplot(gs3[2:5, 8:10])
        ax.hist(zip(*ellip_rjct)[2], histtype='step', lw=2.)
        ax.set_title('e < max distribution ({} stars)'.format(
            len(ellip_rjct)), fontsize=11)
        ax.tick_params(axis='x', which='major', labelsize=8)
        ax.tick_params(axis='y', which='major', labelleft='off')

    ax = plt.subplot(gs3[2:5, 5:8], aspect=1)
    ymax, xmax = np.shape(hdu_data)
    ax.set_xlim(0., xmax)
    ax.set_ylim(0., ymax)
    ax.grid(b=True, which='major', color='grey', linestyle=':', lw=.5,
            zorder=1)
    ax.set_title('Spatial distribution of FWHM stars',
                 fontsize=9)
    ax.tick_params(axis='both', which='major', labelsize=8)
    f = min(xmax, ymax) * .05
    if fwhm_accptd:
        ax.scatter(*positions, s=0.5 * np.exp(zip(*fwhm_accptd)[2]), lw=0.7,
                   facecolors='none', edgecolors='r', zorder=4)
    if fwhm_outl:
        positions = (zip(*fwhm_outl)[0], zip(*fwhm_outl)[1])
        sz = 0.01 * np.exp(zip(*fwhm_outl)[2])
        sz[sz > max(ymax, xmax)] = max(ymax, xmax)
        ax.scatter(*positions, s=sz, lw=0.7,
                   facecolors='none', edgecolors='k', zorder=4)
    if fwhm_min_rjct:
        fwhm_min_pos = (zip(*fwhm_min_rjct)[0], zip(*fwhm_min_rjct)[1])
        r = 10. * np.array(zip(*fwhm_min_rjct)[2])
        r[r > 2. * f] = 2. * f
        ax.scatter(*fwhm_min_pos, s=r, lw=0.7, facecolors='none',
                   edgecolors='g', zorder=4)
    if ellip_rjct:
        fwhm_log = np.log(zip(*ellip_rjct)[2])
        fwhm_log[fwhm_log < 0.] = 1.
        r = f * fwhm_log
        r[r > 2. * f] = 2. * f
        for i, st in enumerate(ellip_rjct):
            e = Ellipse(
                xy=(st[0], st[1]), width=r[i],
                height=r[i] * np.sqrt(1 - st[3] ** 2),
                angle=0., lw=0.7, fc='None', edgecolor='b', zorder=4)
            ax.add_artist(e)

    if getattr(stars, 'size', len(stars)):
        stars = np.array(stars, dtype='float32')
        gs1 = gridspec.GridSpec(10, 12)
        gs1.update(wspace=0., hspace=0.25, top=0.85)
        max_s = min(25, len(stars))
        for i in range(max_s):
            if i < 5:
                pl_n = 60 + i
            elif 5 <= i < 10:
                pl_n = (60 + 12 - 5) + i
            elif 10 <= i < 15:
                pl_n = (60 + 2 * 12 - 10) + i
            elif 15 <= i < 20:
                pl_n = (60 + 3 * 12 - 15) + i
            elif 20 <= i:
                pl_n = (60 + 4 * 12 - 20) + i
            ax = fig.add_subplot(gs1[pl_n])
            ax.imshow(
                stars[i], cmap=plt.cm.viridis, interpolation='nearest',
                origin='lower', vmin=0.)
            ax.set_title("({:.0f}, {:.0f})".format(
                fwhm_accptd[i][0], fwhm_accptd[i][1]), fontsize=8, y=-0.23)
            ax.set_axis_off()

    if getattr(psf_avrg, 'size', len(psf_avrg)):
        gs2 = gridspec.GridSpec(10, 12)
        gs2.update(hspace=0.05, wspace=0.05, left=0.22, right=0.85, top=0.82)
        ax0 = plt.subplot(gs2[6:10, 5:9])
        ax0.imshow(psf_avrg, cmap=plt.cm.hot, interpolation='nearest',
                   origin='lower', vmin=0., aspect='auto')
        # x,y center of the average PSF.
        combined_star_argmax = np.unravel_index(
            np.argmax(psf_avrg), psf_avrg.shape)
        ax0.axhline(combined_star_argmax[0], ls='--', lw=2, c='w')
        ax0.axvline(combined_star_argmax[1], ls='--', lw=2, c='w')
        cent = int(crop_side / 2.) - 1
        apertures = CircularAperture((cent, cent), r=.5 * fwhm_mean)
        apertures.plot(color='b', lw=0.75)

        # x,y histograms of average PSF.
        hx, hy = psf_avrg.sum(axis=0), psf_avrg.sum(axis=1)
        # Top histogram
        axx = plt.subplot(gs2[5:6, 5:9])
        axx.plot(hx)
        axx.set_xlim(ax0.get_xlim())
        axx.set_title("Average PSF model (FWHM={:.2f})".format(fwhm_mean),
                      fontsize=9)
        # Right histogram
        axy = plt.subplot(gs2[6:10, 9:10])
        axy.plot(hy, range(len(hy)))
        axy.set_ylim(ax0.get_ylim())

        # Remove tick labels
        nullfmt = NullFormatter()
        ax0.xaxis.set_major_formatter(nullfmt)
        ax0.yaxis.set_major_formatter(nullfmt)
        axx.xaxis.set_major_formatter(nullfmt)
        axx.yaxis.set_major_formatter(nullfmt)
        axy.xaxis.set_major_formatter(nullfmt)
        axy.yaxis.set_major_formatter(nullfmt)

    if os.path.isfile(out_path):
        fig_name = join(out_path.replace(".fits", ""))
    else:
        fig_name = join(
            out_path, 'filt_' + filter_val,
            imname.split('/')[-1].replace(".fits", ""))
    plt.savefig(fig_name + '.png', dpi=150, bbox_inches='tight')

    # Close to release memory.
    plt.clf()
    plt.close()


def storeCooFile(out_path, imname, filter_val, fwhm_accptd):
    """
    Write .coo file.
    """
    # Save coordinates data to file.
    fn = join(out_path, 'filt_' + filter_val,
              imname.split('/')[-1].replace('.fits', '.coo'))
    ascii.write(
        zip(*fwhm_accptd), fn,
        names=['x', 'y', 'FWHM', 'Ellip', 'Mag'],
        format='commented_header', overwrite=True, delimiter=' ')


def storeOutFile(out_path, out_data):
    """
    Write output file.
    """
    # Sort list by filter names and then by exposure time.
    data_sort = sorted(out_data, key=itemgetter(1, 3))
    out_data_file = join(out_path, "fitstats.dat")
    data = Table(zip(*data_sort), names=(
        '# image         ', 'filter', 'airmass', 'exposure', 'Sky_mean',
        'Sky_STDDEV', 'FWHM_(N_stars)', 'FWHM_(mean)', 'FWHM_(std)'))
    ascii.write(
        data, out_data_file, overwrite=True, formats={
            '# image         ': '%-16s', 'FWHM_(N_stars)': '%14.0f',
            'Sky_mean': '%10.2f', 'Sky_STDDEV': '%10.2f',
            'FWHM_(mean)': '%10.2f', 'FWHM_(std)': '%10.2f'},
        format='fixed_width', delimiter=None)


def main():
    """
    Get FWHM, sky mean, and sky standard deviation from a single .fits file or
    a folder containing .fits files.
    """
    pars, fits_list, out_path = read_params()

    # Generate output dir/subdir if it doesn't exist.
    if not exists(out_path):
        os.makedirs(out_path)

    # HARDCODED: length used to isolate stars used to obtain the average FWHM.
    crop_side = 20

    out_data = []
    # For each .fits image in the root folder.
    for imname in fits_list:
        print("\n* File: {}".format(
            imname.replace(pars['mypath'].replace('tasks', 'input'), "")))

        # Load .fits file.
        hdulist = fits.open(imname)
        # Extract header.
        hdr = hdulist[0].header

        # Check header keys.
        header_flag = headerCheck(
            hdr, pars['gain_key'], pars['rdnoise_key'], pars['filter_key'],
            pars['exposure_key'], pars['airmass_key'])
        if not header_flag:
            print("  ERROR: Missing header information. Skipping this file.")
            break

        # Extract data.
        hdu_data = hdulist[0].data

        # Background estimation.
        sky_mean, sky_median, sky_std = bckg_data(
            hdr, hdu_data, pars['gain_key'], pars['rdnoise_key'],
            pars['sky_method'])

        # Stars selection.
        psf_select, all_sources, n_not_satur = st_fwhm_select(
            float(pars['dmax']), int(pars['max_stars']),
            float(pars['thresh_fit']), float(pars['fwhm_init']),
            sky_std, hdu_data)

        psf_data = []
        if psf_select:
            # IRAF 'psfmeasure' task.
            psf_data = psfmeasure.main(
                float(pars['dmax']), psf_select, imname, hdu_data)

        if psf_data:
            # Filter out stars.
            fwhm_min_rjct, ellip_rjct, fwhm_accptd, fwhm_mean, fwhm_std,\
                fwhm_outl = psf_filter(
                    float(pars['fwhm_min']), float(pars['ellip_max']),
                    psf_data)

            # Create output folder if it does not exist.
            filt_folder = join(out_path, 'filt_' + hdr[pars['filter_key']])
            if not exists(filt_folder):
                os.makedirs(filt_folder)

            stars, psf_avrg = [], []
            if fwhm_accptd:
                # Create .coo file.
                storeCooFile(out_path, imname, hdr[pars['filter_key']],
                             fwhm_accptd)
                # Isolated stars, and average PSF.
                stars, psf_avrg = isolate_stars(
                    hdu_data, fwhm_accptd, crop_side)

            if pars['do_plots_A'] is 'y':
                make_plots(
                    hdu_data, crop_side, int(pars['max_stars']), all_sources,
                    n_not_satur, fwhm_min_rjct, ellip_rjct, fwhm_accptd,
                    fwhm_outl, fwhm_mean, fwhm_std, stars, psf_avrg, out_path,
                    imname, hdr[pars['filter_key']])

        else:
            fwhm_accptd, fwhm_mean, fwhm_std = [], np.nan, np.nan
            print("  ERROR: no stars returned from 'psfmeasure' task.")

        # Store data to write to output file.
        im_name = imname.split('/')[-1]
        out_data.append(
            [im_name, hdr[pars['filter_key']], hdr[pars['airmass_key']],
             hdr[pars['exposure_key']], sky_mean, sky_std,
             len(fwhm_accptd), fwhm_mean, fwhm_std])

    if out_data:
        storeOutFile(out_path, out_data)

    print("\nFinished.")


if __name__ == "__main__":
    main()
