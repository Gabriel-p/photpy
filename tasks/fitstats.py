
import os
from os.path import exists, join, realpath, dirname
import sys
import gc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
from matplotlib.ticker import NullFormatter

from astropy.io import ascii
from astropy.table import Table
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval

from photutils import DAOStarFinder
from photutils import CircularAperture
from photutils.utils import cutout_footprint

import psfmeasure


def create_pars_file(pars_f, pars_list=None):
    """
    Default values for parameters file.
    """
    if pars_list is None:
        pars_list = ['None', 'y', '60000.', '5.', '3.', '100', '0.15',
                     '1.5', 'EGAIN', 'ENOISE', 'FILTER', 'EXPTIME']
    with open(pars_f, 'w') as f:
        f.write(
            "# Default parameters for the 'fitstats' script\n#\n"
            "ff_proc {}\ndo_plots {}\ndmax {}\nthresh_level {}\nfwhm_init {}"
            "\nmax_stars {}\nellip_max {}\nfwhm_min {}\ngain_key {}"
            "\nrdnoise_key {}\nfilter_key {}\nexp_key {}\n".format(
                *pars_list))
    return


def read_params():
    """
    Read parameter values from .pars file.
    """
    pars = {}
    mypath = realpath(join(os.getcwd(), dirname(__file__)))
    pars_f = join(mypath, 'fitstats.pars')
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
    r_path = join(mypath.replace('tasks', 'input'), fname)
    fits_list = []
    if os.path.isdir(r_path):
        for subdir, dirs, files in os.walk(r_path):
            for file in files:
                if file.endswith('.fits'):
                    fits_list.append(os.path.join(subdir, file))
    elif os.path.isfile(r_path):
        print("  Single .fits file.")
        fits_list.append(r_path)
    else:
        print("{}\nis neither a folder nor a file. Exit.".format(r_path))
        sys.exit()
    pars['ff_proc'] = fname
    pars_list.append(pars['ff_proc'])

    answ = raw_input("Create plots? (y/n) ({}): ".format(pars['do_plots']))
    pars['do_plots'] = str(answ) if str(answ) is not '' else\
        str(pars['do_plots'])
    pars['do_plots'] = 'n' if pars['do_plots'] in ('n', 'N') else 'y'
    pars_list.append(pars['do_plots'])

    answ = raw_input("Maximum flux ({}): ".format(pars['dmax']))
    pars['dmax'] = float(answ) if answ is not '' else float(pars['dmax'])
    pars_list.append(pars['dmax'])

    answ = raw_input("Threshold level above STDDEV ({}): ".format(
        pars['thresh_level']))
    pars['thresh_level'] = float(answ) if answ is not '' else\
        float(pars['thresh_level'])
    pars_list.append(pars['thresh_level'])

    answ = raw_input("Initial FWHM ({}): ".format(pars['fwhm_init']))
    pars['fwhm_init'] = float(answ) if answ is not '' else\
        float(pars['fwhm_init'])
    pars_list.append(pars['fwhm_init'])

    answ = raw_input("Max number of stars used to obtain "
                     "FWHM ({}): ".format(pars['max_stars']))
    pars['max_stars'] = int(answ) if answ is not '' else\
        int(pars['max_stars'])
    pars_list.append(pars['max_stars'])

    answ = raw_input("Maximum ellipticity ({}): ".format(pars['ellip_max']))
    pars['ellip_max'] = float(answ) if answ is not '' else\
        float(pars['ellip_max'])
    pars_list.append(pars['ellip_max'])

    answ = raw_input("Minimum FWHM ({}): ".format(pars['fwhm_min']))
    pars['fwhm_min'] = float(answ) if answ is not '' else\
        float(pars['fwhm_min'])
    pars_list.append(pars['fwhm_min'])

    answ = raw_input("GAIN keyword ({}): ".format(pars['gain_key']))
    pars['gain_key'] = str(answ) if answ is not '' else pars['gain_key']
    pars_list.append(pars['gain_key'])

    answ = raw_input("RDNOISE keyword ({}): ".format(pars['rdnoise_key']))
    pars['rdnoise_key'] = str(answ) if answ is not '' else pars['rdnoise_key']
    pars_list.append(pars['rdnoise_key'])

    answ = raw_input("FILTER keyword ({}): ".format(pars['filter_key']))
    pars['filter_key'] = str(answ) if answ is not '' else pars['filter_key']
    pars_list.append(pars['filter_key'])

    answ = raw_input("EXPOSURE keyword ({}): ".format(pars['exp_key']))
    pars['exp_key'] = str(answ) if answ is not '' else pars['exp_key']
    pars_list.append(pars['exp_key'])

    # Write values to file.
    create_pars_file(pars_f, pars_list)

    # Create path to output folder
    out_path = r_path.replace('input', 'output')
    if os.path.isfile(r_path):
        out_path = out_path.replace(out_path.split('/')[-1], '')

    return out_path, fits_list, pars


def headerCheck(hdr, gain_key, rdnoise_key, filter_key, exp_key):
    """
    Check that all necessary header information is present.
    """
    header_flag = True
    passed_keys = [gain_key, rdnoise_key, filter_key, exp_key]
    i = 0
    try:
        hdr[gain_key]
        i += 1
        hdr[rdnoise_key]
        i += 1
        hdr[filter_key]
        i += 1
        hdr[exp_key]
    except KeyError:
        _key = ['Gain', 'Rdnoise', 'Filter', 'Exposure']
        print("{} key '{}' is not present in header.".format(
            _key[i], passed_keys[i]))
        header_flag = False

    return header_flag


def bckg_data(hdr, hdu_data, gain_key, rdnoise_key):
    """
    Estimate sky's mean, median, and standard deviation.
    """
    gain = hdr[gain_key]
    rdnoise = hdr[rdnoise_key]
    print("Gain, Rdnoise: {}, {}".format(gain, rdnoise))

    print("\nEstimating background mean, median, and STDDEV.")
    sky_mean, sky_median, sky_std = sigma_clipped_stats(
        hdu_data, sigma=3.0, iters=5)
    print("Mean, median, and STDDEV: {:.2f}, {:.2f}, {:.2f}".format(
        sky_mean, sky_median, sky_std))
    print("STDDEV estimated from GAIN, RDNOISE, and median: {:.2f}".format(
        np.sqrt(sky_median * gain + rdnoise ** 2) / gain))

    return sky_mean, sky_std


def st_fwhm_select(
        dmax, max_stars, thresh_level, fwhm_init, std, hdu_data):
    """
    Find stars in image with DAOStarFinder, filter out saturated stars (with
    peak values greater than the maximum accepted in 'dmax'), order by the
    largest (most negative) magnitude, and select a maximum of 'max_stars'.
    """
    print("\nFinding stars.")
    thresh = thresh_level * std
    print('Threshold, initial FWHM: {:.1f}, {:.1f}'.format(thresh, fwhm_init))
    stfind = DAOStarFinder(threshold=thresh, fwhm=fwhm_init)
    sources = stfind(hdu_data)
    print("Sources found: {}".format(len(sources)))
    mask = sources['peak'] < dmax
    sour_no_satur = sources[mask]
    print("Non-saturated sources found: {}".format(len(sour_no_satur)))
    sour_no_satur.sort('mag')
    if max_stars == 'max':
        max_stars = len(sour_no_satur)
    psf_select = sour_no_satur[:int(max_stars)]

    return psf_select, len(sources), len(sour_no_satur)


def rm_outliers(fwhm_estim, out_f=2.):
    """
    Remove stars with large FWHMs.
    """
    fwhm_median = np.median(zip(*fwhm_estim)[2])
    fwhm_no_outl, fwhm_outl = [], []
    for st in fwhm_estim:
        if st[2] <= out_f * fwhm_median:
            fwhm_no_outl.append(st)
        else:
            fwhm_outl.append(st)
    print("Outliers with FWHM>{:.2f} rejected: {}".format(
        out_f * fwhm_median, len(fwhm_outl)))

    if fwhm_no_outl:
        print("\nFinal number of accepted stars: {}".format(len(fwhm_no_outl)))
        fwhm_mean, fwhm_std = np.mean(zip(*fwhm_no_outl)[2]),\
            np.std(zip(*fwhm_no_outl)[2])
        print("Mean FWHM +/- std: {:.2f}, {:.2f}".format(fwhm_mean, fwhm_std))
    else:
        print("  WARNING: all selected stars rejected as outliers.\n"
              "  Could not obtain a mean FWHM.")
        fwhm_mean, fwhm_std = 0., 0.

    return fwhm_no_outl, fwhm_mean, fwhm_std, fwhm_outl


def isolate_stars(hdu_data, fwhm_no_outl, crop_side):
    """
    Create crop centered at each star used to obtain the FWHM.
    """
    stars = []
    for st in fwhm_no_outl:
        # Crop image
        crop = cutout_footprint(hdu_data, (st[0], st[1]), crop_side)
        # Skip stars placed on borders.
        if len(crop[0]) == crop_side:
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
    hdu_data, max_stars, crop_side, all_sources, n_not_satur,
        fwhm_min_rjct, ellip_rjct, fwhm_no_outl, fwhm_outl, fwhm_mean,
        fwhm_std, stars, psf_avrg, out_path, imname):
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
    if fwhm_no_outl:
        positions = (zip(*fwhm_no_outl)[0], zip(*fwhm_no_outl)[1])
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
    y = np.array([len(fwhm_min_rjct) / float(max_stars),
                  len(ellip_rjct) / float(max_stars),
                  len(fwhm_outl) / float(max_stars)])
    ax.bar(x, y, align='center', width=0.2, color='g')
    ax.set_xticks(x)
    ax.set_xticklabels(['FWHM min rjct', 'Ellip max rjct', 'Outliers rjct'],
                       fontsize=9)
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.set_title("As % of max FWHM stars ({})".format(max_stars),
                 fontsize=9)

    gs3 = gridspec.GridSpec(10, 12)
    gs3.update(hspace=0.35, wspace=0.1)
    if fwhm_no_outl:
        ax = plt.subplot(gs3[2:5, 8:10])
        ax.hist(
            zip(*fwhm_no_outl)[2], histtype='step', lw=2.,
            range=[fwhm_mean - 2. * fwhm_std, fwhm_mean + 2. * fwhm_std],
            bins=max(1, int(len(fwhm_no_outl) * 0.3)))
        ax.axvline(fwhm_mean, ls='-', lw=3.5, c='g')
        ax.axvline(fwhm_mean + fwhm_std, ls='--', lw=1.5, c='r')
        ax.axvline(fwhm_mean - fwhm_std, ls='--', lw=1.5, c='r')
        ax.set_title('FWHM distribution ({} stars)'.format(
            len(fwhm_no_outl)), fontsize=9)
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
    f = min(xmax, ymax) * .05
    if fwhm_no_outl:
        ax.scatter(*positions, s=0.5 * np.exp(zip(*fwhm_no_outl)[2]), lw=0.7,
                   facecolors='none', edgecolors='r')
    if fwhm_min_rjct:
        fwhm_min_pos = (zip(*fwhm_min_rjct)[0], zip(*fwhm_min_rjct)[1])
        r = 10. * np.array(zip(*fwhm_min_rjct)[2])
        r[r > 2. * f] = 2. * f
        ax.scatter(*fwhm_min_pos, s=r, lw=0.7, facecolors='none',
                   edgecolors='g')
    if ellip_rjct:
        fwhm_log = np.log(zip(*ellip_rjct)[2])
        fwhm_log[fwhm_log < 0.] = 1.
        r = f * fwhm_log
        r[r > 2. * f] = 2. * f
        for i, st in enumerate(ellip_rjct):
            e = Ellipse(
                xy=(st[0], st[1]), width=r[i],
                height=r[i] * np.sqrt(1 - st[3] ** 2),
                angle=0., lw=0.7, fc='None', edgecolor='b')
            ax.add_artist(e)
    ax.set_title('Spatial distribution of (exponential) FWHM sizes',
                 fontsize=9)
    ax.tick_params(axis='both', which='major', labelsize=8)

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
                fwhm_no_outl[i][0], fwhm_no_outl[i][1]), fontsize=8, y=-0.23)
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
        apertures = CircularAperture((cent, cent), r=fwhm_mean)
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
            out_path, imname.split('/')[-1].replace(".fits", ""))
    plt.savefig(fig_name + '.png', dpi=150, bbox_inches='tight')

    # Close to release memory.
    plt.clf()
    plt.close()
    # Force the Garbage Collector to release unreferenced memory.
    gc.collect()


def storeCooFile(out_path, imname, fwhm_no_outl):
    """
    Write .coo file.
    """
    # Save coordinates data to file.
    fn = join(out_path, imname.split('/')[-1].replace('.fits', '.coo'))
    ascii.write(
        zip(*fwhm_no_outl), fn,
        names=['x', 'y', 'FWHM', 'Ellip', 'Mag'],
        format='commented_header', overwrite=True, delimiter=' ')


def storeOutFile(out_path, out_data):
    """
    Write output file.
    """
    out_data_file = join(out_path, "getdata.dat")
    data = Table(zip(*out_data), names=(
        '# image         ', 'filter', 'exposure', 'Sky_mean', 'Sky_STDDEV',
        'FWHM_(N_stars)', 'FWHM_(mean)', 'FWHM_(std)'))
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
    mypath, pars_f, pars = read_params()
    out_path, fits_list, pars = get_params(mypath, pars_f, pars)

    # Generate output dir/subdir if it doesn't exist.
    if not exists(out_path):
        os.makedirs(out_path)

    # HARDCODED: length used to isolate stars used to obtain the average FWHM.
    crop_side = 20

    out_data = []
    # For each .fits image in the root folder.
    for imname in fits_list:
        print("\nFile: {}".format(
            imname.replace(mypath.replace('tasks', 'input'), "")))

        # Load .fits file.
        hdulist = fits.open(imname)
        # Extract header.
        hdr = hdulist[0].header

        # Check header keys.
        header_flag = headerCheck(
            hdr, pars['gain_key'], pars['rdnoise_key'], pars['filter_key'],
            pars['exp_key'])
        if not header_flag:
            print("  ERROR: Missing header information. Skipping this file.")
            break

        # Extract data.
        hdu_data = hdulist[0].data

        # Background estimation.
        sky_mean, sky_std = bckg_data(
            hdr, hdu_data, pars['gain_key'], pars['rdnoise_key'])

        # Stars selection.
        psf_select, all_sources, n_not_satur = st_fwhm_select(
            pars['dmax'], pars['max_stars'], pars['thresh_level'],
            pars['fwhm_init'], sky_std, hdu_data)

        # FWHM selection.
        fwhm_estim, fwhm_min_rjct, ellip_rjct = psfmeasure.main(
            pars['dmax'], pars['ellip_max'], pars['fwhm_min'],
            psf_select, imname, hdu_data)

        if fwhm_estim:
            # Separate outliers.
            fwhm_no_outl, fwhm_mean, fwhm_std, fwhm_outl = rm_outliers(
                fwhm_estim)

            stars, psf_avrg = [], []
            if fwhm_no_outl:
                # Create .coo file.
                storeCooFile(out_path, imname, fwhm_no_outl)
                # Isolated stars, and average PSF.
                stars, psf_avrg = isolate_stars(
                    hdu_data, fwhm_no_outl, crop_side)
            else:
                print("  Could not obtain the average PSF.")
        else:
            fwhm_no_outl, fwhm_outl, fwhm_mean, fwhm_std, stars, psf_avrg =\
                [], [], -1., -1., [], []
            print("  WARNING: no stars left after rejecting\n  by min FWHM"
                  "and max ellipticity.")

        if pars['do_plots'] is 'y':
            make_plots(
                hdu_data, pars['max_stars'], crop_side,
                all_sources, n_not_satur, fwhm_min_rjct, ellip_rjct,
                fwhm_no_outl, fwhm_outl, fwhm_mean, fwhm_std, stars, psf_avrg,
                out_path, imname)

        # Store data to write to output file.
        im_name = imname.split('/')[-1]
        out_data.append(
            [im_name, hdr[pars['filter_key']], hdr[pars['exp_key']], sky_mean,
             sky_std, len(fwhm_no_outl), fwhm_mean, fwhm_std])

    if out_data:
        storeOutFile(out_path, out_data)

    print("\nFinished.")


if __name__ == "__main__":
    main()
