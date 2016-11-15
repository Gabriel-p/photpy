
import os
import sys
from os.path import join, realpath, dirname
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from pyraf import iraf

from astropy.io import fits
from astropy.io import ascii
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval

from photutils import DAOStarFinder
from photutils import CircularAperture
from photutils.utils import cutout_footprint


def bckg_data(hdulist, hdu_data, gain_key, rdnoise_key):
    """
    Estimate sky's mean, median, and standard deviation.
    """
    hdr = hdulist[0].header
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
    return sky_median, sky_std


def st_fwhm_select(dmax, max_psf_stars, thresh_level, std, fwhm_init,
                   hdu_data):
    """
    Find stars in image with DAOStarFinder, filter out saturated stars (with
    peak values greater than the maximum accepted in 'dmax'), order by the
    largest (most negative) magnitude, and select a maximum of 'max_psf_stars'.
    """
    print("\nFinding stars.")
    thresh = thresh_level * std
    print('Threshold, initial FWHM: {:.1f}, {:.1f}'.format(thresh, fwhm_init))
    stfind = DAOStarFinder(threshold=thresh, fwhm=fwhm_init)
    sources = stfind(hdu_data)
    print("Sources found: {}".format(len(sources)))
    print("Filter out saturated stars.")
    mask = sources['peak'] < dmax
    sour_no_satur = sources[mask]
    print("Non-saturated sources found: {}".format(len(sour_no_satur)))
    sour_no_satur.sort('mag')
    psf_select = sour_no_satur[:max_psf_stars]

    return psf_select, len(sour_no_satur)


def psfmeasure(dmax, ellip_max, fwhm_min, psf_select, imname, hdu_data):
    """
    Use the IRAF task 'psfmeasure' to estimate the FWHM and ellipticity of
    all the stars in the 'psf_select' list.
    Reject those with a large ellipticity (> ellip_max), and a very low
    FWHM (<fwhm_min).
    """
    print("\nRun 'psfmeasure' task to estimate the FWHMs.")
    ascii.write(
        psf_select, output='positions',
        include_names=['xcentroid', 'ycentroid'], format='fast_no_header')
    with open('cursor', 'w') as f:
        f.write('q\n')
    iraf.noao()
    iraf.obsutil(Stdout="/dev/null")
    iraf.psfmeasure(
        coords="mark1", wcs="logical", display='no', frame=1, level=0.5,
        size="FWHM", beta='INDEF', scale=1., radius=15, sbuffer=5, swidth=5,
        saturation=dmax, ignore_sat='yes', iterations=5, xcenter='INDEF',
        ycenter='INDEF', logfile="psfmeasure", graphcur="cursor",
        images=imname, imagecur="positions", Stdout="/dev/null")

    # from imexam.imexamine import Imexamine
    # from imexam.math_helper import gfwhm
    # plots = Imexamine()
    # plots.set_data(hdu_data)

    fwhm_estim, psfmeasure_estim, st_rjct = [], 0, 0
    with open("psfmeasure", 'r') as f:
        for i, line in enumerate(f):
            data = line.split()
            if data:
                if data[0] != 'Average':
                    # First line (star) of output file.
                    if i == 3:
                        if float(data[4]) > fwhm_min and\
                                float(data[5]) <= ellip_max:
                            fwhm_estim.append(
                                map(float, [data[1], data[2], data[4],
                                            data[5]]))
                        else:
                            st_rjct += 1
                    # Rest of the lines.
                    elif i > 3:
                        if float(data[3]) > fwhm_min and\
                                float(data[4]) <= ellip_max:
                            fwhm_estim.append(
                                map(float, [data[0], data[1], data[3],
                                            data[4]]))
                            # sys.stdout = open(os.devnull, "w")
                            # gauss_x = plots.line_fit(
                            #     float(data[0]), float(data[1]), genplot=False)
                            # gauss_y = plots.column_fit(
                            #     float(data[0]), float(data[1]), genplot=False)
                            # sys.stdout = sys.__stdout__
                            # print(float(data[3]), float(data[4]),
                            #       gfwhm(gauss_x.stddev)[0],
                            #       gfwhm(gauss_y.stddev)[0])
                        else:
                            st_rjct += 1
                else:
                    # Averaged FWHM by the 'psfmeasure' task.
                    psfmeasure_estim = float(data[-1])

    os.remove('positions')
    os.remove('cursor')
    os.remove('psfmeasure')

    return fwhm_estim, psfmeasure_estim, st_rjct


def rm_outliers(fwhm_estim):
    """
    Remove outlier starts with large FWHMs.
    """
    fwhm_median = np.median(zip(*fwhm_estim)[2])
    fwhm_no_outl = []
    for st in fwhm_estim:
        if st[2] <= 2. * fwhm_median:
            fwhm_no_outl.append(st)

    return fwhm_median, fwhm_no_outl


def isolate_stars(hdu_data, fwhm_no_outl, crop_side):
    """
    Create crop centered at each star used to obtain the FWHM.
    """
    print("\nCrop regions around selected stars.")
    stars = []
    for st in fwhm_no_outl:
        # Crop image
        crop = cutout_footprint(
            hdu_data, (st[0], st[1]), (crop_side, crop_side))
        # Skip stars placed on borders.
        if len(crop[0]) == crop_side:
            stars.append(crop[0])
        else:
            print("Skip star on border (xc, yc)=({}, {})".format(st[0], st[1]))

    return np.array(stars, dtype='float32')


def make_plots(r_path, imname, hdu_data, sky_median, sky_std, fwhm_no_outl,
               fwhm_mean, fwhm_std, stars, crop_side):
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
    positions = (zip(*fwhm_no_outl)[0], zip(*fwhm_no_outl)[1])
    apertures = CircularAperture(positions, r=4.)
    apertures.plot(color='red', lw=1.)
    ax.set_title('Observed image')
    # plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)

    ax = plt.subplot(gs[1:3, 8:10])
    ax.hist(
        zip(*fwhm_no_outl)[2],
        range=[fwhm_mean - 2. * fwhm_std, fwhm_mean + 2. * fwhm_std],
        bins=int(len(fwhm_no_outl) * 0.4))
    ax.axvline(fwhm_mean, ls='--', lw=2, c='r')
    ax.set_title('FWHM distribution ({} stars)'.format(
        len(fwhm_no_outl)), fontsize=11)
    ax.tick_params(axis='both', which='major', labelsize=9)

    ax = plt.subplot(gs[1:4, 5:8])
    ymax, xmax = np.shape(hdu_data)
    ax.set_xlim(0., xmax)
    ax.set_ylim(0., ymax)
    ax.scatter(*positions, s=np.exp(zip(*fwhm_no_outl)[2]), lw=0.7,
               facecolors='none', edgecolors='r')
    ax.set_title('Spatial distribution of (exponential) FWHM sizes',
                 fontsize=11)
    ax.tick_params(axis='both', which='major', labelsize=9)

    gs1 = gridspec.GridSpec(10, 12)
    gs1.update(wspace=0.025, hspace=0.25)
    max_s = min(25, len(stars))
    for i in range(max_s):
        print("(xc, yc), FWHM, ellipticity: ({:.2f}, {:.2f}), {}, {}".format(
            fwhm_no_outl[i][0], fwhm_no_outl[i][1], fwhm_no_outl[i][2],
            fwhm_no_outl[i][3]))
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
            fwhm_no_outl[i][0], fwhm_no_outl[i][1]), fontsize=10, y=-0.2)
        ax.set_axis_off()

    # Create a meshgrid of coordinates (0,1,...,N) times (0,1,...,N)
    y, x = np.mgrid[:len(stars[0, :, 0]), :len(stars[0, 0, :])]
    # duplicating the grids
    xcoord, ycoord = np.array([x] * len(stars)), np.array([y] * len(stars))
    # compute histogram with coordinates as x,y
    h, xe, ye = np.histogram2d(
        xcoord.ravel(), ycoord.ravel(),
        bins=[len(stars[0, 0, :]), len(stars[0, :, 0])], weights=stars.ravel())
    # The argmax of the combined stars.
    combined_star_argmax = np.unravel_index(np.argmax(h), h.shape)

    ax = plt.subplot(gs[5:10, 5:10])
    ax.imshow(h, cmap=plt.cm.viridis, interpolation='nearest',
              origin='lower', vmin=0.)
    # cent = int(crop_side / 2.) - 1
    ax.axhline(combined_star_argmax[0], ls='--', lw=2, c='w')
    ax.axvline(combined_star_argmax[1], ls='--', lw=2, c='w')
    apertures = CircularAperture(combined_star_argmax, r=fwhm_mean)
    apertures.plot(color='r', lw=0.75)
    ax.set_title("Average PSF model (FWHM={:.2f})".format(fwhm_mean))

    fig.tight_layout()
    fig_name = imname.split('/')[-1].replace(".fits", "")
    plt.savefig(join(r_path + '/' + fig_name + '.png'), dpi=300,
                bbox_inches='tight')


def main():
    """
    Get FWHM, sky mean, and sky standard deviation from a .fits file.
    """
    mypath = realpath(join(os.getcwd(), dirname(__file__)))
    fname = raw_input("Root dir where the .fits files are stored: ")
    fname = fname[1:] if fname.startswith('/') else fname
    r_path = join(mypath, fname)
    fits_list = []
    if os.path.isdir(r_path):
        for subdir, dirs, files in os.walk(r_path):
            for file in files:
                if file.endswith('.fits'):
                    fits_list.append(os.path.join(subdir, file))
    else:
        print("Not a folder. Exit.")
        sys.exit()

    do_plots = raw_input("Show plot for FWHM selected stars? (y/n): ")
    do_plots = True if do_plots in ('y', 'Y', 'yes', 'YES') else False

    # Some header info, and hard-coded criteria for selection of stars.
    gain_key = 'EGAIN'
    rdnoise_key = 'ENOISE'
    filter_key = 'FILTER'
    exp_key = 'EXPTIME'
    print("\nKeywords: {}, {}, {}, {}".format(
        gain_key, rdnoise_key, filter_key, exp_key))
    thresh_level = 5.
    fwhm_init = 3.
    print("Threshold level above STDDEV, initial FWHM: {}, {}".format(
        thresh_level, fwhm_init))
    dmax = 60000.
    ellip_max = 0.15
    fwhm_min = 1.5
    max_psf_stars = 100
    crop_side = 20
    print("Dmax, ellip_max, fwhm_min, max_psf_stars, crop side: "
          "{}, {}, {}, {}, {}".format(
              dmax, ellip_max, fwhm_min, max_psf_stars, crop_side))

    with open(join(r_path, "fwhm_final.dat"), 'w') as f:
        f.write("# image      filter   exposure    N (not saturated)    "
                "(% rejected)   FWHM (N stars)     FWHM (mean)    "
                "FWHM (std)\n")

    # For each .fits image in the root folder.
    for imname in fits_list:
        print("\n.fits file: {}".format(imname.replace(mypath, "")))

        # Extract data.
        hdulist = fits.open(imname)
        hdu_data = hdulist[0].data

        # Background estimation.
        sky_median, sky_std = bckg_data(
            hdulist, hdu_data, gain_key, rdnoise_key)

        # Stars selection.
        psf_select, n_not_satur = st_fwhm_select(
            dmax, max_psf_stars, thresh_level, sky_std, fwhm_init, hdu_data)

        # FWHM selection.
        fwhm_estim, psfmeasure_estim, st_rjct = psfmeasure(
            dmax, ellip_max, fwhm_min, psf_select, imname, hdu_data)

        # FWHM median an list with no outliers.
        fwhm_median, fwhm_no_outl = rm_outliers(fwhm_estim)

        print("\nMedian FWHM: {}".format(fwhm_median))
        fwhm_mean, fwhm_std = np.mean(zip(*fwhm_no_outl)[2]),\
            np.std(zip(*fwhm_no_outl)[2])
        print("Mean FWHM +/- std (no outliers, "
              "{} stars): {:.2f} +- {:.2f}".format(
                  len(fwhm_no_outl), fwhm_mean, fwhm_std))
        print("psfmeasure FWHM: {}".format(psfmeasure_estim))

        stars = isolate_stars(hdu_data, fwhm_no_outl, crop_side)

        if do_plots:
            make_plots(
                r_path, imname, hdu_data, sky_median, sky_std, fwhm_no_outl,
                fwhm_mean, fwhm_std, stars, crop_side)

        with open(join(r_path, "fwhm_final.dat"), 'a') as f:
            im_name = imname.split('/')[-1]
            filt = hdulist[0].header[filter_key]
            exptime = hdulist[0].header[exp_key]
            f.write("{}    {}    {}    {}    {:.1f}    {}    {:.2f}    "
                    "{:.3f}\n".format(
                        im_name, filt, exptime, n_not_satur,
                        100. * (st_rjct / float(max_psf_stars)),
                        len(fwhm_no_outl), fwhm_mean, fwhm_std))


if __name__ == "__main__":
    main()
