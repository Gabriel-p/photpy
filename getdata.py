
import os
import sys
from os.path import join, realpath, dirname
import numpy as np
import matplotlib.pyplot as plt
from pyraf import iraf

from astropy.io import fits
from astropy.io import ascii
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder
from photutils import CircularAperture


def bckg_data(hdulist, hdu_data, gain_key, rdnoise_key):
    """
    Estimate sky's mean, median, and standard deviation.
    """
    hdr = hdulist[0].header
    gain = hdr[gain_key]
    rdnoise = hdr[rdnoise_key]
    print("Gain, Rdnoise: {}, {}".format(gain, rdnoise))

    print("\nEstimating background mean, median, and STDDEV.")
    mean, median, std = sigma_clipped_stats(hdu_data, sigma=3.0, iters=5)
    print("Mean, median, and STDDEV: {}, {}, {}".format(mean, median, std))
    print("STDDEV estimated from GAIN, RDNOISE, and median: {}".format(
        np.sqrt(median * gain + rdnoise ** 2) / gain))
    return std


def st_fwhm_select(dmax, max_psf_stars, thresh_level, std, fwhm_init,
                   hdu_data):
    """
    Find stars in image with DAOStarFinder, filter out saturated stars (with
    peak values greater than the maximum accepted in 'dmax'), order by the
    largest (most negative) magnitude, and select a maximum of 'max_psf_stars'.
    """
    print("\nFinding stars.")
    thresh = thresh_level * std
    print('Threshold, initial FWHM: {}, {}'.format(thresh, fwhm_init))
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


def psfmeasure(dmax, ellip_max, psf_select, imname, hdu_data):
    """
    Use the IRAF task 'psfmeasure' to estimate the FWHM and ellipticity of
    all the stars in the 'psf_select' list.
    Reject those with a large ellipticity (> ellip_max).
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
        coords="markall", wcs="logical", display='no', frame=1, level=0.5,
        size="FWHM", beta='INDEF', scale=1., radius=15, sbuffer=5, swidth=5,
        saturation=dmax, ignore_sat='yes', iterations=5, xcenter='INDEF',
        ycenter='INDEF', logfile="psfmeasure", graphcur="cursor",
        images=imname, imagecur="positions", Stdout="/dev/null")

    # from imexam.imexamine import Imexamine
    # from imexam.math_helper import gfwhm
    # plots = Imexamine()
    # plots.set_data(hdu_data)

    fwhm_estim, psfmeasure_estim, ellip_rjct = [], 0, 0
    with open("psfmeasure", 'r') as f:
        for i, line in enumerate(f):
            data = line.split()
            if data:
                if data[0] != 'Average':
                    # First line (star) of output file.
                    if i == 3:
                        if float(data[5]) <= ellip_max:
                            fwhm_estim.append(
                                map(float, [data[1], data[2], data[4],
                                            data[5]]))
                        else:
                            ellip_rjct += 1
                    # Rest of the lines.
                    elif i > 3:
                        if float(data[4]) <= ellip_max:
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
                            ellip_rjct += 1
                else:
                    # Averaged FWHM by the 'psfmeasure' task.
                    psfmeasure_estim = float(data[-1])

    os.remove('positions')
    os.remove('cursor')
    os.remove('psfmeasure')

    return fwhm_estim, psfmeasure_estim, ellip_rjct


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


def make_plots(hdu_data, fwhm_no_outl, fwhm_mean, fwhm_std):
    """
    Make plots.
    """
    fig = plt.figure()
    ax0 = fig.add_subplot(121)
    median, std = np.median(hdu_data), np.std(hdu_data)
    plt.imshow(hdu_data, cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower', vmin=0., vmax=median + std)
    positions = (zip(*fwhm_no_outl)[0], zip(*fwhm_no_outl)[1])
    apertures = CircularAperture(positions, r=4.)
    apertures.plot(color='red', lw=1.5)
    ax0.set_title('Observed image with FWHM selected stars.')
    plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
    #
    ax1 = fig.add_subplot(122)
    ax1.hist(
        zip(*fwhm_no_outl)[2],
        range=[fwhm_mean - 2. * fwhm_std, fwhm_mean + 2. * fwhm_std],
        bins=int(len(fwhm_no_outl) * 0.2))
    ax1.set_title('FWHM distribution for {} stars.'.format(
        len(fwhm_no_outl)))
    plt.show()


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
    print("Keywords: {}, {}, {}, {}".format(
        gain_key, rdnoise_key, filter_key, exp_key))
    thresh_level = 5.
    fwhm_init = 3.
    print("Threshold level above STDDEV, initial FWHM: {}, {}".format(
        thresh_level, fwhm_init))
    dmax = 60000.
    ellip_max = 0.15
    max_psf_stars = 100
    print("Dmax, ellip_max, max_psf_stars: {}, {}".format(
        dmax, ellip_max, max_psf_stars))

    with open(join(r_path, "fwhm_final.dat"), 'w') as f:
        f.write("# image      filter   exposure    N (not saturated)    "
                "N (ellipticity>{})   FWHM (N stars)     FWHM (mean)    "
                "FWHM (std)\n".format(ellip_max))

    # For each .fits image in the root folder.
    for imname in fits_list:
        print("\n.fits file: {}".format(imname.replace(mypath, "")))

        # Extract data.
        hdulist = fits.open(imname)
        hdu_data = hdulist[0].data

        # Background estimation.
        std = bckg_data(hdulist, hdu_data, gain_key, rdnoise_key)

        # Stars selection.
        psf_select, n_not_satur = st_fwhm_select(
            dmax, max_psf_stars, thresh_level, std, fwhm_init, hdu_data)

        # FWHM selection.
        fwhm_estim, psfmeasure_estim, ellip_rjct = psfmeasure(
            dmax, ellip_max, psf_select, imname, hdu_data)

        # FWHM median an list with no outliers.
        fwhm_median, fwhm_no_outl = rm_outliers(fwhm_estim)

        print("\nMedian FWHM: {}".format(fwhm_median))
        fwhm_mean, fwhm_std = np.mean(zip(*fwhm_no_outl)[2]),\
            np.std(zip(*fwhm_no_outl)[2])
        print("Mean FWHM +/- std (no outliers, "
              "{} stars): {:.2f} +- {:.2f}".format(
                  len(fwhm_no_outl), fwhm_mean, fwhm_std))
        print("psfmeasure FWHM: {}".format(psfmeasure_estim))

        if do_plots:
            make_plots(hdu_data, fwhm_no_outl, fwhm_mean, fwhm_std)

        with open(join(r_path, "fwhm_final.dat"), 'a') as f:
            im_name = imname.split('/')[-1]
            filt = hdulist[0].header[filter_key]
            exptime = hdulist[0].header[exp_key]
            f.write("{}    {}    {}    {}    {}    {}    {:.2f}    "
                    "{:.3f}\n".format(
                        im_name, filt, exptime, n_not_satur, ellip_rjct,
                        len(fwhm_no_outl), fwhm_mean, fwhm_std))


if __name__ == "__main__":
    main()
