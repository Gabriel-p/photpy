
"""
Estimate the FWHM, sky mean and STDDEV from an image or a group of images.
"""

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

    return psf_select


def psfmeasure(dmax, ellip_max, psf_select, imname):
    """
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

    fwhm_estim, psfmeasure_estim = [], 0
    with open("psfmeasure", 'r') as f:
        for i, line in enumerate(f):
            data = line.split()
            if data:
                if data[0] != 'Average':
                    if i == 3:
                        if float(data[5]) <= ellip_max:
                            fwhm_estim.append(
                                map(float, [data[1], data[2], data[4],
                                            data[5]]))
                    elif i > 3:
                        if float(data[4]) <= ellip_max:
                            fwhm_estim.append(
                                map(float, [data[0], data[1], data[3],
                                            data[4]]))
                else:
                    psfmeasure_estim = float(data[-1])

    os.remove('positions')
    os.remove('cursor')
    os.remove('psfmeasure')

    return fwhm_estim, psfmeasure_estim


def rm_outliers(fwhm_estim):
    """
    Remove outliers.
    """
    fwhm_median = np.median(zip(*fwhm_estim)[2])
    fwhm_no_outl = []
    for st in fwhm_estim:
        if st[2] <= 2. * fwhm_median:
            fwhm_no_outl.append(st)

    return fwhm_median, fwhm_no_outl


def make_plots(hdu_data, fwhm_no_outl, fwhm_mean, fwhm_std):
    """
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
    fname = raw_input("Root dir where .fits files are stored: ")
    r_path = join(mypath, fname)
    #
    gain_key = 'EGAIN'
    rdnoise_key = 'ENOISE'
    thresh_level = 5.
    fwhm_init = 3.
    dmax = 60000.
    ellip_max = 0.15
    max_psf_stars = 100
    print("Dmax, ellip_max, max_psf_stars: {}, {}".format(
        dmax, ellip_max, max_psf_stars))

    fits_list = []
    if os.path.isdir(r_path):
        for subdir, dirs, files in os.walk(r_path):
            for file in files:
                if file.endswith('.fits'):
                    fits_list.append(os.path.join(subdir, file))
    else:
        print("Not a folder. Exit.")
        sys.exit()

    for imname in fits_list:
        print("\n.fits file: {}".format(imname))

        # Extract data.
        hdulist = fits.open(imname)
        hdu_data = hdulist[0].data

        std = bckg_data(hdulist, hdu_data, gain_key, rdnoise_key)

        psf_select = st_fwhm_select(dmax, max_psf_stars, thresh_level, std,
                                    fwhm_init, hdu_data)

        fwhm_estim, psfmeasure_estim = psfmeasure(
            dmax, ellip_max, psf_select, imname)

        fwhm_median, fwhm_no_outl = rm_outliers(fwhm_estim)

        print("\nMedian FWHM: {}".format(fwhm_median))
        fwhm_mean, fwhm_std = np.mean(zip(*fwhm_no_outl)[2]),\
            np.std(zip(*fwhm_no_outl)[2])
        print("Mean FWHM +/- std (no outliers, "
              "{} stars): {:.2f} +- {:.2f}".format(
                  len(fwhm_no_outl), fwhm_mean, fwhm_std))
        print("psfmeasure FWHM: {}".format(psfmeasure_estim))

        make_plots(hdu_data, fwhm_no_outl, fwhm_mean, fwhm_std)


if __name__ == "__main__":
    main()
