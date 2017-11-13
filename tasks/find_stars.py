
import read_pars_file as rpf

import os
from os.path import join, isfile
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# from photutils import CircularAperture

from astropy.io import ascii, fits
from astropy.visualization import ZScaleInterval
# from photutils.background import MADStdBackgroundRMS
from photutils.detection import IRAFStarFinder
from photutils import DAOStarFinder


def in_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    in_path = join(pars['mypath'].replace('tasks', 'input'), pars['fits_find'])
    out_path = join(
        pars['mypath'].replace('tasks', 'output'), pars['fits_find'])

    fits_list = []
    for file in os.listdir(in_path):
        f = join(in_path, file)
        if isfile(f):
            fits_list.append(f)

    if not fits_list:
        print("No '*_crop.fits' files found in 'output/standards' folder."
              " Exit.")
        sys.exit()

    return pars, in_path, out_path, fits_list


def readData(pars, in_path, imname):
    """
    """
    # Load .fits file.
    hdulist = fits.open(join(in_path, imname))
    # Extract header and data.
    hdr, hdu_data = hdulist[0].header, hdulist[0].data
    filt, exp_time, airmass = hdr[pars['filter_key']],\
        hdr[pars['exposure_key']], hdr[pars['airmass_key']]
    print("Filter {}, Exp time {}, Airmass {}".format(
        filt, exp_time, airmass))

    return hdr, hdu_data


def readStats(out_path, imname):
    """
    """

    fitdata = ascii.read(join(out_path, 'fitstats.dat'))
    im_found = False
    for r in fitdata:
        if r['image'] == imname:
            fwhm, sky_mean, sky_std = r['FWHM_(mean)'], r['Sky_mean'],\
                r['Sky_STDDEV']
            im_found = True

    if not im_found:
        print("{} not found in 'fitstats.dat'  file.".format(imname))
        sys.exit()

    return fwhm, sky_mean, sky_std


def strFind(hdu_data, fwhm, sky_mean, sky_std, find_method, thrsh, round_max):
    """
    """
    if find_method == 'IRAF':
        finder = IRAFStarFinder(threshold=thrsh * sky_std, fwhm=fwhm)
    elif find_method == 'DAO':
        finder = DAOStarFinder(threshold=thrsh * sky_std, fwhm=fwhm)
    sources = finder(hdu_data - sky_mean)
    print("Sources found ({}).".format(len(sources)))

    mask1 = (sources['roundness1'] > -round_max) &\
        (sources['roundness1'] < round_max)
    mask2 = (sources['roundness2'] > -round_max) &\
        (sources['roundness2'] < round_max)
    mask = mask1 | mask2
    sources = sources[mask]
    print("Sources filtered ({}).".format(len(sources)))

    return sources


def writeSources(out_path, imname, filt_val, sources):
    """
    """
    out_file = join(
        out_path, 'filt_' + filt_val,
        imname.split('/')[-1].replace('fits', 'src'))
    ascii.write(
        sources['xcentroid', 'ycentroid'], out_file, format='fixed_width',
        delimiter=' ', formats={'xcentroid': '%10.4f', 'ycentroid': '%10.4f'},
        fill_values=[(ascii.masked, 'nan')], overwrite=True)


def makePlot(
        out_path, imname, filter_val, hdu_data, sky_mean, sky_std, sources):
    """
    """
    print("Plotting.")
    fig = plt.figure(figsize=(20, 10))
    gs = gridspec.GridSpec(10, 20)

    plt.subplot(gs[0:10, 0:10])
    plt.scatter(sources['xcentroid'], sources['ycentroid'],
                marker='.', c='r', s=1, lw=.5)
    interval = ZScaleInterval()
    zmin, zmax = interval.get_limits(hdu_data)
    plt.imshow(hdu_data, cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower', vmin=zmin, vmax=zmax)

    plt.subplot(gs[0:5, 10:15])
    plt.hist(sources['roundness1'], bins=30, alpha=.5, label='round1')
    plt.hist(sources['roundness2'], bins=30, alpha=.5, label='round2')
    plt.legend()

    plt.subplot(gs[0:5, 15:20])
    plt.hist(sources['sharpness'], bins=30, alpha=.5, label='sharpness')
    plt.legend()

    plt.subplot(gs[5:10, 10:15])
    plt.xlabel("flux")
    plt.ylabel("peak")
    plt.scatter(sources['flux'], sources['peak'])

    # xmax = np.median(sources['flux']) + 3. * np.std(sources['flux'])
    # print("flux", xmax)
    # plt.xlim(0., xmax)
    # mask = sources['flux'] < xmax
    # plt.hist(sources['flux'][mask], alpha=.5, label='flux')
    # plt.legend()

    # plt.subplot(gs[5:10, 15:20])
    # xmax = np.median(sources['peak']) + 3. * np.std(sources['peak'])
    # print("peak", xmax)
    # plt.xlim(0., xmax)
    # mask = sources['peak'] < xmax
    # plt.hist(sources['peak'][mask], alpha=.5, label='peak')
    # plt.legend()

    fig.tight_layout()
    fig_name = join(
        out_path, 'filt_' + filter_val,
        imname.split('/')[-1].replace(".fits", "_find"))
    plt.savefig(fig_name + '.png', dpi=150, bbox_inches='tight')
    # Close to release memory.
    plt.clf()
    plt.close("all")


def main():
    """
    """
    pars, in_path, out_path, fits_list = in_params()

    # For each .fits image in the root folder.
    for imname in fits_list:
        print("\n* File: {}".format(
            imname.replace(pars['mypath'].replace('tasks', 'input'), "")))

        hdr, hdu_data = readData(pars, in_path, imname)
        if hdr[pars['filter_key']] == pars['filter_proc']:

            fwhm, sky_mean, sky_std = readStats(
                out_path, imname.split('/')[-1])
            print("Sky mean ({:.2f}) & STDDEV ({:.2f}) estimated.".format(
                sky_mean, sky_std))

            sources = strFind(
                hdu_data, fwhm, sky_mean, sky_std, pars['find_method'],
                float(pars['thresh_find']), float(pars['round_max']))

            writeSources(out_path, imname, hdr[pars['filter_key']], sources)

            if pars['do_plots_F'] is 'y':
                makePlot(
                    out_path, imname, hdr[pars['filter_key']], hdu_data,
                    sky_mean, sky_std, sources)


if __name__ == '__main__':
    main()
