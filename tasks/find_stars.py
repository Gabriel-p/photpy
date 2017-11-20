
import read_pars_file as rpf

import os
from os.path import exists, join, isfile
import sys

import numpy as np
# For server
import matplotlib
matplotlib.use('Agg')
# For server
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
        print("No '*.fits' files found in '{}' folder.".format(in_path))
        sys.exit()

    return pars, in_path, out_path, fits_list


def readData(pars, in_path, imname):
    """
    Read data and header from input .fits file.
    """
    # Load .fits file.
    hdulist = fits.open(join(in_path, imname))
    # Extract header and data.
    hdr, hdu_data = hdulist[0].header, hdulist[0].data

    return hdr, hdu_data


def readStats(out_path, imname):
    """
    Read stats (FWHM, sky mean and STDDEV) for the processed frame, from the
    *existing* 'fitstats.dat' file.
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


def strFind(
        hdu_data, dmax, fwhm, sky_mean, sky_std, find_method, thrsh,
        round_method, round_max):
    """
    Detect sources according to the parameters given.

    IRAF: IRAFStarFinder
    DAO:  DAOStarFinder
    """
    if round_method == 'AND':

        if find_method == 'IRAF':
            finder = IRAFStarFinder(
                threshold=thrsh * sky_std, fwhm=fwhm, roundlo=-round_max,
                roundhi=round_max)
        elif find_method == 'DAO':
            finder = DAOStarFinder(
                threshold=thrsh * sky_std, fwhm=fwhm, roundlo=-round_max,
                roundhi=round_max)
        sources = finder(hdu_data - sky_mean)
        print("Sources after roundness '{}' filter: {}".format(
            round_method, len(sources)))

        mask = sources['peak'] < dmax
        sources = sources[mask]
        print("Sources after 'dmax' filter: {}".format(len(sources)))

    elif round_method == 'OR':

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
        print("Sources after roundness '{}' filter: {}".format(
            round_method, len(sources)))

        mask = sources['peak'] < dmax
        sources = sources[mask]
        print("Sources after 'dmax' filter: {}".format(len(sources)))

    return sources


def writeSources(out_path, imname, filt_val, sources):
    """
    """
    out_file = join(
        out_path, 'filt_' + filt_val,
        imname.split('/')[-1].replace('fits', 'src'))
    ascii.write(
        sources['xcentroid', 'ycentroid', 'flux'],
        out_file, names=['x_0', 'y_0', 'flux_0'], format='fixed_width',
        delimiter=' ',
        formats={'x_0': '%10.4f', 'y_0': '%10.4f', 'flux_0': '%10.4f'},
        fill_values=[(ascii.masked, 'nan')], overwrite=True)


def makePlot(
        out_path, imname, filter_val, hdu_data, sky_mean, sky_std, thresh_find,
        sources):
    """
    """
    print("Plotting.")
    fig = plt.figure(figsize=(20, 10))
    gs = gridspec.GridSpec(10, 20)

    plt.subplot(gs[0:10, 0:10])
    plt.title("Sources detected: {} (thresh = {:.1f})".format(
        len(sources), thresh_find))
    plt.scatter(sources['xcentroid'], sources['ycentroid'],
                marker='.', c='r', s=1, lw=.5)
    interval = ZScaleInterval()
    zmin, zmax = interval.get_limits(hdu_data)
    plt.imshow(hdu_data, cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower', vmin=zmin, vmax=zmax)

    plt.subplot(gs[0:5, 10:15])
    plt.scatter(sources['mag'], sources['roundness1'], s=5, label='round1')
    plt.scatter(sources['mag'], sources['roundness2'], s=5, label='round2')
    plt.xlabel("mag")
    plt.legend()

    plt.subplot(gs[0:5, 15:20])
    plt.hist(sources['sharpness'], bins=30, alpha=.5, label='sharpness')
    # plt.xlabel("mag")
    plt.legend()

    plt.subplot(gs[5:10, 10:15])
    plt.title("Center zoom (5% of length)")
    plt.scatter(sources['xcentroid'], sources['ycentroid'],
                marker='.', c='r', s=1)
    interval = ZScaleInterval()
    zmin, zmax = interval.get_limits(hdu_data)
    plt.imshow(hdu_data, cmap='viridis', interpolation='nearest',
               origin='lower', vmin=zmin, vmax=zmax)
    # Center of image
    xcent, ycent = np.array(hdu_data.shape) / 2.
    # x,y 5% lengths
    xl, yl = np.array(hdu_data.shape) * .05
    plt.xlim(xcent - xl, xcent + xl)
    plt.ylim(ycent - yl, ycent + yl)

    plt.subplot(gs[5:10, 15:20])
    plt.xlabel("flux")
    plt.ylabel("peak")
    plt.scatter(sources['flux'], sources['peak'])

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
        hdr, hdu_data = readData(pars, in_path, imname)

        if hdr[pars['filter_key']] == pars['filter_proc']:

            print("\n* File: {}".format(
                imname.replace(pars['mypath'].replace('tasks', 'input'), "")))
            print("Filter {}, Exp time {}, Airmass {}".format(
                hdr[pars['filter_key']], hdr[pars['exposure_key']],
                hdr[pars['airmass_key']]))

            fwhm, sky_mean, sky_std = readStats(
                out_path, imname.split('/')[-1])
            print("FWHM ({:.2f}), Sky mean ({:.2f}), Sky std ({:.2f}).".format(
                fwhm, sky_mean, sky_std))

            print("Params: {}, {}, {}, {}".format(
                pars['find_method'], float(pars['thresh_find']),
                pars['round_method'], float(pars['round_max'])))
            sources = strFind(
                hdu_data, float(pars['dmax']), fwhm, sky_mean, sky_std,
                pars['find_method'], float(pars['thresh_find']),
                pars['round_method'], float(pars['round_max']))

            # Generate output dir/subdir if it doesn't exist.
            if not exists(join(out_path, 'filt_' + hdr[pars['filter_key']])):
                os.makedirs(join(out_path, 'filt_' + hdr[pars['filter_key']]))

            writeSources(out_path, imname, hdr[pars['filter_key']], sources)

            if pars['do_plots_F'] is 'y':
                makePlot(
                    out_path, imname, hdr[pars['filter_key']], hdu_data,
                    sky_mean, sky_std, float(pars['thresh_find']), sources)


if __name__ == '__main__':
    main()
