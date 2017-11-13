
import read_pars_file as rpf

# import os
from os.path import join  # , isfile
# import sys

# import numpy as np
# import matplotlib.pyplot as plt
# from photutils import CircularAperture

from astropy.io import ascii, fits
# from photutils.background import MADStdBackgroundRMS
from photutils.detection import IRAFStarFinder
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats


def in_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    in_path = pars['mypath'].replace('tasks', 'input')
    out_path = pars['mypath'].replace('tasks', 'output')

    # fits_list = []
    # for file in os.listdir(in_out_path):
    #     f = join(in_out_path, file)
    #     if isfile(f):
    #         if f.endswith('_crop.fits'):
    #             fits_list.append(f)

    # if not fits_list:
    #     print("No '*_crop.fits' files found in 'output/standards' folder."
    #           " Exit.")
    #     sys.exit()

    return pars, in_path, out_path


def readData(pars, in_path, imname):
    """
    """
    # Load .fits file.
    hdulist = fits.open(join(in_path, imname))
    # Extract header and data.
    hdr, hdu_data = hdulist[0].header, hdulist[0].data
    filt, exp_time, airmass = hdr[pars['filter_key']],\
        hdr[pars['exposure_key']], hdr[pars['airmass_key']]
    print("{}: Filter {}, Exp time {}, Airmass {}".format(
        imname, filt, exp_time, airmass))

    return hdu_data


def strFind(hdu_data, find_method, thrsh, fwhm):
    """
    """
    # bkgrms = MADStdBackgroundRMS()
    # std = bkgrms(hdu_data)
    mean, median, std = sigma_clipped_stats(hdu_data, sigma=3.0, iters=5)
    print(mean, median, std)
    print("Sky median & STDDEV estimated.")

    if find_method == 'IRAF':
        finder = IRAFStarFinder(threshold=thrsh * std, fwhm=fwhm)
    elif find_method == 'DAO':
        finder = DAOStarFinder(threshold=thrsh * std, fwhm=fwhm)
    sources = finder(hdu_data - median)
    print("Sources found.")

    mask1 = (sources['roundness1'] > -.2) & (sources['roundness1'] < .2)
    mask2 = (sources['roundness2'] > -.2) & (sources['roundness2'] < .2)
    mask = mask1 | mask2

    # plt.scatter(sources[mask]['xcentroid'], sources[mask]['ycentroid'],
    #             marker='.', c='r', s=2)
    # plt.imshow(hdu_data, origin='lower', cmap='Greys_r', vmin=0.,
    #            vmax=np.median(hdu_data) + np.std(hdu_data))
    # plt.show()

    sources[mask]
    return sources


def writeSources(out_path, imname, sources):
    """
    """
    # tables = []
    # for v in filters.values():
    #     tables.append(Table(zip(*v)))
    # aper_phot = Table(
    #     vstack(tables),
    #     names=('Filt', 'Stnd_field', 'ID', 'file', 'exp_t', 'A', 'ZA_mag',
    #            'Col', 'Mag'))

    out_file = join(out_path, imname.replace('fits', 'src'))
    ascii.write(
        sources, out_file, format='fixed_width', delimiter=' ',
        fill_values=[(ascii.masked, 'nan')], overwrite=True)
    # formats={'ZA_mag': '%10.4f'}


def main():
    """
    """
    pars, in_path, out_path = in_params()
    hdu_data = readData(pars, in_path, pars['fits_find'])

    fwhm = 3.
    sources = strFind(
        hdu_data, pars['find_method'], float(pars['threshold']), fwhm)

    writeSources(out_path, pars['fits_find'], sources)


if __name__ == '__main__':
    main()

# Source Detection
# https://photutils.readthedocs.io/en/stable/detection.html