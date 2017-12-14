
import read_pars_file as rpf
from find_stars import readStats

import os
from os.path import exists, join, isfile
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import ascii, fits
from photutils.psf import IntegratedGaussianPRF
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import ZScaleInterval


def in_params():
    """
    Read and prepare input parameter values.
    """
    pars = rpf.main()

    in_path = join(pars['mypath'].replace('tasks', 'input'), pars['fits_psf'])

    if isfile(in_path):
        fits_list = [in_path]
        out_path = join(
            pars['mypath'].replace('tasks', 'output'),
            pars['fits_psf'].split('/')[0])
    else:
        fits_list = []
        for file in os.listdir(in_path):
            f = join(in_path, file)
            if isfile(f):
                fits_list.append(f)

        if not fits_list:
            print("No '*.fits' files found in '{}' folder.".format(in_path))
            sys.exit()

        out_path = join(
            pars['mypath'].replace('tasks', 'output'), pars['fits_psf'])

    return pars, out_path, fits_list


def readData(pars, imname):
    """
    """
    # Load .fits file.
    hdulist = fits.open(imname)
    # Extract header and data.
    hdr, hdu_data = hdulist[0].header, hdulist[0].data

    return hdr, hdu_data


def getPRF(out_path, imname, filt, hdu_data):
    """
    """
    from photutils.psf.sandbox import DiscretePRF
    psf_file = ascii.read(
        join(out_path, 'filt_' + filt,
             imname.split('/')[-1].replace('.fits', '.coo')))
    psf_coords = psf_file['x', 'y']
    psf_coords.rename_column('x', 'x_0')
    psf_coords.rename_column('y', 'y_0')
    sub_s = 2
    prf_discrete = DiscretePRF.create_from_image(
        hdu_data, psf_coords, size=25, subsampling=sub_s)
    print("  PRF created.")

    fig, axes = plt.subplots(nrows=5, ncols=5)
    fig.set_size_inches(12, 9)
    # Plot kernels
    for i in range(sub_s):
        for j in range(sub_s):
            prf_image = prf_discrete._prf_array[i, j]
            im = axes[i, j].imshow(prf_image, interpolation='None')
    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    plt.colorbar(im, cax=cax)
    plt.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
    plt.show()

    return prf_discrete


def readSources(out_path, imname, filt_val):
    """
    """

    src_f = join(
        out_path, 'filt_' + filt_val,
        imname.split('/')[-1].replace('fits', 'src'))

    sources = ascii.read(src_f)

    return sources


def psfPhot(
        fwhm, niters, fitshape, psf_thresh, group_sep, hdu_data, sources,
        prf_discrete):
    """
    Class to calculate the background in an array...

    MMMBackground            using the DAOPHOT MMM algorithm.
    MeanBackground           as the (sigma-clipped) mean.
    MedianBackground         as the (sigma-clipped) median.
    ModeEstimatorBackground  using a mode estimator of the form
                             (median_factor * median) - (mean_factor * mean).
    SExtractorBackground     using the SExtractor algorithm.

    """

    # from photutils.background import MMMBackground
    # bkg = MMMBackground()

    # from astropy.stats import SigmaClip
    # from photutils import MedianBackground
    # sigma_clip = SigmaClip(sigma=3.)
    # bkg = MedianBackground(sigma_clip)

    # from photutils.detection import IRAFStarFinder
    # finder = IRAFStarFinder(
    #     threshold=psf_thresh * std, fwhm=fwhm, minsep_fwhm=0.01, roundhi=5.0,
    #     roundlo=-5.0, sharplo=0.0, sharphi=2.0)

    # from photutils.psf import DAOGroup
    # daogroup = DAOGroup(group_sep * fwhm)

    psf_model = IntegratedGaussianPRF(sigma=fwhm * gaussian_fwhm_to_sigma)

    # psf_model = prf_discrete

    # from photutils.psf import IterativelySubtractedPSFPhotometry
    # photometry = IterativelySubtractedPSFPhotometry(
    #     group_maker=daogroup, bkg_estimator=bkg, finder=finder,
    #     psf_model=psf_model, fitter=LevMarLSQFitter(), niters=niters,
    #     fitshape=fitshape)

    # from photutils.psf import BasicPSFPhotometry
    # photometry = BasicPSFPhotometry(
    #     group_maker=daogroup, bkg_estimator=bkg, finder=finder,
    #     psf_model=psf_model, fitter=LevMarLSQFitter(),
    #     fitshape=fitshape)

    # psf_model.x_0.fixed = True
    # psf_model.y_0.fixed = True

    from photutils.psf import DAOPhotPSFPhotometry
    photometry = DAOPhotPSFPhotometry(
        crit_separation=group_sep * fwhm, threshold=psf_thresh, fwhm=fwhm,
        psf_model=psf_model, fitshape=fitshape, fitter=LevMarLSQFitter(),
        aperture_radius=fwhm, niters=niters)

    result_tab = photometry(image=hdu_data, init_guesses=sources)
    print("PSF performed.")
    residual_image = photometry.get_residual_image()

    return result_tab, residual_image


def apertCorrect(result_tab):
    """
    """

    return result_tab


def writePSF(out_path, imname, filt_val, results):
    """
    """
    out_file = join(
        out_path, 'filt_' + filt_val,
        imname.split('/')[-1].replace('.fits', '_psf.out'))
    ascii.write(
        results, out_file, format='fixed_width', delimiter=' ',
        formats={
            'x_fit': '%10.4f', 'y_fit': '%10.4f', 'flux_fit': '%10.4f',
            'flux_unc': '%10.4f'},
        fill_values=[(ascii.masked, 'nan')], overwrite=True)


def makePlot(out_path, imname, filter_val, hdu_data, sources, residual_image):
    """
    """
    print("Plotting.")
    interval = ZScaleInterval()

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(10, 10)
    plt.subplot(gs[0:10, 0:10])
    zmin, zmax = interval.get_limits(hdu_data)
    plt.imshow(hdu_data, cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower', vmin=zmin, vmax=zmax,
               extent=[0., hdu_data.shape[1], 0., hdu_data.shape[0]])
    plt.scatter(sources['x_0'], sources['y_0'], marker='.', c='r', s=1, lw=.0)
    plt.xlim(0., hdu_data.shape[1])
    plt.ylim(0., hdu_data.shape[0])
    #
    fig.tight_layout()
    fig_name = join(
        out_path, 'filt_' + filter_val,
        imname.split('/')[-1].replace(".fits", "_orig"))
    plt.savefig(fig_name + '.png', dpi=250, bbox_inches='tight')
    # Close to release memory.
    plt.clf()
    plt.close("all")

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(10, 10)
    plt.subplot(gs[0:10, 0:10])
    # plt.title('Residual Image')
    zmin, zmax = interval.get_limits(residual_image)
    plt.imshow(
        residual_image, cmap='viridis', aspect=1, interpolation='nearest',
        origin='lower', vmin=zmin, vmax=zmax,
        extent=[0., hdu_data.shape[1], 0., hdu_data.shape[0]])
    plt.xlim(0., hdu_data.shape[1])
    plt.ylim(0., hdu_data.shape[0])
    #
    fig.tight_layout()
    fig_name = join(
        out_path, 'filt_' + filter_val,
        imname.split('/')[-1].replace(".fits", "_psf"))
    plt.savefig(fig_name + '.png', dpi=150, bbox_inches='tight')
    # Close to release memory.
    plt.clf()
    plt.close("all")


def main():
    """
    """
    pars, out_path, fits_list = in_params()

    # For each .fits image in the root folder.
    for imname in fits_list:
        hdr, hdu_data = readData(pars, imname)

        if hdr[pars['filter_key']] == pars['filter_psf']:
            import time as t
            s = t.clock()

            print("\n* File: {}".format(
                imname.replace(pars['mypath'].replace('tasks', 'input'), "")))
            filt, exp_time, airmass = hdr[pars['filter_key']],\
                hdr[pars['exposure_key']], hdr[pars['airmass_key']]
            print("  Filter {}, Exp time {}, Airmass {}".format(
                filt, exp_time, airmass))

            fwhm, sky_mean, sky_std = readStats(
                out_path, imname.split('/')[-1])
            print("  FWHM: {:.2f}".format(fwhm))

            # TODO finish PRF creation
            # prf_discrete = getPRF(out_path, imname, filt, hdu_data)
            prf_discrete = []

            sources = readSources(out_path, imname, filt)
            print("  Sources read: {}".format(len(sources)))

            # Odd int
            c_fwhm = np.ceil(float(pars['fitshape']) * fwhm)
            fitshape = int(c_fwhm + 1 if (c_fwhm % 2) < 1. else c_fwhm)
            niters = int(pars['niters'])
            print("  fitshape, group_sep: {:.2f}, {:.2f}".format(
                fitshape, float(pars['group_sep'])))

            result_tab, residual_image = psfPhot(
                fwhm, niters, fitshape, float(pars['psf_thresh']),
                float(pars['group_sep']), hdu_data, sources, prf_discrete)

            print("Perform aperture correction.")
            result_tab = apertCorrect(result_tab)

            # Generate output dir/subdir if it doesn't exist.
            if not exists(join(out_path, 'filt_' + filt)):
                os.makedirs(join(out_path, 'filt_' + filt))

            writePSF(out_path, imname, filt, result_tab)

            if pars['do_plots_G'] is 'y':
                makePlot(
                    out_path, imname, filt, hdu_data, sources, residual_image)

            print(t.clock() - s)


if __name__ == '__main__':
    main()
