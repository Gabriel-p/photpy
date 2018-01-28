
import os
from os.path import join, realpath, dirname
import numpy as np
import matplotlib.pyplot as plt

from pyraf import iraf
import astropy.units as u
from astropy.table import Table
from astropy.io import fits
from astropy.io import ascii
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
from photutils.utils import calc_total_error

from scipy.stats import sigmaclip


def load_data(imname):
    """
    """
    mypath = realpath(join(os.getcwd(), dirname(__file__)))
    image_file = mypath + '/input/' + imname + '.fits'
    coo_file = mypath + '/output/' + imname + '.coo.1'
    print("Fits: {}".format(image_file))
    print("Coordinates: {}".format(coo_file))
    # Load .fits file.
    hdulist = fits.open(image_file)
    hdu_data, hdr = hdulist[0].data, hdulist[0].header
    exp_time, gain, rdnoise = hdr['EXPTIME'], hdr['EGAIN'], hdr['ENOISE']
    print("Itime, Gain, Rdnoise: {}, {}, {}".format(exp_time, gain, rdnoise))

    return image_file, coo_file, hdu_data, exp_time, gain, rdnoise


def iraf_phot(image_file, coo_file, dmin, dmax, gain, rdnoise, exp_time,
              aper_rad):
    """
    """
    try:
        os.remove('iraf_phot.mag')
    except OSError:
        pass

    iraf.daophot()
    iraf.datapars.datamin = dmin
    iraf.datapars.datamax = dmax
    iraf.datapars.epadu = gain
    iraf.datapars.readnoise = rdnoise
    iraf.datapars.itime = exp_time
    iraf.datapars.airmass = "AIRMASS"
    iraf.datapars.filter = "FILTER"

    iraf.centerpars.calgorithm = "centroid"
    iraf.centerpars.cbox = 8.

    iraf.fitskypars.annulus = aper_rad + 5.
    iraf.fitskypars.dannulus = 5.

    iraf.photpars.apertures = aper_rad
    iraf.photpars.zmag = 25.

    iraf.phot(
        image=image_file, coords=coo_file, output="iraf_phot.mag", verify="no",
        mode="h")

    # iraf.epar('phot')

    data = ascii.read('iraf_phot.mag', format='daophot')

    return data


def calibrate_magnitudes(tab, itime=1., zmag=25.):
    tab['cal_mags'] = (zmag - 2.5 * np.log10(tab['flux_fit'] / itime)) * u.mag
    return tab


def calc_aperture_mmm(data, mask, sigma_clip):
    """Helper function to actually calculate the stats for pixels
        falling within some Photutils aperture mask on some array
        of data.
    """
    cutout = mask.cutout(data, fill_value=np.nan)
    if cutout is None:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
    else:
        values = cutout * mask.data / mask.data
        values = values[~np.isnan(values)]
        if sigma_clip:
            values, clow, chigh = sigmaclip(values, low=3, high=3)

        mean = np.mean(values)
        median = np.median(values)
        std = np.std(values)

        mode = 3 * median - 2 * mean
        actual_area = (~np.isnan(values)).sum()
        return (mean, median, mode, std, actual_area)


def aperture_stats_tbl(data, apertures,
                       method='exact', sigma_clip=True):
    """Computes mean/median/mode/std in Photutils apertures.
    Compute statistics for custom local background methods.
    This is primarily intended for estimating backgrounds
    via annulus apertures.  The intent is that this falls easily
    into other code to provide background measurements.
    Parameters
    ----------
    data : array
        The data for the image to be measured.
    apertures : photutils PixelAperture object (or subclass)
        The phoutils aperture object to measure the stats in.
        i.e. the object returned via CirularAperture,
        CircularAnnulus, or RectangularAperture etc.
    method: str
        The method by which to handle the pixel overlap.
        Defaults to computing the exact area.
        NOTE: Currently, this will actually fully include a
        pixel where the aperture has ANY overlap, as a median
        is also being performed.  If the method is set to 'center'
        the pixels will only be included if the pixel's center
        falls within the aperture.
    sigma_clip: bool
        Flag to activate sigma clipping of background pixels
    Returns
    -------
    stats_tbl : astropy.table.Table
        An astropy Table with the colums X, Y, aperture_mean,
        aperture_median, aperture_mode, aperture_std, aperture_area
        and a row for each of the positions of the apertures.
    """

    # Get the masks that will be used to identify our desired pixels.
    masks = apertures.to_mask(method=method)

    # Compute the stats of pixels within the masks
    aperture_stats = [calc_aperture_mmm(data, mask, sigma_clip)
                      for mask in masks]

    aperture_stats = np.array(aperture_stats)

    # Place the array of the x y positions alongside the stats
    stacked = np.hstack([apertures.positions, aperture_stats])
    # Name the columns
    names = ['X','Y','aperture_mean','aperture_median','aperture_mode',
            'aperture_std', 'aperture_area']
    # Make the table
    stats_tbl = Table(data=stacked, names=names)


    return stats_tbl


def photutils_phot(coo_file, hdu_data, gain, exp_time, aper_rad):
    """
    """
    # Coordinates from observed frame.
    coo_t = Table.read(coo_file, format='ascii')
    positions = zip(*[coo_t['XCENTER'], coo_t['YCENTER']])
    apertures = CircularAperture(positions, r=aper_rad)
    annulus_apertures = CircularAnnulus(
        positions, r_in=aper_rad + 5, r_out=aper_rad + 10.)

    apers = [apertures, annulus_apertures]

    # from photutils import Background2D
    # # from astropy.stats import SigmaClip
    # # from photutils import MedianBackground
    # # sigma_clip = SigmaClip(sigma=3., iters=10)
    # # bkg_estimator = MedianBackground()

    # # Selecting the box size requires some care by the user. The box size
    # # should generally be larger than the typical size of sources in the
    # # image, but small enough to encapsulate any background variations. For
    # # best results, the box size should also be chosen so that the data are
    # # covered by an integer number of boxes in both dimensions.
    # bl = 10
    # box_xy = (bl, bl)
    # bkg = Background2D(hdu_data, box_xy)
    # print("background estimated")
    # error = calc_total_error(hdu_data, bkg.background, gain)

    # # from astropy.stats import biweight_midvariance
    # # sky_std = biweight_midvariance(hdu_data)
    # # sky_median = np.median(hdu_data)
    # # error = calc_total_error(hdu_data, sky_median, gain)

    error = None
    phot_table = aperture_photometry(hdu_data, apers, error=error)
    bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area()
    bkg_sum = bkg_mean * apertures.area()
    phot_table['flux_fit'] = phot_table['aperture_sum_0'] - bkg_sum
    phot_table = calibrate_magnitudes(phot_table, itime=exp_time)

    # phot_table['merr'] = 1.0857 *\
    #     phot_table['aperture_sum_err_0'] / phot_table['flux_fit']

    #
    # From https://github.com/astropy/photutils/issues/629#issuecomment-359595642
    error_array = None
    phot_table2 = aperture_photometry(hdu_data, apertures, error=error_array)
    bg_phot = aperture_stats_tbl(hdu_data, annulus_apertures, sigma_clip=False)

    ap_area = apertures.area()
    bg_method = 'mean'
    bg_method_name = 'aperture_{}'.format(bg_method)

    phot_table2['flux_fit'] = phot_table2['aperture_sum'] -\
        bg_phot[bg_method_name] * ap_area

    # Shouldn't 'bkg_mean' and 'bg_phot[bg_method_name]' be equivalent?
    # https://github.com/spacetelescope/wfc3_photometry/issues/2
    import pdb; pdb.set_trace()  # breakpoint e8e7e7d3 //


    # def compute_phot_error(
    #         flux_variance, bg_phot, bg_method, ap_area, epadu=1.0):
    #     """Computes the flux errors using the DAOPHOT style computation"""
    #     bg_variance_terms = (ap_area * bg_phot['aperture_std'] ** 2.) *\
    #         (1. + ap_area / bg_phot['aperture_area'])
    #     variance = flux_variance / epadu + bg_variance_terms
    #     flux_error = variance ** .5
    #     return flux_error

    # flux_error = compute_phot_error(
    #     phot_table['flux_fit'], bg_phot, 'median', ap_area, gain)

    # mag_err = 1.0857 * flux_error / phot_table['flux_fit']
    # mag_err[mag_err < 0] = np.nan
    # phot_table['merr'] = mag_err
    # phot_table = calibrate_magnitudes(phot_table, itime=exp_time)

    return phot_table


def write_data(iraf_data, photu_data):
    """
    """
    tab_comb = Table(iraf_data['XINIT', 'YINIT', 'FLUX', 'MAG', 'MERR'])
    t2 = Table(
        photu_data['xcenter', 'ycenter', 'flux_fit', 'cal_mags', 'merr'])
    tab_comb.add_columns(t2.columns.values())

    flux_diff = iraf_data['FLUX'] - photu_data['flux_fit']
    tab_comb['flux_diff'] = flux_diff
    mag_diff = iraf_data['MAG'] - photu_data['cal_mags']
    tab_comb['mag_diff'] = mag_diff

    tab_comb['MAG'].fill_value = -99.9
    tab_comb['mag_diff'].fill_value = -99.9
    tab_comb['flux_diff'].fill_value = -99.9
    tt = tab_comb.filled()

    tt.sort('mag_diff')
    tt.reverse()

    ascii.write(
        tt, 'iraf_photut_phot.dat', format='fixed_width', delimiter=' ',
        formats={
            'flux_fit': '%10.2f', 'cal_mags': '%10.3f', 'merr': '%10.3f',
            'flux_diff': '%10.3f', 'mag_diff': '%10.6f'}, overwrite=True)


def main():
    """
    Compare results from aperture photometry using IRAF's 'phot' and photutil's
    'CircularAperture'.
    """
    imname = 'field/filt_V/stk1111'
    image_file, coo_file, hdu_data, exp_time, gain, rdnoise = load_data(imname)

    aper_rad, dmax = 15., 60000.
    print("Aperture radius: {}".format(aper_rad))
    dmin = -3. * (rdnoise / gain)
    print("dmin = {:.2f}".format(dmin))

    iraf_data = iraf_phot(image_file, coo_file, dmin, dmax, gain, rdnoise,
                          exp_time, aper_rad)

    photu_data = photutils_phot(coo_file, hdu_data, gain, exp_time, aper_rad)

    # write_data(iraf_data, photu_data)

    #
    plt.subplot(211)
    plt.grid()
    plt.scatter(iraf_data['MAG'], iraf_data['MAG'] - photu_data['cal_mags'])
    plt.xlabel("IRAF mags")
    plt.ylabel("(IRAF - photutil) mags")
    plt.subplot(212)
    plt.grid()
    plt.scatter(
        iraf_data['MERR'], (iraf_data['MERR'] - photu_data['merr']))
    plt.xlabel("IRAF merr")
    plt.ylabel("(IRAF - photutil) merr")
    plt.show()


if __name__ == '__main__':
    main()
