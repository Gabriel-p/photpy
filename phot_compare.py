
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


def load_data(imname):
    """
    """
    mypath = realpath(join(os.getcwd(), dirname(__file__)))
    image_file = mypath + '/input/' + imname + '.fits'
    coo_file = mypath + '/output/' + imname + '.coo'
    print(image_file)
    # Load .fits file.
    hdulist = fits.open(image_file)
    hdu_data, hdr = hdulist[0].data, hdulist[0].header
    exp_time, gain, rdnoise = hdr['EXPTIME'], hdr['EGAIN'], hdr['ENOISE']

    return image_file, coo_file, hdu_data, exp_time, gain, rdnoise


def iraf_phot(image_file, coo_file, dmin, dmax, gain, rdnoise, exp_time,
              aper_rad):
    """
    """
    try:
        os.remove('iraf_phot')
    except OSError:
        pass
    iraf.datapars(
        datamin=dmin, datamax=dmax, epadu=gain, readnoise=rdnoise,
        itime=exp_time, airmass="AIRMASS", filter="FILTER",
        mode="h")
    iraf.centerpars(calgorithm="centroid", cbox=8., mode="h")
    iraf.fitskypars(annulus=aper_rad + 5., dannulus=5., mode="h")
    iraf.photpars(apertures=aper_rad, zmag=25., mode="h")
    iraf.daophot()
    iraf.phot(
        image=image_file, coords=coo_file, output="iraf_phot", verify="no",
        mode="h")

    data = ascii.read('iraf_phot', format='daophot')
    os.remove('iraf_phot')

    return data


def calibrate_magnitudes(tab, itime=1., zmag=25.):
    tab['cal_mags'] = (zmag - 2.5 * np.log10(tab['flux_fit'] / itime)) * u.mag
    return tab


def photutils_phot(coo_file, hdu_data, exp_time, aper_rad):
    """
    """
    # Coordinates from observed frame.
    coo_t = Table.read(coo_file, format='ascii')
    positions = zip(*[coo_t['x'], coo_t['y']])
    apertures = CircularAperture(positions, r=aper_rad)
    annulus_apertures = CircularAnnulus(
        positions, r_in=aper_rad + 5, r_out=aper_rad + 10.)

    apers = [apertures, annulus_apertures]
    phot_table = aperture_photometry(hdu_data, apers)
    bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area()
    bkg_sum = bkg_mean * apertures.area()
    phot_table['flux_fit'] = phot_table['aperture_sum_0'] - bkg_sum
    phot_table = calibrate_magnitudes(phot_table, itime=exp_time)

    return phot_table


def write_data(iraf_data, photu_data):
    """
    """
    tab_comb = Table(iraf_data['XINIT', 'YINIT', 'FLUX', 'MAG'])
    t2 = Table(photu_data['xcenter', 'ycenter', 'flux_fit', 'cal_mags'])
    tab_comb.add_columns(t2.columns.values())

    mag_diff = iraf_data['MAG'] - photu_data['cal_mags']
    tab_comb['mag_diff'] = mag_diff

    tab_comb['MAG'].fill_value = -99.9
    tab_comb['mag_diff'].fill_value = -99.9
    tt = tab_comb.filled()

    tt.sort('mag_diff')
    tt.reverse()

    ascii.write(tt, 'iraf_photut_phot.dat', overwrite=True,
                format='fixed_width', delimiter=' ')


def main():
    """
    Compare results from aperture photometry using IRAF's 'phot' and photutil's
    'CircularAperture'.
    """
    imname = 'standards/stk_2150'
    image_file, coo_file, hdu_data, exp_time, gain, rdnoise = load_data(imname)

    aper_rad, dmax = 15., 60000.
    dmin = -3. * (rdnoise / gain)
    print("dmin = {:.2f}".format(dmin))

    iraf_data = iraf_phot(image_file, coo_file, dmin, dmax, gain, rdnoise,
                          exp_time, aper_rad)

    photu_data = photutils_phot(coo_file, hdu_data, exp_time, aper_rad)

    write_data(iraf_data, photu_data)

    #
    plt.grid()
    plt.scatter(iraf_data['MAG'], iraf_data['MAG'] - photu_data['cal_mags'])
    plt.show()


if __name__ == '__main__':
    main()
