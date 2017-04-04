
import os
from astropy.io import ascii
import numpy as np
from pyraf import iraf


def main(dmax, ellip_max, fwhm_min, psf_select, imname, hdu_data):
    """
    Use the IRAF task 'psfmeasure' to estimate the FWHM and ellipticity of
    all the stars in the 'psf_select' list.
    Reject those with a large ellipticity (> ellip_max), and a very low
    FWHM (<fwhm_min).
    """
    print("\nRun 'psfmeasure' task to estimate the FWHMs.")
    try:
        os.remove('positions')
        os.remove('cursor')
        os.remove('psfmeasure')
    except OSError:
        pass

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

    # Read PSFMEASURE task output, leaving out the last line with the average
    # FWHM.
    psf_data = ascii.read(
        "psfmeasure", format='fixed_width', header_start=1, data_end=-1,
        col_starts=(15, 23, 32, 40, 48, 56))
    psf_data = psf_data['Column', 'Line', 'FWHM', 'Ellip']

    # Extract data
    fwhm_min_rjct = psf_data[psf_data['FWHM'] <= fwhm_min]
    fwhm_min_accpt = psf_data[psf_data['FWHM'] > fwhm_min]
    fwhm_estim = fwhm_min_accpt[fwhm_min_accpt['Ellip'] <= ellip_max]
    ellip_rjct = fwhm_min_accpt[fwhm_min_accpt['Ellip'] > ellip_max]

    # Remove duplicates, if any.
    fwhm_estim = list(set(np.array(fwhm_estim).tolist()))

    os.remove('positions')
    os.remove('cursor')
    os.remove('psfmeasure')

    return fwhm_estim, fwhm_min_rjct, ellip_rjct


if __name__ == "__main__":
    main()
