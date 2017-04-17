
import os
from astropy.io import ascii
from pyraf import iraf


def main(dmax, psf_select, imname, hdu_data):
    """
    Use the IRAF task 'psfmeasure' to estimate the FWHM and ellipticity of
    all the stars in the 'psf_select' list.
    """
    print("\nRun 'psfmeasure' task to estimate the FWHMs.")
    try:
        os.remove('positions')
        os.remove('cursor')
        os.remove('psfmeasure')
    except OSError:
        pass

    print("Total number of analyzed stars: {}".format(len(psf_select)))
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
    psf_data = psf_data['Column', 'Line', 'FWHM', 'Ellip', 'Mag']
    print("Stars rejected by 'psfmeasure': {}".format(
        len(psf_select) - len(psf_data)))

    os.remove('positions')
    os.remove('cursor')
    os.remove('psfmeasure')

    return psf_data


if __name__ == "__main__":
    main()
