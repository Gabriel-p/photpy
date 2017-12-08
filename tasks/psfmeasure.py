
import os
from pyraf import iraf


def main(dmax, imname):
    """
    Use the IRAF task 'psfmeasure' to estimate the FWHM and ellipticity of
    all the stars in the 'psf_select' list.
    """
    dmax = float(dmax)

    print("\nRun 'psfmeasure' task to estimate the FWHMs.")
    try:
        os.remove('cursor')
        os.remove('psfmeasure')
    except OSError:
        pass

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

    os.remove('positions')
    os.remove('cursor')

    return


if __name__ == "__main__":
    import sys
    dmax = sys.argv[1]
    imname = sys.argv[2]

    main(dmax, imname)
