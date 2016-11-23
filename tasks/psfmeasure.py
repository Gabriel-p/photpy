
import os
from astropy.io import ascii
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

    fwhm_estim, fwhm_min_rjct, ellip_rjct, psfmeasure_estim = [], [], [], 0
    with open("psfmeasure", 'r') as f:
        for i, line in enumerate(f):
            data = line.split()
            if data:
                if data[0] != 'Average':
                    # First line (star) of output file.
                    if i == 3:
                        if float(data[4]) > fwhm_min:
                            if float(data[5]) <= ellip_max:
                                fwhm_estim.append(
                                    map(float, [data[1], data[2], data[4],
                                                data[5]]))
                            else:
                                ellip_rjct.append(
                                    map(float, [data[1], data[2], data[4],
                                                data[5]]))
                        else:
                            fwhm_min_rjct.append(
                                map(float, [data[1], data[2], data[4],
                                            data[5]]))
                    # Rest of the lines.
                    elif i > 3:
                        if float(data[3]) > fwhm_min:
                            if float(data[4]) <= ellip_max:
                                fwhm_estim.append(
                                    map(float, [data[0], data[1], data[3],
                                                data[4]]))
                            else:
                                ellip_rjct.append(
                                    map(float, [data[0], data[1], data[3],
                                                data[4]]))
                            # sys.stdout = open(os.devnull, "w")
                            # gauss_x = plots.line_fit(
                            #     float(data[0]), float(data[1]), genplot=False)
                            # gauss_y = plots.column_fit(
                            #     float(data[0]), float(data[1]), genplot=False)
                            # sys.stdout = sys.__stdout__
                            # print(float(data[3]), float(data[4]),
                            #       gfwhm(gauss_x.stddev)[0],
                            #       gfwhm(gauss_y.stddev)[0])
                        else:
                            fwhm_min_rjct.append(
                                map(float, [data[0], data[1], data[3],
                                            data[4]]))
                else:
                    # Averaged FWHM by the 'psfmeasure' task.
                    psfmeasure_estim = float(data[-1])

    os.remove('positions')
    os.remove('cursor')
    os.remove('psfmeasure')

    return fwhm_estim, psfmeasure_estim, fwhm_min_rjct, ellip_rjct


if __name__ == "__main__":
    main()
