
# This takes a range of numbers from the command line, constructs
# the filename by appending each number to the root "ca" and appending
# ".fits[0]", and splots the files one at a time.  It's intended as a
# little pyraf demo script.

import numpy as np
from pyraf import iraf

# This function is invoked when running from the shell.


def main():
    imname = '/home/gabriel/Github/photom/standards/filt_U/stk_2085.fits'
    dmax = 60000.

    gain = iraf.hselect(images=imname, fields='EGAIN', expr='yes', mode='hl')
    rdnoise = iraf.hselect(images=imname, fields='ENOISE', expr='yes',
                           mode='hl')
    import pdb; pdb.set_trace()  # breakpoint 1dec3c10 //
    
    # ira.hselect()

    # Initial FWHM value
    fitrad = 3.
    # Initial SKY MEAN value
    smean = 1.
    # Initial SKY STANDARD DEVIATION value
    sstd = 1.
    sigma = np.sqrt(smean * gain + rdnoise * rdnoise) / gain

    iraf.datapars.fwhmpsf = fitrad
    iraf.datapars.sigma = sstd
    iraf.datapars.datamin = smean - 3 * sigma
    iraf.datapars.datamax = dmax

    # Use a high threshold so only the brighter stars will be found
    i = 5
    iraf.findpars.threshold = i * 3.5 * sigma
    print('\n---------------------------------------------------')
    print(" Daofind task (i): {}".format(i))
    print(" Threshold (i*3.5*sigma): {}".format(i * 3.5 * sigma))
    print('---------------------------------------------------')

    iraf.daofind.verif = 'no'
    iraf.daofind.verb = 'yes'
    iraf.daofind.interactive = 'no'
    iraf.daofind.verbose = 'no'
    iraf.daofind.mode = 'hl'
    iraf.daofind(imname, (imname + '.coo.psf.1'))


if __name__ == "__main__":
    main()
