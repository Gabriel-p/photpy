
import os
import sys
from os.path import join, realpath, dirname
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval
from photutils import CircularAperture
from scipy.spatial import distance
from getdata import st_fwhm_select
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import psfmeasure


def create_pars_file(pars_f, pars_list=None):
    """
    """
    # Default values.
    if pars_list is None:
        pars_list = ['None', 'True', '5.', '3.', '60000.', '0.15', '1.5',
                     '100', '0.', '0.', '-1.', '0.05', 'EGAIN', 'ENOISE',
                     'FILTER', 'EXPTIME']
    with open(pars_f, 'w') as f:
        f.write(
            "# Default parameters for the align_crop.py script\n#\n"
            "ff_proc {}\ndo_plots {}\nthresh_level {}\nfwhm_init {}\ndmax {}\n"
            "ellip_max {}\nfwhm_min {}\nmax_shift_stars {}\nx_init_shift {}\n"
            "y_init_shift {}\nmax_shift {}\ntol {}\ngain_key {}\n"
            "rdnoise_key {}\nfilter_key {}\nexp_key {}\n".format(*pars_list))
    return


def read_params():
    """
    """
    pars = {}
    mypath = realpath(join(os.getcwd(), dirname(__file__)))
    pars_f = join(mypath, '.align_crop.pars')
    if not os.path.isfile(pars_f):
        print("Parameters file missing. Create it.")
        create_pars_file(pars_f)

    with open(pars_f, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                key, value = line.replace('\n', '').split(' ')
                pars[key] = value

    return mypath, pars_f, pars


def get_params(mypath, pars_f, pars):
    """
    """
    pars_list = []

    answ = raw_input("\nDir or file to process ({}): ".format(
        pars['ff_proc']))
    fname = str(answ) if answ is not '' else pars['ff_proc']
    fname = fname[1:] if fname.startswith('/') else fname
    r_path = join(mypath.replace('/tasks', ''), fname)
    fits_list = []
    if os.path.isdir(r_path):
        for subdir, dirs, files in os.walk(r_path):
            for file in files:
                if file.endswith('.fits'):
                    fits_list.append(os.path.join(subdir, file))
    elif os.path.isfile(r_path):
        print("  Single .fits file found in folder.")
        sys.exit()
    else:
        print("{}\nis neither a folder nor a file. Exit.".format(str(r_path)))
        sys.exit()
    pars['ff_proc'] = fname
    pars_list.append(pars['ff_proc'])

    answ = raw_input("Create plots? (y/n) ({}): ".format(
        pars['do_plots']))
    pars['do_plots'] = False if answ in ('n', 'N', 'no', 'NO') else True
    pars_list.append(pars['do_plots'])

    answ = raw_input("Threshold level above STDDEV ({}): ".format(
        pars['thresh_level']))
    pars['thresh_level'] = float(answ) if answ is not '' else\
        float(pars['thresh_level'])
    pars_list.append(pars['thresh_level'])

    answ = raw_input("Initial FWHM ({}): ".format(pars['fwhm_init']))
    pars['fwhm_init'] = float(answ) if answ is not '' else\
        float(pars['fwhm_init'])
    pars_list.append(pars['fwhm_init'])

    answ = raw_input("Maximum flux ({}): ".format(pars['dmax']))
    pars['dmax'] = float(answ) if answ is not '' else float(pars['dmax'])
    pars_list.append(pars['dmax'])

    answ = raw_input("Maximum ellipticity ({}): ".format(pars['ellip_max']))
    pars['ellip_max'] = float(answ) if answ is not '' else\
        float(pars['ellip_max'])
    pars_list.append(pars['ellip_max'])

    answ = raw_input("Minimum FWHM ({}): ".format(pars['fwhm_min']))
    pars['fwhm_min'] = float(answ) if answ is not '' else\
        float(pars['fwhm_min'])
    pars_list.append(pars['fwhm_min'])

    answ = raw_input("Max number of stars used to obtain "
                     "shift ({}): ".format(pars['max_shift_stars']))
    pars['max_shift_stars'] = int(answ) if answ is not '' else\
        int(pars['max_shift_stars'])
    pars_list.append(pars['max_shift_stars'])

    answ = raw_input("Initial estimated shift in x ({}): ".format(
        pars['x_init_shift']))
    pars['x_init_shift'] = float(answ) if answ is not '' else\
        float(pars['x_init_shift'])
    pars_list.append(pars['x_init_shift'])

    answ = raw_input("Initial estimated shift in y ({}): ".format(
        pars['y_init_shift']))
    pars['y_init_shift'] = float(answ) if answ is not '' else\
        float(pars['y_init_shift'])
    pars_list.append(pars['y_init_shift'])

    answ = raw_input("Maximum estimated shift ({}): ".format(
        pars['max_shift']))
    pars['max_shift'] = float(answ) if answ is not '' else\
        float(pars['max_shift'])
    pars_list.append(pars['max_shift'])

    answ = raw_input("Matching tolerance ({}): ".format(pars['tol']))
    pars['tol'] = float(answ) if answ is not '' else float(pars['tol'])
    pars_list.append(pars['tol'])

    answ = raw_input("GAIN keyword ({}): ".format(pars['gain_key']))
    pars['gain_key'] = str(answ) if answ is not '' else pars['gain_key']
    pars_list.append(pars['gain_key'])

    answ = raw_input("RDNOISE keyword ({}): ".format(pars['rdnoise_key']))
    pars['rdnoise_key'] = str(answ) if answ is not '' else pars['rdnoise_key']
    pars_list.append(pars['rdnoise_key'])

    answ = raw_input("FILTER keyword ({}): ".format(pars['filter_key']))
    pars['filter_key'] = str(answ) if answ is not '' else pars['filter_key']
    pars_list.append(pars['filter_key'])

    answ = raw_input("EXPOSURE keyword ({}): ".format(pars['exp_key']))
    pars['exp_key'] = str(answ) if answ is not '' else pars['exp_key']
    pars_list.append(pars['exp_key'])

    # Write values to file.
    create_pars_file(pars_f, pars_list)

    return r_path, fits_list, pars


def avrg_dist(init_shift, max_shift, tol, ref, f):
    """
    Average minimal (Euclidean) distance from points in 'f' to points in
    'ref', until the stopping condition.

    Source: http://codereview.stackexchange.com/a/134918/35351
    """
    # Invert.
    x0, y0 = -1. * init_shift[0], -1. * init_shift[1]
    x0, y0 = init_shift[0], init_shift[1]
    if max_shift < 0.:
        # Use the full length of both sides.
        max_shift = max(max(zip(*ref)[0]), max(zip(*ref)[1]))

    while max_shift > 1.:
        print(x0, y0, max_shift)
        dists = []
        li = 10
        # Shift in x
        for sx in np.linspace(-1. * max_shift, max_shift, li):
            # Shift in y
            for sy in np.linspace(-1. * max_shift, max_shift, li):
                # Apply possible x,y translation
                sf = [zip(*f)[0] + sx + x0, zip(*f)[1] + sy + y0]
                # Average minimal distance.
                d = np.mean(np.min(distance.cdist(zip(*sf), ref), axis=1))
                # Store x,y shift values, and the average minimal distance.
                dists.append([sx + x0, sy + y0, d])

        # Index of the x,y shifts that resulted in the average minimal
        # distance.
        min_idx = np.array(zip(*dists)[2]).argmin()
        x0, y0 = dists[min_idx][0], dists[min_idx][1]
        # Decrease max shift by X%
        max_shift = (1. - tol) * max_shift
        # print(dists[min_idx], max_shift)
        # Break condition: minimum accuracy reached.
        if max_shift < 1.:
            if dists[min_idx][2] > 5.:
                print("  WARNING: match is probably wrong.")

    return dists[min_idx]


def make_plots(hdu, f_list, fig_name):
    """
    Make plots.
    """
    print("\nPlotting.")
    fig = plt.figure(figsize=(20, 15))
    gs = gridspec.GridSpec(10, 12)

    ax = plt.subplot(gs[0:5, 0:5])
    interval = ZScaleInterval()
    zmin, zmax = interval.get_limits(hdu[0])
    plt.imshow(hdu[0], cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower', vmin=zmin, vmax=zmax)
    positions = (zip(*f_list[0][1])[0], zip(*f_list[0][1])[1])
    apertures = CircularAperture(positions, r=10.)
    apertures.plot(color='red', lw=1.)
    ax.tick_params(axis='both', which='major', labelsize=8)
    # plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)

    ax = plt.subplot(gs[0:5, 5:10])
    interval = ZScaleInterval()
    zmin, zmax = interval.get_limits(hdu[1])
    plt.imshow(hdu[1], cmap='viridis', aspect=1, interpolation='nearest',
               origin='lower', vmin=zmin, vmax=zmax)
    positions = (zip(*f_list[1][1])[0], zip(*f_list[1][1])[1])
    apertures = CircularAperture(positions, r=10.)
    apertures.plot(color='red', lw=1.)
    ax.tick_params(axis='both', which='major', labelsize=8)

    fig.tight_layout()
    plt.savefig(fig_name + '.png', dpi=150, bbox_inches='tight')
    # Close to release memory.
    plt.clf()
    plt.close()


def main():
    """
    """

    mypath, pars_f, pars = read_params()
    r_path, fits_list, pars = get_params(mypath, pars_f, pars)

    f_list, l_f = [], []
    hdu = []
    # For each .fits image in the root folder.
    for i, imname in enumerate(fits_list):
        print("\nFile: {}".format(imname.replace(mypath, "")))

        # Extract data.
        hdulist = fits.open(imname)
        hdu_data = hdulist[0].data
        hdu.append(hdu_data)

        # Background estimation.
        print("\nEstimating background mean, median, and STDDEV.")
        sky_mean, sky_median, sky_std = sigma_clipped_stats(
            hdu_data, sigma=3.0, iters=2)

        # Stars selection.
        psf_select, all_sources, n_not_satur = st_fwhm_select(
            pars['dmax'], pars['max_shift_stars'], pars['thresh_level'],
            pars['fwhm_init'], sky_std, hdu_data)

        # FWHM selection.
        fwhm_estim, psfmeasure_estim, fwhm_min_rjct, ellip_rjct =\
            psfmeasure.main(
                pars['dmax'], pars['ellip_max'], pars['fwhm_min'],
                psf_select, imname, hdu_data)

        l_f.append(all_sources)
        f_list.append(
            [imname, zip(*[zip(*fwhm_estim)[0], zip(*fwhm_estim)[1]])])

    ref_i = l_f.index(max(l_f))
    ref = f_list[ref_i][1]
    print('Reference image: {}'.format(f_list[ref_i][0]))
    # For each x,y coords set.
    for i, imname in enumerate(f_list):
        if i != ref_i:
            print("\nFile: {}".format(imname[0].replace(mypath, "")))
            # Obtain shifts
            f = imname[1]
            d = avrg_dist((pars['x_init_shift'], pars['y_init_shift']),
                          pars['max_shift'], pars['tol'], ref, f)
            print(d)

    fig_name = "test"
    if pars['do_plots']:
        make_plots(hdu, f_list, fig_name)


if __name__ == "__main__":
    main()
