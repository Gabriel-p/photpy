
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import ascii
from astropy.table import Table


def main():
    """
    Filter photometry for color ranges given. Also remove stars with all
    'nan' values in their photometry.
    """
    f_id = 'rup152_final.dat'

    # Define acceptable color ranges for this data.
    BV_min, BV_max = 0., 2.
    VI_min, VI_max = 0.5, 2.
    UB_min, UB_max = -.5, 1.5

    # Load final photometry file.
    phot = photLoad(f_id)
    print("Stars in file: {}".format(len(phot)))

    plotCMDs(f_id, phot, 'all', BV_min, BV_max, UB_min, UB_max, VI_min, VI_max)

    phot, phot_rjct = filterPhot(
        phot, BV_min, BV_max, VI_min, VI_max, UB_min, UB_max)
    print("Stars in cleaned file: {}".format(len(phot)))

    plotCMDs(
        f_id, phot_rjct, 'rjct', BV_min, BV_max, UB_min, UB_max, VI_min,
        VI_max)
    plotCMDs(
        f_id, phot, 'accpt', BV_min, BV_max, UB_min, UB_max, VI_min, VI_max)
    print("Plots created")

    fileClean(f_id, phot)


def photLoad(f_id):
    """
    """
    phot = ascii.read('input/' + f_id, fill_values=('INDEF', np.nan))
    phot = Table(
        phot, names=('id', 'x', 'y', 'V', 'eV', 'BV', 'eBV', 'UB', 'eUB', 'VI',
                     'eVI'))
    # To bring closer to APASS
    phot['V'] = phot['V'] + 0.019
    phot['BV'] = phot['BV'] + 0.049
    phot['UB'] = phot['UB'] - 0.068
    phot['VI'] = phot['VI'] + 0.019

    return phot


def filterPhot(phot, BV_min, BV_max, VI_min, VI_max, UB_min, UB_max):
    """
    """
    # Filter colors by range, accepting 'nan' values.
    m1 = np.logical_or(phot['BV'] < BV_max, np.isnan(phot['BV']))
    m2 = np.logical_or(phot['BV'] > BV_min, np.isnan(phot['BV']))
    m3 = np.logical_or(phot['VI'] < VI_max, np.isnan(phot['VI']))
    m4 = np.logical_or(phot['VI'] > VI_min, np.isnan(phot['VI']))
    m5 = np.logical_or(phot['UB'] < UB_max, np.isnan(phot['UB']))
    m6 = np.logical_or(phot['UB'] > UB_min, np.isnan(phot['UB']))

    mask = [m1, m2, m3, m4, m5, m6]
    total_mask = reduce(np.logical_and, mask)

    # Save stars outside color ranges.
    phot_rjct = phot[~total_mask]
    # Save stars within color ranges.
    phot = phot[total_mask]

    # Remove stars with 'nan' values in *all* colors and magnitude.
    nan_msk = [
        ~phot['V'].mask, ~phot['BV'].mask, ~phot['VI'].mask, ~phot['UB'].mask]
    total_mask = reduce(np.logical_or, nan_msk)
    phot = phot[total_mask]

    return phot, phot_rjct


def plotCMDs(
    f_id, phot, acpt_rjct_ID, BV_min, BV_max, UB_min, UB_max, VI_min, VI_max):
    """
    Plot photometry diagrams.
    """
    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(24, 8))
    gs = gridspec.GridSpec(1, 3)

    fig.add_subplot(gs[0])
    plt.title("N={}".format(np.count_nonzero(~np.isnan(phot['BV']))))
    plt.xlabel("(B-V)")
    plt.ylabel("V")
    plt.scatter(phot['BV'], phot['V'], s=5)
    plt.axvline(x=BV_max, c='r')
    plt.axvline(x=BV_min, c='r')
    plt.gca().invert_yaxis()

    fig.add_subplot(gs[1])
    plt.title("N={}".format(np.count_nonzero(~np.isnan(phot['VI']))))
    plt.xlabel("(V-I)")
    plt.ylabel("V")
    plt.scatter(phot['VI'], phot['V'], s=5)
    plt.axvline(x=VI_max, c='r')
    plt.axvline(x=VI_min, c='r')
    plt.gca().invert_yaxis()

    fig.add_subplot(gs[2])
    plt.title("N={}".format(np.count_nonzero(~np.isnan(phot['UB']))))
    plt.xlabel("(B-V)")
    plt.ylabel("U-B")
    plt.scatter(phot['BV'], phot['UB'], s=5)
    plt.axhline(y=UB_max, c='r')
    plt.axhline(y=UB_min, c='r')
    plt.gca().invert_yaxis()

    fig.tight_layout()
    plt.savefig(
        'output/' + f_id.split('.')[0] + '_' + acpt_rjct_ID + '.png', dpi=300,
        bbox_inches='tight')


def fileClean(f_id, phot):
    """
    Create clean photometry file.
    """
    ascii.write(
        phot, 'output/' + f_id, format='fixed_width_no_header',
        delimiter=' ', fill_values=[(ascii.masked, '99.999')], overwrite=True)


if __name__ == '__main__':
    main()
