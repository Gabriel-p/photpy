
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Column


def makeHist(x, y):
    # Estimate the 2D histogram
    nbins = 25
    H, xedges, yedges = np.histogram2d(x, y, bins=nbins)
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    # Mask pixels with a value of zero
    Hmasked = np.ma.masked_where(H == 0, H)
    plt.pcolormesh(xedges, yedges, Hmasked, cmap='viridis')


def magAnalysis(f_id):
    """
    Analyze .mag files
    """
    # Read DAOMASTER .mag data.
    print("U mag")
    ufilt = ascii.read('input/' + f_id + '/ufilter.mag')
    xu, yu = ufilt['col2'], ufilt['col3']
    print("B mag")
    bfilt = ascii.read('input/' + f_id + '/bfilter.mag')
    xb, yb = bfilt['col2'], bfilt['col3']
    print("V mag")
    vfilt = ascii.read('input/' + f_id + '/vfilter.mag')
    xv, yv = vfilt['col2'], vfilt['col3']
    print("I mag")
    ifilt = ascii.read('input/' + f_id + '/ifilter.mag')
    xi, yi = ifilt['col2'], ifilt['col3']

    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(10, 10))

    # xy density histograms of DAOMASTER's .mag matched filter files.
    plt.subplot(221)
    plt.title("Filter U")
    makeHist(xu, yu)
    plt.subplot(222)
    plt.title("Filter B")
    makeHist(xb, yb)
    plt.subplot(223)
    plt.title("Filter V")
    makeHist(xv, yv)
    plt.subplot(224)
    plt.title("Filter I")
    makeHist(xi, yi)

    fig.tight_layout()
    plt.savefig(
        'output/' + f_id + '_filters_dens.png', dpi=300, bbox_inches='tight')


def createObsfile(f_id, AU, AB, AV, AI):
    """
    Create 'daom.obs' file from 'daom.raw' file and longest exposure airmasses.
    """
    data = ascii.read('input/' + f_id + '/daom.raw')
    XV = Column(data=[AV] * len(data), name='v')
    XB = Column(data=[AB] * len(data), name='b')
    XU = Column(data=[AU] * len(data), name='u')
    XI = Column(data=[AI] * len(data), name='i')
    data.add_columns([XV, XB, XU, XI])
    ascii.write(
        data, 'output/daom_' + f_id + '.obs', format='fixed_width_no_header',
        delimiter=' ', overwrite=True)


def plotRaw(f_id):
    """
    Plot 'daom.raw' photometry diagrams.
    """
    data = ascii.read(
        'input/' + f_id + '/daom.raw',
        fill_values=[('INDEF', np.nan), ('99.9999', np.nan)])
    v, b, u, i = data['col4'], data['col6'], data['col8'], data['col10']
    v, b, u, i = np.array(v), np.array(b), np.array(u), np.array(i)

    N_vb = np.count_nonzero(~np.isnan(b - v))
    N_vi = np.count_nonzero(~np.isnan(v - i))
    N_ub = np.count_nonzero(~np.isnan(u - b))
    print("Non zero BV: {}".format(N_vb))
    print("Non zero VI: {}".format(N_vi))
    print("Non zero UB: {}".format(N_ub))

    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(30, 10))

    plt.subplot(131)
    plt.xlabel("(B-V)")
    plt.ylabel("V")
    plt.scatter(b - v, v, s=5)
    plt.gca().invert_yaxis()
    plt.title("N={}".format(N_vb))

    plt.subplot(132)
    plt.xlabel("(V-I)")
    plt.ylabel("V")
    plt.scatter(v - i, v, s=5)
    plt.gca().invert_yaxis()
    plt.title("N={}".format(N_vi))

    plt.subplot(133)
    plt.xlabel("(B-V)")
    plt.ylabel("U-B")
    plt.scatter(b - v, u - b, s=5)
    plt.gca().invert_yaxis()
    plt.title("N={}".format(N_ub))

    # plt.show()
    fig.tight_layout()
    plt.savefig('output/' + f_id + '_daom.png', dpi=300, bbox_inches='tight')


def main():
    """
    """
    # ID for this observation.
    f_id = 'bh73'
    magAnalysis(f_id)

    # Airmasses of longest exposures
    # Rup42
    # AU, AB, AV, AI = 1.003, 1.002, 1.057, 1.01
    # Haffner 14
    # AU, AB, AV, AI = 1.006, 1.001, 1.001, 1.022
    # Ruprecht 152
    # AU, AB, AV, AI = 1.019, 1.03, 1.013, 1.015
    # Ruprecht 41
    # AU, AB, AV, AI = 1.001, 1.008, 1.002, 1.008
    # vdB-Hagen 73 (longexp)
    # AU, AB, AV, AI = 1.092, 1.072, 1.075, 1.164
    # vdB-Hagen 73 (large FWHM)
    AU, AB, AV, AI = 1.201, 1.213, 1.192, 1.164

    createObsfile(f_id, AU, AB, AV, AI)

    plotRaw(f_id)


if __name__ == '__main__':
    main()
