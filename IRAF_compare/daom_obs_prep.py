
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
# from astropy.table import Column


def makeHist(x, y):
    # Estimate the 2D histogram
    nbins = 25
    H, xedges, yedges = np.histogram2d(x, y, bins=nbins)
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    # Mask pixels with a value of zero
    Hmasked = np.ma.masked_where(H == 0, H)
    plt.pcolormesh(xedges, yedges, Hmasked)  # , vmin=0, vmax=150)


# xy density histograms of DAOMASTER's .mag matched filter files.
ufilt = ascii.read('3_DAOM/ufilter.mag')
xu, yu = ufilt['col2'], ufilt['col3']
bfilt = ascii.read('3_DAOM/bfilter.mag')
xb, yb = bfilt['col2'], bfilt['col3']
vfilt = ascii.read('3_DAOM/vfilter.mag')
xv, yv = vfilt['col2'], vfilt['col3']
ifilt = ascii.read('3_DAOM/ifilter.mag')
xi, yi = ifilt['col2'], ifilt['col3']

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
plt.show()

WT = '4'
print("WT={}".format(WT))

# Create 'daom.obs' file.
# data = ascii.read('daom_' + WT + '.raw')
# XV = Column(data=[1.075] * len(data), name='v')
# XB = Column(data=[1.072] * len(data), name='b')
# XU = Column(data=[1.092] * len(data), name='u')
# XI = Column(data=[1.164] * len(data), name='i')
# data.add_columns([XV, XB, XU, XI])
# ascii.write(data, 'daom.obs', format='fixed_width', delimiter=' ',
#             overwrite=True)

# Explore 'daom.raw' photometry
data = ascii.read('3_DAOM/daom_' + WT + '.raw', fill_values=('INDEF', np.nan))
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
plt.xlim(-2., 4.)
plt.ylim(25., 10.)
plt.scatter(b - v, v, s=5)
plt.title("N={}".format(N_vb))

plt.subplot(132)
plt.xlim(0., 4.)
plt.ylim(25., 10.)
plt.scatter(v - i, v, s=5)
plt.title("N={}".format(N_vi))

plt.subplot(133)
plt.xlim(.5, 2.5)
plt.ylim(4., 2.)
plt.scatter(b - v, u - b, s=5)
plt.title("N={}".format(N_ub))

# plt.show()
fig.tight_layout()
plt.savefig('output/daomaster_wt_' + WT + '.png', dpi=300, bbox_inches='tight')
