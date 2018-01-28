
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import ascii


# Load IRAF final photometry load

# phot = ascii.read(
#     'input/phot_compare/BH73_IRAF_2coeffs.dat', fill_values=('INDEF', np.nan))
# id_i, x_i, y_i, v_i, ev_i, bv_i, ebv_i, ub_i, eub_i, vi_i, evi_i =\
#     phot['ID'], phot['x'], phot['y'], phot['V'], phot['eV'],\
#     phot['BV'], phot['eBV'], phot['UB'], phot['eUB'],\
#     phot['VI'], phot['eVI']

# phot = ascii.read(
#     'input/phot_compare/BH73_IRAF_4coeffs.dat', fill_values=('INDEF', np.nan))
# # id_i, x_i, y_i, v_i, ev_i, bv_i, ebv_i, ub_i, eub_i, vi_i, evi_i
# id_i, x_i, y_i, v_i, ev_i, bv_i, ebv_i, ub_i, eub_i, vi_i, evi_i =\
#     phot['ID'], phot['x'], phot['y'], phot['V'], phot['eV'],\
#     phot['BV'], phot['eBV'], phot['UB'], phot['eUB'],\
#     phot['VI'], phot['eVI']

phot = ascii.read(
    'input/phot_compare/BH73_IRAF_4_2_coeffs.dat',
    fill_values=('INDEF', np.nan))
id_i, x_i, y_i, v_i, ev_i, bv_i, ebv_i, ub_i, eub_i, vi_i, evi_i =\
    phot['ID'], phot['x'], phot['y'], phot['V'], phot['eV'],\
    phot['BV'], phot['eBV'], phot['UB'], phot['eUB'],\
    phot['VI'], phot['eVI']

# Load Python final photometry

phot = ascii.read(
    'input/phot_compare/BH73_phot.dat', fill_values=('INDEF', np.nan))
id_p, x_p, y_p, v_p, ev_p, bv_p, ebv_p, vi_p, evi_p, ub_p, eub_p =\
    phot['ID'], phot['x'], phot['y'], phot['V'], phot['eV'],\
    phot['BV'], phot['eBV'], phot['VI'], phot['eVI'], phot['UB'],\
    phot['eUB']

print("Data read")

v_i, bv_i, vi_i, ub_i, v_p, bv_p, vi_p, ub_p =\
    [np.array(_) for _ in [v_i, bv_i, vi_i, ub_i, v_p, bv_p, vi_p, ub_p]]

plt.style.use('seaborn-darkgrid')
fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(4, 1)

ax = fig.add_subplot(gs[0])
med = np.nanmedian(v_i - v_p)
plt.title("V ({:.3f})".format(med))
plt.ylim(-.05, .05)
plt.xlabel(r"$V_{IRAF}$")
plt.ylabel(r"$V_{IRAF}-V_{Python}$")
plt.scatter(v_i, v_i - v_p, s=4)
plt.axhline(y=med, c='r')

ax = fig.add_subplot(gs[1])
# plt.ylim(-.05, .05)
plt.xlabel(r"$V_{IRAF}$")
plt.ylabel(r"$BV_{IRAF}-BV_{Python}$")
med = np.nanmedian(bv_i - bv_p)
plt.title("BV ({:.3f})".format(med))
plt.scatter(v_i, bv_i - bv_p, s=4)
plt.axhline(y=med, c='r')

ax = fig.add_subplot(gs[2])
plt.ylim(-.05, .05)
plt.xlabel(r"$V_{IRAF}$")
plt.ylabel(r"$VI_{IRAF}-VI_{Python}$")
med = np.nanmedian(vi_i - vi_p)
plt.title("VI ({:.3f})".format(med))
plt.scatter(v_i, vi_i - vi_p, s=4)
plt.axhline(y=med, c='r')

ax = fig.add_subplot(gs[3])
# plt.ylim(0., -.15)
plt.xlabel(r"$V_{IRAF}$")
plt.ylabel(r"$UB_{IRAF}-UB_{Python}$")
med = np.nanmedian(ub_i - ub_p)
plt.title("UB ({:.3f})".format(med))
plt.scatter(v_i, ub_i - ub_p, s=4)
plt.axhline(y=med, c='r')

fig.tight_layout()
plt.savefig('output/phot_compare.png', dpi=300, bbox_inches='tight')
