
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.io import ascii
from astropy.table import Table


phot = ascii.read('3_DAOM/daom.obs', fill_values=('INDEF', np.nan))
id_, x, y, v, ev, b, eb, u, eu, i, ei, XV, XB, XU, XI =\
    phot['col1'], phot['col2'], phot['col3'], phot['col4'], phot['col5'],\
    phot['col6'], phot['col7'], phot['col8'], phot['col9'], phot['col10'],\
    phot['col11'], phot['col14'], phot['col15'], phot['col16'], phot['col17']

# print(len(x))
# plt.subplot(121)
# plt.scatter(x, y, s=5)
# msk = ~np.isnan(v)
# x, y = x[msk], y[msk]
# print(len(x))
# plt.subplot(122)
# plt.scatter(x, y, s=5)
# plt.show()

ext_coeffs = {'U': .442, 'B': .232, 'V': .136, 'I': .048}
transf_coeffs = ascii.read('fit_coeffs.dat')

# Transform to zero airmass instrumental magnitudes.
V_za = v - ext_coeffs['V'] * XV[0]
B_za = b - ext_coeffs['B'] * XB[0]
U_za = u - ext_coeffs['U'] * XU[0]
I_za = i - ext_coeffs['I'] * XI[0]

# Generate instrumental colors.
BV_za, UB_za, VI_za = B_za - V_za, U_za - B_za, V_za - I_za

# Transform to calibrated standard system.
BV = transf_coeffs[1][3] * BV_za + transf_coeffs[1][4]
UB = transf_coeffs[2][3] * UB_za + transf_coeffs[2][4]
VI = transf_coeffs[3][3] * VI_za + transf_coeffs[3][4]
V = V_za + transf_coeffs[0][3] * BV + transf_coeffs[0][4]

# Errors in calibrated system.
eBV = (eb + ev) * transf_coeffs[1][3]
eVI = (ev + ei) * transf_coeffs[2][3]
eUB = (eu + eb) * transf_coeffs[3][3]
eV = ev

plt.style.use('seaborn-darkgrid')
fig = plt.figure(figsize=(30, 10))
gs = gridspec.GridSpec(1, 3)

ax = fig.add_subplot(gs[0])
plt.xlim(-1., 4.5)
plt.ylim(22.5, 10.)
plt.title("BV vs V")
plt.scatter(BV, V, s=3)

ax = fig.add_subplot(gs[1])
plt.xlim(-1., 4.5)
plt.ylim(22.5, 10.)
plt.title("VI vs V")
plt.scatter(VI, V, s=3)

ax = fig.add_subplot(gs[2])
plt.xlim(0., 3.)
plt.ylim(1.7, -1.)
plt.title("BV vs UB")
plt.scatter(BV, UB, s=3)

fig.tight_layout()
plt.savefig('output/final_phot.png', dpi=300, bbox_inches='tight')

tab = Table([id_, x, y, V, eV, BV, eBV, VI, eVI, UB, eUB],
            names=('ID', 'x', 'y', 'V', 'eV', 'BV', 'eBV', 'VI', 'eVI',
                   'UB', 'eUB'))
ascii.write(
    tab, 'input/phot_compare/final_phot.dat', format='fixed_width',
    delimiter=' ', formats={_: '%10.4f' for _ in tab.keys()[1:]},
    fill_values=[(ascii.masked, 'nan')], overwrite=True)
