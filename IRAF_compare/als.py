
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii


# .fits file to be read
N_als = '42'

# IRAF OLD photometry
phot_old = ascii.read(
    '2_PSF/test/' + N_als + '/stk_fcd00' + N_als + '.als_OLD.1',
    format='daophot')
id_old, x_old, y_old, v_old = phot_old['ID'], phot_old['XCENTER'],\
    phot_old['YCENTER'], phot_old['MAG']

# # IRAF NEW1 photometry
# phot_new = ascii.read(
#     '2_PSF/test/' + N_als + '/stk_fcd00' + N_als + '.als_NEW1.1',
#     format='daophot')
# id_new, x_new, y_new, v_new = phot_new['ID'], phot_new['XCENTER'],\
#     phot_new['YCENTER'], phot_new['MAG']
# # id_old, x_old, y_old, v_old = phot_new['ID'], phot_new['XCENTER'],\
# #     phot_new['YCENTER'], phot_new['MAG']

# IRAF NEW2 photometry
phot_new = ascii.read(
    '2_PSF/test/' + N_als + '/stk_fcd00' + N_als + '.als_NEW2.1',
    format='daophot')
id_new, x_new, y_new, v_new = phot_new['ID'], phot_new['XCENTER'],\
    phot_new['YCENTER'], phot_new['MAG']

print("Data read")

# IDs in common
ids = np.intersect1d(id_old, id_new)
# Index of shared IDs for both arrays
m_old = np.nonzero(ids[:, None] == id_old)[1]
m_new = np.nonzero(ids[:, None] == id_new)[1]

# Compare matching stars.
v_new, v_old = v_new[m_new], v_old[m_old]
x_new, y_new = x_new[m_new], y_new[m_new]
x_old, y_old = x_old[m_old], y_old[m_old]

plt.style.use('seaborn-darkgrid')
# fig = plt.figure(figsize=(10, 10))
# gs = gridspec.GridSpec(4, 1)

# ax = fig.add_subplot(gs[0])
med = np.nanmedian(v_old - v_new)
plt.title("V ({:.3f})".format(med))
# plt.ylim(-.05, .05)
plt.xlabel(r"$V_{OLD}$")
plt.ylabel(r"$V_{OLD}-V_{NEW}$")
plt.scatter(v_old, v_old - v_new, s=4)
plt.axhline(y=med, c='r')
plt.axhline(y=0., c='g')

plt.show()
# fig.tight_layout()
# plt.savefig('phot_compare.png', dpi=300, bbox_inches='tight')
