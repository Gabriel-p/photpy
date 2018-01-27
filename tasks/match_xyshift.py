
import numpy as np
from match_funcs import xyTrans, autoSrcDetect


def main(pars, hdulist, xy_ref, mags_ref):
    """
    """

    _, xy_dtct, mags_dtct = autoSrcDetect(pars, hdulist)

    xmax, ymax = np.array(xy_dtct).max(axis=0)

    xtr, ytr = pars['xtr_min-xtr_max'][0], pars['ytr_min-ytr_max'][0]
    none_list = ['none', 'None', 'NONE']
    xmi = float(xtr[0]) if xtr[0] not in none_list else -xmax
    xma = float(xtr[1]) if xtr[1] not in none_list else xmax
    ymi = float(ytr[0]) if ytr[0] not in none_list else -ymax
    yma = float(ytr[1]) if ytr[1] not in none_list else ymax
    max_shift = ((xmi, xma), (ymi, yma))
    print(" xy trans limits: ({:.1f}, {:.1f}), ({:.1f}, {:.1f})".format(
        xmi, xma, ymi, yma))

    mtoler = float(pars['match_toler'])
    xy_shift = xyTrans(
        max_shift, xy_ref, np.array(mags_ref), xy_dtct, np.array(mags_dtct),
        mtoler)

    scale, rot_angle = np.nan, np.nan

    return scale, rot_angle, xy_shift, xy_dtct
