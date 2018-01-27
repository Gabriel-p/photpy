
import numpy as np
from match_funcs import xyTrans, autoSrcDetect


def main(pars, hdulist, xy_ref, mags_ref):
    """
    """

    _, xy_dtct, mags_dtct = autoSrcDetect(pars, hdulist)

    xmi, xma = map(float, pars['xtr_min-xtr_max'][0])
    ymi, yma = map(float, pars['ytr_min-ytr_max'][0])
    max_shift = [[xmi, xma], [ymi, yma]]
    mtoler = float(pars['match_toler'])
    xy_shift = xyTrans(
        max_shift, xy_ref, np.array(mags_ref), xy_dtct, np.array(mags_dtct),
        mtoler)

    scale, rot_angle = np.nan, np.nan

    return scale, rot_angle, xy_shift, xy_dtct
