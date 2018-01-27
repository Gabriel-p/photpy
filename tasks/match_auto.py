
import numpy as np
from match_funcs import triangleMatch, reCenter, autoSrcDetect


def main(pars, xy_ref, hdulist):
    """
    """
    hdu_data = hdulist[0].data

    # Scale and rotation ranges, and match tolerance.
    scale_range = (float(pars['scale_min']), float(pars['scale_max']))
    rot_range = (float(pars['rot_min']), float(pars['rot_max']))

    xtr, ytr = pars['xtr_min-xtr_max'][0], pars['ytr_min-ytr_max'][0]
    none_list = ['none', 'None', 'NONE']
    tr_x_min = float(xtr[0]) if xtr[0] not in none_list else -np.inf
    tr_x_max = float(xtr[1]) if xtr[1] not in none_list else np.inf
    tr_y_min = float(ytr[0]) if ytr[0] not in none_list else -np.inf
    tr_y_max = float(ytr[1]) if ytr[1] not in none_list else np.inf
    trans_range = ((tr_x_min, tr_x_max), (tr_y_min, tr_y_max))

    mtoler = float(pars['match_toler'])

    _, xy_choice, _ = autoSrcDetect(pars, hdulist)

    print("\nFinding scale, translation, and rotation.")
    std_tr_match, obs_tr_match, scale, rot_angle, xy_shift, xy_transf =\
        triangleMatch(
            xy_ref, xy_choice, mtoler, scale_range, rot_range, trans_range,
            hdu_data)

    print("Re-center final coordinates.")
    xy_transf = reCenter(hdu_data, xy_transf, side=40)

    return std_tr_match, obs_tr_match, scale, rot_angle, xy_shift, xy_transf
