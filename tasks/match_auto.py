
from match_funcs import triangleMatch, reCenter, autoSrcDetect


def main(pars, xy_ref, hdulist):
    """
    """
    hdu_data = hdulist[0].data

    # Scale and rotation ranges, and match tolerance.
    scale_range = (float(pars['scale_min']), float(pars['scale_max']))
    rot_range = (float(pars['rot_min']), float(pars['rot_max']))
    trans_range = (
        map(float, pars['xtr_min-xtr_max'][0]),
        map(float, pars['ytr_min-ytr_max'][0]))
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
