
import numpy as np
from match_funcs import reCenter, getTriangles, findRotAngle, standard2observed
from plot_events import refDisplay, refStrSlct


def srcSelect(
    f_name, hdu_data, ref_field_img, id_ref, xy_ref, id_selec, xy_selec):
    """
    Manually mark sources on 'fr' fits file.
    """
    answ = 'n'
    if id_selec:

        for i, xy in enumerate(xy_selec):
            # Re-center coordinates.
            xy = reCenter(hdu_data, xy, side=40).tolist()
            # Display frames.
            refDisplay(
                ref_field_img, f_name, hdu_data, id_ref, xy_ref, id_selec[i],
                xy)

            while True:
                answ = raw_input(
                    "Use coordinates ({})? (y/n): ".format(i)).strip()
                if answ == 'y':
                    idx_c = i
                    # Re-write with re-centered coordinates.
                    xy_selec[idx_c] = xy
                    break
                elif answ == 'n':
                    break
                else:
                    print("Wrong answer. Try again.")

    # Select new coordinates
    if answ != 'y':
        while True:
            # Manually pick three reference stars.
            id_pick, xy_pick = refStrSlct(
                ref_field_img, f_name, hdu_data, id_ref, xy_ref)
            if len(id_pick) < 1:
                print("ERROR: at least one star must be selected.")
            elif '' in id_pick:
                print("ERROR: all stars must be identified. Try again.")
            else:
                break

        # Re-center coordinates.
        xy_cent = reCenter(hdu_data, xy_pick, side=20).tolist()

        # Select latest coordinates.
        idx_c = -1
        id_selec.append(id_pick)
        xy_selec.append(xy_cent)

        print(id_pick)
        print(xy_cent)

    print(" Stars selected for match: {}".format(len(xy_selec[idx_c])))

    return id_selec, xy_selec, idx_c


def main(f_name, hdu_data, ref_field_img, id_ref, xy_ref, id_selec, xy_selec):
    """
    """
    # Coordinates from observed frame.
    id_selec, xy_selec, idx_c = srcSelect(
        f_name, hdu_data, ref_field_img, id_ref, xy_ref, id_selec, xy_selec)

    # Selected IDs and coordinates.
    id_choice, xy_choice = id_selec[idx_c], xy_selec[idx_c]

    # Match selected reference stars to standard stars by
    # the IDs given by the user.
    xy_ref_sel = []
    for r_id in id_choice:
        i = id_ref.index(r_id)
        xy_ref_sel.append(xy_ref[i])

    # Matched stars in ref and obs
    std_tr_match, obs_tr_match = xy_ref_sel, xy_choice

    if len(id_choice) < 3:
        print("  WARNING: Less than 3 stars selected,"
              " no match can be performed.")
        xy_transf = []
        for st_id in id_ref:
            if st_id in id_choice:
                i = id_choice.index(st_id)
                xy_transf.append(
                    [xy_choice[i][0], xy_choice[i][1]])
            else:
                xy_transf.append([np.nan, np.nan])

        scale, rot_angle, xy_shift = np.nan, np.nan,\
            [np.nan, np.nan]

    else:
        _, A_tr_not_scaled = getTriangles(xy_ref_sel, [(0, 1, 2)])
        _, B_tr_not_scaled = getTriangles(xy_choice, [(0, 1, 2)])

        # Obtain scale.
        scale = np.mean(
            np.array(B_tr_not_scaled[0]) / A_tr_not_scaled[0])
        # Rotation angle between triangles.
        rot_angle, ref_cent, obs_cent = findRotAngle(xy_ref_sel, xy_choice)

        # Apply translation, scaling, and rotation to reference
        # coordinates.
        xy_transf, rot_angle = standard2observed(
            xy_ref, scale, rot_angle, ref_cent, obs_cent, obs_tr_match)

        xy_shift = ref_cent - obs_cent

    print("Re-center final coordinates.")
    xy_transf = reCenter(hdu_data, xy_transf, side=40)

    return std_tr_match, obs_tr_match, scale, rot_angle, xy_shift, xy_transf
