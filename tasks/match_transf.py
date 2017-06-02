
import operator
import numpy as np
import random
from scipy.spatial.distance import cdist

import timeit
import matplotlib.pyplot as plt


def genData():
    """
    Generate random (x,y) coordinates
    """
    random.seed(9001)
    np.random.seed(117)
    frames = {'U': {'60': [], '100': [], '250': []},
              'B': {'30': [], '70': [], '200': []},
              'V': {'30': [], '70': [], '100': []},
              'I': {'30': [], '50': [], '100': []}}
    # Initial positions
    x = np.random.uniform(0., 4000., 50000)
    y = np.random.uniform(0., 4000., 50000)
    # Factor that determines the number of stars per filter.
    fact = {'U': [5, .5], 'B': [7, .7], 'V': [8, .8], 'I': [10, 1.]}
    for filt, expDict in frames.iteritems():
        for exps in expDict.keys():
            N = 10 * int(float(exps) * np.random.randint(4, fact[filt][0]) *
                    fact[filt][1])
            # Sigma and mu for normal distribution
            sigma, mu = 1., 1.
            # Generate normally distributed noise.
            noise = sigma * np.random.randn(len(x)) + mu
            x_n, y_n = x + noise, y + noise
            x_n, y_n = x_n[:N], y_n[:N]
            xy_n = list(zip(*[x_n, y_n]))
            random.shuffle(xy_n)
            x_n, y_n = np.asarray(zip(*xy_n))
            mag = np.random.uniform(10., 24., N)
            e_mag = np.random.uniform(.01, .2, N)
            # plt.scatter(x, y)
            # plt.scatter(x_n, y_n)
            # plt.show()

            # This is the information stored for each star.
            expDict[exps] = [x_n, y_n, mag, e_mag]

    return frames


def framesOrder(frames):
    """
    Assign reference frame as the one with the largest number of detected
    stars. Order the remaining frames with the largest one on top.

    Parameters
    ----------
    frames : dictionary
        Observed photometry. Contains sub-dictionaries with the exposure
        times as keys, and the photometry data as values.

    Returns
    -------
    refFrameInfo : list
        Information about the reference frame. Each sub-list contains the
        filter name, exposure time, and number of stars added to the list
        by appending this frame.
    refFrame : list
        Photometry of the reference frame, selected as the frame with the
        largest number of stars.
    framesOrdered : list
        Each sub-list is an observed frame containing its filter name,
        exposure time, and photometry. Ordered putting the frame with the
        largest number of stars on top.

    """
    flen = []
    for filt, f in frames.iteritems():
        for expTime, fexpT in f.iteritems():
            N = len(fexpT[0])
            flen.append([filt, expTime, N])
    # Sort putting the frames with the largest number of stars first.
    sortIdx = sorted(flen, key=operator.itemgetter(2), reverse=True)

    # Isolate reference frame.
    filt, expTime = sortIdx[0][0], sortIdx[0][1]
    # Store each star separately.
    stars = list(zip(*frames[filt][expTime]))
    # Convert each float to list, to prepare for later appending of data.
    # refFrame = [star1, star2, ...starN]
    # starX = [x, y, mag, e_mag]
    # x = [x1, x2, ...] ; y = [y1, y2, ...] ; mag = [mag1, mag2, ...]
    refFrame = [list(list([_]) for _ in st) for st in stars]
    # Store here names of filters, exposure times, and number of stars added
    # when processing each frame.
    refFrameInfo = [[filt, expTime, len(refFrame)]]

    framesOrdered = []
    for filt, expTime, dummy in sortIdx[1:]:
        framesOrdered.append([filt, expTime, frames[filt][expTime]])

    return refFrameInfo, refFrame, framesOrdered


def closestStar(x_fr1, y_fr1, x_fr2, y_fr2):
    """
    For every star in fr1, find the closest star in fr2.

    Parameters
    ----------
    x_fr1 : list
       x coordinates for stars in the reference frame.
    y_fr1 : list
       y coordinates for stars in the reference frame.
    x_fr2 : list
       x coordinates for stars in the processed frame.
    y_fr2 : list
       y coordinates for stars in the processed frame.

    Returns
    -------
    min_dist_idx : numpy array
        Index to the processed star closest to the reference star, for each
        reference star:
        * fr2[min_dist_idx[i]]: closest star in fr2 to this star in fr1.
        Also the index of the minimum distance in dist[i], i.e.: distance to
        the closest processed star to the ith reference star:
        * dist[i][min_dist_idx[i]]: distance between these two stars.

    Notes
    -----
    len(fr1) = len(dist) = len(min_dist_idx)


    """
    fr1 = np.array(zip(*[x_fr1, y_fr1]))
    fr2 = np.array(zip(*[x_fr2, y_fr2]))
    # Distance to all stars in fr2, for each star in fr1.
    dist = cdist(fr1, fr2, 'euclidean')
    # Indexes to minimum distances.
    min_dist_idx = np.argmin(dist, axis=1)
    # Store only the minimum distances for each star in fr1, to a star in fr2.
    min_dists = [dist[i][md_idx] for i, md_idx in enumerate(min_dist_idx)]

    return min_dist_idx, min_dists


def starMatch(c1_ids, c2_ids, d2d, maxrad):
    """
    Reject duplicated matches and matched stars with distances beyond
    the maximum separation defined.

    Parameters
    ----------
    c1_ids : list
        IDs of stars in the reference frame.
    c2_ids : list
        IDs of stars in the processed frame, closest to each star in the
        reference frame.
    d2d : list
        Distances between each star in the reference frame, and the closest
        star in the processed frame.
    maxrad : float
        Maximum allowed distance (radius) for a match to be valid.

    Returns
    -------
    match_c1_ids : list
        IDs in reference frame for unique matches between frames.
    match_c2_ids : list
        IDs in the processed frame for unique matches between frames.
    no_match_c1 : list
        IDs of reference stars with no match found within maxrad.
    dupl_c1_ids : list
        IDs of reference stars that were matched to the same processed star as
        another reference star, and had a larger distance. These stars will
        be re-processed.

    """

    match_c1_ids, dupl_c1_ids, match_c2_ids, match_d = [], [], [], []
    # Indexes of matched star in both frames and distance between them.
    for c1_i, c2_i, d in zip(*[c1_ids, c2_ids, d2d]):
        # Filter by maximum allowed match distance.
        if d <= maxrad:
            # Check if this processed star was already stored as a match with
            # another reference star.
            if c2_i in match_c2_ids:
                # Index of this stored c2 star.
                j = match_c2_ids.index(c2_i)
                # If the previous match had a larger distance than this match,
                # replace with this match.
                if match_d[j] > d:
                    # Store replaced reference star here.
                    dupl_c1_ids.append(match_c1_ids[j])
                    # Now replace this star.
                    match_c1_ids[j] = c1_i
                    match_d[j] = d
                else:
                    dupl_c1_ids.append(c1_i)
            else:
                # Store IDs of both matched stars, and their distance.
                match_c1_ids.append(c1_i)
                match_c2_ids.append(c2_i)
                match_d.append(d)
        # else:
            # This reference star has no processed star closer than the max
            # distance allowed.

    return dupl_c1_ids, match_c1_ids, match_c2_ids, match_d


def frameCoordsUpdt(x_fr, y_fr, match_fr2_ids):
    """
    Identify stars in the processed frame that where not matched to any star
    in the reference frame.
    To avoid messing with the indexes, change the coordinates of already
    matched 'frame' stars so that they will not be matched again.

    Parameters
    ----------
    x_fr : list
        Original x coordinates of the stars in the processed frame.
    y_fr : list
        Original y coordinates of the stars in the processed frame.
    match_fr2_ids : list
        IDs of stars in the processed frame that were matched to a star in
        the reference frame.

    Returns
    -------
    x_fr_updt, y_fr_updt : list, list
        Coordinates of frame stars with those identified as matches changed
        so that they will not be matched again.

    """
    # Modify coordinates of matched stars. Use copy of arrays to avoid
    # overwriting the original coordinate values in 'frame'.
    x_fr_updt, y_fr_updt = np.copy(x_fr), np.copy(y_fr)
    x_fr_updt[match_fr2_ids] = -1000.
    y_fr_updt[match_fr2_ids] = -1000.

    return x_fr_updt, y_fr_updt


def UpdtRefFrame(refFrameInfo, refFrame, frame, match_fr1_ids_all,
                 match_fr2_ids_all):
    """
    Update the reference frame adding the stars in the processed frame that
    were assigned as matches to each reference star. If a reference star was
    not assigned any match from the processed frame, add a Nan value.

    Also add to the end of the list (thereby increasing the length of the
    reference frame) those processed stars that could not be matched to any
    reference star.

    Parameters
    ----------
    """
    # Extract processed frame data.
    fr_filt, fr_expTime = frame[:2]
    x_fr, y_fr, mag_fr, emag_fr = frame[2]

    start_time = timeit.default_timer()
    # for each reference frame star.
    for ref_id, ref_st in enumerate(refFrame):
        # Check if this reference star was uniquely associated with a
        # processed star.
        if ref_id in match_fr1_ids_all:
            # Index of the associated processed star.
            j = match_fr1_ids_all.index(ref_id)
            fr_id = match_fr2_ids_all[j]
            #
            ref_st[0].append(x_fr[fr_id])
            ref_st[1].append(y_fr[fr_id])
            ref_st[2].append(mag_fr[fr_id])
            ref_st[3].append(emag_fr[fr_id])
        else:
            # If this reference star could not be matched to any processed star
            # within the maximum match radius defined, add a NaN value
            # to mean that no match in the processed frame was found
            # for this reference star.
            ref_st[0].append(np.nan)
            ref_st[1].append(np.nan)
            ref_st[2].append(np.nan)
            ref_st[3].append(np.nan)
    print("C1", timeit.default_timer() - start_time)

    # Number of frames processed this far including the reference frame, but
    # excluding this one.
    N_fr = len(refFrame[0][0]) - 1

    start_time = timeit.default_timer()
    # For each processed frame star.
    fr_st_no_match = 0
    for fr_id, fr_st in enumerate(zip(*[x_fr, y_fr, mag_fr, emag_fr])):
        if fr_id not in match_fr2_ids_all:
            # This frame star was not matched to any reference star.
            x = [np.nan for _ in range(N_fr)] + [fr_st[0]]
            y = [np.nan for _ in range(N_fr)] + [fr_st[1]]
            mag = [np.nan for _ in range(N_fr)] + [fr_st[2]]
            emag = [np.nan for _ in range(N_fr)] + [fr_st[3]]
            refFrame.append([x, y, mag, emag])
            fr_st_no_match += 1
    print("C2", timeit.default_timer() - start_time)

    # Update the information stored on the frames processed.
    refFrameInfo.append([fr_filt, fr_expTime, fr_st_no_match])

    return refFrameInfo, refFrame


def main():
    """
    """
    maxrad = 5.

    frames = genData()
    # Read data in the correct order.
    refFrameInfo, refFrame, framesOrdered = framesOrder(frames)

    for frame in framesOrdered:
        print("\n--------------------------------------")
        print("Reference frame (N={}), composed of:".format(
            np.sum(zip(*refFrameInfo)[2])))
        for _ in refFrameInfo:
            print("{}, {} (N={})".format(*_))

        # Extract (x,y) coordinates, averaging the values assigned to the
        # same star.
        x_ref, y_ref = [], []
        for st in refFrame:
            x_ref.append(np.mean(st[0]))
            y_ref.append(np.mean(st[1]))

        # Extract filter name and exposure time of the processed frame.
        fr_filt, fr_expTime = frame[:2]
        x_fr, y_fr = frame[2][:2]
        print('\nProcessing frame: {}, {} (N={})'.format(
            fr_filt, fr_expTime, len(x_fr)))

        # Initial full list of IDs for the reference and processed frame.
        fr1_ids = np.arange(len(refFrame)).tolist()
        fr2_ids = np.arange(len(x_fr)).tolist()

        match_fr1_ids_all, match_fr2_ids_all, match_d_all = [], [], []
        # Continue until no more duplicate matches exist.
        counter = 1
        while fr1_ids:

            start_time = timeit.default_timer()
            # Find closest stars between reference and processed frame.
            fr2_ids_dup, fr1fr2_d2d = closestStar(x_ref, y_ref, x_fr, y_fr)
            print("A", timeit.default_timer() - start_time)

            start_time = timeit.default_timer()
            # Match reference and processed frame.
            fr1_ids, match_fr1_ids, match_fr2_ids, match_d =\
                starMatch(fr1_ids, fr2_ids_dup, fr1fr2_d2d, maxrad)
            # Store unique matches and distances.
            match_fr1_ids_all += match_fr1_ids
            match_fr2_ids_all += match_fr2_ids
            match_d_all += match_d
            print("B", timeit.default_timer() - start_time)

            print("{}.".format(counter))
            counter += 1

            print("Matched reference stars: {}".format(len(match_fr1_ids)))
            if match_fr1_ids:
                print("(Mean match dist: {:.2f} px)".format(
                    np.mean(match_d_all)))
            print("Reference stars w/ no match within maxrad: {}".format(
                len(refFrame) - len(match_fr1_ids_all) - len(fr1_ids)))

            # If there are any stars from the reference frame that had
            # duplicated matches and were stored for re-matching.
            if fr1_ids:
                print("Reference stars for re-match: {}".format(
                    len(fr1_ids)))
                print("Frame stars for re-match: {}".format(
                    len(fr2_ids) - len(match_fr2_ids_all)))
                # Update coordinates of matched stars in processed frame.
                x_fr, y_fr = frameCoordsUpdt(x_fr, y_fr, match_fr2_ids)
            else:
                print("Processed stars w/ no match within maxrad: {}".format(
                    len(fr2_ids) - len(match_fr2_ids_all)))

        # Update reference frame associating all the matches found in the
        # processed frame to a given reference star. Also append those stars
        # from the processed frame with no match to the end of the list.
        refFrameInfo, refFrame = UpdtRefFrame(
            refFrameInfo, refFrame, frame, match_fr1_ids_all,
            match_fr2_ids_all)

    print("\nFinal combined reference frame (N={})".format(
        np.sum(zip(*refFrameInfo)[2])))
    for _ in refFrameInfo:
        print("{}, {} (N={})".format(*_))


if __name__ == '__main__':
    main()
