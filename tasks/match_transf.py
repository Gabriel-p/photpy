
import operator
import numpy as np
from scipy.spatial.distance import cdist

import matplotlib.pyplot as plt


def genData():
    """
    Generate random (x,y) coordinates
    """
    frames = {'U': {'60': [], '100': [], '250': []},
              'B': {'30': [], '70': [], '200': []},
              'V': {'30': [], '70': [], '100': []},
              'I': {'30': [], '50': [], '100': []}}
    # Initial positions
    x = np.random.uniform(0., 4000., 5000)
    y = np.random.uniform(0., 4000., 5000)
    # Factor that determines the number of stars per filter.
    fact = {'U': [5, .5], 'B': [7, .7], 'V': [8, .8], 'I': [10, 1.]}
    for filt, expDict in frames.iteritems():
        for exps in expDict.keys():
            N = int(float(exps) * np.random.randint(4, fact[filt][0]) *
                    fact[filt][1])
            # Sigma and mu for normal distribution
            sigma, mu = 1., 1.
            # Generate normally distributed noise.
            noise = sigma * np.random.randn(len(x)) + mu
            x_n, y_n = x + noise, y + noise
            x_n, y_n = x_n[:N], y_n[:N]
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

    * i: index of ith star in fr1.
    * dist[i]: stores the distances from this star to *all* stars in fr2
    * min_dist_idx[i]: index to the star in fr2 closest to the ith star in fr1
      It is also the index of the minimum distance in dist[i], i.e.: distance
      to the closest star in fr2 to the ith star in fr1.
    * fr2[min_dist[i]]: closest star in fr2 to this star in fr1.
    * dist[i][min_dist_idx[i]]: distance between these two stars.

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

    * match_c1_ids, match_c2_ids: ids for unique matches between frames
    * no_match_c1: ids of reference stars with no match found
    * dupl_c1_ids: ids of reference stars with a duplicated match that will
      be re-processed.
    """

    match_c1_ids, dupl_c1_ids, no_match_c1 = [], [], []
    match_c2_ids, match_d2d, no_match_d2d = [], [], []
    # Indexes of matched star in both frames and distance between them.
    for c1_i, c2_i, d in zip(*[c1_ids, c2_ids, d2d]):
        # Filter by maximum allowed match distance.
        if d <= maxrad:
            # Check if this queried star in c2 was already stored.
            if c2_i in match_c2_ids:
                # Index of this stored c2 star.
                j = match_c2_ids.index(c2_i)
                # Only replace with this match if the previous match had a
                # larger distance.
                if match_d2d[j] > d:
                    # Store rejected replaced star.
                    dupl_c1_ids.append(match_c1_ids[j])
                    # Replace star
                    match_c1_ids[j] = c1_i
                    match_d2d[j] = d
                else:
                    dupl_c1_ids.append(c1_i)
            else:
                match_c1_ids.append(c1_i)
                match_c2_ids.append(c2_i)
                match_d2d.append(d)
        else:
            # If this observed star has no queried star closer than the max
            # distance allowed, discard.
            no_match_c1.append(c1_i)
            no_match_d2d.append(d)

    return match_c1_ids, match_c2_ids, no_match_c1, dupl_c1_ids, match_d2d,\
        no_match_d2d


def main():
    """
    """
    maxrad = 5.

    frames = genData()
    # Read data in the correct order.
    refFrameInfo, refFrame, framesOrdered = framesOrder(frames)

    for frame in framesOrdered:
        print("\nReference frame:")
        for _ in refFrameInfo:
            print("{}, {} (N={})".format(*_))

        # Extract (x,y) coordinates, averaging the values assigned to the
        # same star.
        x_ref, y_ref = [], []
        for st in refFrame:
            x_ref.append(np.mean(st[0]))
            y_ref.append(np.mean(st[1]))

        fr_filt, fr_expTime = frame[:2]
        x_fr_orig, y_fr_orig = frame[2][:2]
        # Copy of original lists, as these will be edited.
        x_fr, y_fr = list(x_fr_orig), list(y_fr_orig)
        print('\nProcessing frame: {}, {} (N={})\n'.format(
            fr_filt, fr_expTime, len(x_fr)))

        # Initial full list of IDs for the reference and processed frame.
        fr1_ids = [_ for _ in range(len(refFrame))]
        fr2_ids = [_ for _ in range(len(x_fr_orig))]

        match_fr1_ids_all, match_fr2_ids_all, no_match_fr1_all, match_d2d_all,\
            no_match_d2d_all = [], [], [], [], []
        # Continue until no more duplicate matches exist.
        i_while = 1
        while fr1_ids:

            # Match reference and processed frame.
            fr2_ids_dup, fr1fr2_d2d = closestStar(x_ref, y_ref, x_fr, y_fr)

            match_fr1_ids, match_fr2_ids, no_match_fr1, fr1_ids, match_d2d,\
                no_match_d2d = starMatch(
                    fr1_ids, fr2_ids_dup, fr1fr2_d2d, maxrad)
            # Store all unique solutions and no match solutions.
            match_fr1_ids_all += match_fr1_ids
            match_fr2_ids_all += match_fr2_ids
            no_match_fr1_all += no_match_fr1
            match_d2d_all += match_d2d
            no_match_d2d_all += no_match_d2d

            print("{}.".format(i_while))
            print("Matched reference stars: {}".format(len(match_fr1_ids)))
            if match_fr1_ids:
                print("(Mean match dist: {:.2f} px)".format(
                    np.mean(match_d2d_all)))
            print("Reference stars w/ no match within maxrad: {}".format(
                len(no_match_fr1)))

            if fr1_ids:
                print("Reference stars for re-match: {}".format(
                    len(fr1_ids)))
                # 'frame' stars that were not matched to an observed star.
                fr2_ids_r = [_ for _ in fr2_ids if _ not in match_fr2_ids_all]
                print("Frame stars for re-match: {}".format(len(fr2_ids_r)))
                # To avoid messing with the indexes, change the coordinates
                # of already matched 'frame' stars so that they can not
                # possibly be matched again.
                x_fr, y_fr = [], []
                for c2_i in fr2_ids:
                    if c2_i in match_fr2_ids_all:
                        x_fr.append(-1000.)
                        y_fr.append(-1000.)
                    else:
                        x_fr.append(x_fr_orig[c2_i])
                        y_fr.append(y_fr_orig[c2_i])
            i_while += 1

        print('\nTotal stars matched: {}'.format(len(match_fr1_ids_all)))
        print('Total stars not matched: {}'.format(len(no_match_fr1_all)))


if __name__ == '__main__':
    main()
