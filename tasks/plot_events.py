
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
from astropy.io import fits
from photutils import CircularAperture
import time


def zoom(event, ax_all):
    """
    Source: https://gist.github.com/tacaswell/3144287
    """
    # get event location
    xdata, ydata = event.xdata, event.ydata
    base_scale = 2.

    for axX in ax_all:
        if event.inaxes == axX:
            ax = axX
            break
    # else:
    #     ax = ax2

    # get the current x and y limits
    cur_xlim = ax.get_xlim()
    cur_ylim = ax.get_ylim()
    # set the range
    cur_xrange = (cur_xlim[1] - cur_xlim[0]) * .5
    cur_yrange = (cur_ylim[1] - cur_ylim[0]) * .5
    if event.button == 'up':
        # deal with zoom in
        scale_factor = 1 / base_scale
    elif event.button == 'down':
        # deal with zoom out
        scale_factor = base_scale
    else:
        # deal with something that should never happen
        scale_factor = 1
        print event.button
    # set new limits
    ax.set_xlim([xdata - cur_xrange * scale_factor,
                 xdata + cur_xrange * scale_factor])
    ax.set_ylim([ydata - cur_yrange * scale_factor,
                 ydata + cur_yrange * scale_factor])
    # force re-draw
    ax.figure.canvas.draw()


def drawCirclesIDs(ax, hdu_data, id_selec, xy):
    """
    """
    apertures = CircularAperture(xy, r=10.)
    apertures.plot(ax=ax, color='red', lw=1.)
    xmax, ymax = hdu_data.shape[1], hdu_data.shape[0]
    for j, txt in enumerate(id_selec):
        xtxt, ytxt = xy[j][0] + .01 * xmax, xy[j][1] + .01 * ymax
        ax.annotate(
            txt, xy=(xy[j][0], xy[j][1]),
            xytext=(xtxt, ytxt),
            fontsize=15, color='r', arrowprops=dict(
                arrowstyle="-", color='b',
                connectionstyle="angle3,angleA=90,angleB=0"))


def refDisplay(
    ref_field_img, f_name, hdu_data, id_ref, xy_ref, id_selec, xy_select):
    """
    Manual selection of reference stars in the observed frame.
    """
    # Manual selection of reference stars in the observed frame.
    fig, (ax1, ax2) = plt.subplots(1, 2)
    interval = ZScaleInterval()

    if ref_field_img.endswith('.fits'):
        # Load .fits file.
        hdulist = fits.open(ref_field_img)
        hdu_data_ref = hdulist[0].data
        zmin, zmax = interval.get_limits(hdu_data_ref)
        ax1.imshow(
            hdu_data_ref, cmap='Greys', aspect=1, interpolation='nearest',
            origin='lower', vmin=zmin, vmax=zmax)
        drawCirclesIDs(ax1, hdu_data_ref, id_ref, xy_ref)
    else:
        ax1.imshow(plt.imread(ref_field_img))

    zmin, zmax = interval.get_limits(hdu_data)
    ax2.imshow(
        hdu_data, cmap='Greys', aspect=1, interpolation='nearest',
        origin='lower', vmin=zmin, vmax=zmax)
    drawCirclesIDs(ax2, hdu_data, id_selec, xy_select)
    plt.title("Coordinates over {}".format(f_name))

    fig.canvas.mpl_connect(
        'scroll_event', lambda event: zoom(event, [ax1, ax2]))
    plt.show()
    plt.close('all')


def refStrSlct(ref_field_img, f_name, hdu_data, id_ref, xy_ref):
    """
    """
    id_pick, xy_pick = [], []
    # Ref stars selection

    # Manual selection of reference stars in the observed frame.
    fig, (ax1, ax2) = plt.subplots(1, 2)
    interval = ZScaleInterval()

    ax1.set_title("Reference image")
    if ref_field_img.endswith('.fits'):
        # Load .fits file.
        hdulist = fits.open(ref_field_img)
        hdu_data_ref = hdulist[0].data
        zmin, zmax = interval.get_limits(hdu_data_ref)
        ax1.imshow(
            hdu_data_ref, cmap='Greys', aspect=1, interpolation='nearest',
            origin='lower', vmin=zmin, vmax=zmax)
        drawCirclesIDs(ax1, hdu_data_ref, id_ref, xy_ref)
    else:
        ax1.imshow(plt.imread(ref_field_img))

    zmin, zmax = interval.get_limits(hdu_data)
    ax2.imshow(
        hdu_data, cmap='Greys', aspect=1, interpolation='nearest',
        origin='lower', vmin=zmin, vmax=zmax)
    plt.title("Observed frame: {}".format(f_name))

    def onclick(event, ax):
        ax.time_onclick = time.time()

    def onrelease(event, ax):
        # Only clicks inside this axis.
        if event.inaxes == ax:
            if event.button == 1 and\
                    ((time.time() - ax.time_onclick) < .1):
                apertures = CircularAperture(
                    (event.xdata, event.ydata), r=10.)
                apertures.plot(color='green', lw=1.)
                ax.figure.canvas.draw()
                ref_id = raw_input(" ID of selected star: ").strip()
                print(" {} added to list: ({:.2f}, {:.2f})".format(
                    ref_id, event.xdata, event.ydata))
                id_pick.append(ref_id)
                xy_pick.append((event.xdata, event.ydata))
                if len(id_pick) == 3:
                    # Exit when 3 stars have been identified
                    plt.close()
            elif event.button == 2:
                print("scroll click")
            elif event.button == 3:
                print("right click")
            else:
                pass

    # Mouse click / scroll zoom events.
    fig.canvas.mpl_connect(
        'scroll_event', lambda event: zoom(event, [ax1, ax2]))
    fig.canvas.mpl_connect(
        'button_press_event', lambda event: onclick(event, ax2))
    fig.canvas.mpl_connect(
        'button_release_event', lambda event: onrelease(event, ax2))

    print(" (Select three reference stars covering\n"
          "  as much frame as possible)")
    plt.show()
    plt.close('all')

    return id_pick, xy_pick
