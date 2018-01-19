

def zoom(event, ax):
    """
    Source: https://gist.github.com/tacaswell/3144287
    """
    base_scale = 2.
    # get the current x and y limits
    cur_xlim = ax.get_xlim()
    cur_ylim = ax.get_ylim()
    # set the range
    cur_xrange = (cur_xlim[1] - cur_xlim[0]) * .5
    cur_yrange = (cur_ylim[1] - cur_ylim[0]) * .5
    # get event location
    xdata, ydata = event.xdata, event.ydata
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
