def center_axis(ax, left_spine=True, bottom_spine=True):
    # ax.set_aspect('equal')
    ax.grid(True, which='both')

    if left_spine:
        # set the x-spine (see below for more info on `set_position`)
        ax.spines['left'].set_position('zero')

    # turn off the right spine/ticks
    ax.spines['right'].set_color('none')
    ax.yaxis.tick_left()

    if bottom_spine:
        # set the y-spine
        ax.spines['bottom'].set_position('zero')

    # turn off the top spine/ticks
    ax.spines['top'].set_color('none')
    ax.xaxis.tick_bottom()


def center_right_axis(ax):
    # ax.set_aspect('equal')
    ax.grid(True, which='both')

    ax.spines['right'].set_position('zero')
    # ax.yaxis.tick_left()

    # set the y-spine
    ax.spines['bottom'].set_position('zero')

    # turn off the top spine/ticks
    ax.spines['top'].set_color('none')
    ax.xaxis.tick_bottom()
