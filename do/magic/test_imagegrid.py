
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np

from pts.magic.tools import plotting
from pts.core.basics.plot import MPLPlot, MPLFigure
from pts.magic.core.frame import Frame
from pts.magic.plot.imagegrid import StandardImageGridPlotter, ResidualImageGridPlotter

# -----------------------------------------------------------------

# Logging
from pts.core.basics.log import setup_log
setup_log("DEBUG")

# -----------------------------------------------------------------

# Set seed
np.random.seed(1)

# -----------------------------------------------------------------

def make_random_frame(nxpixels, nypixels=None):

    """
    This function ...
    :param nxpixels:
    :param nypixels:
    :return:
    """

    # Set shape
    if nypixels is None: nypixels = nxpixels
    shape = (nypixels, nxpixels)

    # Return the frame
    return Frame.random(shape)

# -----------------------------------------------------------------

def make_random_frames(nframes, min_npixels=50, max_npixels=200, xsize=None, ysize=None):

    """
    This function ...
    :param nframes:
    :param min_npixels:
    :param max_npixels:
    :param xsize:
    :param ysize:
    :return:
    """

    frames = []

    # Make nframes frames
    for _ in range(nframes):

        # Determine shape
        if xsize is None: xsize = int(np.random.uniform(min_npixels, max_npixels))
        if ysize is None: ysize = int(np.random.uniform(min_npixels, max_npixels))

        # Create the random frame
        frame = make_random_frame(xsize, ysize)

        # Add the frame
        frames.append(frame)

    # Return the frames
    return frames

# -----------------------------------------------------------------

def test_direct(figsize=(4,4), nframes=5):

    """
    This function ...
    :param figsize:
    :param nframes:
    :return:
    """

    # Get the frames
    frames = make_random_frames(nframes)

    # Setup the figure
    figure = plt.figure(figsize=figsize)
    plt.clf()
    ax = figure.gca()

    # No axes for the main figure
    ax.set_axis_off()

    ncols = 2
    nrows = 6

    wspace = 0.0
    hspace = 0.0

    cbar_mode = "single"
    axes_pad = (wspace, hspace)

    colorbar_relsize=0.05
    cbar_size = str(colorbar_relsize*100) + "%"

    axes_class = None
    label_mode = "L"

    # Create image grid
    grid = ImageGrid(figure, 111,  # similar to subplot(111)
                     nrows_ncols=(nrows, ncols),  # creates 2x2 grid of axes
                     axes_pad=axes_pad,  # pad between axes in inch.
                     aspect=True,
                     cbar_mode=cbar_mode, add_all=True, cbar_set_cax=False, cbar_size=cbar_size, axes_class=axes_class,
                     label_mode=label_mode)

    # Initialize structure to contain the plots
    plots = [[None for i in range(ncols)] for j in range(nrows)]

    # Loop over the images
    index = 0
    for row in range(nrows):
        for col in range(ncols):
            # Get axes, create subplot?
            ax = grid[index]
            plot = ax
            # Create plot
            plot = MPLPlot(plot=plot)
            # Add the plot
            plots[row][col] = plot
            index += 1

    index = 0
    for i in range(nrows):
        for j in range(ncols):

            # Get the plot
            plot = plots[i][j]

            # Color spines
            plot.axes.spines['bottom'].set_color("white")
            plot.axes.spines['top'].set_color("white")
            plot.axes.spines['left'].set_color("white")
            plot.axes.spines['right'].set_color("white")

            # Color ticks
            # plot.axes.xaxis.label.set_color("white")
            # plot.axes.yaxis.label.set_color("white")
            plot.axes.tick_params(axis='x', colors="white", direction="inout")
            plot.axes.tick_params(axis='y', colors="white", direction="inout")

            #plot.axes.set_facecolor("black")
            #plot.axes.set_adjustable('box-forced')

            #im = np.arange(100)

            #im.shape = xsize, ysize

            # Get the next frame
            frame = frames[index]

            #grid[i].imshow(im)  # The AxesGrid object work as a list of axes.
            plotting.plot_frame(frame, axes=plot.axes)

            # Add the label
            plot.axes.text(0.95, 0.95, "text", color='white', transform=plot.axes.transAxes, fontsize=10, va="top", ha="right")

    plt.show()
    plt.close()

# -----------------------------------------------------------------

def test_standard(nframes=5):

    """
    This function ...
    :param nframes:
    :return:
    """

    # Get the frames
    frames = make_random_frames(nframes)

    # Initialize the plotter
    plotter = StandardImageGridPlotter()

    # Loop over the frames
    for index, frame in enumerate(frames):

        # Add the frame
        plotter.add_frame(frame, str(index))

    # Run the plotter
    plotter.run()

# -----------------------------------------------------------------

def test_residual(nframes=5, ngrids=2, max_nrows=3, add_small=False, small_size=6, small_where="last",
                  share_scale=True, scale_reference=None,
                  share_scale_residuals=False, scale_residuals_reference=None, shape=None, adjust_grid=None,
                  relative=True, absolute=False, distributions=False):

    """
    This function ...
    :param nframes:
    :param ngrids:
    :param max_nrows:
    :param add_small:
    :param small_size:
    :param small_where:
    :param share_scale:
    :param scale_reference:
    :param share_scale_residuals:
    :param scale_residuals_reference:
    :param shape:
    :param adjust_grid:
    :param relative:
    :param absolute:
    :param distributions:
    :return:
    """

    # Same shape
    if shape is not None: frames = make_random_frames(nframes, xsize=shape[1], ysize=shape[0])

    # Get frames
    elif add_small:

        if small_where == "last":

            frames = make_random_frames(nframes-1)
            small_frame = make_random_frame(small_size)
            frames.append(small_frame)

        elif small_where == "first":

            small_frame = make_random_frame(small_size)
            frames = [small_frame]
            frames.extend(make_random_frames(nframes-1))

        else: raise ValueError("Invalid option for 'small_where'")

    # Make all random frames
    else: frames = make_random_frames(nframes)

    # Initialize the plotter
    plotter = ResidualImageGridPlotter()
    plotter.config.distributions = distributions
    plotter.config.max_nrows = max_nrows
    plotter.config.ngrids = ngrids

    # Set scale references
    plotter.config.share_scale = share_scale
    plotter.config.scale_reference = scale_reference
    plotter.config.share_scale_residuals = share_scale_residuals
    plotter.config.scale_residuals_reference = scale_residuals_reference

    plotter.config.adjust_grid = adjust_grid

    plotter.config.relative = relative
    plotter.config.absolute = absolute

    # Loop over the frames
    for index, frame in enumerate(frames):

        name = str(index)

        # Show the shape of the image
        #print(name, frame.xsize, frame.ysize)

        # Make observation and model frame
        observation = frame
        model = frame + Frame.random_normal(frame.shape, mean=0.0, sigma=0.5)

        # Add row
        plotter.add_row(observation, model, name, with_residuals=True)

    # Run the plotter
    plotter.run()

# -----------------------------------------------------------------

#test_residual(add_small=True, small_where="first")
#test_residual(add_small=True, small_where="last", share_scale_residuals=True, scale_residuals_reference="4", adjust_grid=True)
#test_residual(add_small=True, small_where="last", adjust_grid=True)
#test_residual(add_small=True, shape=(100,100), adjust_grid=True)

test_residual(add_small=True, shape=(100,100), adjust_grid=True, relative=False, distributions=True)
#test_residual(add_small=True, shape=(100,100), adjust_grid=True, relative=False)

# -----------------------------------------------------------------
