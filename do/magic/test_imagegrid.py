
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

def make_random_frames(nframes, min_npixels=50, max_npixels=200):

    """
    This function ...
    :param nframes:
    :param min_npixels:
    :param max_npixels:
    :return:
    """

    frames = []

    for _ in range(nframes):

        xsize = int(np.random.uniform(min_npixels, max_npixels))
        ysize = int(np.random.uniform(min_npixels, max_npixels))
        shape = (ysize, xsize)

        print(xsize, ysize)

        # Create the random frame
        frame = Frame.random(shape)

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

            #plot.axes.set_axis_bgcolor("black")
            #plot.axes.set_adjustable('box-forced')

            #im = np.arange(100)

            #im.shape = xsize, ysize

            # Get the next frame
            frame = frames[index]

            #grid[i].imshow(im)  # The AxesGrid object work as a list of axes.
            plotting.plot_frame(frame, axes=plot.axes, return_image=True, return_normalization=True)

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

def test_residual(nframes=5, ngrids=2, max_nrows=3):

    """
    This function ...
    :param nframes:
    :param ngrids:
    :param max_nrows:
    :return:
    """

    # Get frames
    frames = make_random_frames(nframes)

    # Initialize the plotter
    plotter = ResidualImageGridPlotter()
    #plotter.config.distributions = True
    plotter.config.max_nrows = max_nrows
    plotter.config.ngrids = ngrids

    # Loop over the frames
    for index, frame in enumerate(frames):

        name = str(index)

        # Make observation and model frame
        observation = frame
        model = frame * Frame.random_normal(frame.shape, mean=0.0, sigma=0.5)

        # Add row
        plotter.add_row(observation, model, name, with_residuals=True)

    # Run the plotter
    plotter.run()

# -----------------------------------------------------------------

test_residual()

# -----------------------------------------------------------------
