
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np

from pts.magic.tools import plotting
from pts.core.basics.plot import MPLPlot, MPLFigure
from pts.magic.core.frame import Frame
from pts.magic.plot.imagegrid import StandardImageGridPlotter, ResidualImageGridPlotter

# Logging
from pts.core.basics.log import setup_log
setup_log("DEBUG")

size = (4,4)

# Setup the figure
figure = plt.figure(figsize=size)
plt.clf()
ax = figure.gca()

# No axes for the main figure
ax.set_axis_off()

#im = np.arange(100)
#im.shape = 10, 10

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

frames = []

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

        plot.axes.set_adjustable('box-forced')

        #im = np.arange(100)
        xsize = int(np.random.uniform(50,200))
        ysize = int(np.random.uniform(50,200))
        #im.shape = xsize, ysize
        print(xsize, ysize)
        im = np.random.rand(xsize, ysize)
        frame = Frame(im)
        frames.append(frame)

        #grid[i].imshow(im)  # The AxesGrid object work as a list of axes.
        plotting.plot_frame(frame, axes=plot.axes, return_image=True, return_normalization=True)

        # Add the label
        plot.axes.text(0.95, 0.95, "text", color='white', transform=plot.axes.transAxes, fontsize=10, va="top", ha="right")

plt.show()
plt.close()

#plotter = StandardImageGridPlotter()
#for index, frame in enumerate(frames):
#    plotter.add_frame(frame, str(index))
#plotter.run()


plotter = ResidualImageGridPlotter()
#plotter.config.distributions = True
plotter.config.max_nrows = 3
plotter.config.ngrids = 2

# Have 5 frames
frames = frames[:5]

for index, frame in enumerate(frames):
    observation = frame
    model = frame * (1+np.random.rand(observation.ysize, observation.xsize))
    name = str(index)
    # Add row
    plotter.add_row(observation, model, name, with_residuals=True)
    #plotter.add_row(None, None, name, with_residuals=True)

plotter.run()

