#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.transmission Contains the TransmissionPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from textwrap import wrap
import matplotlib.pyplot as plt
from collections import OrderedDict
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.configurable import Configurable
from ..data.transmission import TransmissionCurve
from ..basics.emissionlines import EmissionLines
from ..filter.filter import parse_filter
from ..filter.broad import BroadBandFilter
from ..filter.narrow import NarrowBandFilter
from ..basics.range import RealRange
from ..basics.plot import MPLFigure
from ..tools import filesystem as fs
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

rc('text', usetex=True)

# -----------------------------------------------------------------

line_styles = ['-', '--', '-.', ':']
filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd']
pretty_colors = ["r", "dodgerblue", "purple", "darkorange", "lawngreen", "yellow", "darkblue", "teal", "darkgreen", "lightcoral", "crimson", "saddlebrown"]

# -----------------------------------------------------------------

class TransmissionPlotter(Configurable):
    
    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(TransmissionPlotter, self).__init__(*args, **kwargs)

        # The title
        self.title = None

        # Output
        self.output = None

        # The different curves
        self.curves = OrderedDict()

        # The wavelengths
        self.wavelengths = []
        self.wavelength_labels = dict()

        # The axis limits
        self._min_wavelength = None
        self._max_wavelength = None
        self._min_transmission = None
        self._max_transmission = None

        # The figure
        self._figure = None

        # Properties
        self.colormap = "rainbow"  # or "nipy_spectral"
        self.format = None
        self.transparent = False

        # The plot
        self.plt = None

    # -----------------------------------------------------------------

    def set_title(self, title):

        """
        This function ...
        :param title:
        :return:
        """

        self.title = title

    # -----------------------------------------------------------------

    def add_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Parse the filter
        fltr = parse_filter(fltr)

        # Broad band filter
        if isinstance(fltr, BroadBandFilter):

            # Create transmission curve
            curve = TransmissionCurve.from_filter(fltr)

            # Add the transmission curve
            self.add_transmission_curve(curve, fltr.description())

        # Narrow band filter
        elif isinstance(fltr, NarrowBandFilter):

            # Add wavelength
            self.add_wavelength(fltr.wavelength, fltr.description())

        # Invalid
        else: raise ValueError("A Filter object must be passed")

    # -----------------------------------------------------------------

    def add_transmission_curve(self, transmission_curve, label):

        """
        This function ...
        :param transmission_curve:
        :param label:
        :return:
        """

        # Add the transmission curve
        self.curves[label] = transmission_curve

    # -----------------------------------------------------------------

    @property
    def ncurves(self):

        """
        This function ...
        :return: 
        """

        return len(self.curves)

    # -----------------------------------------------------------------

    @property
    def has_curves(self):

        """
        This function ...
        :return: 
        """

        return self.ncurves > 0

    # -----------------------------------------------------------------

    def add_wavelength(self, wavelength, label=None):

        """
        This function ...
        :param wavelength:
        :param label:
        :return:
        """

        # Add the wavelength
        self.wavelengths.append(wavelength)

        # Set the label if a label is given
        if label is not None: self.wavelength_labels[len(self.wavelengths)-1] = label

    # -----------------------------------------------------------------

    @property
    def nwavelengths(self):

        """
        This function ...
        :return: 
        """

        return len(self.wavelengths)

    # -----------------------------------------------------------------

    @property
    def has_wavelengths(self):

        """
        This function ...
        :return: 
        """

        return self.nwavelengths > 0

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Write
        if self.config.write: self.write()

        # 2. Make the plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TransmissionPlotter, self).setup(**kwargs)

        # Check kwargs input
        if "curves" in kwargs: self.curves = kwargs.pop("curves")
        if "filters" in kwargs:
            for fltr in kwargs.pop("filters"): self.add_filter(fltr)
        elif not self.has_curves and self.config.filters is not None:
            for fltr in self.config.filters: self.add_filter(fltr)

        # Add emission lines
        if self.config.emission:
            emission_lines = EmissionLines()
            for line in emission_lines:
                center = line.center * Unit("micron")
                self.add_wavelength(center, label=line.label)

        # Add wavelengths from kwargs
        if "wavelengths" in kwargs:
            wavelengths = kwargs.pop("wavelengths")
            if isinstance(wavelengths, dict):
                for label in wavelengths: self.add_wavelength(wavelengths[label], label=label)
            else:
                for wavelength in wavelengths: self.add_wavelength(wavelength)

        # Wavelengths are specified in the configuration
        elif not self.has_wavelengths and self.config.wavelengths is not None:
            for wavelength in self.config.wavelengths:
                self.add_wavelength(wavelength)

        # Set properties
        if "title" in kwargs: self.title = kwargs.pop("title")
        else: self.title = self.config.title

        # Set limits
        self._min_wavelength = kwargs.pop("min_wavelength", None)
        self._max_wavelength = kwargs.pop("max_wavelength", None)
        self._min_transmission = kwargs.pop("min_transmission", None)
        self._max_transmission = kwargs.pop("max_transmission", None)

        # Set output
        if "output" in kwargs: self.output = kwargs.pop("output")
        else: self.output = self.config.output

        # Initialize the plot
        self.plt = MPLFigure(size=self.figsize)

    # -----------------------------------------------------------------

    @property
    def figsize(self):

        """
        This function ...
        :return:
        """

        return (self.config.plot.xsize, self.config.plot.ysize)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the transmission plotter ...")

        # Set default values for all attributes
        self.title = None
        self.curves = OrderedDict()
        self.wavelengths = []
        self._min_wavelength = None
        self._max_wavelength = None
        self._min_transmission = None
        self._max_transmission = None
        self._figure = None
        self.colormap = "rainbow"
        self.format = None
        self.transparent = False

    # -----------------------------------------------------------------

    @property
    def min_wavelength(self):

        """
        This function ...
        :return:
        """

        if self._min_wavelength is not None: return self._min_wavelength

        min_wavelength = None

        # Check the curves
        for label in self.curves:
            if min_wavelength is None or self.curves[label].min_wavelength < min_wavelength: min_wavelength = self.curves[label].min_wavelength

        # Check the wavelengths
        for wavelength in self.wavelengths:
            if min_wavelength is None or wavelength.to("micron").value < min_wavelength: min_wavelength = wavelength.to("micron").value

        # Return the minimum wavelength
        return min_wavelength

    # -----------------------------------------------------------------

    @property
    def max_wavelength(self):

        """
        This function ...
        :return:
        """

        if self._max_wavelength is not None: return self._max_wavelength

        max_wavelength = None
        for label in self.curves:
            if max_wavelength is None or self.curves[label].max_wavelength > max_wavelength: max_wavelength = self.curves[label].max_wavelength

        # CHeck the wavelengths
        for wavelength in self.wavelengths:
            if max_wavelength is None or wavelength.to("micron").value > max_wavelength: max_wavelength = wavelength.to("micron").value

        # Return the maximum wavelength
        return max_wavelength

    # -----------------------------------------------------------------

    @property
    def min_transmission(self):

        """
        This function ...
        :return:
        """

        if self._min_transmission is not None: return self._min_transmission

        min_transmission = None
        for label in self.curves:
            if min_transmission is None or self.curves[label].min_transmission < min_transmission: min_transmission = self.curves[label].min_transmission

        return min_transmission

    # -----------------------------------------------------------------

    @property
    def max_transmission(self):

        """
        This function ...
        :return:
        """

        if self._max_transmission is not None: return self._max_transmission

        max_transmission = None
        for label in self.curves:
            if max_transmission is None or self.curves[label].max_transmission > max_transmission: max_transmission = self.curves[label].max_transmission

        return max_transmission

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

        # Write the transmission data
        self.write_data()

    # -----------------------------------------------------------------

    @lazyproperty
    def sorted_labels(self):

        """
        This function ...
        :return: 
        """

        return sorted(self.curves.keys(), key=lambda label: self.curves[label].peak_wavelength)

    # -----------------------------------------------------------------

    def write_data(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the transmission data ...")

        # Determine the path for the data file
        path = fs.join(self.config.path, "transmission.dat")

        # Initialize list to contain the data (titles and columns)
        data = []

        # Loop over the labels in a sorted order
        for label in self.sorted_labels:

            curve = self.curves[label]

            wavelengths = curve.wavelengths(unit="micron", add_unit=False)
            transmissions = curve.transmissions()

            # Add the title
            data.append(label)

            # Add the columns
            data.append(wavelengths)
            data.append(transmissions)

        # Write the data
        fs.write_multi_data(path, *data)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the transmission plot ...")

        # Set axes limits
        plt.xlim(self.min_wavelength, self.max_wavelength)
        plt.ylim(self.min_transmission, self.max_transmission*1.05)

        plt.xscale('log')

        plt.tick_params(labelsize=17)

        # Get the color map
        cm = plt.get_cmap(self.colormap)

        cNorm = colors.Normalize(vmin=0, vmax=len(self.curves) - 1.)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

        min_wavelength = float("inf")
        max_wavelength = 0.0

        # Plot the transmission curves
        counter = 0
        for label in self.sorted_labels:

            curve = self.curves[label]

            wavelengths = curve.wavelengths(unit="micron", add_unit=False)
            transmissions = curve.transmissions()

            colorVal = scalarMap.to_rgba(counter)

            # Plot the curve
            plt.fill(wavelengths, transmissions, label=label, color=colorVal, alpha=0.5, linewidth=self.config.plot.linewidth)

            minw = np.min(wavelengths)
            maxw = np.max(wavelengths)
            if minw < min_wavelength: min_wavelength = minw
            if maxw > max_wavelength: max_wavelength = maxw

            counter += 1

        # Plot the wavelengths
        for index, wavelength in enumerate(self.wavelengths):

            label = self.wavelength_labels[index] if index in self.wavelength_labels else None
            plt.axvline(wavelength.to("micron").value, color="0.8", label=label)

        # Get axes
        ax = plt.gca()

        # Set background
        #set_background(plt.gcf(), plt.gca(), self.config.plot) # it crashes with this

        # Set the plot title
        if self.config.plot.add_titles and self.title is not None:
            #title = "Transmission curves"
            title = "\n".join(wrap(self.title, 60))
            plt.suptitle(self.title, fontsize=self.config.plot.title_fontsize)

        # Set the ticks
        set_ticks(ax, self.config.plot, RealRange(min_wavelength, max_wavelength), len(self.curves) * 2)

        # Set grid
        set_grid(self.config.plot)

        # Set the borders
        set_borders(ax, self.config.plot)

        # Set the labels
        self.plt.set_labels('$\lambda [\mu m]$', '$T_\lambda$', fontsize=self.config.plot.label_fontsize)

        # Add the legend
        add_legend(ax, self.config.plot, "Filters")

        # Finish
        self.plt.finish(self.output)

# -----------------------------------------------------------------

def set_borders(ax, config):

    """
    This function ...
    :param ax:
    :param config:
    :return:
    """

    # Set border width
    if config.add_borders: [i.set_linewidth(config.borderwidth) for i in ax.spines.itervalues()]

    # Remove borders
    else:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis=u'both', which=u'both', length=0)

# -----------------------------------------------------------------

def add_legend(ax, config, legend_title=None):

    """
    This function ...
    :param ax:
    :param config:
    :param legend_title:
    :return:
    """

    #if nlegends > 1: percentage = 25.
    #else: percentage = 10.
    percentage = 20.

    # Shrink current axis's height by a certain percentage on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * percentage/100., box.width, box.height * (100-percentage)/100.])

    # Plot legend
    legend_title = r"\underline{" + legend_title + "}"
    legend = plt.legend(loc="upper center", title=legend_title, bbox_to_anchor=(0.5, -0.25), fancybox=False, shadow=False, ncol=4)

    frame = legend.get_frame()

    # Set legend frame color and line width
    if config.add_legend_border:
        frame.set_linewidth(config.legend_borderwidth)
    else:
        frame.set_linewidth(0)
    frame.set_edgecolor(config.legend_bordercolor)

    # Set background color
    frame.set_facecolor('0.85')
    legend.legendPatch.set_alpha(0.75)

    index = 0

    # Move to foreground
    legend.set_zorder(100 + index)

    # Set fontsize
    plt.setp(legend.get_title(), fontsize=str(config.legend_title_fontsize))

# -----------------------------------------------------------------

def set_ticks(ax, config, x_range, nxticks):

    """
    This function ...
    :param config:
    :param x_range:
    :param nxticks:
    :return:
    """

    # Format the axis ticks and create a grid
    ticks = x_range.log(nxticks)
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks)

    # Tick formatter
    ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))

    # Set ticks fontsize
    plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=config.ticks_fontsize)
    plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=config.ticks_fontsize)

# -----------------------------------------------------------------

def set_grid(config):

    """
    This function ...
    :param config:
    :return:
    """

    if config.add_grid: plt.grid(linewidth=config.grid_linewidth, linestyle=config.grid_linestyle, color=config.grid_color)

# -----------------------------------------------------------------

def set_background(fig, ax, config):

    """
    This function ...
    :param fig:
    :param ax:
    :param config:
    :return:
    """

    if config.transparent_background:

        # Set transparent background
        for item in [fig, ax]:
            item.patch.set_visible(False)

# -----------------------------------------------------------------
