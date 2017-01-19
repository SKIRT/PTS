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
from ..tools.logging import log
from ..basics.configurable import Configurable
from ..data.transmission import TransmissionCurve
from ...modeling.core.emissionlines import EmissionLines
from ..basics.filter import identifiers, Filter
from ..basics.range import RealRange

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

    def __init__(self, config=None):

        """
        This function ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(TransmissionPlotter, self).__init__(config)

        # The title
        self.title = None

        # Output
        self.output = None

        # The different curves
        self.curves = OrderedDict()

        # The wavelengths
        self.wavelengths = []

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

    # -----------------------------------------------------------------

    def set_title(self, title):

        """
        This function ...
        :param title:
        :return:
        """

        self.title = title

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

    def add_wavelength(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        # Add the wavelength
        self.wavelengths.append(wavelength)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Make the plot
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

        # Load the images (from config or input kwargs)
        if "curves" in kwargs: self.curves = kwargs.pop("curves")
        if "filters" in kwargs:
            for fltr in kwargs.pop("filters"):
                curve = TransmissionCurve.from_filter(fltr)
                self.add_transmission_curve(curve, fltr.description())
        elif self.config.filters is not None:
            for fltr in self.config.filters:
                curve = TransmissionCurve.from_filter(fltr)
                self.add_transmission_curve(curve, fltr.description())

        # Add emission lines
        if self.config.emission:
            emission_lines = EmissionLines()
            for line in emission_lines:
                center = line.center * Unit("micron")
                self.add_wavelength(center)

        # No curves
        if len(self.curves) == 0:
            for spec in identifiers:
                fltr = Filter(spec)
                curve = TransmissionCurve.from_filter(fltr)
                self.add_transmission_curve(curve, fltr.description())

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
        for label in self.curves:
            if min_wavelength is None or self.curves[label].min_wavelength < min_wavelength: min_wavelength = self.curves[label].min_wavelength

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

    def plot(self):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Plotting the transmission plot ...")

        # Create the figure
        self._figure = plt.figure(figsize=self.config.plot.figsize)
        plt.clf()

        # Set axes limits
        plt.xlim(self.min_wavelength, self.max_wavelength)
        plt.ylim(self.min_transmission, self.max_transmission*1.05)

        plt.xscale('log')

        plt.tick_params(labelsize=17)

        # Get the color map
        cm = plt.get_cmap(self.colormap)

        cNorm = colors.Normalize(vmin=0, vmax=len(self.curves) - 1.)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

        # Sort the labels based on the peak wavelength
        sorted_labels = sorted(self.curves.keys(), key=lambda label: self.curves[label].peak_wavelength)

        min_wavelength = float("inf")
        max_wavelength = 0.0

        # Plot the transmission curves
        counter = 0
        for label in sorted_labels:

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
        for wavelength in self.wavelengths:
            plt.axvline(wavelength.to("micron").value, color="0.8")

        # Get axes
        ax = plt.gca()

        # Set background
        #set_background(plt.gcf(), plt.gca(), self.config.plot) # it crashes with this

        # Set the plot title
        if self.config.plot.add_titles and self.title is not None:
            #title = "Transmission curves"
            title = "\n".join(wrap(self.title, 60))
            plt.suptitle(self.title, fontsize=self.config.plot.title_fontsize)

        # Format the axis ticks and create a grid
        ticks = RealRange(min_wavelength, max_wavelength).log(len(self.curves)*2, fancy=True)
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks)

        # Tick formatter
        #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))

        # Set ticks fontsize
        plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)
        plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=self.config.plot.ticks_fontsize)

        # Set grid
        set_grid(self.config.plot)

        # Set the borders
        set_borders(ax, self.config.plot)

        # Labels
        plt.xlabel('$\lambda [\mu m]$', fontsize=self.config.plot.label_fontsize)
        plt.ylabel('$T_\lambda$', fontsize=self.config.plot.label_fontsize)

        # Add the legend
        add_legend(ax, self.config.plot, "Filters")

        # Finish the plot
        self.finish_plot()

    # -----------------------------------------------------------------

    def finish_plot(self):

        """
        This function ...
        :return:
        """

        # plt.tight_layout()

        # Debugging
        if type(self.output).__name__ == "BytesIO": log.debug("Saving the SED plot to a buffer ...")
        elif self.output is None: log.debug("Showing the transmission plot ...")
        else: log.debug("Saving the transmission plot to " + str(self.output) + " ...")

        # Save the figure
        if self.output is not None: plt.savefig(self.output, bbox_inches='tight', pad_inches=0.25, transparent=self.transparent, format=self.format)
        else: plt.show()
        plt.close()

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
