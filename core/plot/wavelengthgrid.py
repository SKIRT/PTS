#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.wavelengthgrid Contains the WavelengthGridPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import math
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.legend import Legend
from collections import defaultdict
import matplotlib.patches as patches
from scipy.interpolate import interp1d

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.configurable import Configurable
from ..tools import filesystem as fs
from ..tools.serialization import load_list, load_dict
from ..basics.plot import MPLFigure, BokehFigure, mpl, bokeh
from ..tools.utils import lazyproperty
from ..basics.map import Map
from ..filter.broad import BroadBandFilter
from ..filter.narrow import NarrowBandFilter
from ..basics import containers
from ..basics.range import RealRange, QuantityRange
from ..basics.emissionlines import EmissionLine
from ...magic.tools.wavelengths import find_single_wavelength_range
from ..basics.plot import dark_pretty_colors
from ..tools import sequences
from ..filter.filter import parse_filter
from ..tools import types
from ..tools import nr
from ..tools.stringify import tostr

# -----------------------------------------------------------------

line_styles = ['-', '--', '-.', ':']
filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd']

# -----------------------------------------------------------------

def get_sed_template(name, **kwargs):

    """
    This function ...
    :param name:
    :return:
    """

    from ...modeling.core.mappings import create_mappings_sed
    from ...modeling.core.bruzualcharlot import create_bruzual_charlot_sed

    # Mappings
    if name == "mappings": return create_mappings_sed(**kwargs)

    # Bruzual-Charlot
    elif name == "bruzual_charlot": return create_bruzual_charlot_sed(**kwargs)

    # Not recognized
    else: raise ValueError("Template not recognized")

# -----------------------------------------------------------------

def plot_wavelength_grid(wavelength_grid, label, **kwargs):

    """
    This function ...
    :param wavelength_grid:
    :param label:
    :param kwargs:
    :return:
    """

    # Create the wavelength plotter
    plotter = WavelengthGridPlotter()

    # Add the wavelength grid
    plotter.add_wavelength_grid(wavelength_grid, label, **kwargs)

    # Plot
    plotter.run()

# -----------------------------------------------------------------

def plot_wavelength_grids(grids, **kwargs):

    """
    This function ...
    :param grids:
    :param kwargs:
    :return:
    """

    # Create the wavelength grid plotter
    plotter = WavelengthGridPlotter()

    # Loop over the grids
    for label in grids: plotter.add_wavelength_grid(grids[label], label, **kwargs)

    # Plot
    plotter.run()

# -----------------------------------------------------------------

class WavelengthGridPlotter(Configurable):
    
    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(WavelengthGridPlotter, self).__init__(*args, **kwargs)

        # The figure
        self.figure = None

        # The title
        self.title = None

        # The output path
        self.out_path = None

        # Reference wavelength grids
        self.references = OrderedDict()

        # The wavelength grids
        self.grids = OrderedDict()

        # SEDs to plot together with the wavelength grids
        self.seds = OrderedDict()

        # Emission lines
        self.emission_lines = []

        # Single wavelengths
        self.wavelengths = []

        # Single wavelengths to be removed from the complete grid
        self.remove = []

        # Wavelength ranges
        self.wavelength_ranges = []

        # Filters
        self.filters = OrderedDict()

        # The axis limits
        self.min_wavelength = None
        self.max_wavelength = None

        # Subplots
        self.main = None
        self.separate = []
        self.groups = OrderedDict()
        self.lines = None
        self.individual = None
        self.complete = None
        self.differences = None
        self.refs = []
        self.delta = None
        self.residuals = None

        # Legends of the main plot
        self.grids_legend = None
        self.filters_legend = None
        self.seds_legend = None

        # The grid scatter handles
        self.grid_handles = []

        # Shared labels with their colors
        self.shared_labels = dict()

        # Grid colour iterator
        self.grids_colors = None

        # Have filters been added from config?
        self._filters_added = False

        # Backgrounds
        self.backgrounds = []

    # -----------------------------------------------------------------

    def set_title(self, title):

        """
        This function ...
        :param title:
        :return:
        """

        self.title = title

    # -----------------------------------------------------------------

    def add_reference_grid(self, grid, label, pointsize=None, linewidth=None, linealpha=None, color=None, in_legend=False):

        """
        This function ...
        :param grid:
        :param label:
        :param pointsize:
        :param linewidth:
        :param linealpha:
        :param color:
        :param in_legend:
        :return:
        """

        # Properties
        props = Map()
        props.grid = grid
        props.pointsize = pointsize
        props.linewidth = linewidth
        props.linealpha = linealpha
        props.color = color
        props.in_legend = in_legend

        # Add reference grid
        self.references[label] = props

    # -----------------------------------------------------------------

    def add_wavelength_grid(self, grid, label, pointsize=None, linewidth=None, linealpha=None, color=None,
                            in_legend=True, y_value=None, copy_grid=True, shared_label=None, add_tag=False,
                            separate=None, plot_on_filter=None, group=None):

        """
        This function ...
        :param grid:
        :param label:
        :param pointsize:
        :param linewidth:
        :param linealpha:
        :param color:
        :param in_legend:
        :param y_value:
        :param copy_grid:
        :param shared_label:
        :param add_tag:
        :param separate:
        :param plot_on_filter:
        :param group:
        :return:
        """

        # Make copy?
        if copy_grid: grid = copy.deepcopy(grid)

        # Properties
        props = Map()
        props.grid = grid
        props.pointsize = pointsize
        props.linewidth = linewidth
        props.linealpha = linealpha
        props.color = color
        props.in_legend = in_legend
        props.y_value = y_value
        props.shared_label = shared_label
        props.add_tag = add_tag
        props.separate = separate
        props.plot_on_filter = plot_on_filter
        props.group = group

        # Add grid
        self.grids[label] = props

    # -----------------------------------------------------------------

    def add_grid_from_file(self, path, label=None, pointsize=None, linewidth=None, linealpha=None, color=None,
                           in_legend=True, y_value=None, shared_label=None, add_tag=False, separate=None,
                           plot_on_filter=None, group=None):

        """
        Thisf unction ...
        :param path:
        :param label:
        :param pointsize:
        :param linewidth:
        :param linealpha:
        :param color:
        :param in_legend:
        :param y_value:
        :param shared_label:
        :param add_tag:
        :param separate:
        :param plot_on_filter:
        :param group:
        :return:
        """

        from ..simulation.wavelengthgrid import WavelengthGrid
        if label is None: label = fs.strip_extension(fs.name(path))
        grid = WavelengthGrid.from_file(path)
        self.add_wavelength_grid(grid, label, pointsize=pointsize, linewidth=linewidth, linealpha=linealpha,
                                 color=color, in_legend=in_legend, y_value=y_value, copy_grid=False,
                                 shared_label=shared_label, add_tag=add_tag, separate=separate,
                                 plot_on_filter=plot_on_filter, group=group)

    # -----------------------------------------------------------------

    def add_wavelengths(self, wavelengths, label, unit=None, pointsize=None, linewidth=None, linealpha=None,
                        color=None, in_legend=True, y_value=None, shared_label=None, add_tag=False, separate=None,
                        plot_on_filter=None, group=None):

        """
        This function ...
        :param wavelengths:
        :param label:
        :param unit:
        :param pointsize:
        :param linewidth:
        :param linealpha:
        :param color:
        :param in_legend:
        :param y_value:
        :param shared_label:
        :param add_tag:
        :param separate:
        :param plot_on_filter:
        :param group:
        :return:
        """

        from ..simulation.wavelengthgrid import WavelengthGrid
        grid = WavelengthGrid.from_wavelengths(wavelengths, unit=unit, sort=True)
        self.add_wavelength_grid(grid, label, pointsize=pointsize, linewidth=linewidth, linealpha=linealpha,
                                 color=color, in_legend=in_legend, y_value=y_value, copy_grid=False,
                                 shared_label=shared_label, add_tag=add_tag, separate=separate,
                                 plot_on_filter=plot_on_filter, group=group)

    # -----------------------------------------------------------------

    def add_wavelength(self, wavelength, label=None, colour=None, linewidth=None, group=None, connect=False,
                       connect_in_style=False):

        """
        This function ...
        :param wavelength:
        :param label:
        :param colour:
        :param linewidth:
        :param group:
        :param connect:
        :param connect_in_style:
        :return:
        """

        # Set properties
        props = Map()
        props.wavelength = wavelength
        props.label = label
        props.color = colour
        props.linewidth = linewidth
        props.group = group
        props.connect = connect
        props.connect_in_style = connect_in_style

        # Add
        self.wavelengths.append(props)

    # -----------------------------------------------------------------

    def remove_wavelength(self, wavelength):

        """
        This function ...
        :param wavelength:
        :return:
        """

        # Remove
        self.remove.append(wavelength)

    # -----------------------------------------------------------------

    def remove_wavelengths(self, wavelengths):

        """
        This function ...
        :param wavelengths:
        :return:
        """

        for wavelength in wavelengths: self.remove_wavelength(wavelength)

    # -----------------------------------------------------------------

    def remove_wavelengths_in_range(self, wavelength_range, except_grids=None):

        """
        This function ...
        :param wavelength_range:
        :param except_grids:
        :return:
        """

        self.remove_wavelengths_between(wavelength_range.min, wavelength_range.max, except_grids=except_grids)

    # -----------------------------------------------------------------

    def remove_wavelengths_between(self, min_wavelength, max_wavelength, except_grids=None):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param except_grids:
        :return:
        """

        # Loop over the wavelength grids
        for label in self.grids:

            #print(label)

            # Except?
            if except_grids is not None and label in except_grids: continue

            # Get the wavelength points
            grid = self.grids[label].grid

            # Get the wavelengths between the min and max
            wavelengths = grid.wavelengths(self.config.wavelength_unit, min_wavelength=min_wavelength, max_wavelength=max_wavelength, inclusive=False)

            # Remove wavelengths
            self.remove_wavelengths(wavelengths)

    # -----------------------------------------------------------------

    def add_wavelength_point(self, label, wavelength):

        """
        This function ...
        :param label:
        :param wavelength:
        :return:
        """

        self.grids[label].add_point(wavelength)

    # -----------------------------------------------------------------

    def add_sed(self, sed, label):

        """
        This function ...
        :param sed:
        :param label:
        :return:
        """

        self.seds[label] = sed

    # -----------------------------------------------------------------

    def add_seds(self, seds):

        """
        This function ...
        :param seds:
        :return:
        """

        for label in seds: self.add_sed(seds[label], label)

    # -----------------------------------------------------------------

    def add_emission_line(self, line):

        """
        This function ...
        :param line:
        :return:
        """

        self.emission_lines.append(line)

    # -----------------------------------------------------------------

    def add_filter(self, fltr, label=None):

        """
        This function ...
        :param fltr:
        :param label:
        """

        if label is None: label = str(fltr)
        self.filters[label] = fltr

    # -----------------------------------------------------------------

    def add_wavelength_range(self, wavelength_range, label=None):

        """
        This function ...
        :param wavelength_range:
        :param label:
        :return:
        """

        props = Map()
        props.range = wavelength_range
        props.label = label

        # Add
        self.wavelength_ranges.append(props)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Make the plot
        self.plot()

        # 3. Write
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(WavelengthGridPlotter, self).setup(**kwargs)

        # Set the title
        self.title = kwargs.pop("title", None)

        # Set the output path
        self.out_path = kwargs.pop("output", None)
        self.out_path = kwargs.pop("output_path", self.out_path)
        if self.out_path is None and self.config.output is not None: self.out_path = fs.join(self.config.output, self.config.filename + "." + self.config.format)
        #print("OUT PATH:", self.out_path, self.config.output)

        # Set the axis limits
        self.min_wavelength = kwargs.pop("min_wavelength", self.config.min_wavelength)
        self.max_wavelength = kwargs.pop("max_wavelength", self.config.max_wavelength)

        # Load grids
        if self.config.grids is not None:
            for path in self.config.grids:
                self.add_grid_from_file(path)

        # Add filters: IMPORTANT, BEFORE SUBGRIDS (for add_filter_wavelengths)
        #if self.config.add_filters: self.add_filters()
        # NOW ADDING THE FILTERS IS ADDED TO 'sorted_filter_labels' and 'categorized_filter_labels' methods,
        # TO MAKE SURE THE FILTERS ARE LOADED BEFORE E.G. LOAD_ELEMENTS IS CALLED FROM EXTERNAL
        # NEW: added if not yet added
        if not self._filters_added and self.config.add_filters: self.add_filters()

        # Load grids from subgrid generation
        if self.config.load_subgrids: self.load_subgrids()

        # Generate subgrids
        elif self.config.create_subgrids: self.create_subgrids()

        # Are there grids
        if not self.has_grids: raise RuntimeError("No wavelength grids are added")

        # Add emission lines
        if self.config.add_lines: self.create_lines()

        # Add SEDs
        if self.config.add_seds: self.create_seds()

        # Add SEDs from input
        if kwargs.get("seds", None) is not None: self.add_seds(kwargs.pop("seds"))

        # Add regimes
        if self.config.add_regimes: self.add_regimes()

        # Create the plot
        if self.config.library == mpl: self.figure = MPLFigure(size=self.figsize)
        elif self.config.library == bokeh: self.figure = BokehFigure()
        else: raise ValueError("Invalid libary: " + self.config.library)

        # Set the 'show' flag
        if self.config.show is None:
            if self.out_path is not None: self.config.show = False
            else: self.config.show = True

        # Initialize plot
        self.initialize_plot()

        # Check
        if self.config.lines_in_group is not None and self.config.lines_in_group not in self.group_names: raise ValueError("Group '" + self.config.lines_in_group + "' does not exist")

        # Check
        if self.config.lines_in_group is not None and self.config.separate_lines: raise ValueError("Cannot plot lines in group and plot lines separately")

        # Check
        if self.config.plot_differences and not self.has_single_reference_grid: raise ValueError("Has to have one reference grid to plot differences")

    # -----------------------------------------------------------------

    @property
    def figsize(self):

        """
        This function ...
        :return:
        """

        return (self.config.plot.xsize, self.config.plot.ysize)

    # -----------------------------------------------------------------

    @property
    def min_y(self):

        """
        This function ...
        :return:
        """

        return 0

    # -----------------------------------------------------------------

    @property
    def max_y(self):

        """
        This function ...
        :return:
        """

        return 3

    # -----------------------------------------------------------------

    @lazyproperty
    def x_limits(self):

        """
        This function ...
        :return:
        """

        return [self.min_wavelength.to(self.config.wavelength_unit).value, self.max_wavelength.to(self.config.wavelength_unit).value]

    # -----------------------------------------------------------------

    @lazyproperty
    def main_y_limits(self):

        """
        This function ...
        :return:
        """

        return (self.min_y, self.max_y)

    # -----------------------------------------------------------------

    @lazyproperty
    def main_y_span(self):

        """
        This function ...
        :return:
        """

        return self.max_y - self.min_y

    # -----------------------------------------------------------------

    @lazyproperty
    def min_y_delta(self):

        """
        This function ...
        :return:
        """

        radius = 0.5 * self.complete_grid_logdelta_span
        center = self.complete_grid_center_logdelta
        min_logdelta = center - 1.1 * radius
        return min_logdelta

    # -----------------------------------------------------------------

    @lazyproperty
    def max_y_delta(self):

        """
        This function ...
        :return:
        """

        radius = 0.5 * self.complete_grid_logdelta_span
        center = self.complete_grid_center_logdelta
        max_logdelta = center + 1.1 * radius
        return max_logdelta

    # -----------------------------------------------------------------

    @lazyproperty
    def delta_y_limits(self):

        """
        Thisfunction ...
        :return:
        """

        return (self.min_y_delta, self.max_y_delta)

    # -----------------------------------------------------------------

    @lazyproperty
    def delta_y_span(self):

        """
        Thisfunction ...
        :return:
        """

        return self.max_y_delta - self.min_y_delta

    # -----------------------------------------------------------------

    @property
    def min_y_residuals(self):

        """
        This function ...
        :return:
        """

        return -100.

    # -----------------------------------------------------------------

    @property
    def max_y_residuals(self):

        """
        This function ...
        :return:
        """

        return 100.

    # -----------------------------------------------------------------

    @lazyproperty
    def residuals_y_limits(self):

        """
        This function ...
        :return:
        """

        return (self.min_y_residuals, self.max_y_residuals)

    # -----------------------------------------------------------------

    @lazyproperty
    def residuals_y_span(self):

        """
        This function ...
        :return:
        """

        return self.max_y_residuals - self.min_y_residuals

    # -----------------------------------------------------------------

    @property
    def min_y_complete_grid(self):

        """
        This function ...
        :return:
        """

        return 0.

    # -----------------------------------------------------------------

    @property
    def max_y_complete_grid(self):

        """
        This function ...
        :return:
        """

        return 1.

    # -----------------------------------------------------------------

    @lazyproperty
    def complete_y_limits(self):

        """
        This function ...
        :return:
        """

        return [self.min_y_complete_grid, self.max_y_complete_grid]

    # -----------------------------------------------------------------

    @property
    def nrows(self):

        """
        This function ...
        :return:
        """

        # Main and delta
        total = 2

        # Complete grid plot
        if self.config.separate_complete_grid: total += 1

        # Reference grid plots
        if self.has_reference_grids: total += self.nreference_grids

        # Separate grid plots
        if self.has_separate_grids: total += self.nseparate_grids

        # Group plots
        if self.has_groups: total += self.ngroups

        # Individual wavelengths plot
        if self.has_wavelengths_not_in_group and self.config.group_wavelengths: total += 1

        # Lines
        if self.has_lines and self.config.separate_lines: total += 1

        # Differences
        if self.config.plot_differences: total += 1

        # Residuals
        if self.has_seds and self.config.plot_residuals: total += 1

        # Return the number of rows
        return total

    # -----------------------------------------------------------------

    @property
    def min_y_references(self):

        """
        This function ...
        :return:
        """

        return 0.

    # -----------------------------------------------------------------

    @property
    def max_y_references(self):

        """
        This function ...
        :return:
        """

        return 1.

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_y_limits(self):

        """
        This function ...
        :return:
        """

        return [self.min_y_references, self.max_y_references]

    # -----------------------------------------------------------------

    @property
    def min_y_separate(self):

        """
        This function ...
        :return:
        """

        return 0.

    # -----------------------------------------------------------------

    @property
    def max_y_separate(self):

        """
        This function ...
        :return:
        """

        return 1.

    # -----------------------------------------------------------------

    @lazyproperty
    def separate_y_limits(self):

        """
        This function ...
        :return:
        """

        return [self.min_y_separate, self.max_y_separate]

    # -----------------------------------------------------------------

    @property
    def min_y_groups(self):

        """
        This function ...
        :return:
        """

        return 0.

    # -----------------------------------------------------------------

    @property
    def max_y_groups(self):

        """
        This function ...
        :return:
        """

        return 1.

    # -----------------------------------------------------------------

    @property
    def groups_y_limits(self):

        """
        This function ...
        :return:
        """

        return [self.min_y_groups, self.max_y_groups]

    # -----------------------------------------------------------------

    @property
    def min_y_individual(self):

        """
        This function ...
        :return:
        """

        return 0.

    # -----------------------------------------------------------------

    @property
    def max_y_individual(self):

        """
        This function ...
        :return:
        """

        return 1.

    # -----------------------------------------------------------------

    @property
    def individual_y_limits(self):

        """
        This function ...
        :return:
        """

        return [self.min_y_individual, self.max_y_individual]

    # -----------------------------------------------------------------

    @property
    def min_y_separate_lines(self):

        """
        This function ...
        :return:
        """

        return 0.

    # -----------------------------------------------------------------

    @property
    def max_y_separate_lines(self):

        """
        This function ...
        :return:
        """

        return 1.

    # -----------------------------------------------------------------

    @property
    def separate_lines_y_limits(self):

        """
        This function ...
        :return:
        """

        return [self.min_y_separate_lines, self.max_y_separate_lines]

    # -----------------------------------------------------------------

    @property
    def min_y_differences(self):

        """
        This function ...
        :return:
        """

        return 0.

    # -----------------------------------------------------------------

    @property
    def max_y_differences(self):

        """
        This function ...
        :return:
        """

        return 1.

    # -----------------------------------------------------------------

    @property
    def differences_y_limits(self):

        """
        This function ...
        :return:
        """

        return [self.min_y_differences, self.max_y_differences]

    # -----------------------------------------------------------------

    @lazyproperty
    def y_limits(self):


        """
        This function ...
        :return:
        """

        # Main plot
        limits = [self.main_y_limits]

        # Separate grid plots
        if self.has_separate_grids:
            for _ in range(self.nseparate_grids): limits.append(self.separate_y_limits)

        # Group plots
        if self.has_groups:
            for _ in range(self.ngroups): limits.append(self.groups_y_limits)

        # Lines
        if self.has_lines and self.config.separate_lines: limits.append(self.separate_lines_y_limits)

        # Individual wavelengths plot
        if self.has_wavelengths_not_in_group and self.config.group_wavelengths: limits.append(self.individual_y_limits)

        # Complete grid plot
        if self.config.separate_complete_grid: limits.append(self.complete_y_limits)

        # Reference grid plots
        if self.has_reference_grids:
            for _ in range(self.nreference_grids): limits.append(self.reference_y_limits)

        # Differences plot
        if self.config.plot_differences: limits.append(self.differences_y_limits)

        # Delta plot
        limits.append(self.delta_y_limits)

        # Residuals
        if self.has_seds and self.config.plot_residuals: limits.append(self.residuals_y_limits)

        # Return the limits
        return limits

    # -----------------------------------------------------------------

    @lazyproperty
    def height_ratios(self):

        """
        This function ...
        :return:
        """

        # Main plot
        ratios = [4]

        # Separate grid plots
        if self.has_separate_grids:
            for _ in range(self.nseparate_grids): ratios.append(0.2)

        # Group plots
        if self.has_groups:
            for _ in range(self.ngroups): ratios.append(0.4)

        # Lines
        if self.has_lines and self.config.separate_lines: ratios.append(0.4)

        # Individual wavelengths plot
        if self.has_wavelengths_not_in_group and self.config.group_wavelengths: ratios.append(0.3)

        # Complete grid plot
        if self.config.separate_complete_grid: ratios.append(0.3)

        # Reference grid plots
        if self.has_reference_grids:
            for _ in range(self.nreference_grids): ratios.append(0.3)

        # Differences plot
        if self.config.plot_differences: ratios.append(0.2)

        # Delta plot
        ratios.append(1)

        # Residuals plot
        if self.has_seds and self.config.plot_residuals: ratios.append(1)

        # Return the list of ratios
        return ratios

    # -----------------------------------------------------------------

    @property
    def x_label(self):

        """
        This function ...
        :return:
        """

        return '$\lambda/\mu$m'

    # -----------------------------------------------------------------

    @lazyproperty
    def y_labels(self):

        """
        This function ...
        :return:
        """

        labels = [None] * self.nrows

        # Delta plot
        labels[self.delta_plot_index] = r"$\Delta\lambda\,(\mathrm{dex})$"

        # Residuals
        if self.has_seds and self.config.plot_residuals: labels[self.residuals_plot_index] = r"$\epsilon\,(\%)$"

        # Return the labels
        return labels

    # -----------------------------------------------------------------

    @property
    def x_scale(self):

        """
        This function ...
        :return:
        """

        return "log"

    # -----------------------------------------------------------------

    @property
    def nseparate_grids(self):

        """
        Thisn function ...
        :return:
        """

        total = 0
        for label in self.grids:
            props = self.grids[label]
            separate = props.separate

            # Don't count as separate when in group
            group = props.group
            if group is not None: continue

            if separate: total += 1
            elif separate is None and self.config.separate_grids: total += 1

        return total

    # -----------------------------------------------------------------

    @property
    def has_separate_grids(self):

        """
        Thisn function ...
        :return:
        """

        return self.nseparate_grids > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def group_names(self):

        """
        This function ...
        :return:
        """

        names = set()

        # Loop over the grid props
        for label in self.grids:

            props = self.grids[label]
            group = props.group
            if group is None: continue
            names.add(group)

        # Loop over the individual wavelength props
        for props in self.wavelengths:
            group = props.group
            if group is None: continue
            names.add(group)

        # Return the group names
        return list(names)

    # -----------------------------------------------------------------

    @property
    def ngroups(self):

        """
        The number of groups of grids that go on the same row
        :return:
        """

        return len(self.group_names)

    # -----------------------------------------------------------------

    @property
    def has_groups(self):

        """
        This function ...
        :return:
        """

        return self.ngroups > 0

    # -----------------------------------------------------------------

    def initialize_plot(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Initializing the plot ...")

        # Create subplots
        plots = self.figure.create_column(self.nrows, share_axis=True, height_ratios=self.height_ratios, x_label=self.x_label, y_labels=self.y_labels, x_limits=self.x_limits, x_scale=self.x_scale, y_limits=self.y_limits, x_log_scalar=True)

        # Set main plot
        self.main = plots[0]

        # Set separate grid plots
        if self.has_separate_grids: self.separate = plots[1:self.nseparate_grids + 1]

        # Set group plots
        if self.has_groups:
            start = self.nseparate_grids + 1
            group_plots = plots[start:start+self.ngroups]
            for name, plot in zip(self.group_names, group_plots): self.groups[name] = plot

        # Set lines plot
        if self.has_lines and self.config.separate_lines: self.lines = plots[self.nseparate_grids + 1 + self.ngroups]

        # Set individual wavelengths plot
        if self.has_wavelengths_not_in_group and self.config.group_wavelengths:
            if self.has_lines and self.config.separate_lines: index = self.nseparate_grids + 2 + self.ngroups
            else: index = self.nseparate_grids + 1 + self.ngroups
            self.individual = plots[index]

        # Set complete grid plot
        if self.config.separate_complete_grid:
            if self.has_seds and self.config.plot_residuals: offset = -1
            else: offset = 0
            if self.has_reference_grids and self.config.plot_differences: index = -2 - self.nreference_grids - 1 + offset
            else: index = -1 - self.nreference_grids - 1 + offset
            self.complete = plots[index]

        # Set reference grid plots
        if self.has_reference_grids:
            if self.has_seds and self.config.plot_residuals: offset = -1
            else: offset = 0
            if self.has_reference_grids and self.config.plot_differences:
                start = -2 - self.nreference_grids + offset
                end = -2 + offset
            else:
                start = -1-self.nreference_grids + offset
                end = -1 + offset
            self.refs = plots[start:end]

        # Set differences plot
        if self.has_reference_grids and self.config.plot_differences: self.differences = plots[self.differences_plot_index]

        # Set delta plot
        self.delta = plots[self.delta_plot_index]

        # Set residuals plot
        if self.has_seds and self.config.plot_residuals: self.residuals = plots[self.residuals_plot_index]

        # Set background
        self.set_backgrounds()

    # -----------------------------------------------------------------

    @property
    def differences_plot_index(self):

        """
        This function ...
        :return:
        """

        if self.has_seds and self.config.plot_residuals: return -3
        else: return -2

    # -----------------------------------------------------------------

    @property
    def delta_plot_index(self):

        """
        This function ...
        :return:
        """

        if self.has_seds and self.config.plot_residuals: return -2
        else: return -1

    # -----------------------------------------------------------------

    @property
    def residuals_plot_index(self):

        """
        This function ..
        :return:
        """

        return -1

    # -----------------------------------------------------------------

    def set_backgrounds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting backgrounds ...")

        # Fill backgrounds
        for props in self.backgrounds:

            min_wavelength = props.min_wavelength
            max_wavelength = props.max_wavelength
            group = props.group
            color = props.color
            alpha = props.alpha
            hatch = props.hatch
            fill = props.fill

            #x = (min_wavelength.to(self.config.wavelength_unit).value, max_wavelength.to(self.config.wavelength_unit).value)
            #print(x, self.min_y_groups, self.max_y_groups, alpha, color)

            width = max_wavelength.to(self.config.wavelength_unit).value - min_wavelength.to(self.config.wavelength_unit).value

            if group is not None:
                xy = (min_wavelength.to(self.config.wavelength_unit).value, self.min_y_groups)
                height = self.max_y_groups - self.min_y_groups
                self.groups[group].axes.add_patch(patches.Rectangle(xy, width, height, alpha=alpha, facecolor=color, linewidth=0, hatch=hatch, fill=fill))
            else:
                xy = (min_wavelength.to(self.config.wavelength_unit).value, self.min_y)
                height = self.max_y - self.min_y
                self.main_axes.add_patch(patches.Rectangle(xy, width, height, alpha=alpha, facecolor=color, linewidth=0, hatch=hatch, fill=fill))

    # -----------------------------------------------------------------

    def load_subgrids(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading from subgrids wavelength grid generation ...")

        # Set paths
        if self.config.subgrids_path is not None:
            fixed_path = fs.join(self.config.subgrids_path, "fixed_grid.dat")
            new_path = fs.join(self.config.subgrids_path, "new.dat")
            replaced_path = fs.join(self.config.subgrids_path, "replaced.dat")
            filter_wavelengths_path = fs.join(self.config.subgrids_path, "filter_wavelengths.dat")
            line_wavelengths_path = fs.join(self.config.subgrids_path, "line_wavelengths.dat")
        else:
            fixed_path = self.input_path_file("fixed_grid.dat")
            new_path = self.input_path_file("new.dat")
            replaced_path = self.input_path_file("replaced.dat")
            filter_wavelengths_path = self.input_path_file("filter_wavelengths.dat")
            line_wavelengths_path = self.input_path_file("line_wavelengths.dat")

        # Add the subgrid files
        if self.config.subgrids_path is not None: self.add_subgrids_in_path(self.config.subgrids_path)
        else: self.add_subgrids_in_path(self.input_path)

        # Check which files are present
        has_fixed = fs.is_file(fixed_path)
        has_new = fs.is_file(new_path)
        has_replaced = fs.is_file(replaced_path)
        has_filters = fs.is_file(filter_wavelengths_path)
        has_lines = fs.is_file(line_wavelengths_path)

        # Fixed?
        if has_fixed: self.add_fixed_grid_from_file(fixed_path)

        # New?
        if has_new: self.add_new_wavelengths_from_file(new_path)

        # Filter
        if has_filters: filter_wavelengths = self.add_filter_wavelengths_from_file(filter_wavelengths_path)
        else: filter_wavelengths = None

        # Replaced?
        if has_replaced: self.add_replaced_wavelengths_from_file(replaced_path, connect=has_filters, filter_wavelengths=filter_wavelengths)

        # Emission lines?
        if has_lines: self.add_line_wavelengths_from_file(line_wavelengths_path)

    # -----------------------------------------------------------------

    def create_subgrids(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Generating a subgrids wavelength grid ...")

        from ...modeling.build.wavelengthgrid import WavelengthGridBuilder

        # Create wavelength grid builder
        builder = WavelengthGridBuilder()

        # Set options
        builder.config.add_emission_lines = self.config.wg.add_emission_lines
        builder.config.emission_lines = self.config.wg.emission_lines
        builder.config.min_wavelength = self.config.wg.range.min
        builder.config.max_wavelength = self.config.wg.range.max
        builder.config.check_filters = self.config.wg.check_filters
        builder.config.npoints = self.config.wg.npoints
        builder.config.adjust_to = self.config.wg.adjust_to
        builder.config.fixed = self.config.wg.fixed

        # Don't write or plot anything
        builder.config.write = False
        builder.config.plot = False

        # Run the wavelength grid builder
        builder.run()

        # Add the elements from the wavelength grid builder
        self.add_elements(builder.subgrids, fixed=builder.fixed, filter_wavelengths=builder.filter_wavelengths,
                          replaced=builder.replaced, new=builder.new, line_wavelengths=builder.line_wavelengths)

    # -----------------------------------------------------------------

    def add_elements(self, subgrids, fixed=None, fixed_grid=None, filter_wavelengths=None, replaced=None, new=None,
                     line_wavelengths=None):

        """
        This function ...
        :param subgrids:
        :param fixed:
        :param fixed_grid:
        :param filter_wavelengths:
        :param replaced:
        :param new:
        :param line_wavelengths:
        :return:
        """

        # Add the subgrids
        self.add_subgrids(subgrids)

        # Check what is present
        has_fixed = fixed is not None and len(fixed) > 0
        has_fixed_grid = fixed_grid is not None and len(fixed_grid) > 0
        has_filters = filter_wavelengths is not None and len(filter_wavelengths) > 0
        has_replaced = replaced is not None and len(replaced) > 0
        has_new = new is not None and len(new) > 0
        has_line_wavelengths = line_wavelengths is not None and len(line_wavelengths) > 0

        # Grid of fixed wavelengths
        if has_fixed: self.add_fixed_wavelengths(fixed)
        if has_fixed_grid:
            if has_fixed: raise ValueError("Cannot pass fixed wavelengths and fixed grid")
            self.add_fixed_grid(fixed_grid)

        # Get filter wavelengths
        if has_filters: self.add_filter_wavelengths(filter_wavelengths)

        # Get replaced
        if has_replaced: self.add_replaced_wavelengths(replaced, connect=has_filters, filter_wavelengths=filter_wavelengths)

        # Get new
        if has_new: self.add_new_wavelengths(new)

        # Get line wavelengths
        if has_line_wavelengths: self.add_line_wavelengths(line_wavelengths)

    # -----------------------------------------------------------------

    def add_subgrids(self, subgrids):

        """
        This function ...
        :param subgrids:
        :return:
        """

        # Loop over the subgrids
        for name in subgrids:

            # Debugging
            log.debug("Adding '" + name + "' subgrid ...")

            # Add grid
            subgrid = subgrids[name]
            self.add_wavelength_grid(subgrid, label=name, group="subgrids")

    # -----------------------------------------------------------------

    def add_subgrid_from_file(self, name, path):

        """
        This function ...
        :param name:
        :param path:
        :return:
        """

        # Debugging
        log.debug("Adding '" + name + "' subgrid ...")

        # Add grid
        self.add_grid_from_file(path, label=name, group="subgrids")

    # -----------------------------------------------------------------

    def add_subgrids_in_path(self, input_path):

        """
        This function ...
        :param input_path:
        :return:
        """

        # Loop over the subgrid files
        for path, name in fs.files_in_path(input_path, extension="dat", startswith="subgrid_", returns=["path", "name"]):

            # Get subgrid label
            subgrid = name.split("subgrid_")[1]

            # Check
            if subgrid not in self.config.subgrids: continue

            # Add subgrid
            self.add_subgrid_from_file(subgrid, path)

    # -----------------------------------------------------------------

    def add_fixed_wavelengths(self, fixed):

        """
        This function ...
        :param fixed:
        :return:
        """

        from ..simulation.wavelengthgrid import WavelengthGrid

        # Debugging
        log.debug("Adding fixed wavelengths ...")

        # Add
        fixed_grid = WavelengthGrid.from_wavelengths(fixed)
        self.add_wavelength_grid(fixed_grid, label="fixed", linewidth=1., pointsize=20)

    # -----------------------------------------------------------------

    def add_fixed_grid(self, fixed_grid):

        """
        This function ...
        :param fixed_grid:
        :return:
        """

        # Debugging
        log.debug("Adding fixed wavelengths ...")

        # Add
        self.add_wavelength_grid(fixed_grid, label="fixed", linewidth=1., pointsize=20)

    # -----------------------------------------------------------------

    def add_fixed_grid_from_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        from ..simulation.wavelengthgrid import WavelengthGrid

        # Load
        fixed_grid = WavelengthGrid.from_file(path)
        if len(fixed_grid) == 0:
            log.warning("No points in the fixed wavelength grid: skipping ...")
            return

        # Add
        self.add_fixed_grid(fixed_grid)

    # -----------------------------------------------------------------

    def add_filter_wavelengths(self, filter_wavelengths):

        """
        This function ...
        :param filter_wavelengths:
        :return:
        """

        # Debugging
        log.debug("Adding filter wavelengths ...")

        # Sort the filters
        sorted_filters = sorted(filter_wavelengths.keys(), key=lambda fltr: fltr.wavelength.to("micron").value)
        nfilters = len(filter_wavelengths)

        # Set colours if needed
        cm = self.filters_colormap
        ncolours = nfilters
        cNorm = colors.Normalize(vmin=0, vmax=ncolours - 1.)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
        filter_colors = [scalarMap.to_rgba(i) for i in range(ncolours)]
        filter_colors = iter(filter_colors) # make iterable

        # Get the names of all filter wavelength grids
        filter_wavelength_grid_names = [str(fltr) for fltr in sorted_filters]

        remove_betweens = []

        # Loop over the filters
        for fltr in sorted_filters:

            # Get the wavelengths
            wavelengths = filter_wavelengths[fltr]
            nwavelengths = len(wavelengths)
            filter_name = str(fltr)

            add_tag = False

            # Set color for filter
            if self.has_colour_for_filter(fltr): color = self.get_filter_colour(fltr)
            elif self.has_filter_colours:
                color = "k" # black: there are colours for other filters
                add_tag = True
            else:
                color = next(filter_colors)
                add_tag = True

            # Remove wavelengths, if more than one
            if nwavelengths > 1:
                min_wavelength = min(wavelengths)
                max_wavelength = max(wavelengths)
                remove_betweens.append((min_wavelength, max_wavelength))

            # Add the wavelengths as a grid
            self.add_wavelengths(wavelengths, label=filter_name, in_legend=False, color=color, add_tag=add_tag,
                                 y_value=0.7, separate=False, plot_on_filter=fltr)

        # Remove wavelengths of subgrids between filter wavelengths
        hatch = '///'
        for min_wavelength, max_wavelength in remove_betweens:
            # Remove wavelengths between
            self.remove_wavelengths_between(min_wavelength, max_wavelength, except_grids=filter_wavelength_grid_names)

            # Colour background
            self.colour_background_in_group(min_wavelength, max_wavelength, "subgrids", "grey", alpha=0.2, hatch=hatch, fill=False)
            if hatch == '///': hatch = "\\\\\\"
            elif hatch == "\\\\\\": hatch = "///"
            else: raise ValueError("Something went wrong")

    # -----------------------------------------------------------------

    def colour_background_in_group(self, min_wavelength, max_wavelength, group, color, alpha=None, hatch=None, fill=True):

        """
        This function ...
        :param min_wavelength:
        :param max_wavelength:
        :param group:
        :param color:
        :param alpha:
        :param hatch:
        :param fill:
        :return:
        """

        props = Map()
        props.min_wavelength = min_wavelength
        props.max_wavelength = max_wavelength
        props.group = group
        props.color = color
        props.alpha = alpha
        props.hatch = hatch
        props.fill = fill

        #print(props)

        # Add
        self.backgrounds.append(props)

    # -----------------------------------------------------------------

    def add_filter_wavelengths_from_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Load and add
        filter_wavelengths = load_dict(path)
        self.add_filter_wavelengths(filter_wavelengths)
        return filter_wavelengths

    # -----------------------------------------------------------------

    def add_replaced_wavelengths(self, replaced, connect=False, filter_wavelengths=None):

        """
        This function ...
        :param replaced:
        :param connect:
        :param filter_wavelengths:
        :return:
        """

        # Debugging
        log.debug("Adding replaced wavelengths ...")

        # Set wavelength range for each filter
        if filter_wavelengths is not None:
            filter_ranges = dict()
            for fltr in filter_wavelengths:
                wavelengths = filter_wavelengths[fltr]
                nwavelengths = len(wavelengths)
                if nwavelengths == 1: continue
                #print(wavelengths)
                wavelength_range = QuantityRange.limits(wavelengths)
                filter_ranges[fltr] = wavelength_range
        else: filter_ranges = None

        # Add
        for old, replacement in replaced:

            if connect:
                if filter_ranges is None: connect_line = True
                elif in_some_range(old, filter_ranges): connect_line = True
                else: connect_line = False
            else: connect_line = False

            # Add the old wavelength and the replacement wavelength
            self.add_wavelength(old, colour="red", linewidth=0.5, group="adjusted", connect=connect_line, connect_in_style=True)
            self.add_wavelength(replacement, colour="green", linewidth=0.5, group="adjusted")

            # Remove the old wavelength from the complete grid
            self.remove_wavelength(old)

    # -----------------------------------------------------------------

    def add_replaced_wavelengths_from_file(self, path, connect=False, filter_wavelengths=None):

        """
        This function ...
        :param path:
        :param connect:
        :param filter_wavelengths:
        :return:
        """

        # Load and add
        replaced = load_list(path)
        self.add_replaced_wavelengths(replaced, connect=connect, filter_wavelengths=filter_wavelengths)

    # -----------------------------------------------------------------

    def add_new_wavelengths(self, new):

        """
        This function ...
        :param new:
        :return:
        """

        # Debugging
        log.debug("Adding new wavelengths ...")

        # Add wavelengths
        self.add_wavelengths(new, label="new", color="b", linewidth=0.5, group="adjusted")

    # -----------------------------------------------------------------

    def add_new_wavelengths_from_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Load and add
        new = load_list(path)
        self.add_new_wavelengths(new)

    # -----------------------------------------------------------------

    def add_line_wavelengths(self, line_wavelengths):

        """
        This function ...
        :param line_wavelengths: 
        :return: 
        """

        # Debugging
        log.debug("Adding emission line wavelengths ...")

        # Loop over the lines
        for identifier in line_wavelengths:

            wavelengths = line_wavelengths[identifier]
            if identifier[0] is None: continue
            nwavelengths = len(wavelengths)
            label = identifier[0] + identifier[1]

            # Remove wavelengths of subgrids between line wavelengths
            if nwavelengths > 1:

                min_wavelength = min(wavelengths)
                max_wavelength = max(wavelengths)

                #print(label, min_wavelength, max_wavelength)

                # Remove wavelengths between min and max of line wavelengths
                self.remove_wavelengths_between(min_wavelength, max_wavelength)

                # Colour background
                self.colour_background_in_group(min_wavelength, max_wavelength, "subgrids", "grey", alpha=0.5)

            # Add line wavelengths
            self.add_wavelengths(wavelengths, label=label, color="grey", y_value=0.6, shared_label="lines", separate=False, group="lines", linewidth=0.3)

    # -----------------------------------------------------------------

    def add_line_wavelengths_from_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Load and add
        line_wavelengths = load_dict(path)
        self.add_line_wavelengths(line_wavelengths)

    # -----------------------------------------------------------------

    def add_filters(self):

        """
        This function ....
        :return:
        """

        # Inform the user
        log.info("Adding filters ...")

        self._filters_added = True

        # Loop over the filters
        for fltr in self.config.filters:

            # Debugging
            log.debug("Adding the '" + str(fltr) + "' filter ...")

            # Add
            self.add_filter(fltr)

    # -----------------------------------------------------------------

    def create_lines(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Adding emission lines ...")

        # Add
        for line_id in self.config.lines:

            # Create line
            line = EmissionLine.from_string(line_id)

            # Add line
            self.add_emission_line(line)

    # -----------------------------------------------------------------

    def create_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the suer
        log.info("Creating SED templates ...")

        # Loop over the template names
        for name in self.config.seds:

            # MAPPINGS
            if name == "mappings":

                # Debugging
                log.debug("Creating MAPPINGS SED template ...")

                properties = dict()
                properties["metallicity"] = self.config.metallicity
                properties["compactness"] = self.config.compactness
                properties["pressure"] = self.config.pressure
                properties["covering_factor"] = self.config.covering_factor

                # Set label
                label = "MAPPINGS"

                sed = get_sed_template(name, **properties)
                self.add_sed(sed, label=label)

            # Stellar Bruzual Charlot
            elif name == "bruzual_charlot":

                # Debugging
                log.debug("Creating Bruzual-Charlot SED templates ...")

                # Loop over the ages
                for age in self.config.ages:

                    properties = dict()
                    properties["metallicity"] = self.config.metallicity
                    properties["age"] = age

                    #label = name + "_" + str(age).replace(" ", "")
                    label = "Bruzual-Charlot " + str(age)
                    sed = get_sed_template(name, **properties)
                    self.add_sed(sed, label=label)

            # Invalid
            else: raise ValueError("Invalid SED template name")

    # -----------------------------------------------------------------

    def add_regimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding wavelength regimes ...")

        # Loop over the regimes
        for regime in self.config.regimes:

            # Debugging
            log.debug("Adding '" + regime + "' wavelength regime ...")

            # Determine wavelength range
            wavelength_range, string = find_single_wavelength_range(regime, return_string=True, only_subregime_string=self.config.only_subregimes)
            wavelength_range.adjust_inwards(self.grid_wavelength_range_with_unit)
            self.add_wavelength_range(wavelength_range, label=string)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the wavelength grid plotter ...")

        # Set default values for all attributes
        self.figure = None
        self.title = None
        self.grids = OrderedDict()
        self.seds = OrderedDict()
        self.emission_lines = []
        self.filters = OrderedDict()
        self.min_wavelength = None
        self.max_wavelength = None
        #self._figure = None
        #self.colormap = "rainbow"
        #self.format = None
        #self.transparent = False
        self.main = None
        self.delta = None

    # -----------------------------------------------------------------

    @property
    def nremove(self):

        """
        This function ...
        :return:
        """

        return len(self.remove)

    # -----------------------------------------------------------------

    @property
    def has_remove(self):

        """
        This function ...
        :return:
        """

        return self.nremove > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def remove_no_unit(self):

        """
        This function ...
        :return:
        """

        return [wav.to(self.config.wavelength_unit).value for wav in self.remove]

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelengths_no_unit(self):

        """
        This function ...
        :return:
        """

        return [props.wavelength.to(self.config.wavelength_unit).value for props in self.wavelengths]

    # -----------------------------------------------------------------

    @lazyproperty
    def all_grid_wavelengths(self):

        """
        This function ...
        :return:
        """

        all_grid_wavelengths = []

        # Loop over the wavelength grids
        for label in self.grids:

            grid = self.grids[label].grid
            wavelengths = grid.wavelengths(self.config.wavelength_unit, asarray=True)

            # Add wavelengths but not those to remove from complete grid
            for wavelength in wavelengths:
                all_grid_wavelengths.append(wavelength)

        # Add single wavelengths
        if self.has_wavelengths:
            # Loop over the wavelengths
            for wavelength in self.wavelengths_no_unit: all_grid_wavelengths.append(wavelength)

        # Remove wavelengths?
        # SOME WAVELENGTHS IN SELF.WAVELENGTHS CAN ALSO BE ADDED TO REMOVE
        if self.has_remove:
            remove_indices = []
            for index in range(len(all_grid_wavelengths)):
                wavelength = all_grid_wavelengths[index]
                if wavelength in self.remove_no_unit: remove_indices.append(index)
            sequences.remove_indices(all_grid_wavelengths, remove_indices)

        # Sort all wavelengths and create array
        all_grid_wavelengths = np.array(sorted(all_grid_wavelengths))
        return all_grid_wavelengths

    # -----------------------------------------------------------------

    @lazyproperty
    def log_all_grid_wavelengths(self):

        """
        This function ...
        :return:
        """

        return np.log10(self.all_grid_wavelengths)

    # -----------------------------------------------------------------

    @lazyproperty
    def complete_grid(self):

        """
        This function ...
        :return:
        """

        from ..simulation.wavelengthgrid import WavelengthGrid
        return WavelengthGrid.from_wavelengths(self.all_grid_wavelengths, unit=self.config.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def complete_grid_deltas(self):

        """
        This function ...
        :return:
        """

        return self.complete_grid.deltas(unit=self.config.wavelength_unit, asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def complete_grid_logdeltas(self):

        """
        This function ...
        :return:
        """

        return self.complete_grid.logdeltas(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def complete_grid_min_logdelta(self):

        """
        This function ...
        :return:
        """

        return 0.

    # -----------------------------------------------------------------

    @lazyproperty
    def complete_grid_max_logdelta(self):

        """
        This function ...
        :return:
        """

        return np.max(self.complete_grid_logdeltas)

    # -----------------------------------------------------------------

    @lazyproperty
    def complete_grid_logdelta_span(self):

        """
        This function ...
        :return:
        """

        return self.complete_grid_max_logdelta - self.complete_grid_min_logdelta

    # -----------------------------------------------------------------

    @lazyproperty
    def complete_grid_center_logdelta(self):

        """
        Thisn function ...
        :return:
        """

        return 0.5 * (self.complete_grid_min_logdelta + self.complete_grid_max_logdelta)

    # -----------------------------------------------------------------

    @lazyproperty
    def all_grid_nwavelengths(self):

        """
        Thisfunction ...
        :return:
        """

        return len(self.all_grid_wavelengths)

    # -----------------------------------------------------------------

    @lazyproperty
    def grid_min_wavelength(self):

        """
        This function ...
        :return:
        """

        return min(self.all_grid_wavelengths)

    # -----------------------------------------------------------------

    @lazyproperty
    def grid_max_wavelength(self):

        """
        This function ...
        :return:
        """

        return max(self.all_grid_wavelengths)

    # -----------------------------------------------------------------

    @lazyproperty
    def grid_wavelength_range(self):

        """
        This function ...
        :return:
        """

        return RealRange(self.grid_min_wavelength, self.grid_max_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def grid_wavelength_range_with_unit(self):

        """
        This function ...
        :return:
        """

        return QuantityRange(self.grid_min_wavelength, self.grid_max_wavelength, unit=self.config.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_wavelengths(self):

        """
        This function ...
        :return:
        """

        wavelengths = dict()

        # Loop over the SEDs
        for label in self.seds:

            # Get the wavelengths array
            warray = self.seds[label].wavelengths(unit=self.config.wavelength_unit, asarray=True)

            # Add to the dictionary
            wavelengths[label] = warray

        # Return the dictionary
        return wavelengths

    # -----------------------------------------------------------------

    @lazyproperty
    def log_sed_wavelengths(self):

        """
        This function ...
        :return:
        """

        log_wavelengths = dict()

        # Loop over the SEDs
        for label in self.seds: log_wavelengths[label] = np.log10(self.sed_wavelengths[label])

        # Return
        return log_wavelengths

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_fluxes(self):

        """
        This function ...
        :return:
        """

        fluxes = dict()

        # Loop over the SEDs
        for label in self.seds:

            # Get the fluxes
            farray = self.seds[label].normalized_photometry(method="max", asarray=True)

            # Add to the dictionary
            fluxes[label] = farray

        # Return the dictionary
        return fluxes

    # -----------------------------------------------------------------

    @lazyproperty
    def log_sed_fluxes(self):

        """
        This function ...
        :return:
        """

        log_fluxes = dict()

        # Loop over the SEDs
        for label in self.seds: log_fluxes[label] = np.log10(self.sed_fluxes[label])

        # Return
        return log_fluxes

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_grid_fluxes(self):

        """
        This function ...
        :return:
        """

        fluxes = dict()

        # Loop over the SEDs
        for label in self.seds:

            # Get the fluxes
            farray = self.seds[label].normalized_photometry(method="max", asarray=True)

            # Calculate interpolated (resampled) values
            interpolated = nr.resample_log_log(self.all_grid_wavelengths, self.sed_wavelengths[label], farray)

            # Add to the dictionary
            fluxes[label] = interpolated

        # Return the dictionary
        return fluxes

    # -----------------------------------------------------------------

    @lazyproperty
    def log_sed_grid_fluxes(self):

        """
        This function ...
        :return:
        """

        log_fluxes = dict()

        # Loop over the SEDs
        for label in self.seds: log_fluxes[label] = np.log10(self.sed_grid_fluxes[label])

        # Return
        return log_fluxes

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_interpolated_fluxes(self):

        """
        This function ...
        :return:
        """

        fluxes = dict()

        # Loop over the SEDs
        for label in self.seds:

            # Create interpolation function from grid fluxes in log space
            f2 = interp1d(self.log_all_grid_wavelengths, self.log_sed_grid_fluxes[label], kind=self.config.interpolation_method, bounds_error=False, fill_value=float("NaN"))

            # Evaluate the function at the original wavelengths of the SED
            #wavelengths = self.sed_wavelengths[label]
            log_wavelengths = self.log_sed_wavelengths[label]
            log_interpolated = f2(log_wavelengths)
            interpolated = 10**log_interpolated

            # Set the interpolated fluxes
            fluxes[label] = interpolated

        # Return the interpolated fluxes
        return fluxes

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_residuals(self):

        """
        This function ...
        :return:
        """

        residuals = dict()

        # Loop over the SEDs
        for label in self.seds:

            # Get the original fluxes and original wavelengths
            #wavelengths = self.sed_wavelengths[label]
            fluxes = self.sed_fluxes[label]

            # Get the interpolated fluxes
            interpolated = self.sed_interpolated_fluxes[label]

            # Calculate residuals
            res = (interpolated - fluxes) / fluxes

            # Add residuals
            residuals[label] = res

        # Return
        return residuals

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_residuals_percentual(self):

        """
        This function ...
        :return:
        """

        residuals = dict()
        for label in self.seds: residuals[label] = self.sed_residuals[label] * 100
        return residuals

    # -----------------------------------------------------------------

    @lazyproperty
    def min_sed_residuals(self):

        """
        This function ...
        :return:
        """

        minima = dict()

        # Loop over the SEDs
        for label in self.seds: minima[label] = np.nanmin(self.sed_residuals[label])

        # Return the minima
        return minima

    # -----------------------------------------------------------------

    @lazyproperty
    def min_sed_residuals_percentual(self):

        """
        This function ...
        :return:
        """

        minima = dict()
        for label in self.seds: minima[label] = self.min_sed_residuals[label] * 100.
        return minima

    # -----------------------------------------------------------------

    @lazyproperty
    def max_sed_residuals(self):

        """
        This function ...
        :return:
        """

        maxima = dict()

        # Loop over the SEDs
        for label in self.seds: maxima[label] = np.nanmax(self.sed_residuals[label])

        # Return the maxima
        return maxima

    # -----------------------------------------------------------------

    @lazyproperty
    def max_sed_residuals_percentual(self):

        """
        This function ...
        :return:
        """

        maxima = dict()
        for label in self.seds: maxima[label] = self.max_sed_residuals[label] * 100.
        return maxima

    # -----------------------------------------------------------------

    @lazyproperty
    def min_residual(self):

        """
        This function ...
        :return:
        """

        minimum = None

        # Loop over the SEDs
        for label in self.seds:
            if minimum is None or self.min_sed_residuals[label] < minimum:
                minimum = self.min_sed_residuals[label]

        return minimum

    # -----------------------------------------------------------------

    @lazyproperty
    def max_residual(self):

        """
        This function ...
        :return:
        """

        maximum = None

        # Loop over the SEDs
        for label in self.seds:
            if maximum is None or self.max_sed_residuals[label] > maximum:
                maximum = self.max_sed_residuals[label]

        return maximum

    # -----------------------------------------------------------------

    @property
    def min_residual_percentual(self):

        """
        This function ...
        :return:
        """

        return self.min_residual * 100.

    # -----------------------------------------------------------------

    @property
    def max_residual_percentual(self):

        """
        This function ...
        :return:
        """

        return self.max_residual * 100.

    # -----------------------------------------------------------------

    # @lazyproperty
    # def delta_wavelengths(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     return self.all_grid_wavelengths[1:]
    #
    # # -----------------------------------------------------------------
    #
    # @lazyproperty
    # def deltas(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     return np.log10(self.all_grid_wavelengths[1:]) - np.log10(self.all_grid_wavelengths[:-1])
    #
    # # -----------------------------------------------------------------
    #
    # @lazyproperty
    # def max_delta(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     return np.max(self.deltas)
    #
    # # -----------------------------------------------------------------
    #
    # @lazyproperty
    # def min_delta(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     return np.min(self.deltas)

    # -----------------------------------------------------------------

    @property
    def nseds(self):

        """
        Thisn function ...
        :return:
        """

        return len(self.seds)

    # -----------------------------------------------------------------

    @property
    def has_seds(self):

        """
        Thisn function ...
        :return:
        """

        return self.nseds > 0

    # -----------------------------------------------------------------

    @property
    def nfilters(self):

        """
        This function ...
        :return:
        """

        return len(self.filters)

    # -----------------------------------------------------------------

    @property
    def has_filters(self):

        """
        This function ...
        :return:
        """

        return self.nfilters > 0

    # -----------------------------------------------------------------

    @property
    def has_lines(self):

        """
        This function ...
        :return:
        """

        #return self.nlines_with_label > 0
        return self.nlines > 0 # otherwise no -nolabel- warnings

    # -----------------------------------------------------------------

    @property
    def nwavelengths(self):

        """
        This function ...
        :return:
        """

        return len(self.wavelengths)

    # -----------------------------------------------------------------

    @lazyproperty
    def nwavelengths_not_in_group(self):

        """
        This function ...
        :return:
        """

        total = 0

        # Loop over the wavelengths
        for properties in self.wavelengths:

            # Get properties
            group = properties.group
            if group is not None: continue
            total += 1

        # Return the number
        return total

    # -----------------------------------------------------------------

    @lazyproperty
    def nwavelengths_in_group(self):

        """
        This function ...
        :return:
        """

        total = 0

        # Loop over the wavelengths
        for properties in self.wavelengths:

            # Get properties
            group = properties.group
            if group is None: continue
            total += 1

        # Return the number
        return total

    # -----------------------------------------------------------------

    @property
    def has_wavelengths(self):

        """
        This function ...
        :return:
        """

        return self.nwavelengths > 0

    # -----------------------------------------------------------------

    @property
    def has_wavelengths_not_in_group(self):

        """
        This function ...
        :return:
        """

        return self.nwavelengths_not_in_group > 0

    # -----------------------------------------------------------------

    @property
    def has_wavelengths_in_group(self):

        """
        This function ...
        :return:
        """

        return self.nwavelengths_in_group > 0

    # -----------------------------------------------------------------

    @property
    def nranges(self):

        """
        This fnction ...
        :return:
        """

        return len(self.wavelength_ranges)

    # -----------------------------------------------------------------

    @property
    def has_ranges(self):

        """
        This function ...
        :return:
        """

        return self.nranges > 0

    # -----------------------------------------------------------------

    @property
    def nreference_grids(self):

        """
        This function ...
        :return:
        """

        return len(self.references)

    # -----------------------------------------------------------------

    @property
    def has_reference_grids(self):

        """
        This function ...
        :return:
        """

        return self.nreference_grids > 0

    # -----------------------------------------------------------------

    @property
    def has_single_reference_grid(self):

        """
        This function ...
        :return:
        """

        return self.nreference_grids == 1

    # -----------------------------------------------------------------

    @property
    def single_reference_grid_label(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_reference_grid: raise ValueError("Not a single reference grid")
        return self.references.keys()[0]

    # -----------------------------------------------------------------

    @property
    def single_reference_grid(self):

        """
        This function ...
        :return:
        """

        return self.references[self.single_reference_grid_label].grid

    # -----------------------------------------------------------------

    @property
    def single_reference_grid_wavelengths(self):

        """
        This function ...
        :return:
        """

        #print(self.single_reference_grid, type(self.single_reference_grid))
        return self.single_reference_grid.wavelengths(self.config.wavelength_unit, asarray=True)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelength grids ...")

        # 1. SEDs
        if self.has_seds: self.plot_seds()

        # 2. Filters
        if self.has_filters: self.plot_filters()

        # 3. Grids
        self.plot_grids()

        # 4. Plot complete grid
        if self.config.plot_complete_grid: self.plot_complete_grid()

        # 5. Plot reference grids
        if self.has_reference_grids: self.plot_reference_grids()

        # Plot difference grids
        if self.config.plot_differences: self.plot_differences()

        # 6. Plot legend
        self.plot_legend()

        # 7. Plot complete grid deltas
        self.plot_deltas()

        # 8. Emission lines
        if self.has_lines: self.plot_lines()

        # 9. Single wavelengths
        if self.has_wavelengths: self.plot_wavelengths()

        # 10. Wavelength ranges
        if self.has_ranges: self.plot_ranges()

        # 11. Residuals
        if self.has_seds and self.config.plot_residuals: self.plot_residuals()

        # 12. Finish
        self.finish_plot()

    # -----------------------------------------------------------------

    @property
    def grid_points_y(self):

        """
        This function ...
        :return:
        """

        return 0.5

    # -----------------------------------------------------------------

    @property
    def ngrids(self):

        """
        This function ..
        :return:
        """

        return len(self.grids)

    # -----------------------------------------------------------------

    @property
    def has_grids(self):

        """
        This function ...
        :return:
        """

        return self.ngrids > 0

    # -----------------------------------------------------------------

    @property
    def grid_labels_with_fixed_y(self):

        """
        This function ...
        :return:
        """

        return [label for label in self.grids if self.grids[label].y_value is not None]

    # -----------------------------------------------------------------

    @property
    def grid_labels_without_fixed_y(self):

        """
        Thisn function ...
        :return:
        """

        return [label for label in self.grids if self.grids[label].y_value is None]

    # -----------------------------------------------------------------

    @property
    def ngrids_with_fixed_y(self):

        """
        Thisfunction ...
        :return:
        """

        return len(self.grid_labels_with_fixed_y)

    # -----------------------------------------------------------------

    @property
    def ngrids_without_fixed_y(self):

        """
        This function ...
        :return:
        """

        return len(self.grid_labels_without_fixed_y)

    # -----------------------------------------------------------------

    @property
    def grid_points_y_range(self):

        """
        This function ...
        :return:
        """

        return RealRange(0.3, 0.5, inclusive=True)

    # -----------------------------------------------------------------

    @property
    def grid_y_values(self):

        """
        This function ...
        :return:
        """

        return self.grid_points_y_range.linear(self.ngrids_without_fixed_y)

    # -----------------------------------------------------------------

    @property
    def complete_grid_y_value(self):

        """
        This function ...
        :return:
        """

        return 0.2

    # -----------------------------------------------------------------

    @property
    def complete_grid_color(self):

        """
        This function ...
        :return:
        """

        return "black"

    # -----------------------------------------------------------------

    @property
    def complete_grid_pointsize(self):

        """
        This function ...
        :return:
        """

        return 20

    # -----------------------------------------------------------------

    @property
    def complete_grid_linewidth(self):

        """
        This function ...
        :return:
        """

        return 0.8

    # -----------------------------------------------------------------

    @property
    def complete_grid_linealpha(self):

        """
        This function ...
        :return:
        """

        return 0.8

    # -----------------------------------------------------------------

    @property
    def main_axes(self):

        """
        This function ...
        :return:
        """

        return self.main.axes

    # -----------------------------------------------------------------

    @property
    def delta_axes(self):

        """
        This function ...
        :return:
        """

        return self.delta.axes

    # -----------------------------------------------------------------

    @property
    def grids_legend_properties(self):

        """
        This function ...
        :return:
        """

        legend_properties = dict()

        legend_properties["loc"] = 1
        legend_properties["numpoints"] = 4
        legend_properties["scatterpoints"] = 4
        legend_properties["ncol"] = 2
        legend_properties["shadow"] = False
        legend_properties["frameon"] = True
        legend_properties["facecolor"] = None
        legend_properties["fontsize"] = "smaller"

        return legend_properties

    # -----------------------------------------------------------------

    @property
    def lines_legend_properties(self):

        """
        This function ...
        :return:
        """

        legend_properties = dict()

        legend_properties["loc"] = 1
        legend_properties["ncol"] = self.nlines_with_label
        legend_properties["shadow"] = False
        legend_properties["frameon"] = True
        legend_properties["facecolor"] = None
        legend_properties["fontsize"] = 7

        return legend_properties

    # -----------------------------------------------------------------

    @property
    def delta_linewidth(self):

        """
        Thisn function ...
        :return:
        """

        return 0.6

    # -----------------------------------------------------------------

    @property
    def tag_label_alpha(self):

        """
        This function ...
        :return:
        """

        return 0.7

    # -----------------------------------------------------------------

    def plot_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the grids ...")

        # Make iterable from distinct colors
        self.grids_colors = iter(dark_pretty_colors)

        # Get y values iterator
        y_values = iter(self.grid_y_values)

        tag_position = "above"

        # Loop over the wavelength grids
        separate_index = 0
        for label in self.grids:

            # Debugging
            log.debug("Adding '" + label + "' wavelength grid to the plot ...")

            # Get the wavelength points
            grid = self.grids[label].grid
            wavelengths = grid.wavelengths(self.config.wavelength_unit, asarray=True)
            nwavelengths = grid.nwavelengths

            # Get properties
            pointsize = self.grids[label].pointsize
            linewidth = self.grids[label].linewidth
            linealpha = self.grids[label].linealpha
            color = self.grids[label].color
            in_legend = self.grids[label].in_legend
            y_value = self.grids[label].y_value
            shared_label = self.grids[label].shared_label
            add_tag = self.grids[label].add_tag
            separate = self.grids[label].separate
            plot_on_filter = self.grids[label].plot_on_filter
            group = self.grids[label].group

            # Check
            if shared_label is not None and not in_legend: raise ValueError("in_legend is disabled but shared label is specified")

            # Set shared label
            add_label = True
            if shared_label is not None:
                if shared_label in self.shared_labels:
                    if color is None: color = self.shared_labels[shared_label]
                    elif color != self.shared_labels[shared_label]: raise ValueError("Shared label but colours are not the same")
                    add_label = False # not a new (shared) label
                else:
                    if color is None: color = next(self.grids_colors)
                    self.shared_labels[shared_label] = color
                    add_label = True # a new (shared) label

            # Get next color, if needed
            if color is None: color = next(self.grids_colors)

            # Set line width and line alpha
            if pointsize is None: pointsize = self.config.pointsize
            if linewidth is None: linewidth = self.config.linewidth
            if linealpha is None: linealpha = self.config.linealpha

            # Get next y value, if needed
            if y_value is None: y_value = next(y_values)

            # Plot points
            if add_label:
                if shared_label is not None: legend_label = shared_label
                else: legend_label = label + " (" + str(nwavelengths) + " points)"
            else: legend_label = None

            # Set separate flag
            if separate is None and self.config.separate_grids: separate = True
            used_separate = False  # can stay False when group is specified

            # Group?
            if group is not None:

                if plot_on_filter is not None: raise ValueError("Cannot plot on filter when plotting in group")

                y = [0.5 for _ in wavelengths]

                # Plot
                sc = self.groups[group].scatter(wavelengths, y, s=pointsize, marker='.', color=color, linewidths=0, label=legend_label)

            # Separate?
            elif separate:

                if plot_on_filter is not None: raise ValueError("Cannot plot on filter when plotting separately")

                y = [0.5 for _ in wavelengths]

                # Plot
                sc = self.separate[separate_index].scatter(wavelengths, y, s=pointsize, marker='.', color=color, linewidths=0, label=legend_label)

                # Set flag
                used_separate = True

            # Not separate: on main plot
            else:

                # Plot on filter transmission curve?
                if plot_on_filter is not None:

                    from ..data.transmission import TransmissionCurve

                    # Create transmission curve
                    curve = TransmissionCurve.from_filter(plot_on_filter)
                    curve.normalize(value=self.max_y_filters, method="max")

                    #print(plot_on_filter.min, plot_on_filter.max, curve.min_wavelength, curve.max_wavelength, wavelengths)

                    # Find interpolated normalized transmission value for the wavelength
                    #y = [curve.transmission_at(wavelength * self.config.wavelength_unit) for wavelength in wavelengths]
                    y = []
                    for wavelength_scalar in wavelengths:
                        wavelength = wavelength_scalar * self.config.wavelength_unit
                        if not curve.in_range(wavelength):
                            log.warning("The wavelength " + tostr(wavelength) + " is not in the range of the transmission curve of the '" + tostr(plot_on_filter) + "' filter")
                            y.append(0.)
                        else: y.append(curve.transmission_at(wavelength))

                # Set same y value for each point
                else: y = [y_value for _ in wavelengths]

                # Plot
                sc = self.main.scatter(wavelengths, y, s=pointsize, marker='.', color=color, linewidths=0, label=legend_label)

            # Plot a vertical line for each grid point
            for w in wavelengths:

                alpha = linealpha
                linestyle = "solid"
                if self.config.mark_removed and self.has_remove:
                    if w in self.remove_no_unit:
                        linestyle = "dashed"
                        alpha = 0.5 * linealpha

                if group is not None: self.groups[group].vlines(w, self.min_y_groups, self.max_y_groups, color=color, lw=linewidth, alpha=alpha, linestyles=linestyle)
                elif separate: self.separate[separate_index].vlines(w, self.min_y_separate, self.max_y_separate, color=color, lw=linewidth, alpha=alpha, linestyles=linestyle)
                else: self.main.vlines(w, self.min_y, self.max_y, color=color, lw=linewidth, alpha=alpha, linestyles=linestyle)

            # Add text
            if add_tag:

                # Determine text position
                center = grid.geometric_mean_wavelength.to(self.config.wavelength_unit).value
                if tag_position == "above":
                    if separate or group is not None: tag_y_value = 0.75
                    else: tag_y_value = y_value + 0.1
                elif tag_position == "below":
                    if separate or group is not None: tag_y_value = 0.15
                    else: tag_y_value = y_value - 0.2
                else: raise ValueError("Invalid tag position")

                # Plot text
                if group is not None: t = self.groups[group].text(center, tag_y_value, label, horizontalalignment='center', fontsize='xx-small', color=color, backgroundcolor='w', zorder=200)
                elif separate: t = self.separate[separate_index].text(center, tag_y_value, label, horizontalalignment='center', fontsize='xx-small', color=color, backgroundcolor='w', zorder=200)
                else: t = self.main.text(center, tag_y_value, label, horizontalalignment='center', fontsize='xx-small', color=color, backgroundcolor='w', zorder=200)

                # Set text box properties
                t.set_bbox(dict(color='w', alpha=self.tag_label_alpha, edgecolor='w'))

            # Add the handle
            if in_legend and add_label: self.grid_handles.append(sc)

            # Switch position
            if tag_position == "above": tag_position = "below"
            elif tag_position == "below": tag_position = "above"
            else: raise ValueError("Invalid tag position")

            # Increment
            if separate and used_separate: separate_index += 1

    # -----------------------------------------------------------------

    def plot_complete_grid(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the complete grid ...")

        # Set legend label
        legend_label = "all (" + str(self.all_grid_nwavelengths) + " points)"

        # Plot the complete grid in a separate subplot
        if self.config.separate_complete_grid:

            linewidth = self.complete_grid_linewidth
            linealpha = self.complete_grid_linealpha

            # Plot points
            complete_y = [0.5 for _ in self.all_grid_wavelengths]
            sc = self.complete.scatter(self.all_grid_wavelengths, complete_y, s=self.complete_grid_pointsize, marker='.', color=self.complete_grid_color, linewidths=0, label=legend_label)
            self.grid_handles.append(sc)

            # Add grid lines
            for w in self.all_grid_wavelengths: self.complete.vlines(w, self.min_y_complete_grid, self.max_y_complete_grid, color=self.complete_grid_color, lw=linewidth, alpha=linealpha)

        # Plot the complete grid in the main plot
        else:

            # Plot as one complete grid
            complete_y = [self.complete_grid_y_value for _ in self.all_grid_wavelengths]
            sc = self.main.scatter(self.all_grid_wavelengths, complete_y, s=self.complete_grid_pointsize, marker='.', color=self.complete_grid_color, linewidths=0, label=legend_label)
            self.grid_handles.append(sc)

    # -----------------------------------------------------------------

    def plot_reference_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the reference grids ...")

        # Loop over the wavelength grids
        index = 0
        for label in self.references:

            # Debugging
            log.debug("Adding '" + label + "' reference grid to the plot ...")

            # Get the wavelength points
            grid = self.references[label].grid
            wavelengths = grid.wavelengths(self.config.wavelength_unit, asarray=True)
            nwavelengths = grid.nwavelengths

            # Get properties
            pointsize = self.references[label].pointsize
            linewidth = self.references[label].linewidth
            linealpha = self.references[label].linealpha
            color = self.references[label].color
            in_legend = self.references[label].in_legend

            # Set color
            if color is None: color = next(self.grids_colors) #color = "black"

            # Set line width and line alpha
            if pointsize is None: pointsize = self.config.pointsize
            if linewidth is None: linewidth = self.config.linewidth
            if linealpha is None: linealpha = self.config.linealpha

            # Plot points
            #if add_label: legend_label = label + " (" + str(nwavelengths) + " points)"
            #else: legend_label = None
            legend_label = label + " (" + str(nwavelengths) + " points)"

            #print(wavelengths)
            #print(pointsize)

            # Plot points
            complete_y = [0.5 for _ in wavelengths]
            #print(complete_y)
            sc = self.refs[index].scatter(wavelengths, complete_y, s=pointsize, marker='.', color=color, linewidths=0, label=legend_label)
            if in_legend: self.grid_handles.append(sc)

            # Add grid lines
            for w in wavelengths: self.refs[index].vlines(w, self.min_y_references, self.max_y_references,
                                                         color=color, lw=linewidth,
                                                         alpha=linealpha)

            index += 1

    # -----------------------------------------------------------------

    def plot_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting differences between complete grid and reference grid ...")

        # Get reference grid wavelengths
        #reference = self.single_reference_grid
        reference_wavelengths = self.single_reference_grid_wavelengths
        nreference_wavelengths = len(reference_wavelengths)

        linewidth = 0.4

        # Add grid lines
        matching_indices = []
        for w in self.all_grid_wavelengths:

            #
            index = sequences.find_index(reference_wavelengths, w)
            if index is not None: matching_indices.append(index)
            else: # not found in reference
                self.differences.vlines(w, self.min_y_differences, self.max_y_differences, color="green", lw=linewidth, alpha=self.config.linealpha)

        # Loop over the reference wavelengths
        for index in range(nreference_wavelengths):
            if index in matching_indices: continue
            wavelength = reference_wavelengths[index]
            self.differences.vlines(wavelength, self.min_y_differences, self.max_y_differences, color="red", lw=linewidth, alpha=self.config.linealpha)

    # -----------------------------------------------------------------

    def plot_legend(self):

        """
        This function ...
        :return:
        """

        # Create the legend
        labels = [handle.get_label() for handle in self.grid_handles]
        self.grids_legend = Legend(self.main_axes, self.grid_handles, labels, **self.grids_legend_properties)

        # Add the legend
        self.main_axes.add_artist(self.grids_legend)

    # -----------------------------------------------------------------

    def plot_deltas(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the grid deltas ...")

        # Plot the deltas
        #self.delta.plot(self.delta_wavelengths, self.deltas, color="m", linestyle="solid", lw=self.delta_linewidth)
        self.delta.plot(self.all_grid_wavelengths, self.complete_grid_logdeltas, color="m", linestyle="solid", lw=self.delta_linewidth) # (correct) log delta calculation with minimas and maximas like how SKIRT calculates them

    # -----------------------------------------------------------------

    @property
    def seds_legend_properties(self):

        """
        This function ...
        :return:
        """

        legend_properties = dict()

        legend_properties["loc"] = 9
        legend_properties["ncol"] = 1
        legend_properties["shadow"] = False
        legend_properties["frameon"] = True
        legend_properties["facecolor"] = None
        legend_properties["fontsize"] = "smaller"

        return legend_properties

    # -----------------------------------------------------------------

    @property
    def seds_min_y(self):

        """
        This function ...
        :return:
        """

        return 1.

    # -----------------------------------------------------------------

    @property
    def seds_max_y(self):

        """
        This function ...
        :return:
        """

        return 2.5

    # -----------------------------------------------------------------

    @property
    def sed_linewidth(self):

        """
        This function ...
        :return:
        """

        #return 0.3
        #return 1.
        return 1.

    # -----------------------------------------------------------------

    @property
    def sed_linealpha(self):

        """
        This function ...
        :return:
        """

        return 0.7

    # -----------------------------------------------------------------

    @lazyproperty
    def seds_colours(self):

        """
        This function ...
        :return:
        """

        return dark_pretty_colors[-4:]

    # -----------------------------------------------------------------

    @property
    def seds_min_wavelength_scaling(self):

        """
        This function ...
        :return:
        """

        return 0.9

    # -----------------------------------------------------------------

    @property
    def seds_max_wavelength_scaling(self):

        """
        This function ...
        :return:
        """

        return 1.1

    # -----------------------------------------------------------------

    @property
    def seds_min_wavelength(self):

        """
        This function ...
        :return:
        """

        return self.seds_min_wavelength_scaling * self.grid_min_wavelength

    # -----------------------------------------------------------------

    @property
    def seds_max_wavelength(self):

        """
        This function ...
        :return:
        """

        return self.seds_max_wavelength_scaling * self.grid_max_wavelength

    # -----------------------------------------------------------------

    def plot_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting SEDs ...")

        handles = []

        #line_styles_iterator = iter(line_styles)
        colours = iter(self.seds_colours)

        # Loop over the SEDs
        for label in self.seds:

            # Debugging
            log.debug("Adding '" + label + "' SED to the plot ...")

            # Get the wavelengths and fluxes array
            #wavelengths = self.seds[label].wavelengths(unit=self.config.wavelength_unit, asarray=True)
            #fluxes = self.seds[label].normalized_photometry(method="max")
            wavelengths = self.sed_wavelengths[label]
            fluxes = self.sed_fluxes[label]

            nonzero = fluxes != 0
            wavelengths = wavelengths[nonzero]
            fluxes = fluxes[nonzero]

            # Calculate minlogflux here, before removing outside min and max wavelength
            logfluxes = np.log10(fluxes)
            #minlogflux = np.mean(logfluxes)
            maxlogflux = np.max(logfluxes)
            minlogflux = maxlogflux - 3 # set zero point to be 3 magitudes lower than the maximum

            # Remove below min and above max
            below_min = wavelengths < self.seds_min_wavelength
            above_max = wavelengths > self.seds_max_wavelength
            outside = below_min + above_max
            not_outside = np.logical_not(outside)
            wavelengths = wavelengths[not_outside]
            fluxes = fluxes[not_outside]

            # Normalize in log space
            logfluxes = np.log10(fluxes)
            #minlogflux = np.mean(logfluxes)
            logfluxes = logfluxes - minlogflux

            maxlogflux = np.max(logfluxes)
            logfluxes = logfluxes / maxlogflux

            log_fluxes = self.seds_min_y + logfluxes  # normalized
            #log_fluxes = logfluxes

            # Get line style
            #line_style = next(line_styles_iterator)
            line_style = "-"

            #color = "c"
            #color = "grey"
            color = next(colours)

            # Plot the SED
            h = self.main.plot(wavelengths, log_fluxes, color=color, lw=self.sed_linewidth, label=label, ls=line_style, alpha=self.sed_linealpha)
            handle = h[0]

            # Add handle
            handles.append(handle)

            # Plot resampled
            if self.config.plot_resampled:

                resampled = self.sed_grid_fluxes[label]
                logresampled = np.log10(resampled)
                logresampled = logresampled - minlogflux
                logresampled = logresampled / maxlogflux
                logresampled = self.seds_min_y + logresampled

                pointsize = 15
                pointcolor = "black"

                # wavelengths, y, s=pointsize, marker='.', color=color, linewidths=0, label=legend_label
                self.main.scatter(self.all_grid_wavelengths, logresampled, s=pointsize, marker=".", color=pointcolor, linewidths=0, zorder=100) # make sure is above all lines

            # Plot interpolated
            if self.config.plot_interpolated:

                # Get resampled, mask out
                resampled = self.sed_interpolated_fluxes[label][nonzero][not_outside]

                # Normalize in the same way
                logresampled = np.log10(resampled)
                logresampled = logresampled - minlogflux
                logresampled = logresampled / maxlogflux
                logresampled = self.seds_min_y + logresampled

                #line_style = '--'
                linewidth = 0.8 * self.sed_linewidth
                color = "grey"

                # Plot
                self.main.plot(wavelengths, logresampled, color=color, lw=linewidth, ls=line_style, alpha=self.sed_linealpha)

        # Create the legend
        labels = [handle.get_label() for handle in handles]
        self.seds_legend = Legend(self.main_axes, handles, labels, **self.seds_legend_properties)

        # Add the legend
        self.main_axes.add_artist(self.seds_legend)

    # -----------------------------------------------------------------

    @lazyproperty
    def sorted_filter_labels(self):

        """
        This function ...
        :return:
        """

        # NEW
        if self.config.add_filters: self.add_filters()

        return list(sorted(self.filters.keys(), key=lambda label: self.filters[label].wavelength.to("micron").value))

    # -----------------------------------------------------------------

    @lazyproperty
    def categorized_filter_labels(self):

        """
        This function ...
        :return:
        """

        # NEW
        if self.config.add_filters: self.add_filters()

        categorized = defaultdict(list)

        for label in self.filters:
            fltr = self.filters[label]
            categorized[fltr.instrument].append(label)

        #print(self.filters)
        #print("categorized", categorized)

        sort_key = lambda labels: min([self.filters[label].wavelength.to("micron").value for label in labels])
        return containers.ordered_by_value(categorized, key=sort_key)

    # -----------------------------------------------------------------

    @property
    def max_y_filters(self):

        """
        This function ...
        :return:
        """

        #return 1.5
        #return 1.1
        return 0.8

    # -----------------------------------------------------------------

    @property
    def filters_legend_properties(self):

        """
        This function ...
        :return:
        """

        legend_properties = dict()

        legend_properties["loc"] = 2
        #legend_properties["numpoints"] = 4
        #legend_properties["scatterpoints"] = 4
        legend_properties["ncol"] = 4
        legend_properties["shadow"] = False
        legend_properties["frameon"] = True
        legend_properties["facecolor"] = None
        legend_properties["fontsize"] = "smaller"

        return legend_properties

    # -----------------------------------------------------------------

    @property
    def filters_colormap(self):

        """
        This function ...
        :return:
        """

        name = "rainbow"
        return plt.get_cmap(name)

    # -----------------------------------------------------------------

    @lazyproperty
    def filters_hierarchy(self):

        """
        This function ...
        :return:
        """

        # Loop over the filters
        if self.config.categorize_filters: categorized = self.categorized_filter_labels
        else: categorized = {"dummy": self.sorted_filter_labels}
        return categorized

    # -----------------------------------------------------------------

    @lazyproperty
    def filter_colours(self):

        """
        This function ...
        :return:
        """

        colours = dict()

        cm = self.filters_colormap

        if self.config.categorize_filters: ncolours = len(self.filters_hierarchy)
        else: ncolours = self.nfilters

        cNorm = colors.Normalize(vmin=0, vmax=ncolours - 1.)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

        # Loop over the instrument
        counter = 0
        instrument_counter = 0
        for instrument in self.filters_hierarchy:

            if self.config.categorize_filters: colorVal = scalarMap.to_rgba(instrument_counter)
            else: colorVal = None

            # Loop over the instruments
            for label in self.filters_hierarchy[instrument]:

                # Get the filter
                fltr = self.filters[label]

                # Get color for this filter
                if not self.config.categorize_filters: colorVal = scalarMap.to_rgba(counter)

                # Set the colour
                colours[fltr] = colorVal

                # Increment counter
                counter += 1

            # Increment instrument counter
            instrument_counter += 1

        # Return the colours
        return colours

    # -----------------------------------------------------------------

    @property
    def nfilter_colours(self):

        """
        This function ...
        :return:
        """

        return len(self.filter_colours)

    # -----------------------------------------------------------------

    @property
    def has_filter_colours(self):

        """
        This function ....
        :return:
        """

        return self.nfilter_colours > 0

    # -----------------------------------------------------------------

    def has_colour_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        if not self.has_filter_colours: return None
        else: return fltr in self.filter_colours

    # -----------------------------------------------------------------

    def get_filter_colour(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        if not self.has_filter_colours: raise ValueError("No filter colours")
        if types.is_string_type(fltr): fltr = parse_filter(fltr)
        return self.filter_colours[fltr]

    # -----------------------------------------------------------------

    def plot_filters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting filters ...")

        min_wavelength = float("inf")
        max_wavelength = 0.0

        alpha = 0.2
        linewidth = 1.
        linealpha = 0.5

        handles = []

        #print(self.filters_hierarchy)

        # Loop over the instrument
        for instrument in self.filters_hierarchy:

            first_for_instrument = True

            # Loop over the instruments
            for label in self.filters_hierarchy[instrument]:

                # Get the filter
                fltr = self.filters[label]
                description = fltr.description()

                # Get color for this filter
                colorVal = self.filter_colours[fltr]

                # Determine legend label
                if self.config.categorize_filters:
                    if first_for_instrument: legend_label = instrument
                    else: legend_label = None
                else: legend_label = label

                # Broad band filter
                if isinstance(fltr, BroadBandFilter):

                    from ..data.transmission import TransmissionCurve

                    # Create transmission curve
                    curve = TransmissionCurve.from_filter(fltr)
                    curve.normalize(value=self.max_y_filters, method="max")

                    wavelengths = curve.wavelengths(unit=self.config.wavelength_unit, add_unit=False)
                    transmissions = curve.transmissions()

                    # Plot the curve
                    h = self.main.fill(wavelengths, transmissions, label=legend_label, color=colorVal, alpha=alpha, linewidth=linewidth)
                    handle = h[0]

                    minw = np.min(wavelengths)
                    maxw = np.max(wavelengths)
                    if minw < min_wavelength: min_wavelength = minw
                    if maxw > max_wavelength: max_wavelength = maxw

                # Narrow band filter
                elif isinstance(fltr, NarrowBandFilter):

                    wavelength = fltr.wavelength.to(self.config.wavelength_unit).value
                    handle = self.main.vlines(wavelength, self.min_y, self.max_y_filters, color=colorVal, lw=linewidth, alpha=linealpha, label=legend_label)

                # Invalid
                else: raise ValueError("Unrecognized filter object")

                # Add the handle
                handles.append(handle)

                # Set flag
                first_for_instrument = False

        # Create legend
        valid_handles = [handle for handle in handles if handle.get_label() is not None]
        labels = [handle.get_label() for handle in valid_handles]
        self.filters_legend = Legend(self.main_axes, valid_handles, labels, **self.filters_legend_properties)

        # Add legend
        self.main_axes.add_artist(self.filters_legend)

    # -----------------------------------------------------------------

    @property
    def max_y_lines(self):

        """
        This function ...
        :return:
        """

        return 0.9

    # -----------------------------------------------------------------

    @property
    def lines_linewidth(self):

        """
        This function ...
        :return:
        """

        return 0.5
        #return 0.3

    # -----------------------------------------------------------------

    @property
    def lines_label_alpha(self):

        """
        This function ...
        :return:
        """

        return 0.7

    # -----------------------------------------------------------------

    @property
    def lines_alpha(self):

        """
        This function ...
        :return:
        """

        return 0.9

    # -----------------------------------------------------------------

    @lazyproperty
    def lines_colours(self):

        """
        This function ...
        :return:
        """

        #return ['m', 'r', 'y', 'g']
        return dark_pretty_colors[-4:]

    # -----------------------------------------------------------------

    @property
    def nline_colours(self):

        """
        This function ...
        :return:
        """

        return len(self.lines_colours)

    # -----------------------------------------------------------------

    @property
    def separate_line_colours(self):

        """
        This function ...
        :return:
        """

        return dark_pretty_colors[::-1]

    # -----------------------------------------------------------------

    @property
    def nseparate_line_colours(self):

        """
        This function ...
        :return:
        """

        return len(self.separate_line_colours)

    # -----------------------------------------------------------------

    @lazyproperty
    def lines_y_range(self):

        """
        This function ...
        :return:
        """

        return RealRange(1., 1.5, inclusive=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def lines_y_values(self):

        """
        This function ...
        :return:
        """

        return self.lines_y_range.linear(self.lines_nsteps)

    # -----------------------------------------------------------------

    @property
    def lines_nsteps(self):

        """
        This function ...
        :return:
        """

        return 4

    # -----------------------------------------------------------------

    @property
    def nlines(self):

        """
        This function ...
        :return:
        """

        return len(self.emission_lines)

    # -----------------------------------------------------------------

    @lazyproperty
    def nlines_with_label(self):

        """
        This function ...
        :return:
        """

        return len([line for line in self.emission_lines if len(line.label) > 0])

    # -----------------------------------------------------------------

    def plot_lines(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting emission lines ...")

        # Plot a labeled vertical line for each emission line
        absolute_index = 0
        index = 0

        handles = []

        # Set line colours and styles list
        if self.config.separate_lines or self.config.lines_in_group:

            #print(self.nseparate_line_colours, self.nlines)

            if self.nlines > self.nseparate_line_colours:

                factor_more = float(self.nlines) / float(self.nseparate_line_colours)
                #print(factor_more)
                int_factor_more = int(math.ceil(factor_more))
                #print(int_factor_more)
                line_colours = self.separate_line_colours * int_factor_more
                linestyles = sequences.repeat(line_styles, self.nseparate_line_colours)

            else:
                line_colours = self.separate_line_colours
                linestyles = ["solid"] * self.nlines

        else:
            line_colours = self.lines_colours
            linestyles = ["solid"] * self.nline_colours

        # Loop over the emission lines
        for line in self.emission_lines:

            # Get center wavelength and label
            center = line.center.to(self.config.wavelength_unit).value
            label = line.label

            if len(label) == 0:
                log.warning("Emission line without label at " + str(center) + " " + str(self.config.wavelength_unit) + ": not plotting ...")
                continue

            # Plot lines separately
            if self.config.separate_lines:

                # Get line colour and line style
                colour = line_colours[absolute_index]
                linestyle = linestyles[absolute_index]

                # Plot line
                handle = self.lines.vlines(center, self.min_y_separate_lines, self.max_y_separate_lines, color=colour, lw=self.lines_linewidth, alpha=self.lines_alpha, label=label, linestyles=linestyle)
                handles.append(handle)

                t = None

            # Plot lines in group
            elif self.config.lines_in_group:

                # Get line colour
                colour = line_colours[absolute_index]
                linestyle = linestyles[absolute_index]

                # Plot line
                handle = self.groups[self.config.lines_in_group].vlines(center, self.min_y_groups, self.max_y_groups, color=colour, lw=self.lines_linewidth, alpha=self.lines_alpha, label=label, linestyles=linestyle)
                handles.append(handle)

                # Add text
                #t = self.groups[self.config.lines_in_group].text(center, 0.75, label, horizontalalignment='center', fontsize='xx-small', color=colour, backgroundcolor='w')
                t = None

            # Plot lines on main row
            else:

                # Get colour, rotating
                colour = line_colours[index]
                linestyle = linestyles[index]

                # Plot line
                self.main.vlines(center, self.min_y, self.max_y_lines, color=colour, lw=self.lines_linewidth, alpha=self.lines_alpha, linestyles=linestyle)

                # Add text
                t = self.main.text(center, self.lines_y_values[index], label, horizontalalignment='center', fontsize='xx-small', color=colour, backgroundcolor='w')

            # Set text box options
            if t is not None: t.set_bbox(dict(color='w', alpha=self.lines_label_alpha, edgecolor='w'))

            # Adapt indices
            index = (index + 1) % self.lines_nsteps
            absolute_index += 1

        # Lines legend on separate row
        if self.config.separate_lines:

            # Create the legend
            labels = [handle.get_label() for handle in handles]
            lines_legend = Legend(self.lines.axes, handles, labels, **self.lines_legend_properties)

            # Add the legend to the row
            self.lines.add_artist(lines_legend)

        # Lines legend on group row
        elif self.config.lines_in_group:

            # Create the legend
            labels = [handle.get_label() for handle in handles]
            lines_legend = Legend(self.groups[self.config.lines_in_group].axes, handles, labels, **self.lines_legend_properties)

            # Add the legend to the group plot
            self.groups[self.config.lines_in_group].add_artist(lines_legend)

    # -----------------------------------------------------------------

    @property
    def max_y_wavelengths(self):

        """
        This function ...
        :return:
        """

        return 2.

    # -----------------------------------------------------------------

    @property
    def wavelengths_linewidth(self):

        """
        Thi function ...
        :return:
        """

        return 1.

    # -----------------------------------------------------------------

    @property
    def wavelengths_label_alpha(self):

        """
        This function ...
        :return:
        """

        return 0.7

    # -----------------------------------------------------------------

    @property
    def wavelengths_alpha(self):

        """
        This function ...
        :return:
        """

        return 0.8

    # -----------------------------------------------------------------

    def plot_wavelengths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting individual wavelengths ...")

        # Loop over the wavelengths
        for properties in self.wavelengths:

            # Get properties
            wavelength = properties.wavelength.to(self.config.wavelength_unit).value
            label = properties.label
            color = properties.color
            linewidth = properties.linewidth
            group = properties.group
            connect = properties.connect
            connect_in_style = properties.connect_in_style

            # Determine linewidth
            if linewidth is None: linewidth = self.wavelengths_linewidth

            # Plot line
            if group is not None: self.groups[group].vlines(wavelength, self.min_y_groups, self.max_y_groups, color=color, lw=linewidth, alpha=self.wavelengths_alpha)
            elif self.config.group_wavelengths: self.individual.vlines(wavelength, self.min_y_individual, self.max_y_individual, color=color, lw=linewidth, alpha=self.wavelengths_alpha)
            else: self.main.vlines(wavelength, self.min_y, self.max_y_wavelengths, color=color, lw=linewidth, alpha=self.wavelengths_alpha)

            # Add text
            if label is not None:
                if group is not None: t = self.groups[group].text(wavelength, 0.7, label, horizontalalignment='center', fontsize='xx-small', color=color, backgroundcolor='w')
                elif self.config.group_wavelengths: t = self.individual.text(wavelength, 0.7, label, horizontalalignment='center', fontsize='xx-small', color=color, backgroundcolor='w')
                else: t = self.main.text(wavelength, self.max_y_wavelengths, label, horizontalalignment='center', fontsize='xx-small', color=color, backgroundcolor='w')
                t.set_bbox(dict(color='w', alpha=self.wavelengths_label_alpha, edgecolor='w'))

            # Connect to main plot?
            if connect:

                if connect_in_style:

                    connect_color = color
                    connect_linewidth = linewidth
                    connect_linealpha = self.wavelengths_alpha
                    connect_linestyle = "solid"

                else:

                    connect_color = "grey"
                    connect_linewidth = 0.7
                    connect_linealpha = 0.7
                    connect_linestyle = "dotted"

                # Order: main, separate grids, groups, lines, individual wavelengths grouped

                if group is not None:

                    # All separate grids
                    for index in range(self.nseparate_grids):
                        self.separate[index].vlines(wavelength, self.min_y_separate, self.max_y_separate, color=connect_color, lw=connect_linewidth, alpha=connect_linealpha, linestyles=connect_linestyle)

                    # Groups before
                    groups_before = sequences.before(self.groups.keys(), group)
                    for name in groups_before:
                        self.groups[name].vlines(wavelength, self.min_y_groups, self.max_y_groups, color=connect_color, lw=connect_linewidth, alpha=connect_linealpha, linestyles=connect_linestyle)

                elif self.config.group_wavelengths:

                    # All separate grids
                    for index in range(self.nseparate_grids):
                        self.separate[index].vlines(wavelength, self.min_y_separate, self.max_y_separate, color=connect_color, lw=connect_linewidth, alpha=connect_linealpha, linestyles=connect_linestyle)

                    # All groups
                    for name in self.group_names:
                        self.groups[name].vlines(wavelength, self.min_y_groups, self.max_y_groups, color=connect_color, lw=connect_linewidth, alpha=connect_linealpha, linestyles=connect_linestyle)

                    # Lines
                    if self.has_lines and self.config.separate_lines:

                        self.lines.vlines(wavelength, self.min_y_separate_lines, self.max_y_separate_lines, color=connect_color, lw=connect_linewidth, alpha=connect_linealpha, linestyles=connect_linestyle)

                else: raise ValueError("Cannot connect when plotted on main row")

    # -----------------------------------------------------------------

    @property
    def ranges_y_range(self):

        """
        This function ...
        :return:
        """

        #return RealRange(1.8, 2.2, inclusive=True) # above emission line labels (to 1.6)

        #print(self.max_y_delta)
        if self.config.ranges_on_deltas: return RealRange(0.7 * self.max_y_delta, 0.8 * self.max_y_delta)
        else: return RealRange(1.8, 2., inclusive=True)

    # -----------------------------------------------------------------

    @property
    def ranges_nsteps(self):

        """
        This function ...
        :return:
        """

        return self.config.ranges_nsteps

    # -----------------------------------------------------------------

    @property
    def ranges_y_values(self):

        """
        This function ...
        :return:
        """

        return self.ranges_y_range.linear(self.ranges_nsteps)

    # -----------------------------------------------------------------

    def plot_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting wavelength ranges ...")

        index = 0
        label_position = self.config.ranges_label_position

        # Loop over the wavelength ranges
        for properties in self.wavelength_ranges:

            # Get properties
            wavelength_range = properties.range.to(self.config.wavelength_unit).value
            label = properties.label
            color = properties.color

            # Get next y value
            #y = next(y_values)
            y = self.ranges_y_values[index]

            # Plot arrow
            if self.config.ranges_on_deltas: self.delta.horizontal_arrow(wavelength_range.min, wavelength_range.max, y)
            else: self.main.horizontal_arrow(wavelength_range.min, wavelength_range.max, y)

            #average_wavelength = wavelength_range.mean
            average_wavelength = wavelength_range.geometric_mean

            if color is None: color = "black"

            # Plot label
            if label is not None:

                if self.config.ranges_on_deltas:

                    if label_position == "above": label_y = y + 0.05 * self.delta_y_span
                    elif label_position == "below": label_y = y - 0.15 * self.delta_y_span
                    else: raise ValueError("Invalid label position")
                    self.delta.text(average_wavelength, label_y, label, horizontalalignment='center', fontsize='xx-small', color=color)

                else:

                    if label_position == "above": label_y = y + 0.05 * self.main_y_span
                    elif label_position == "below": label_y = y - 0.15 * self.main_y_span
                    else: raise ValueError("Invalid label position")
                    self.main.text(average_wavelength, label_y, label, horizontalalignment='center', fontsize='xx-small', color=color)

            # Increment or reset the index
            index = (index + 1) % self.ranges_nsteps

            # Flip position
            if self.config.ranges_alternate_labels:

                if label_position == "above": label_position = "below"
                elif label_position == "below": label_position = "above"
                else: raise ValueError("Invalid label position")

    # -----------------------------------------------------------------

    def plot_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting SED residuals ...")

        # Get colours for the SEDs
        colours = iter(self.seds_colours)

        # Loop over the SEDs
        for label in self.seds:

            # Debugging
            log.debug("Adding '" + label + "' SED residuals to the plot ...")

            # Get the wavelengths and fluxes array
            #wavelengths = self.seds[label].wavelengths(unit=self.config.wavelength_unit, asarray=True)
            #fluxes = self.seds[label].normalized_photometry(method="max")

            # Get the wavelengths and residuals
            wavelengths = self.sed_wavelengths[label]
            residuals = self.sed_residuals_percentual[label]

            # Get colour and line style
            color = next(colours)
            line_style = "-"

            # Plot the residuals
            h = self.residuals.plot(wavelengths, residuals, color=color, lw=self.sed_linewidth, label=label, ls=line_style, alpha=self.sed_linealpha)
            #handle = h[0]

    # -----------------------------------------------------------------

    def finish_plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finishing the plot ...")

        # Set the title
        if self.title is not None: self.figure.set_title(self.title, width=60)

        # Plot
        self.figure.finish(out=self.out_path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the complete grid
        self.write_complete()

        # Write the individual grids
        self.write_grids()

    # -----------------------------------------------------------------

    def write_complete(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the complete grid ...")

        # Determine the path
        path = self.output_path_file("complete_grid.dat")

        # Write the grid
        self.complete_grid.saveto(path)

    # -----------------------------------------------------------------

    def write_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the individual grids ...")

        # Loop over the grids
        for label in self.grids:

            # Debugging
            log.debug("Writing the '" + label + "' wavelength grid ...")

            # Get the wavelength points
            grid = self.grids[label].grid

            # Determine the path
            path = self.output_path_file(label + ".dat")

            # Save the grid
            grid.saveto(path)

# -----------------------------------------------------------------

def in_some_range(wavelength, filter_ranges):

    """
    This function ...
    :param wavelength:
    :param filter_ranges:
    :return:
    """

    for fltr in filter_ranges:
        wavelength_range = filter_ranges[fltr]
        if wavelength in wavelength_range: return True
    return False

# -----------------------------------------------------------------
