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
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.legend import Legend
from collections import defaultdict

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
        self.complete = None
        self.delta = None

        # Legends of the main plot
        self.grids_legend = None
        self.filters_legend = None
        self.seds_legend = None

        # The grid scatter handles
        self.grid_handles = []

        # Shared labels with their colors
        self.shared_labels = dict()

    # -----------------------------------------------------------------

    def set_title(self, title):

        """
        This function ...
        :param title:
        :return:
        """

        self.title = title

    # -----------------------------------------------------------------

    def add_wavelength_grid(self, grid, label, pointsize=None, linewidth=None, linealpha=None, color=None,
                            in_legend=True, y_value=None, copy_grid=True, shared_label=None):

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

        # Add grid
        self.grids[label] = props

    # -----------------------------------------------------------------

    def add_grid_from_file(self, path, label=None, pointsize=None, linewidth=None, linealpha=None, color=None,
                           in_legend=True, y_value=None, shared_label=None):

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
        :return:
        """

        from ..simulation.wavelengthgrid import WavelengthGrid
        if label is None: label = fs.strip_extension(fs.name(path))
        grid = WavelengthGrid.from_file(path)
        self.add_wavelength_grid(grid, label, pointsize=pointsize, linewidth=linewidth, linealpha=linealpha,
                                 color=color, in_legend=in_legend, y_value=y_value, copy_grid=False,
                                 shared_label=shared_label)

    # -----------------------------------------------------------------

    def add_wavelengths(self, wavelengths, label, unit=None, pointsize=None, linewidth=None, linealpha=None,
                        color=None, in_legend=True, y_value=None, shared_label=None):

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
        :return:
        """

        from ..simulation.wavelengthgrid import WavelengthGrid
        grid = WavelengthGrid.from_wavelengths(wavelengths, unit=unit)
        self.add_wavelength_grid(grid, label, pointsize=pointsize, linewidth=linewidth, linealpha=linealpha,
                                 color=color, in_legend=in_legend, y_value=y_value, copy_grid=False,
                                 shared_label=shared_label)

    # -----------------------------------------------------------------

    def add_wavelength(self, wavelength, label=None, colour=None):

        """
        This function ...
        :param wavelength:
        :param label:
        :param colour:
        :return:
        """

        # Set properties
        props = Map()
        props.wavelength = wavelength
        props.label = label
        props.color = colour

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

    def run(self, **kwargs):

        """
        This function ...
        :param output_path:
        :param min_wavelength:
        :param max_wavelength:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

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
        super(WavelengthGridPlotter, self).setup(**kwargs)

        # Set the title
        self.title = kwargs.pop("title", None)

        # Set the output path
        self.out_path = kwargs.pop("output", None)
        self.out_path = kwargs.pop("output_path", self.out_path)

        # Set the axis limits
        self.min_wavelength = kwargs.pop("min_wavelength", self.config.min_wavelength)
        self.max_wavelength = kwargs.pop("max_wavelength", self.config.max_wavelength)

        # Load grids
        if self.config.grids is not None:
            for path in self.config.grids:
                self.add_grid_from_file(path)

        # Load grids from subgrid generation
        if self.config.load_subgrids: self.load_subgrids()

        # Generate subgrids
        elif self.config.create_subgrids: self.create_subgrids()

        # Are there grids
        if not self.has_grids: raise RuntimeError("No wavelength grids are added")

        # Add filters
        if self.config.add_filters: self.add_filters()

        # Add emission lines
        if self.config.add_lines: self.create_lines()

        # Add SEDs
        if self.config.add_seds: self.create_seds()

        # Add regimes
        if self.config.add_regimes: self.add_regimes()

        # Get figure size
        figsize = self.config.plot.figsize

        # Create the plot
        if self.config.library == mpl: self.figure = MPLFigure(size=figsize)
        elif self.config.library == bokeh: self.figure = BokehFigure()
        else: raise ValueError("Invalid libary: " + self.config.library)

        # Set the 'show' flag
        if self.config.show is None:
            if self.out_path is not None: self.config.show = False
            else: self.config.show = True

        # Initialize plot
        self.initialize_plot()

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

    @lazyproperty
    def y_limits(self):


        """
        This function ...
        :return:
        """

        if self.config.separate_complete_grid: return [self.main_y_limits, self.complete_y_limits, self.delta_y_limits]
        else: return [self.main_y_limits, self.delta_y_limits]

    # -----------------------------------------------------------------

    def initialize_plot(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Initializing the plot ...")

        # max_delta = 0.1 # was 0.05
        # max_delta = None

        if self.config.separate_complete_grid:

            nplots = 3
            height_ratios = (4, 0.3, 1)
            x_label = '$\lambda/\mu$m'
            y_labels = [None, None, r"$\Delta\lambda\,(\mathrm{dex})$"]
            x_scale = "log"

            # Create subplots
            self.main, self.complete, self.delta = self.figure.create_column(nplots, share_axis=True,
                                                                             height_ratios=height_ratios,
                                                                             x_label=x_label, y_labels=y_labels,
                                                                             x_limits=self.x_limits,
                                                                             x_scale=x_scale, y_limits=self.y_limits,
                                                                             x_log_scalar=True)

        else:

            nplots = 2
            height_ratios = (4, 1)
            x_label = '$\lambda/\mu$m'
            y_labels = [None, r"$\Delta\lambda\,(\mathrm{dex})$"]
            x_scale = "log"

            # Create subplots
            self.main, self.delta = self.figure.create_column(nplots, share_axis=True, height_ratios=height_ratios,
                                                              x_label=x_label, y_labels=y_labels,
                                                              x_limits=self.x_limits,
                                                              x_scale=x_scale, y_limits=self.y_limits,
                                                              x_log_scalar=True)

    # -----------------------------------------------------------------

    def load_subgrids(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading from subgrids wavelength grid generation ...")

        # Set paths
        fixed_path = self.input_path_file("fixed_grid.dat")
        new_path = self.input_path_file("new.dat")
        replaced_path = self.input_path_file("replaced.dat")
        filter_wavelengths_path = self.input_path_file("filter_wavelengths.dat")
        line_wavelengths_path = self.input_path_file("line_wavelengths.dat")

        # Loop over the subgrid files
        for path, name in fs.files_in_path(self.input_path, extension="dat", startswith="subgrid_",
                                           returns=["path", "name"]):
            # Get subgrid label
            subgrid = name.split("subgrid_")[1]

            # Check
            if subgrid not in self.config.subgrids: continue

            # Debugging
            log.debug("Adding '" + subgrid + "' subgrid ...")

            # Add grid
            self.add_grid_from_file(path, label=subgrid)

        # Fixed?
        if fs.is_file(fixed_path):

            # Debugging
            log.debug("Adding fixed wavelengths ...")

            # Add
            self.add_grid_from_file(fixed_path, label="fixed", linewidth=1., pointsize=20)

        # New?
        if fs.is_file(new_path):

            # Debugging
            log.debug("Adding new wavelengths ...")

            # Add wavelengths
            new = load_list(new_path)
            self.add_wavelengths(new, label="new", color="b")

        # Replaced?
        if fs.is_file(replaced_path):

            # Debugging
            log.debug("Adding replaced wavelengths ...")

            # Add
            replaced = load_list(replaced_path)
            for old, replacement in replaced:
                self.add_wavelength(old, colour="red")
                self.add_wavelength(replacement, colour="green")

                # Remove the old wavelength from the complete grid
                self.remove_wavelength(old)

        # Filter
        if fs.is_file(filter_wavelengths_path):

            # Debugging
            log.debug("Adding filter wavelengths ...")

            # Add
            filter_wavelengths = load_dict(filter_wavelengths_path)
            for fltr in filter_wavelengths:
                wavelengths = filter_wavelengths[fltr]
                filter_name = str(fltr)
                self.add_wavelengths(wavelengths, label=filter_name)

        # Emission lines?
        if fs.is_file(line_wavelengths_path):

            # Debugging
            log.debug("Adding emission lines ...")

            line_wavelengths = load_dict(line_wavelengths_path)
            for identifier in line_wavelengths:
                wavelengths = line_wavelengths[identifier]
                self.add_wavelengths(wavelengths, label=identifier, color="grey", shared_label="lines", y_value=0.6)

    # -----------------------------------------------------------------

    def create_subgrids(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Generating a subgrids wavelength grid ...")

        from ...modeling.build.wavelengthgrid import WavelengthGridBuilder
        from ..simulation.wavelengthgrid import WavelengthGrid

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

        # Loop over the subgrids
        for name in builder.subgrids:

            # Debugging
            log.debug("Adding '" + name + "' subgrid ...")

            # Add grid
            subgrid = builder.subgrids[name]
            self.add_wavelength_grid(subgrid, label=name)

        # Grid of fixed wavelengths
        if builder.has_fixed:

            # Debugging
            log.debug("Adding fixed wavelengths ...")

            # Add
            fixed_grid = WavelengthGrid.from_wavelengths(builder.fixed)
            self.add_wavelength_grid(fixed_grid, label="fixed", linewidth=1., pointsize=20)

        # Get filter wavelengths
        if builder.has_filters:

            # Debugging
            log.debug("Adding filter wavelengths ...")

            # Get
            filter_wavelengths = builder.filter_wavelengths

            # Loop over the filters
            for fltr in filter_wavelengths:

                wavelengths = filter_wavelengths[fltr]
                filter_name = str(fltr)
                self.add_wavelengths(wavelengths, label=filter_name)

        # Get replaced
        if builder.has_replaced:

            replaced = builder.replaced

            # Debugging
            log.debug("Adding replaced wavelengths ...")

            # Add
            for old, replacement in replaced:
                self.add_wavelength(old, colour="red")
                self.add_wavelength(replacement, colour="green")

                # Remove the old wavelength from the complete grid
                self.remove_wavelength(old)

        # Get new
        if builder.has_new:

            new = builder.new

            # Debugging
            log.debug("Adding new wavelengths ...")

            # Add wavelengths
            self.add_wavelengths(new, label="new", color="b")

        # Get emission line wavelengths
        if builder.has_line_wavelengths:

            line_wavelengths = builder.line_wavelengths

            # Debugging
            log.debug("Adding emission lines ...")

            # Loop over the lines
            for identifier in line_wavelengths:
                wavelengths = line_wavelengths[identifier]
                if identifier[0] is None: continue
                #print(identifier)
                label = identifier[0] + identifier[1]
                self.add_wavelengths(wavelengths, label=label, color="grey", y_value=0.6, shared_label="lines")

    # -----------------------------------------------------------------

    def add_filters(self):

        """
        This function ....
        :return:
        """

        # Inform the user
        log.info("Adding filters ...")

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

            # Add wavelengths
            #all_grid_wavelengths.extend(wavelengths)

            # Add wavelengths but not those to remove from complete grid
            for wavelength in wavelengths:

                # NO: DON'T REMOVE HERE: CAN BE RE-ADDED IF SOME WAVELENGTHS ARE ADDED AS SINGLE WAVELENGTHS
                # BUT THEN ALSO ADDED TO REMOVE
                #if self.has_remove:
                    #wavelength in self.remove_no_unit: continue
                    #fw = sequences.find_closest_value(self.remove_no_unit, wavelength)
                    #print(wavelength, fw)
                    #if np.isclose(fw, wavelength): print(wavelength in self.remove_no_unit)
                    #if wavelength in self.remove_no_unit:
                    #    print(wavelength)
                    #    continue
                all_grid_wavelengths.append(wavelength)

        # Add single wavelengths
        if self.has_wavelengths:
            # Loop over the wavelengths
            for wavelength in self.wavelengths_no_unit: all_grid_wavelengths.append(wavelength)

        # Remove wavelengths?
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
    def nfilters(self):

        """
        This function ...
        :return:
        """

        return len(self.filters)

    # -----------------------------------------------------------------

    @property
    def has_lines(self):

        """
        This function ...
        :return:
        """

        return self.nfilters > 0

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

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelength grids ...")

        # SEDs
        if self.has_seds: self.plot_seds()

        # Filters
        if self.has_filters: self.plot_filters()

        # Grids
        self.plot_grids()

        # Plot complete grid
        if self.config.plot_complete_grid: self.plot_complete_grid()

        # Plot legend
        self.plot_legend()

        # Plot complete grid deltas
        self.plot_deltas()

        # Emission lines
        if self.has_lines: self.plot_lines()

        # Single wavelengths
        if self.has_wavelengths: self.plot_wavelengths()

        # Wavelength ranges
        if self.has_ranges: self.plot_ranges()

        # Finish
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
    def delta_linewidth(self):

        """
        Thisn function ...
        :return:
        """

        return 0.6

    # -----------------------------------------------------------------

    def plot_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the grids ...")

        # Make iterable from distinct colors
        different_colors = iter(dark_pretty_colors)

        # Get y values iterator
        y_values = iter(self.grid_y_values)

        # Loop over the wavelength grids
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

            # Check
            if shared_label is not None and not in_legend: raise ValueError("in_legend is disabled but shared label is specified")

            add_label = True
            if shared_label is not None:
                if shared_label in self.shared_labels:
                    if color is None: color = self.shared_labels[shared_label]
                    elif color != self.shared_labels[shared_label]: raise ValueError("Shared label but colours are not the same")
                    add_label = False # not a new (shared) label
                else:
                    if color is None: color = next(different_colors)
                    self.shared_labels[shared_label] = color
                    add_label = True # a new (shared) label

            # Get next color, if needed
            if color is None: color = next(different_colors)

            # Set line width and line alpha
            if pointsize is None: pointsize = self.config.pointsize
            if linewidth is None: linewidth = self.config.linewidth
            if linealpha is None: linealpha = self.config.linealpha

            # Get next y value, if needed
            if y_value is None: y_value = next(y_values)

            # Set y for each point
            y = [y_value for _ in wavelengths]

            # Plot points
            if add_label:
                if shared_label is not None: legend_label = shared_label
                else: legend_label = label + " (" + str(nwavelengths) + " points)"
            else: legend_label = None
            sc = self.main.scatter(wavelengths, y, s=pointsize, marker='.', color=color, linewidths=0, label=legend_label)

            # Plot a vertical line for each grid point
            for w in wavelengths: self.main.vlines(w, self.min_y, self.max_y, color=color, lw=linewidth, alpha=linealpha)

            # Add the handle
            if in_legend and add_label: self.grid_handles.append(sc)

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

        # Plot the complete grid in a seperate subplot
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
            wavelengths = self.seds[label].wavelengths(unit=self.config.wavelength_unit, asarray=True)
            fluxes = self.seds[label].normalized_photometry(method="max")

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

            logfluxes = logfluxes / np.max(logfluxes)

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

        return list(sorted(self.filters.keys(), key=lambda label: self.filters[label].wavelength.to("micron").value))

    # -----------------------------------------------------------------

    @lazyproperty
    def categorized_filter_labels(self):

        """
        This function ...
        :return:
        """

        categorized = defaultdict(list)

        for label in self.filters:
            fltr = self.filters[label]
            categorized[fltr.instrument].append(label)

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

    def plot_filters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting filters ...")

        # Get the color map
        colormap = "rainbow"
        cm = plt.get_cmap(colormap)

        counter = 0
        min_wavelength = float("inf")
        max_wavelength = 0.0

        alpha = 0.2
        linewidth = 1.
        linealpha = 0.5

        handles = []

        # Loop over the filters
        if self.config.categorize_filters: categorized = self.categorized_filter_labels
        else: categorized = {"dummy": self.sorted_filter_labels}

        if self.config.categorize_filters: ncolours = len(categorized)
        else: ncolours = self.nfilters

        cNorm = colors.Normalize(vmin=0, vmax=ncolours - 1.)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

        # Loop over the instrument
        instrument_counter = 0
        for instrument in categorized:

            if self.config.categorize_filters: colorVal = scalarMap.to_rgba(instrument_counter)
            else: colorVal = None

            first_for_instrument = True

            # Loop over the instruments
            for label in categorized[instrument]:

                # Get the filter
                fltr = self.filters[label]
                description = fltr.description()

                # Get color for this filter
                if not self.config.categorize_filters: colorVal = scalarMap.to_rgba(counter)

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

                # Increment counter
                counter += 1

                # Set flag
                first_for_instrument = False

            # Increment instrument counter
            instrument_counter += 1

        #print(handles)

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

        return 0.8

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

    @lazyproperty
    def lines_y_range(self):

        """
        This function ...
        :return:
        """

        #return RealRange(1.15, 1.6, inclusive=True)
        return RealRange(1., 1.5, inclusive=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def lines_y_values(self):

        """
        This function ...
        :return:
        """

        #return [1.3, 1.45, 1.6, 1.75]
        #return [1.15, 1.3, 1.45, 1.6]
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

    def plot_lines(self):


        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting emission lines ...")

        # Plot a labeled vertical line for each emission line
        index = 0

        # Loop over the emission lines
        for line in self.emission_lines:

            # Get center wavelength and label
            center = line.center.to(self.config.wavelength_unit).value
            label = line.label

            if len(label) > 0:

                # Get colour
                colour = self.lines_colours[index]

                # Plot line
                self.main.vlines(center, self.min_y, self.max_y_lines, color=colour, lw=self.lines_linewidth, alpha=self.lines_alpha)

                # Add text
                t = self.main.text(center, self.lines_y_values[index], label, horizontalalignment='center', fontsize='xx-small', color=colour, backgroundcolor='w')
                t.set_bbox(dict(color='w', alpha=self.lines_label_alpha, edgecolor='w'))
                index = (index + 1) % self.lines_nsteps

            # Warn
            else: log.warning("Emission line without label at " + str(center) + " " + str(self.config.wavelength_unit))

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

            # Plot line
            self.main.vlines(wavelength, self.min_y, self.max_y_wavelengths, color=color, lw=self.wavelengths_linewidth, alpha=self.wavelengths_alpha)

            # Add text
            if label is not None:
                t = self.main.text(wavelength, self.max_y_wavelengths, label, horizontalalignment='center',
                                   fontsize='xx-small', color=color, backgroundcolor='w')
                t.set_bbox(dict(color='w', alpha=self.wavelengths_label_alpha, edgecolor='w'))

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

    def finish_plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finishing the plot ...")

        # Set the title
        if self.title is not None: self.figure.set_title(self.title, width=60)  # self._figure.suptitle("\n".join(wrap(self.title, 60)))

        # # Save or show the plot
        # # if self.out_path is None: self.figure.show()
        # if self.config.show:
        #     #self.figure.show()
        #     log.debug("Showing the SED plot ...")
        #     plt.show()
        #     plt.close()
        #
        # # Save the figure
        # if self.out_path is not None: self.save_figure()

        self.figure.finish(out=self.out_path)

# -----------------------------------------------------------------
