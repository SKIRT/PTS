#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.wavelengthgrids Contains the WavelengthGridGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import defaultdict, OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ..basics.emissionlines import EmissionLines, EmissionLine
from ..basics.configurable import Configurable
from ..basics.table import SmartTable
from ..basics.range import QuantityRange
from ..units.parsing import parse_unit as u
from ..filter.broad import BroadBandFilter
from ..filter.narrow import NarrowBandFilter
from ..tools import sequences
from ..basics.containers import DefaultOrderedDict
from ..tools.stringify import stringify_list_fancy
from ..tools import parsing
from ..tools.utils import lazyproperty
from ..tools import filesystem as fs
from ..tools.serialization import write_dict, write_list
from ..units.parsing import parse_quantity as q

# -----------------------------------------------------------------

euv = "EUV"
stellar = "stellar"
aromatic = "aromatic"
thermal = "thermal"
microwave = "microwave"

# -----------------------------------------------------------------

# The names of the subgrids
subgrids = [euv, stellar, aromatic, thermal, microwave]

# Define the ranges of the subgrids
ranges = OrderedDict()
ranges[euv] = QuantityRange(0.02, 0.085, unit="micron")
ranges[stellar] = QuantityRange(0.085, 3., unit="micron")
ranges[aromatic] = QuantityRange(3., 27., unit="micron")
ranges[thermal] = QuantityRange(27., 1000., unit="micron")
ranges[microwave] = QuantityRange(1000., 2000, unit="micron")

# Define the relative number of points of the subgrids
relpoints = OrderedDict()
relpoints[euv] = 25./325.           # 25
relpoints[stellar] = 100./325.     # 100
relpoints[aromatic] = 125./325.         # 125
relpoints[thermal] = 50./325.         # 50
relpoints[microwave] = 25./325.    # 25

# -----------------------------------------------------------------

class WavelengthGridsTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Label"] = (str, None, "wavelength grid label")
    _column_info["EUV points"] = (int, None, "number of points in EUV spectrum (range: " + str(ranges[euv]) + ")")
    _column_info["Stellar points"] = (int, None, "number of points in stellar spectrum (range: " + str(ranges[stellar]) + ")")
    _column_info["Aromatic points"] = (int, None, "number of points in aromatic spectrum (range: " + str(ranges[aromatic]) + ")")
    _column_info["Thermal points"] = (int, None, "number of points in thermal spectrum (range: " + str(ranges[thermal]) + ")")
    _column_info["Microwave points"] = (int, None, "number of points in microwave spectrum (range: " + str(ranges[microwave]) + ")")
    _column_info["Broad band filters"] = (str, None, "broad band filters for which the wavelength range was resampled")
    _column_info["Narrow band filters"] = (str, None, "narrow band filters for which the wavelength was added")
    _column_info["Adjusted points"] = (int, None, "number of points that were adjusted")
    _column_info["New points"] = (int, None, "number of new points")
    _column_info["Emission lines"] = (int, None, "number of emission lines")
    _column_info["Fixed points"] = (int, None, "number of fixed points")
    _column_info["Total points"] = (int, None, "total number of points")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(WavelengthGridsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

    # -----------------------------------------------------------------

    def add_grid(self, label, npoints, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added,
                 adjusted_npoints, new_npoints):

        """
        This function ...
        :param label:
        :param npoints:
        :param subgrid_npoints:
        :param emission_npoints:
        :param fixed_npoints:
        :param broad_resampled:
        :param narrow_added:
        :param adjusted_npoints:
        :param new_npoints:
        :return:
        """

        # Get values
        euv_npoints = subgrid_npoints[euv] if euv in subgrid_npoints else 0
        stellar_npoints = subgrid_npoints[stellar] if stellar in subgrid_npoints else 0
        aromatic_npoints = subgrid_npoints[aromatic] if aromatic in subgrid_npoints else 0
        thermal_npoints = subgrid_npoints[thermal] if thermal in subgrid_npoints else 0
        microwave_npoints = subgrid_npoints[microwave] if microwave in subgrid_npoints else 0

        #print("BROAD:", [str(fltr) for fltr in broad_resampled])

        # Create strings from the filter lists
        broad_string = ",".join([str(fltr) for fltr in broad_resampled])
        narrow_string = ",".join([str(fltr) for fltr in narrow_added])

        # Add row
        self.add_row([label, euv_npoints, stellar_npoints, aromatic_npoints, thermal_npoints, microwave_npoints, broad_string,
                      narrow_string, adjusted_npoints, new_npoints, emission_npoints, fixed_npoints, npoints])

    # -----------------------------------------------------------------

    def get_label(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self["Label"][index]

    # -----------------------------------------------------------------

    def get_index_for_label(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        from ..tools import tables
        return tables.find_index(self, label, column_name="Label")

    # -----------------------------------------------------------------

    def get_npoints_euv(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self["EUV points"][index]

    # -----------------------------------------------------------------

    def get_npoints_stellar(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self["Stellar points"][index]

    # -----------------------------------------------------------------

    def get_npoints_aromatic(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self["Aromatic points"][index]

    # -----------------------------------------------------------------

    def get_npoints_thermal(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self["Thermal points"][index]

    # -----------------------------------------------------------------

    def get_npoints_microwave(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self["Microwave points"][index]

    # -----------------------------------------------------------------

    def get_broad_band_filters(self, index):

        """
        Thisn function ...
        :param index:
        :return:
        """

        if self["Broad band filters"].mask[index]: return None
        return parsing.broad_band_filter_list(self["Broad band filters"][index])

    # -----------------------------------------------------------------

    def get_narrow_band_filters(self, index):

        """
        Thisn function ...
        :param index:
        :return:
        """

        if self["Narrow band filters"].mask[index]: return None
        return parsing.narrow_band_filter_list(self["Narrow band filters"][index])

    # -----------------------------------------------------------------

    def get_nemission_lines(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self["Emission lines"][index]

    # -----------------------------------------------------------------

    def get_npoints_fixed(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self["Fixed points"][index]

    # -----------------------------------------------------------------

    def get_npoints(self, index):

        """
        Thisn function ...
        :param index:
        :return:
        """

        return self["Total points"][index]

# -----------------------------------------------------------------

class WavelengthGridGenerator(Configurable):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(WavelengthGridGenerator, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The wavelength grid property table
        self.table = None

        # Lists to contain the elements of each grid
        self.npoints = []
        self.grids = []
        self.subgrids = []
        self.filter_wavelengths = []
        self.replaced = []
        self.new = []
        self.line_wavelengths = []
        self.fixed = []

        # SEDs for plotting
        self.seds = OrderedDict()

        # Out paths and plot paths
        self.out_paths = None
        self.plot_paths = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Generate the grids
        self.generate()

        # Show
        if self.config.show: self.show()

        # Plot
        if self.config.plot: self.plot()

        # Write
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    @property
    def ngrids(self):
        return len(self.grids)

    # -----------------------------------------------------------------

    @property
    def no_grids(self):
        return self.ngrids == 0

    # -----------------------------------------------------------------

    @property
    def single_grid(self):
        if self.no_grids: raise RuntimeError("No grids")
        elif self.ngrids == 1: return self.grids[0]
        else: raise RuntimeError("More than one grid")

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(WavelengthGridGenerator, self).setup(**kwargs)

        # Get the wavelength grids table
        self.table = kwargs.pop("table", None)
        if self.table is None and self.config.table: self.table = WavelengthGridsTable()

        # Get SEDS for plotting
        if kwargs.get("seds", None) is not None:
            seds = kwargs.pop("seds")
            for label in seds: self.add_sed(seds[label], label)

        # Create SEDs
        if self.config.plot_seds: self.create_seds()

        # Get out paths
        self.out_paths = kwargs.pop("out_paths", None)

        # Get plot paths
        self.plot_paths = kwargs.pop("plot_paths", None)

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

        self.seds.update(seds)

    # -----------------------------------------------------------------

    def create_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the suer
        log.info("Creating SED templates ...")

        # Create
        seds = create_template_seds(self.config.seds, ages=self.config.ages, metallicity=self.config.metallicity,
                                    compactness=self.config.compactness, pressure=self.config.pressure,
                                    covering_factor=self.config.covering_factor)

        # Add
        self.add_seds(seds)

    # -----------------------------------------------------------------

    @property
    def has_table(self):
        return self.table is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def min_wavelength(self):
        return get_min_wavelength(self.config.range.min, self.config.check_filters, self.config.adjust_minmax)

    # -----------------------------------------------------------------

    @lazyproperty
    def max_wavelength(self):
        return get_max_wavelength(self.config.range.max, self.config.check_filters, self.config.adjust_minmax)

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :param:
        """

        # Inform the user
        log.info("Generating the wavelength grids ...")

        # Loop over the different number of points
        for index, npoints in enumerate(self.config.npoints_range.linear(self.config.ngrids)):

            # Debugging
            log.debug("Creating a wavelength grid with a target of " + str(npoints) + " points ...")

            # Create the grid
            wavelength_grid, subgrid_wavelengths, filter_wavelengths, replaced, new, line_wavelengths, fixed = create_one_subgrid_wavelength_grid(
                npoints, self.emission_lines, self.config.fixed,
                min_wavelength=self.min_wavelength, max_wavelength=self.max_wavelength,
                filters=self.config.filters, adjust_to=self.config.adjust_to,
                min_wavelengths_in_filter=self.config.min_wavelengths_in_filter,
                min_wavelengths_in_fwhm=self.config.min_wavelengths_in_fwhm,
                return_elements=True)

            #print("Filters:", filter_wavelengths)

            # Add to lists
            self.npoints.append(npoints)
            self.grids.append(wavelength_grid)
            self.subgrids.append(subgrid_wavelengths)
            self.filter_wavelengths.append(filter_wavelengths)
            self.replaced.append(replaced)
            self.new.append(new)
            self.line_wavelengths.append(line_wavelengths)
            self.fixed.append(fixed)

            # Add entry to the table
            if self.has_table: self.add_to_table(index, npoints)

    # -----------------------------------------------------------------

    @lazyproperty
    def emission_lines(self):

        """
        This function ...
        :return:
        """

        # Create the emission lines instance
        if self.config.add_emission_lines:

            # Line IDs are specified
            if self.config.emission_lines is not None:

                lines = []

                # Loop over the line IDS
                for line_id in self.config.emission_lines:

                    # Create line
                    line = EmissionLine.from_string(line_id)

                    # Add line
                    lines.append(line)

                # Return the lines
                return lines

            # No lines are specified: take all
            else: return EmissionLines()

        # No emission lines to be used
        else: return None

    # -----------------------------------------------------------------

    def get_label(self, npoints):

        """
        This function ...
        :param npoints:
        :return:
        """

        # Determine the label
        return self.config.label + str(npoints) if self.config.label is not None else str(npoints)

    # -----------------------------------------------------------------

    def add_to_table(self, index, target_npoints):

        """
        This function ...
        :param index:
        :param target_npoints:
        :return:
        """

        # Debugging
        log.debug("Adding row to the wavelength grids table ...")

        # Get label (with target npoints)
        label = self.get_label(target_npoints)

        # Get properties
        npoints = self.get_npoints(index)
        subgrid_npoints = self.get_subgrid_npoints(index)
        emission_npoints = self.get_emission_npoints(index)
        fixed_npoints = self.get_fixed_npoints(index)
        broad_resampled = self.get_broad_resampled(index)
        narrow_added = self.get_narrow_added(index)
        adjusted_npoints = self.get_replaced_npoints(index)
        new_npoints = self.get_new_npoints(index)

        # Add to table
        self.table.add_grid(label, npoints, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled,
                            narrow_added, adjusted_npoints, new_npoints)

    # -----------------------------------------------------------------

    def get_target_npoints(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.npoints[index]

    # -----------------------------------------------------------------

    def get_grid(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.grids[index]

    # -----------------------------------------------------------------

    def get_npoints(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return len(self.get_grid(index))

    # -----------------------------------------------------------------

    def get_subgrids(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        subgrid_wavelengths = self.subgrids[index]

        grids = OrderedDict()
        for name in subgrid_wavelengths:
            grid = WavelengthGrid.from_wavelengths(subgrid_wavelengths[name], sort=True)
            grids[name] = grid
        return grids

    # -----------------------------------------------------------------

    def get_subgrid_npoints(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        subgrid_wavelengths = self.subgrids[index]

        # Keep track of the number of points per subgrid
        subgrid_npoints = OrderedDict()
        for subgrid in subgrid_wavelengths:
            subgrid_npoints[subgrid] = len(subgrid_wavelengths[subgrid])
        return subgrid_npoints

    # -----------------------------------------------------------------

    def get_emission_npoints(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Set emission npoints
        emission_npoints = 0
        for line_identifier in self.line_wavelengths[index]:
            emission_npoints += len(self.line_wavelengths[index][line_identifier])
        return emission_npoints

    # -----------------------------------------------------------------

    def get_fixed_npoints(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return len(self.fixed[index])

    # -----------------------------------------------------------------

    def has_fixed(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_fixed_npoints(index) > 0

    # -----------------------------------------------------------------

    def get_fixed_grid(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return WavelengthGrid.from_wavelengths(self.fixed[index])

    # -----------------------------------------------------------------

    def get_nfilters(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return len(self.filter_wavelengths[index])

    # -----------------------------------------------------------------

    def has_filters(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_nfilters(index) > 0

    # -----------------------------------------------------------------

    def get_broad_resampled(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Set broad resampled nwavelengths
        broad_resampled = []
        for fltr in self.filter_wavelengths[index]:
            if not isinstance(fltr, BroadBandFilter): continue
            broad_resampled.append(fltr)
        return broad_resampled

    # -----------------------------------------------------------------

    def get_narrow_added(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Set narrow band filter added nwavelengths
        narrow_added = []
        for fltr in self.filter_wavelengths[index]:
            if not isinstance(fltr, NarrowBandFilter): continue
            narrow_added.append(fltr)
        return narrow_added

    # -----------------------------------------------------------------

    def get_nreplaced(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return len(self.replaced[index])

    # -----------------------------------------------------------------

    get_replaced_npoints = get_nreplaced

    # -----------------------------------------------------------------

    def has_replaced(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_nreplaced(index) > 0

    # -----------------------------------------------------------------

    def get_replaced(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.replaced[index]

    # -----------------------------------------------------------------

    def get_nnew(self, index):

        """
        Thisn function ...
        :param index:
        :return:
        """

        return len(self.new[index])

    # -----------------------------------------------------------------

    get_new_npoints = get_nnew

    # -----------------------------------------------------------------

    def has_new(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_nnew(index) > 0

    # -----------------------------------------------------------------

    def get_new(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.new[index]

    # -----------------------------------------------------------------

    def get_nlines(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return len(self.line_wavelengths[index])

    # -----------------------------------------------------------------

    def has_line_wavelengths(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        return self.get_nlines(index) > 0

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the wavelength grids ...")

        # Loop over the grids
        for index in range(self.ngrids):

            # Get properties
            grid = self.get_grid(index)
            npoints = self.get_npoints(index)

            # Show
            print("Wavelength grid with " + str(npoints) + " wavelength points:")
            print("")

            print(" - number of points in subgrids:")
            subgrid_npoints = self.get_subgrid_npoints(index)
            for subgrid in subgrid_npoints: print("    * " + subgrid + ": " + str(subgrid_npoints[subgrid]))
            print(" - number of emission points: " + str(self.get_emission_npoints(index)))
            print(" - number of fixed points: " + str(self.get_fixed_npoints(index)))
            print(" - filters for which extra sampling was performed: " + str(self.get_broad_resampled(index)))
            print(" - narrow band filters for which wavelength was added: " + str(self.get_narrow_added(index)))
            if self.has_replaced(index): print(" - replaced wavelengths:")
            for old, new in self.get_replaced(index): print("    * " + str(old) + " -> " + str(new))
            if self.has_new(index): print(" - new wavelengths:")
            for line in stringify_list_fancy(self.get_new(index))[1].split("\n"): print("    " + line)
            print("")
            print(grid)
            print("")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the wavelength grids ...")

        from ..plot.wavelengthgrid import WavelengthGridPlotter

        # Loop over the grids
        for index in range(self.ngrids):

            # Get label
            target_npoints = self.get_target_npoints(index)
            label = self.get_label(target_npoints)

            # Get the elements of this grid
            grid = self.get_grid(index)
            subgrids = self.get_subgrids(index)
            filter_wavelengths = self.filter_wavelengths[index]
            replaced = self.replaced[index]
            new = self.new[index]
            line_wavelengths = self.line_wavelengths[index]
            fixed = self.fixed[index]

            # Create the plotter
            plotter = WavelengthGridPlotter()

            # Set settings
            plotter.config.add_regimes = self.config.plot_regimes
            plotter.config.regimes = self.config.regimes
            plotter.config.only_subregimes = self.config.only_subregimes

            #print("FILTERS", self.config.filters)
            plotter.config.add_filters = self.config.plot_filters
            if self.config.plotting_filters is not None: plotter.config.filters = self.config.plotting_filters
            else: plotter.config.filters = self.config.filters
            plotter.config.categorize_filters = self.config.categorize_filters

            # Lines
            plotter.config.add_lines = self.config.plot_lines
            plotter.config.lines = self.config.emission_lines

            # Plot separate
            plotter.config.separate_grids = True
            plotter.config.group_wavelengths = True
            #plotter.config.lines_in_group = "lines"
            plotter.config.separate_lines = True
            plotter.config.mark_removed = True

            # Plot resampled and residuals
            plotter.config.plot_resampled = self.config.plot_resampled
            plotter.config.plot_interpolated = self.config.plot_interpolated
            plotter.config.plot_residuals = self.config.plot_residuals

            # Add the elements
            plotter.add_elements(subgrids, fixed=fixed, filter_wavelengths=filter_wavelengths, replaced=replaced, new=new, line_wavelengths=line_wavelengths)

            # Add complete grid
            if self.config.plot_reference:
                plotter.add_reference_grid(grid, label="reference", in_legend=True)
                plotter.config.plot_differences = True

            # Determine plot filepath
            if self.plot_paths is not None and target_npoints in self.plot_paths: plot_filepath = fs.join(self.plot_paths[target_npoints], "grid.pdf")
            elif self.config.plot_path is not None: plot_filepath = fs.join(self.config.plot_path, label + ".pdf")
            else: plot_filepath = None

            # Run the plotter
            plotter.run(output=plot_filepath, seds=self.seds, min_wavelength=self.min_wavelength, max_wavelength=self.max_wavelength)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the grids
        if self.config.write_grids: self.write_grids()

        # Write the subgrids
        if self.config.write_elements: self.write_subgrids()

        # Write the grids of fixed wavelengths
        if self.config.write_elements: self.write_fixed()

        # Write filter wavelengths
        if self.config.write_elements: self.write_filter_wavelengths()

        # Write replaced wavelengths
        if self.config.write_elements: self.write_replaced()

        # Write new wavelengths
        if self.config.write_elements: self.write_new()

        # Write emission line wavelengths
        if self.config.write_elements: self.write_line_wavelengths()

        # Write the table
        if self.has_table and self.config.write_table: self.write_table()

    # -----------------------------------------------------------------

    def write_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the grids ...")

        # Loop over the grids
        for index in range(self.ngrids):

            # Get target npoints
            target_npoints = self.get_target_npoints(index)

            # Debugging
            log.debug("Writing the " + str(target_npoints) + " points wavelength grid ...")

            # Get label (with target npoints)
            label = self.get_label(target_npoints)

            # Get wavelength grid
            wavelength_grid = self.grids[index]

            # Determine filepath
            if self.out_paths is not None and target_npoints in self.out_paths: path = fs.join(self.out_paths[target_npoints], "grid.dat")
            else: path = self.output_path_file(label + ".dat")

            # Write the wavelength grid
            wavelength_grid.saveto(path)

    # -----------------------------------------------------------------

    def write_subgrids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the subgrids ...")

        # Loop over the grids
        for index in range(self.ngrids):

            # Get target npoints
            target_npoints = self.get_target_npoints(index)

            # Debugging
            log.debug("Writing the subgrids of the " + str(target_npoints) + " points wavelength grid ...")

            # Get label (with target npoints)
            label = self.get_label(target_npoints)

            # Get the elements
            subgrids = self.get_subgrids(index)

            # Determine directory path
            if self.out_paths is not None and target_npoints in self.out_paths: elements_path = self.out_paths[target_npoints]
            else:  elements_path = self.output_path_directory(label, create=True)

            # Loop over the subgrids
            for name in subgrids:

                # Debugging
                log.debug("Writing the '" + name + "' subgrid ...")

                # Get path
                path = fs.join(elements_path, "subgrid_" + name + ".dat")

                # Write subgrid
                subgrids[name].saveto(path)

    # -----------------------------------------------------------------

    def write_fixed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing fixed wavelengths ...")

        # Loop over the grids
        for index in range(self.ngrids):

            # Has fixed?
            if not self.has_fixed(index): continue

            # Get target npoints
            target_npoints = self.get_target_npoints(index)

            # Debugging
            log.debug("Writing the fixed wavelengths of the " + str(target_npoints) + " points wavelength grid ...")

            # Get label (with target npoints)
            label = self.get_label(target_npoints)

            # Get the fixed grid
            fixed_grid = self.get_fixed_grid(index)

            # Determine directory path
            if self.out_paths is not None and target_npoints in self.out_paths: elements_path = self.out_paths[target_npoints]
            else: elements_path = self.output_path_directory(label, create=True)

            # Determine filepath
            path = fs.join(elements_path, "fixed_grid.dat")

            # Write the grid
            fixed_grid.saveto(path)

    # -----------------------------------------------------------------

    def write_filter_wavelengths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the filter wavelengths ...")

        # Loop over the grids
        for index in range(self.ngrids):

            # Has filter wavelengths?
            if not self.has_filters(index): continue

            # Get target npoints
            target_npoints = self.get_target_npoints(index)

            # Debugging
            log.debug("Writing the filter wavelengths of the " + str(target_npoints) + " wavelength grid ...")

            # Get label (with target npoints)
            label = self.get_label(target_npoints)

            # Get the elements
            filter_wavelengths = self.filter_wavelengths[index]

            # Determine directory path
            if self.out_paths is not None and target_npoints in self.out_paths: elements_path = self.out_paths[target_npoints]
            else: elements_path = self.output_path_directory(label, create=True)

            # Determine filepath
            path = fs.join(elements_path, "filter_wavelengths.dat")

            # Write
            write_dict(filter_wavelengths, path)

    # -----------------------------------------------------------------

    def write_replaced(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing replaced wavelengths ...")

        # Loop over the grids
        for index in range(self.ngrids):

            # Has replaced?
            if not self.has_replaced(index): continue

            # Get target npoints
            target_npoints = self.get_target_npoints(index)

            # Debugging
            log.debug("Writing replaced wavelengths of the " + str(target_npoints) + " points wavelength grid ...")

            # Get label (with target npoints)
            label = self.get_label(target_npoints)

            # Get the elements
            replaced = self.replaced[index]

            # Determine directory path
            if self.out_paths is not None and target_npoints in self.out_paths: elements_path = self.out_paths[target_npoints]
            else: elements_path = self.output_path_directory(label, create=True)

            # Determine filepath
            path = fs.join(elements_path, "replaced.dat")

            # Write
            write_list(replaced, path)

    # -----------------------------------------------------------------

    def write_new(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing new wavelengths ...")

        # Loop over the grids
        for index in range(self.ngrids):

            # Has new?
            if not self.has_new(index): continue

            # Get target npoints
            target_npoints = self.get_target_npoints(index)

            # Debugging
            log.debug("Writing new wavelengths of the " + str(target_npoints) + " points wavelength grid ...")

            # Get label (with target npoints)
            label = self.get_label(target_npoints)

            # Get the elements
            new = self.new[index]

            # Determine directory path
            if self.out_paths is not None and target_npoints in self.out_paths: elements_path = self.out_paths[target_npoints]
            else: elements_path = self.output_path_directory(label, create=True)

            # Determine filepath
            path = fs.join(elements_path, "new.dat")

            # Write
            write_list(new, path)

    # -----------------------------------------------------------------

    def write_line_wavelengths(self):

        """
        Thins function ...
        :return:
        """

        # Inform the user
        log.info("Writing line wavelengths ...")

        # Loop over the grids
        for index in range(self.ngrids):

            # Has line wavelengths?
            if not self.has_line_wavelengths(index): continue

            # Get target npoints
            target_npoints = self.get_target_npoints(index)

            # Debugging
            log.debug("Writing line wavelengths of the " + str(target_npoints) + " points wavelength grid ...")

            # Get label (with target npoints)
            label = self.get_label(target_npoints)

            # Get the line wavelengths
            line_wavelengths = self.line_wavelengths[index]

            # Determine directory path
            if self.out_paths is not None and target_npoints in self.out_paths: elements_path = self.out_paths[target_npoints]
            else: elements_path = self.output_path_directory(label, create=True)

            # Determine file path
            path = fs.join(elements_path, "line_wavelengths.dat")

            # Write
            write_dict(line_wavelengths, path)

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grids table ...")

        # Determine the path
        path = self.output_path_file("table.dat")

        # Save the table
        self.table.saveto(path)

# -----------------------------------------------------------------

def get_subgrid_wavelengths(npoints, min_wavelength=None, max_wavelength=None):

    """
    This function ...
    :param min_wavelength:
    :param max_wavelength:
    :param npoints:
    :return:
    """

    subgrid_wavelengths = OrderedDict()

    # Loop over the subgrids
    for subgrid in subgrids:

        # Debugging
        log.debug("Checking " + subgrid + " subgrid wavelengths ...")

        # Determine minimum, maximum
        min_lambda = ranges[subgrid].min
        max_lambda = ranges[subgrid].max

        # Skip subgrids out of range
        if min_wavelength is not None and max_lambda < min_wavelength: continue
        if max_wavelength is not None and min_lambda > max_wavelength: continue

        # Determine the normal number of wavelength points for this subgrid
        points = int(round(relpoints[subgrid] * npoints))

        # Generate and add the wavelength points
        wavelengths_subgrid_original = make_grid(min_lambda, max_lambda, points)

        # Filter based on given boundaries
        wavelengths_subgrid = []
        for wav in wavelengths_subgrid_original:
            if min_wavelength is not None and wav < min_wavelength: continue
            if max_wavelength is not None and wav > max_wavelength: continue
            wavelengths_subgrid.append(wav)

        # Add the sequence of wavelengths
        subgrid_wavelengths[subgrid] = wavelengths_subgrid

    # Return
    return subgrid_wavelengths

# -----------------------------------------------------------------

def resample_filter_wavelengths(wavelengths, filters, min_wavelengths_in_filter=5, min_wavelengths_in_fwhm=3):

    """
    This function ...
    :param wavelengths:
    :param filters:
    :param min_wavelengths_in_filter:
    :param min_wavelengths_in_fwhm:
    :return:
    """

    # Initialize dictionary for the wavelengths for each filter
    filter_wavelengths = DefaultOrderedDict(list)

    # Debugging
    log.debug("Adding wavelengths for sampling filter bandpasses ...")

    # Keep track of the wavelengths that need to remain exactly as they are
    exact_wavelengths = defaultdict(list)

    # Loop over the filters
    for fltr in filters:

        # Debugging
        log.debug("Adding wavelength(s) for the " + str(fltr) + " filter ...")

        # Broad band filter: make sure there are at least 10 wavelength points in the range min_wavelength > max_wavelength
        if isinstance(fltr, BroadBandFilter):

            # Get filter minimum and maximum wavelength
            min_wavelength = fltr.min
            max_wavelength = fltr.max

            # Get wavelength of maximum transmission

            # Check that at least 5 wavelength points sample the range of the filter
            # Get the indices of the wavelengths that fall within this range in the current list of wavelengths
            current_indices = [i for i in range(len(wavelengths)) if min_wavelength < wavelengths[i] < max_wavelength]

            # Check if there at least ..
            #if len(current_indices) >= min_wavelengths_in_filter: continue

            if len(current_indices) < min_wavelengths_in_filter:

                # Otherwise, delete the current wavelengths and add 10 new ones
                # NO: DON'T DELETE JUST YET: DON'T DELETE WAVELENGTHS OF OTHER FILTERS
                # FIRST DELETE AT THE END AND THEN ADD THEM ALL
                #for index in sorted(current_indices, reverse=True): del wavelengths[index]

                # Generate new wavelengths for sampling the filter range on a logarithmic grid
                new_wavelengths = fltr.range.log(min_wavelengths_in_filter, as_list=True)

                #print(fltr, fltr.range, [str(wavelength) for wavelength in new_wavelengths])

                # Add the new wavelengths
                # NO: DON'T ADD JUST YET: DON'T DELETE WAVELENGTHS OF OTHER FILTERS
                # FIRST DELETE AT THE END AND THEN ADD THEM ALL
                #wavelengths += new_wavelengths

                # Add the wavelengths
                filter_wavelengths[fltr].extend(new_wavelengths)

            else: new_wavelengths = None

            # If FWHM of the filter is defined
            if fltr.has_fwhm:

                # Check that at least 3 wavelength points sample the inner range of the filter
                current_indices = [i for i in range(len(wavelengths)) if wavelengths[i] in fltr.fwhm_range]
                if new_wavelengths is not None: current_indices_new = [i for i in range(len(new_wavelengths)) if new_wavelengths[i] in fltr.fwhm_range]
                else: current_indices_new = None

                # Check if there are at least 3
                #if len(current_indices) >= min_wavelengths_in_fwhm: continue

                if len(current_indices) < min_wavelengths_in_fwhm:

                    # Otherwise, delete the current wavelengths and add 3 new ones
                    # NO: DON'T DELETE JUST YET: DON'T DELETE WAVELENGTHS OF OTHER FILTERS
                    # FIRST DELETE AT THE END AND THEN ADD THEM ALL
                    #for index in sorted(current_indices, reverse=True): del wavelengths[index]

                    # Generate new wavelengths for sampling the inner filter range on a logarithmic grid
                    new_wavelengths = fltr.fwhm_range.log(min_wavelengths_in_fwhm, as_list=True)

                    #print("(FWHM)", fltr, fltr.fwhm_range, [str(wavelength) for wavelength in new_wavelengths])

                    # Add the new wavelengths
                    # NO: DON'T ADD JUST YET: DON'T DELETE WAVELENGTHS OF OTHER FILTERS
                    # FIRST DELETE AT THE END AND THEN ADD THEM ALL
                    #wavelengths += new_wavelengths

                    # Add to dictionary
                    if current_indices_new is not None:
                        for index in sorted(current_indices_new, reverse=True): del filter_wavelengths[fltr][index]
                    filter_wavelengths[fltr].extend(new_wavelengths)

            #print(str(fltr), fltr.peak)

            # NEW: Make sure the peak wavelength is included
            if fltr.peak is not None:

                # New wavelengths have been added
                filter_wavelengths[fltr].extend([]) # make sure there is a list
                if len(filter_wavelengths[fltr]) > 0:

                    # Find closest
                    index = sequences.find_closest_index(filter_wavelengths[fltr], fltr.peak)
                    #all_index = sequences.find_closest_index(wavelengths, fltr.peak)

                    # Replace
                    filter_wavelengths[fltr][index] = fltr.peak
                    #wavelengths[all_index] = fltr.peak

                # Nothing has been added
                else:

                    #wavelengths.append(fltr.peak)
                    filter_wavelengths[fltr].append(fltr.peak)

                # ADD AS EXACT WAVELENGTH
                exact_wavelengths[fltr].append(fltr.peak)

        # For a narrow band filter, add the exact wavelength of the filter to the wavelength grid
        elif isinstance(fltr, NarrowBandFilter):

            # Add the wavelength
            #wavelengths.append(fltr.wavelength)

            # One more filter for which we have added a wavelength
            #narrow_added.append(str(fltr))

            # Add the wavelength
            filter_wavelengths[fltr].append(fltr.wavelength)

            # Add as exact wavelength
            exact_wavelengths[fltr].append(fltr.wavelength)

        # Unrecognized filter
        else: raise ValueError("Unrecognized filter object: " + str(fltr))

    # DELETE WAVELENGTHS
    for fltr in filter_wavelengths:

        fltr_wavelengths = filter_wavelengths[fltr]
        nwavelengths = len(fltr_wavelengths)
        if nwavelengths == 1: continue

        min_wavelength = min(fltr_wavelengths)
        max_wavelength = max(fltr_wavelengths)

        #print(str(fltr), min_wavelength, max_wavelength)

        # DELETE WAVELENGTHS BETWEEN
        # Get indices of wavelengths between
        indices = [i for i in range(len(wavelengths)) if min_wavelength < wavelengths[i] < max_wavelength]

        # Remove indices
        for index in sorted(indices, reverse=True): del wavelengths[index]

    # ADD FILTER WAVELENGTHS
    for fltr in filter_wavelengths: wavelengths.extend(filter_wavelengths[fltr])

    # Return the filter wavelengths and exact wavelength
    return filter_wavelengths, exact_wavelengths

# -----------------------------------------------------------------

def adjust_to_wavelengths(wavelengths, adjust_to, keep=None):

    """
    This function ...
    :param wavelengths:
    :param adjust_to:
    :param keep:
    :return:
    """

    # Debugging
    log.debug("Adjusting wavelength grid to include exact wavelengths ...")

    # Replace dict
    replace_dict = defaultdict(list)

    # Original wavelengths with replaced
    replaced = []

    # New wavelengths
    new = []

    # Loop over the wavelengths
    for wavelength in adjust_to:

        # Debugging
        log.debug("Adjusting the wavelength grid for " + str(wavelength) + "...")

        # Find closest
        index = sequences.find_closest_index(wavelengths, wavelength)

        # Check if wavelength needs to be kept
        if keep is not None and wavelengths[index] in keep:
            wavelengths.append(wavelength)
            new.append(wavelength)

        # Add wavelength to be replaced
        else: replace_dict[index].append(wavelength)

    # Loop over the indices for which a replacement will be performed
    for index in replace_dict:

        # Only one replacement
        if len(replace_dict[index]) == 1:

            # Get the new wavelength
            wavelength = replace_dict[index][0]

            # Set
            replaced.append((wavelengths[index], wavelength))

            # Replace the wavelength
            wavelengths[index] = wavelength

        # Multiple wavelengths for one index
        else:

            # Get the wavelengths
            adjust_wavelengths = replace_dict[index]

            # Get the closest wavelength to the previous wavelength
            previous_wavelength = wavelengths[index]

            # Find closest
            adjust_index = sequences.find_closest_index(adjust_wavelengths, previous_wavelength)

            # Set
            replaced.append((wavelengths[index], adjust_wavelengths[adjust_index]))

            # Replace closest
            wavelengths[index] = adjust_wavelengths[adjust_index]

            # Add other
            for j in range(len(adjust_wavelengths)):

                # Already added
                if j == adjust_index: continue

                # Add the wavelength
                wavelengths.append(adjust_wavelengths[j])

                # Add to new
                new.append(adjust_wavelengths[j])

    # Return
    return replaced, new

# -----------------------------------------------------------------

def create_one_subgrid_wavelength_grid(npoints, emission_lines=None, fixed=None, min_wavelength=None, max_wavelength=None,
                                       filters=None, min_wavelengths_in_filter=5, min_wavelengths_in_fwhm=3, adjust_to=None,
                                       return_elements=False):

    """
    This function ...
    :param npoints:
    :param emission_lines:
    :param fixed:
    :param min_wavelength:
    :param max_wavelength:
    :param filters:
    :param min_wavelengths_in_filter:
    :param min_wavelengths_in_fwhm:
    :param adjust_to:
    :param return_elements:
    :return:
    """

    # Debugging
    log.debug("Creating wavelength grid with " + str(npoints) + " points ...")

    # A list of the wavelength points
    wavelengths = []

    # Get subgrid wavelength sequences
    subgrid_wavelengths = get_subgrid_wavelengths(npoints, min_wavelength=min_wavelength, max_wavelength=max_wavelength)

    # Loop over the subgrids
    for subgrid in subgrid_wavelengths:

        # Debugging
        log.debug("Adding the " + subgrid + " subgrid ...")

        # Add the wavelength points
        wavelengths += subgrid_wavelengths[subgrid]

    # Loop over the filters
    if filters is not None: filter_wavelengths, exact_filter_wavelengths = resample_filter_wavelengths(wavelengths, filters, min_wavelengths_in_filter=min_wavelengths_in_filter, min_wavelengths_in_fwhm=min_wavelengths_in_fwhm)
    else:
        filter_wavelengths = dict()
        exact_filter_wavelengths = dict()

    # Get a list of all the exact wavelengths that can't be adjusted
    exact_wavelengths = []
    for fltr in exact_filter_wavelengths: exact_wavelengths.extend(exact_filter_wavelengths[fltr])

    # Adjust to passed wavelengths
    if adjust_to is not None: replaced, new = adjust_to_wavelengths(wavelengths, adjust_to, keep=exact_wavelengths)
    else: replaced = new = []

    # Add the emission lines
    if emission_lines is not None:

        # Debugging
        log.debug("Adding the emission lines ...")

        # Add the mission lines
        wavelengths, line_wavelengths = added_emission_lines(wavelengths, emission_lines, min_wavelength, max_wavelength, return_added=True)

    # No emission lines
    else: line_wavelengths = dict()

    # Add fixed wavelength points
    new_fixed = []
    if fixed is not None:

        # Debugging
        log.debug("Adding the fixed points to the grid ...")

        # Loop over the wavelengths
        for wavelength in fixed:

            # Check if not already there (e.g. from adjust_to)
            if wavelength in wavelengths: continue

            # Debugging
            log.debug("Adding fixed wavelength " + str(wavelength) + " ...")

            if min_wavelength is not None and wavelength < min_wavelength: continue
            if max_wavelength is not None and wavelength > max_wavelength: continue
            wavelengths.append(wavelength)

            # Add to fixed
            new_fixed.append(wavelength)

    #print(wavelengths)

    # Sort the wavelength points
    wavelengths = sorted(wavelengths)

    # Create the wavelength grid
    grid = WavelengthGrid.from_wavelengths(wavelengths)

    # Return the grid and some information about the subgrids
    if return_elements: return grid, subgrid_wavelengths, filter_wavelengths, replaced, new, line_wavelengths, new_fixed
    else:

        # Keep track of the number of points per subgrid
        subgrid_npoints = OrderedDict()
        for subgrid in subgrid_wavelengths:
            subgrid_npoints[subgrid] = len(subgrid_wavelengths[subgrid])

        # Set broad resampled nwavelengths
        broad_resampled = []
        for fltr in filter_wavelengths:
            if not isinstance(fltr, BroadBandFilter): continue
            broad_resampled.append(fltr)

        # Set narrow band filter added nwavelengths
        narrow_added = []
        for fltr in filter_wavelengths:
            if not isinstance(fltr, NarrowBandFilter): continue
            narrow_added.append(fltr)

        # Set emission npoints
        emission_npoints = 0
        for line_identifier in line_wavelengths:
            emission_npoints += len(line_wavelengths[line_identifier])

        # Set fixed npoints
        fixed_npoints = len(new_fixed)

        # Return
        return grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added, replaced, new

# -----------------------------------------------------------------

def create_one_logarithmic_wavelength_grid(wrange, npoints, emission_lines=None, fixed=None):

    """
    This function ...
    :param wrange:
    :param npoints:
    :param emission_lines:
    :param fixed:
    :return:
    """

    # Verify the grid parameters
    if npoints < 2: raise ValueError("the number of points in the grid should be at least 2")
    if wrange.min <= 0: raise ValueError("the shortest wavelength should be positive")

    # Calculate log of boundaries
    logmin = np.log10(float(wrange.min))
    logmax = np.log10(float(wrange.max))

    # Calculate the grid points
    wavelengths = np.logspace(logmin, logmax, num=npoints, endpoint=True, base=10)

    # Add the emission lines
    emission_npoints = 0
    if emission_lines is not None:

        # Add the mission lines
        wavelengths = added_emission_lines(wavelengths, emission_lines)
        emission_npoints = len(emission_lines)

    # Add fixed wavelength points
    fixed_npoints = 0
    if fixed is not None:
        fixed_npoints = len(fixed)
        for wavelength in fixed: wavelengths.append(wavelength)

    # Sort the wavelength points
    wavelengths = sorted(wavelengths)

    # Create the wavelength grid
    grid = WavelengthGrid.from_wavelengths(wavelengths)

    # Return the grid
    return grid, emission_npoints, fixed_npoints

# -----------------------------------------------------------------

def create_one_nested_log_wavelength_grid(wrange, npoints, wrange_zoom, npoints_zoom, emission_lines=None, fixed=None):

    """
    This function ...
    :param wrange:
    :param npoints:
    :param wrange_zoom:
    :param npoints_zoom:
    :param emission_lines:
    :param fixed:
    :return:
    """

    # Verify the grid parameters
    if npoints < 2: raise ValueError("the number of points in the low-resolution grid should be at least 2")
    if npoints_zoom < 2: raise ValueError("the number of points in the high-resolution subgrid should be at least 2")
    if wrange.min <= 0: raise ValueError("the shortest wavelength should be positive")
    if (wrange_zoom.min <= wrange.min
        or wrange_zoom.max <= wrange_zoom.min
        or wrange.max <= wrange_zoom.max):
        raise ValueError("the high-resolution subgrid should be properly nested in the low-resolution grid")

    logmin = np.log10(float(wrange.min))
    logmax = np.log10(float(wrange.max))
    logmin_zoom = np.log10(float(wrange_zoom.min))
    logmax_zoom = np.log10(float(wrange_zoom.max))

    # Build the high- and low-resolution grids independently
    base_grid = np.logspace(logmin, logmax, num=npoints, endpoint=True, base=10)
    zoom_grid = np.logspace(logmin_zoom, logmax_zoom, num=npoints_zoom, endpoint=True, base=10)

    # Merge the two grids
    wavelengths = []

    # Add the wavelengths of the low-resolution grid before the first wavelength of the high-resolution grid
    for wavelength in base_grid:
        if wavelength < wrange_zoom.min: wavelengths.append(wavelength)

    # Add the wavelengths of the high-resolution grid
    for wavelength in zoom_grid: wavelengths.append(wavelength)

    # Add the wavelengths of the low-resolution grid after the last wavelength of the high-resolution grid
    for wavelength in base_grid:
        if wavelength > wrange_zoom.max: wavelengths.append(wavelength)

    # Add the emission lines
    emission_npoints = 0
    if emission_lines is not None:

        # Add the mission lines
        wavelengths = added_emission_lines(wavelengths, emission_lines)
        emission_npoints = len(emission_lines)

    # Add fixed wavelength points
    fixed_npoints = 0
    if fixed is not None:
        fixed_npoints = len(fixed)
        for wavelength in fixed: wavelengths.append(wavelength)

    # Sort the wavelength points
    wavelengths = sorted(wavelengths)

    # Create the wavelength grid
    grid = WavelengthGrid.from_wavelengths(wavelengths)

    # Return the grid
    return grid, emission_npoints, fixed_npoints

# -----------------------------------------------------------------

def added_emission_lines(wavelengths, emission_lines, min_wavelength=None, max_wavelength=None, return_added=False):

    """
    This function ...
    :param wavelengths:
    :param emission_lines:
    :param min_wavelength:
    :param max_wavelength:
    :param return_added:
    :return:
    """

    # Initialize dictionary to contain the wavelengths added to the grid for each emission line
    line_wavelengths = DefaultOrderedDict(list)

    # Add emission line grid points
    logdelta = 0.001
    for line in emission_lines:

        center_micron = line.center.to("micron").value
        left_micron = line.left.to("micron").value
        right_micron = line.right.to("micron").value

        if min_wavelength is not None and line.center < min_wavelength: continue
        if max_wavelength is not None and line.center > max_wavelength: continue

        # logcenter = np.log10(center)
        logleft = np.log10(left_micron if left_micron > 0 else center_micron) - logdelta
        logright = np.log10(right_micron if right_micron > 0 else center_micron) + logdelta

        newgrid = []

        for w in wavelengths:

            logw = np.log10(w.to("micron").value)
            if logw < logleft or logw > logright: newgrid.append(w)

        newgrid.append(line.center)
        line_wavelengths[line.identifier].append(line.center)

        if left_micron > 0:
            newgrid.append(line.left)
            line_wavelengths[line.identifier].append(line.left)

        if right_micron > 0:
            newgrid.append(line.right)
            line_wavelengths[line.identifier].append(line.right)

        wavelengths = newgrid

    # Return the new wavelength list
    if return_added: return wavelengths, line_wavelengths
    else: return wavelengths

# -----------------------------------------------------------------

def make_grid(wmin, wmax, N):

    """
    This function returns a wavelength grid (in micron) with a given resolution (nr of points per decade)
    # in the specified range (in micron), aligned with the 10^n grid points.
    """

    result = []

    wmin = wmin.to("micron").value
    wmax = wmax.to("micron").value

    # generate wavelength points p on a logarithmic scale with lambda = 10**p micron
    #  -2 <==> 0.01
    #   4 <==> 10000
    for i in range(-2*N,4*N+1):
        p = float(i)/N
        w = 10.**p
        if wmin <= w < wmax: result.append(w * u("micron"))

    # Return the grid
    return result

# -----------------------------------------------------------------

def get_min_wavelength(minimum, check_filters=None, adjust_minmax=False):

    """
    This function checks and gets the minimum wavelength for a grid
    :param minimum:
    :param check_filters:
    :param adjust_minmax:
    :return:
    """

    # Check filters?
    if check_filters is not None:

        # Get specified minimum wavelength
        min_wavelength = minimum

        # Loop over the filters
        for fltr in check_filters:

            # Check below
            if fltr.wavelength < min_wavelength:

                # Warning
                log.warning("The wavelength range does not contain the wavelength of the '" + str(fltr) + "' filter")

                # Adjust?
                if adjust_minmax:
                    log.debug("Adjusting the minimum wavelength to incorporate the '" + str(fltr) + "' filter")
                    min_wavelength = 0.99 * fltr.wavelength

        # Return the lower wavelength
        return min_wavelength

    # Return the specified minimum wavelength
    else: return minimum

# -----------------------------------------------------------------

def get_max_wavelength(maximum, check_filters=None, adjust_minmax=False):

    """
    This function ...
    :param maximum:
    :param check_filters:
    :param adjust_minmax:
    :return:
    """

    # Check filters?
    if check_filters is not None:

        # Get specified maximum wavelength
        max_wavelength = maximum

        # Loop over the filters
        for fltr in check_filters:

            # Check above
            if fltr.wavelength > max_wavelength:

                # Warning
                log.warning("The wavelength range does not contain the wavelength of the '" + str(fltr) + "' filter")

                # Adjust?
                if adjust_minmax:
                    log.debug("Adjusting the maximum wavelength to incorporate the '" + str(fltr) + "' filter")
                    max_wavelength = 1.01 * fltr.wavelength

        # Return the higher wavelength
        return max_wavelength

    # Return the specified maximum wavelength
    else: return maximum

# -----------------------------------------------------------------

def create_template_seds(names, **kwargs):

    """
    This function ...
    :param names:
    :param kwargs:
    :return:
    """

    from ..plot.wavelengthgrid import get_sed_template

    # Get properties
    ages = kwargs.pop("ages", None)
    metallicity = kwargs.pop("metallicity", None)
    compactness = kwargs.pop("compactness", None)
    pressure = kwargs.pop("pressure", None)
    covering_factor = kwargs.pop("covering_factor", None)

    # Set defaults
    if ages is None: ages = [q("8 Gyr"), q("0.1 Gyr")]
    if metallicity is None: metallicity = 0.02
    if compactness is None: compactness = 5.5
    if pressure is None: pressure = q("1e12 K/m3")
    if covering_factor is None: covering_factor = 0.2

    # Initialize dictionary
    seds = OrderedDict()

    # Loop over the template names
    for name in names:

        # MAPPINGS
        if name == "mappings":

            # Debugging
            log.debug("Creating MAPPINGS SED template ...")

            properties = dict()
            properties["metallicity"] = metallicity
            properties["compactness"] = compactness
            properties["pressure"] = pressure
            properties["covering_factor"] = covering_factor

            # Set label
            label = "MAPPINGS"

            sed = get_sed_template(name, **properties)
            #self.add_sed(sed, label=label)
            seds[label] = sed

        # Stellar Bruzual Charlot
        elif name == "bruzual_charlot":

            # Debugging
            log.debug("Creating Bruzual-Charlot SED templates ...")

            # Loop over the ages
            #if ages is None: raise ValueError("Ages are not specified")
            for age in ages:

                properties = dict()
                properties["metallicity"] = metallicity
                properties["age"] = age

                # label = name + "_" + str(age).replace(" ", "")
                label = "Bruzual-Charlot " + str(age)
                sed = get_sed_template(name, **properties)
                #self.add_sed(sed, label=label)
                seds[label] = sed

        # Invalid
        else: raise ValueError("Invalid SED template name")

    # Return
    return seds

# -----------------------------------------------------------------
