#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.wavelengthgrid Contains the WavelengthGridBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.basics.configurable import Configurable
from ...core.tools.utils import lazyproperty
from ...core.basics.emissionlines import EmissionLines, EmissionLine
from ...core.prep.wavelengthgrids import create_one_subgrid_wavelength_grid, create_template_seds, get_min_wavelength, get_max_wavelength
from ...core.filter.broad import BroadBandFilter
from ...core.filter.narrow import NarrowBandFilter
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.tools.serialization import write_dict, write_list, load_dict, load_list
from ...core.tools.stringify import tostr, stringify_list_fancy, stringify_list
from ...core.tools import filesystem as fs
from ...core.plot.wavelengthgrid import WavelengthGridPlotter

# -----------------------------------------------------------------

class WavelengthGridBuilder(Configurable):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(WavelengthGridBuilder, self).__init__(*args, **kwargs)

        # The grid
        self.wavelength_grid = None

        # Grid elements
        self.subgrid_wavelengths = None
        self.filter_wavelengths = None
        self.replaced = None
        self.new = None
        self.line_wavelengths = None
        self.fixed = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Build
        self.build()

        # Show
        self.show()

        # Writing
        if self.config.write: self.write()

        # Plotting
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(WavelengthGridBuilder, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def min_wavelength(self):

        """
        This function ...
        :return:
        """

        # Check filters?
        if self.config.check_filters is not None:

            # Get specified minimum wavelength
            min_wavelength = self.config.min_wavelength

            # Loop over the filters
            for fltr in self.config.check_filters:

                # Check below
                if fltr.wavelength < min_wavelength:

                    # Warning
                    log.warning("The wavelength range does not contain the wavelength of the '" + str(fltr) + "' filter")

                    # Adjust?
                    if self.config.adjust_minmax:
                        log.debug("Adjusting the minimum wavelength to incorporate the '" + str(fltr) + "' filter")
                        min_wavelength = 0.99 * fltr.wavelength

            # Return the lower wavelength
            return min_wavelength

        # Return the specified minimum wavelength
        else: return self.config.min_wavelength

    # -----------------------------------------------------------------

    @lazyproperty
    def max_wavelength(self):

        """
        This function ...
        :return:
        """

        # Check filters?
        if self.config.check_filters is not None:

            # Get specified maximum wavelength
            max_wavelength = self.config.max_wavelength

            # Loop over the filters
            for fltr in self.config.check_filters:

                # Check above
                if fltr.wavelength > max_wavelength:

                    # Warning
                    log.warning("The wavelength range does not contain the wavelength of the '" + str(fltr) + "' filter")

                    # Adjust?
                    if self.config.adjust_minmax:
                        log.debug("Adjusting the maximum wavelength to incorporate the '" + str(fltr) + "' filter")
                        max_wavelength = 1.01 * fltr.wavelength

            # Return the higher wavelength
            return max_wavelength

        # Return the specified maximum wavelength
        else: return self.config.max_wavelength

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

    def build(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")

        # Create the grid
        self.wavelength_grid, self.subgrid_wavelengths, self.filter_wavelengths, self.replaced, self.new, self.line_wavelengths, self.fixed = create_one_subgrid_wavelength_grid(self.config.npoints, self.emission_lines, self.config.fixed,
                                                                                                                                                                                   min_wavelength=self.min_wavelength, max_wavelength=self.max_wavelength,
                                                                                                                                                                                   filters=self.config.filters,
                                                                                                                                                                                   adjust_to=self.config.adjust_to,
                                                                                                                                                                                 min_wavelengths_in_filter=self.config.min_wavelengths_in_filter,
                                                                                                                                                                                 min_wavelengths_in_fwhm=self.config.min_wavelengths_in_fwhm,
                                                                                                                                                                                 return_elements=True)

    # -----------------------------------------------------------------

    @property
    def nwavelengths(self):
        return len(self.wavelength_grid)

    # -----------------------------------------------------------------

    @property
    def nfixed(self):
        return len(self.fixed)

    # -----------------------------------------------------------------

    @property
    def has_fixed(self):
        return self.nfixed > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def fixed_grid(self):
        return WavelengthGrid.from_wavelengths(self.fixed)

    # -----------------------------------------------------------------

    @property
    def nfilters(self):
        return len(self.filter_wavelengths)

    # -----------------------------------------------------------------

    @property
    def has_filters(self):
        return self.nfilters > 0

    # -----------------------------------------------------------------

    @property
    def nfilter_wavelengths(self):
        nwavelengths = 0
        for fltr in self.filter_wavelengths: nwavelengths += len(self.filter_wavelengths[fltr])
        return nwavelengths

    # -----------------------------------------------------------------

    @property
    def nreplaced(self):
        return len(self.replaced)

    # -----------------------------------------------------------------

    @property
    def has_replaced(self):
        return self.nreplaced > 0

    # -----------------------------------------------------------------

    @property
    def nnew(self):
        return len(self.new)

    # -----------------------------------------------------------------

    @property
    def has_new(self):
        return self.nnew > 0

    # -----------------------------------------------------------------

    @property
    def nlines(self):
        return len(self.line_wavelengths)

    # -----------------------------------------------------------------

    @property
    def has_line_wavelengths(self):
        return self.nlines > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def subgrids(self):

        """
        This function ...
        :return:
        """

        grids = OrderedDict()
        for name in self.subgrid_wavelengths:
            grid = WavelengthGrid.from_wavelengths(self.subgrid_wavelengths[name], sort=True)
            grids[name] = grid
        return grids

    # -----------------------------------------------------------------

    @property
    def subgrid_npoints(self):

        """
        This function ...
        :return:
        """

        # Keep track of the number of points per subgrid
        subgrid_npoints = OrderedDict()
        for subgrid in self.subgrid_wavelengths:
            subgrid_npoints[subgrid] = len(self.subgrid_wavelengths[subgrid])
        return subgrid_npoints

    # -----------------------------------------------------------------

    @property
    def broad_resampled(self):

        """
        This function ...
        :return:
        """

        # Set broad resampled nwavelengths
        broad_resampled = []
        for fltr in self.filter_wavelengths:
            if not isinstance(fltr, BroadBandFilter): continue
            broad_resampled.append(fltr)
        return broad_resampled

    # -----------------------------------------------------------------

    @property
    def narrow_added(self):

        """
        This function ...
        :return:
        """

        # Set narrow band filter added nwavelengths
        narrow_added = []
        for fltr in self.filter_wavelengths:
            if not isinstance(fltr, NarrowBandFilter): continue
            narrow_added.append(fltr)
        return narrow_added

    # -----------------------------------------------------------------

    @property
    def emission_npoints(self):

        """
        This function ...
        :return:
        """

        # Set emission npoints
        emission_npoints = 0
        for line_identifier in self.line_wavelengths:
            emission_npoints += len(self.line_wavelengths[line_identifier])
        return emission_npoints

    # -----------------------------------------------------------------

    @property
    def fixed_npoints(self):
        return len(self.fixed)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Created a wavelength grid with:")
        log.debug("")
        log.debug(" - number of points: " + str(self.nwavelengths))
        log.debug(" - number of points in subgrids:")
        for subgrid in self.subgrid_npoints: log.debug("    * " + subgrid + ": " + str(self.subgrid_npoints[subgrid]))
        log.debug(" - number of emission points: " + str(self.emission_npoints))
        log.debug(" - number of fixed points: " + str(self.fixed_npoints))
        log.debug(" - filters for which extra sampling was performed: " + str(self.broad_resampled))
        log.debug(" - narrow band filters for which wavelength was added: " + str(self.narrow_added))
        if self.has_replaced: log.debug(" - replaced wavelengths:")
        for old, new in self.replaced: log.debug("    * " + str(old) + " -> " + str(new))
        if self.has_new: log.debug(" - new wavelengths:")
        for line in stringify_list_fancy(self.new)[1].split("\n"): log.debug("    " + line)
        log.debug("")

        if log.is_debug:
            print("")
            print(self.wavelength_grid)
            print("")

    # -----------------------------------------------------------------

    def write(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the grid
        if self.config.write_grid: self.write_grid()

        # Write the subgrids
        self.write_subgrids()

        # Write the grid of fixed wavelengths
        if self.has_fixed: self.write_fixed_grid()

        # Write filter wavelengths
        if self.has_filters: self.write_filter_wavelengths()

        # Write replaced wavelengths
        if self.has_replaced: self.write_replaced()

        # Write new wavelengths
        if self.has_new: self.write_new()

        # Write emission line wavelengths
        if self.has_line_wavelengths: self.write_line_wavelengths()

    # -----------------------------------------------------------------

    @property
    def grid_path(self):
        return self.output_path_file("grid.dat")

    # -----------------------------------------------------------------

    def write_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the grid ...")

        # Write the wavelength grid
        self.wavelength_grid.saveto(self.grid_path)

    # -----------------------------------------------------------------

    def subgrid_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.output_path_file("subgrid_" + name + ".dat")

    # -----------------------------------------------------------------

    def write_subgrids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the subgrids ...")

        # Loop over the subgrids
        for name in self.subgrids:

            # Get path
            path = self.subgrid_path(name)

            # Write subgrid
            self.subgrids[name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def fixed_grid_path(self):
        return self.output_path_file("fixed_grid.dat")

    # -----------------------------------------------------------------

    def write_fixed_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the grid of fixed wavelength points ...")

        # Write
        self.fixed_grid.saveto(self.fixed_grid_path)

    # -----------------------------------------------------------------

    @property
    def filter_wavelengths_path(self):
        return self.output_path_file("filter_wavelengths.dat")

    # -----------------------------------------------------------------

    def write_filter_wavelengths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the filter wavelengths ...")

        # Write
        write_dict(self.filter_wavelengths, self.filter_wavelengths_path)

    # -----------------------------------------------------------------

    @property
    def replaced_path(self):
        return self.output_path_file("replaced.dat")

    # -----------------------------------------------------------------

    def write_replaced(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.info("Writing the replaced wavelengths ...")

        # Write
        write_list(self.replaced, self.replaced_path)

    # -----------------------------------------------------------------

    @property
    def new_path(self):
        return self.output_path_file("new.dat")

    # -----------------------------------------------------------------

    def write_new(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing the new wavelengths ...")

        # Write
        write_list(self.new, self.new_path)

    # -----------------------------------------------------------------

    @property
    def line_wavelengths_path(self):
        return self.output_path_file("line_wavelengths.dat")

    # -----------------------------------------------------------------

    def write_line_wavelengths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the emission line wavelengths ...")

        # Write
        write_dict(self.line_wavelengths, self.line_wavelengths_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------

def load_wavelength_grid_data(data_path):

    """
    This function ...
    :param data_path:
    :return:
    """

    # Load grid
    wavelength_grid_path = fs.join(data_path, "grid.dat")
    wavelength_grid = WavelengthGrid.from_file(wavelength_grid_path)

    # Load subgrids
    subgrids = OrderedDict()
    for filepath, filename in fs.files_in_path(data_path, startswith="subgrid_", returns=["path", "name"]):
        name = filename.split("subgrid_")[1].split(".dat")[0]
        grid = WavelengthGrid.from_file(filepath)
        subgrids[name] = grid

    # Load fixed grid
    fixed_grid_path = fs.join(data_path, "fixed_grid.dat")
    if fs.is_file(fixed_grid_path): fixed_grid = WavelengthGrid.from_file(fixed_grid_path)
    else: fixed_grid = None

    # Load filter wavelengths
    filter_wavelengths_path = fs.join(data_path, "filter_wavelengths.dat")
    if fs.is_file(filter_wavelengths_path): filter_wavelengths = load_dict(filter_wavelengths_path)
    else: filter_wavelengths = None

    # Replaced
    replaced_wavelengths_path = fs.join(data_path, "replaced.dat")
    if fs.is_file(replaced_wavelengths_path): replaced_wavelengths = load_list(replaced_wavelengths_path)
    else: replaced_wavelengths = None

    # Load new
    new_wavelengths_path = fs.join(data_path, "new.dat")
    if fs.is_file(new_wavelengths_path): new_wavelengths = load_list(new_wavelengths_path)
    else: new_wavelengths = None

    # Load lines
    line_wavelengths_path = fs.join(data_path, "line_wavelengths.dat")
    if fs.is_file(line_wavelengths_path): line_wavelengths = load_dict(line_wavelengths_path)
    else: line_wavelengths = None

    # Return the data
    return wavelength_grid, subgrids, fixed_grid, filter_wavelengths, replaced_wavelengths, new_wavelengths, line_wavelengths

# -----------------------------------------------------------------

def show_wavelength_grid_data(data_path):

    """
    This function ...
    :param data_path:
    :return:
    """

    # Get data
    wavelength_grid, subgrids, fixed_grid, filter_wavelengths, replaced_wavelengths, new_wavelengths, line_wavelengths = load_wavelength_grid_data(data_path)

    # Get properties
    nwavelengths = len(wavelength_grid)
    if fixed_grid is not None:
        fixed_npoints = len(fixed_grid)
        has_fixed = fixed_npoints > 0
    else: has_fixed, fixed_npoints = False, 0

    # Get emission npoints
    if line_wavelengths is not None:
        emission_npoints = 0
        for line_identifier in line_wavelengths:
            emission_npoints += len(line_wavelengths[line_identifier])
        has_lines = emission_npoints > 0
    else: has_lines, emission_npoints = False, 0

    # Subgrid
    # Keep track of the number of points per subgrid
    subgrid_npoints = OrderedDict()
    for subgrid in subgrids: subgrid_npoints[subgrid] = len(subgrids[subgrid])

    # Get broad resampled nwavelengths
    if filter_wavelengths is not None:
        broad_resampled = []
        for fltr in filter_wavelengths:
            if not isinstance(fltr, BroadBandFilter): continue
            broad_resampled.append(fltr)
        has_broad = len(broad_resampled) > 0
    else: has_broad, broad_resampled = False, None

    # Set narrow band filter added nwavelengths
    if filter_wavelengths is not None:
        narrow_added = []
        for fltr in filter_wavelengths:
            if not isinstance(fltr, NarrowBandFilter): continue
            narrow_added.append(fltr)
        has_narrow = len(narrow_added) > 0
    else: has_narrow, narrow_added = False, None

    # Replaced?
    if replaced_wavelengths is not None:
        nreplaced = len(replaced_wavelengths)
        has_replaced = nreplaced > 0
    else: has_replaced = False

    # New?
    if new_wavelengths is not None:
        nnew = len(new_wavelengths)
        has_new = nnew > 0
    else: has_new = False

    # Show
    print("")
    print(" - number of points: " + str(nwavelengths))
    print(" - number of points in subgrids:")
    for subgrid in subgrid_npoints: print("    * " + subgrid + ": " + str(subgrid_npoints[subgrid]))

    # Lines
    if has_lines:
        print(" - number of emission lines: " + str(len(line_wavelengths)))
        print(" - number of emission points: " + str(emission_npoints))

    # Fixed
    if has_fixed: print(" - number of fixed points: " + str(fixed_npoints))

    # Broad band filters
    if has_broad:
        print(" - broad band filters for which extra sampling was performed (" + str(len(broad_resampled)) + "): " + stringify_list(broad_resampled)[1])
        broad_nwavelengths = 0
        for fltr in broad_resampled:
            wavelengths = filter_wavelengths[fltr]
            print("   * " + str(fltr) + "(" + str(len(wavelengths)) + "): " + stringify_list(wavelengths)[1])
            broad_nwavelengths += len(wavelengths)
        print("  -> total wavelengths for broad band filters: " + str(broad_nwavelengths))

    # Narrow band filters
    if has_narrow:
        print(" - narrow band filters for which wavelength was added (" + str(len(narrow_added)) + "): " + stringify_list(narrow_added)[1])
        narrow_nwavelengths = 0
        for fltr in narrow_added:
            wavelengths = filter_wavelengths[fltr]
            print("   * " + str(fltr) + ": " + stringify_list(wavelengths)[1])
            narrow_nwavelengths += len(wavelengths)
        print("  -> total wavelengths for narrow band filters: " + str(narrow_nwavelengths))

    # Replaced
    if has_replaced:
        print(" - replaced wavelengths:")
        for old, new in replaced_wavelengths: print("    * " + str(old) + " -> " + str(new))

    # New
    if has_new:
        print(" - new wavelengths:")
        for line in stringify_list_fancy(new_wavelengths)[1].split("\n"): print("    " + line)

    # Show the grid
    print("")
    print(" - all wavelengths:")
    print("")
    print(wavelength_grid)
    print("")

# -----------------------------------------------------------------

def plot_wavelength_grid(data_path, wavelength_range=None, filepath=None, add_seds=None, ages=None,
                         metallicity=None, compactness=None, pressure=None, covering_factor=None, plot_reference=False):

    """
    This function ...
    :param data_path:
    :param wavelength_range:
    :param filepath:
    :param add_seds:
    :param ages:
    :param metallicity:
    :param compactness:
    :param pressure:
    :param covering_factor:
    :return:
    """

    # Get the elements of this grid
    wavelength_grid, subgrids, fixed_grid, filter_wavelengths, replaced_wavelengths, new_wavelengths, line_wavelengths = load_wavelength_grid_data(data_path)

    # Get range
    if wavelength_range is None: wavelength_range = wavelength_grid.range

    # Set limits
    check_filters = filter_wavelengths.keys()
    plot_filters = check_filters
    min_wavelength = get_min_wavelength(wavelength_range.min, check_filters, False)
    max_wavelength = get_max_wavelength(wavelength_range.max, check_filters, False)

    # Create seds?
    if add_seds is not None:
        seds = create_template_seds(add_seds, ages=ages, metallicity=metallicity,
                                    compactness=compactness, pressure=pressure,
                                    covering_factor=covering_factor)
    else: seds = {}

    # Create the plotter
    plotter = WavelengthGridPlotter()

    # Set settings
    plotter.config.add_regimes = True
    plotter.config.regimes = ["euv", "fuv", "muv", "nuv", "optical", "nir", "mir", "fir", "submm", "microwave"]
    plotter.config.only_subregimes = True

    # Filters
    plotter.config.add_filters = True
    #if self.config.plotting_filters is not None: plotter.config.filters = self.config.plotting_filters
    #else: plotter.config.filters = self.config.filters
    plotter.config.filters = plot_filters
    plotter.config.categorize_filters = True

    # Lines
    if line_wavelengths is not None:
        plotter.config.add_lines = True
        lines = line_wavelengths.keys()
        plotter.config.lines = lines

    # Plot separate
    plotter.config.separate_grids = True
    plotter.config.group_wavelengths = True
    # plotter.config.lines_in_group = "lines"
    plotter.config.separate_lines = True
    plotter.config.mark_removed = True

    # Plot resampled and residuals
    plotter.config.plot_resampled = True
    plotter.config.plot_interpolated = False
    plotter.config.plot_residuals = True

    # Add the elements
    plotter.add_elements(subgrids, fixed=fixed_grid, filter_wavelengths=filter_wavelengths, replaced=replaced_wavelengths, new=new_wavelengths, line_wavelengths=line_wavelengths)

    # Add complete grid
    if plot_reference:
        plotter.add_reference_grid(wavelength_grid, label="reference", in_legend=True)
        plotter.config.plot_differences = True

    # Run the plotter
    plotter.run(output=filepath, seds=seds, min_wavelength=min_wavelength, max_wavelength=max_wavelength)

# -----------------------------------------------------------------
