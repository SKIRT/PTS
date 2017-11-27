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
from ...core.basics.emissionlines import EmissionLines
from ...core.prep.wavelengthgrids import create_one_subgrid_wavelength_grid
from ...core.filter.broad import BroadBandFilter
from ...core.filter.narrow import NarrowBandFilter
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.tools.serialization import write_dict, write_list
from ...core.tools.stringify import tostr, stringify_list_fancy

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

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Build
        self.build()

        # 3. Show
        self.show()

        # 4. Writing
        if self.config.write: self.write()

        # 5. Plotting
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
        if self.config.add_emission_lines: return EmissionLines()
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
                                                                                                                                                                                   adjust_to=self.config.adjust_to, return_elements=True)

    # -----------------------------------------------------------------

    @property
    def nwavelengths(self):

        """
        Thisf unction ...
        :return:
        """

        return len(self.wavelength_grid)

    # -----------------------------------------------------------------

    @property
    def nfixed(self):

        """
        Thisf unction ...
        :return:
        """

        return len(self.fixed)

    # -----------------------------------------------------------------

    @property
    def has_fixed(self):

        """
        This function ...
        :return:
        """

        return self.nfixed > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def fixed_grid(self):

        """
        This function ...
        :return:
        """

        return WavelengthGrid.from_wavelengths(self.fixed)

    # -----------------------------------------------------------------

    @property
    def nfilters(self):

        """
        This function ...
        :return:
        """

        return len(self.filter_wavelengths)

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
    def nfilter_wavelengths(self):

        """
        This function ...
        :return:
        """

        nwavelengths = 0
        for fltr in self.filter_wavelengths: nwavelengths += len(self.filter_wavelengths[fltr])
        return nwavelengths

    # -----------------------------------------------------------------

    @property
    def nreplaced(self):

        """
        This function ...
        :return:
        """

        return len(self.replaced)

    # -----------------------------------------------------------------

    @property
    def has_replaced(self):

        """
        This function ...
        :return:
        """

        return self.nreplaced > 0

    # -----------------------------------------------------------------

    @property
    def nnew(self):

        """
        Thisj function ...
        :return:
        """

        return len(self.new)

    # -----------------------------------------------------------------

    @property
    def has_new(self):

        """
        This function ...
        :return:
        """

        return self.nnew > 0

    # -----------------------------------------------------------------

    @property
    def nlines(self):

        """
        This function ...
        :return:
        """

        return len(self.line_wavelengths)

    # -----------------------------------------------------------------

    @property
    def has_line_wavelengths(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

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
        log.debug(" - replaced wavelengths:")
        for old, new in self.replaced: log.debug("    * " + str(old) + " -> " + str(new))
        log.debug(" - new wavelengths:")
        for line in stringify_list_fancy(self.new)[1].split("\n"): log.debug("    " + line)
        log.debug("")

        if log.is_debug():
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

        """
        Thisfunction ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

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

        """
        Thisn function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

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
