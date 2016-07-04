#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.dustgridadjuster Contains the DustGridAdjuster class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.simulation.execute import SkirtExec
from ...core.simulation.arguments import SkirtArguments
from ...core.simulation.skifile import SkiFile
from ..basics.grids import load_grid

# -----------------------------------------------------------------

class DustGridAdjuster(FittingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(DustGridAdjuster, self).__init__(config)

        # -- Attributes --

        # The ski file
        self.ski = None

        # The adjusted dust grids
        self.grids = []

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the ski file
        self.load_ski()

        # 3. Adjust
        self.adjust()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DustGridAdjuster, self).setup()

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Open the ski file
        self.ski = SkiFile(self.fit_ski_path)

    # -----------------------------------------------------------------

    def adjust(self):

        """
        This function ...
        :return:
        """

        # Loop over the dust grid files
        for path, name in fs.files_in_path(self.fit_dust_grids_path, extension="dg", sort=int, returns=["path", "name"]):

            # Create corresponding directory
            out_path = fs.create_directory_in(path, name)

            # Load the dust grid instance
            grid = load_grid(path)

            # Generate the grid data
            optical_depth = self.generate_grid(grid, out_path)

            # Debugging
            log.debug("90% of the cells have an optical depth smaller than " + str(optical_depth))

            # Adapt the maximal optical depth criterion
            grid.max_optical_depth = optical_depth

            # Inform the user
            log.info("Generating the high-resolution grid data ...")

            # Rerun the simulation
            optical_depth = self.generate_grid(grid, out_path)

            # Debugging
            log.debug("90% of the cells have an optical depth smaller than " + str(optical_depth))

            # Add the adjusted grid
            self.grids.append(grid)

    # -----------------------------------------------------------------

    def generate_grid(self, grid, output_path):

        """
        This function ...
        :param grid:
        :param output_path:
        :return:
        """

        # Inform the user
        log.info("Running a simulation just to generate the dust grid data files ...")

        # Create a copy of the ski file
        ski = self.ski.copy()

        # Set the dust grid
        ski.set_dust_grid(grid)

        # Convert to oligochromatic simulation
        ski.to_oligochromatic([1. * Unit("micron")])

        # Remove the instrument system
        ski.remove_instrument_system()

        # Set the number of photon packages to zero
        ski.setpackages(0)

        # Disable all writing options, except the one for writing the dust grid and cell properties
        ski.disable_all_writing_options()
        ski.set_write_grid()
        ski.set_write_cell_properties()

        # Write the ski file
        ski_path = fs.join(output_path, self.galaxy_name + ".ski")
        ski.saveto(ski_path)

        # Create the local SKIRT execution context
        skirt = SkirtExec()

        # Create the SKIRT arguments object
        arguments = SkirtArguments()
        arguments.ski_pattern = ski_path
        arguments.input_path = self.fit_in_path
        arguments.output_path = output_path

        # Run SKIRT to generate the dust grid data files
        skirt.run(arguments)

        # Determine the path to the cell properties file
        cellprops_path = fs.join(output_path, self.galaxy_name + "_ds_cellprops.dat")

        # Get the optical depth for which 90% of the cells have a smaller value
        optical_depth = None
        for line in reversed(open(cellprops_path).readlines()):
            if "of the cells have optical depth smaller than" in line:
                optical_depth = float(line.split("than: ")[1])
                break

        # Return the optical depth
        return optical_depth

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        self.write_grids()

    # -----------------------------------------------------------------

    def write_grids(self):

        """
        This function ...
        :return:
        """

        # Loop over the grids
        for index, grid in enumerate(self.grids):

            path = fs.join(self.fit_dust_grids_path, str(index) + ".dg")
            grid.save(path)

# -----------------------------------------------------------------
