#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.dustgridtool Contains the DustGridTool class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.simulation.execute import SkirtExec
from ...core.simulation.arguments import SkirtArguments
from ...core.basics.map import Map
from ..simulation.grids import load_grid
from ..basics.unit import parse_unit as u

# -----------------------------------------------------------------

class DustGridTool(object):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # -- Attributes --

        # Local SKIRT execution environment
        self.skirt = SkirtExec()

    # -----------------------------------------------------------------

    def optimize(self, ski):

        """
        This function ...
        :param ski:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def get_statistics(self, ski, simulation_path, input_path, prefix):

        """
        This function ...
        :return:
        """

        out_path = fs.create_directory_in(simulation_path, "out")

        # Make a copy of the ski file
        ski = ski.copy()

        # Set npackages to zero
        ski.to_oligochromatic(1. * u("micron"))
        ski.setpackages(0)

        # ski.remove_all_stellar_components()
        ski.remove_all_instruments()

        # Disable all writing options, except the one for writing the dust grid and cell properties
        ski.disable_all_writing_options()
        ski.set_write_grid()
        ski.set_write_cell_properties()

        # Save the ski file
        ski_path = fs.join(simulation_path, prefix + ".ski")
        ski.saveto(ski_path)

        # Create arguments
        arguments = SkirtArguments()

        arguments.input_path = input_path
        arguments.output_path = out_path

        arguments.single = True
        arguments.ski_pattern = ski_path
        arguments.logging.verbose = True

        # Run SKIRT
        simulation = self.skirt.run(arguments, progress_bar=True)

        # Get and parse the log file
        log_file = simulation.log_file

        # Check for errors
        if log_file.last_message.startswith("*** Error"):
            for line in log_file.messages: log.error(line)
            raise RuntimeError("The simulation has crashed")

        # Initilize statistics map
        statistics = Map()

        # Get the number of dust cells
        statistics.ncells = log_file.dust_cells_tree

        # Get the number of tree nodes
        statistics.tree_nodes = log_file.tree_nodes

        # Get the tree leaf distribution
        statistics.tree_leaf_distribution = log_file.tree_leaf_distribution

        # Get the number of tree levels
        statistics.tree_levels = log_file.tree_levels

        # Determine the path to the cell properties file
        cellprops_path = fs.join(out_path, prefix + "_ds_cellprops.dat")

        # Get the optical depth for which 90% of the cells have a smaller value
        optical_depth = None
        for line in reversed(open(cellprops_path).readlines()):
            if "of the cells have optical depth smaller than" in line:
                optical_depth = float(line.split("than: ")[1])
                break
        statistics.optical_depth_90 = optical_depth

        # Return the statistics
        return statistics

    # -----------------------------------------------------------------

    def adjust(self):

        """
        This function ...
        :return:
        """

        # Loop over the dust grid files
        for path, name in fs.files_in_path(self.dust_grids_path, extension="dg", sort=int, returns=["path", "name"]):

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
        ski.to_oligochromatic([1. * u("micron")])

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
