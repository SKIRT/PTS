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
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.simulation.execute import SkirtExec
from ...core.simulation.arguments import SkirtArguments
from ...core.basics.map import Map
from ..simulation.grids import load_grid
from ..units.parsing import parse_unit as u
from ...core.simulation.definition import SingleSimulationDefinition

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
            prefix = generate_grid(grid, out_path)
            optical_depth = get_optical_depth_criterium(out_path, prefix)

            # Debugging
            log.debug("90% of the cells have an optical depth smaller than " + str(optical_depth))

            # Adapt the maximal optical depth criterion
            grid.max_optical_depth = optical_depth

            # Inform the user
            log.info("Generating the high-resolution grid data ...")

            # Rerun the simulation
            prefix = generate_grid(grid, out_path)
            optical_depth = get_optical_depth_criterium(out_path, prefix)

            # Debugging
            log.debug("90% of the cells have an optical depth smaller than " + str(optical_depth))

            # Add the adjusted grid
            self.grids.append(grid)

# -----------------------------------------------------------------

def generate_grid(ski, grid, output_path, input_path):

    """
    This function ...
    :param ski:
    :param grid:
    :param output_path:
    :param input_path:
    :return:
    """

    # Inform the user
    log.info("Running a simulation just to generate the dust grid data files ...")

    # Determine ski prefix
    prefix = "dg"

    # Create a copy of the ski file
    ski = ski.copy()

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

    # WRITE THE DUST GRID TREE
    if ski.has_tree_dust_grid: ski.set_write_grid_tree()

    # WRITE OTHER
    ski.set_write_quality()
    ski.set_write_convergence()
    ski.set_write_density()
    ski.set_write_depth_map()
    ski.set_write_cell_properties()

    # OR: ski.enable_all_dust_system_writing_options() ?

    # Write the ski file
    ski_path = fs.join(output_path, prefix + ".ski")
    ski.saveto(ski_path)

    # Create the local SKIRT execution context
    skirt = SkirtExec()

    # Create simulation definition
    definition = SingleSimulationDefinition(ski_path, output_path, input_path, name=None)

    # Run SKIRT to generate the dust grid data files
    skirt.run(definition)

    # Check
    if fs.is_empty(output_path, besides=ski_path): raise RuntimeError("Something went wrong: no output in " + output_path)

    # Return the prefix
    return prefix

# -----------------------------------------------------------------

def get_optical_depth_criterium(output_path, prefix):

    """
    This function ...
    :param output_path: 
    :param prefix: 
    :return: 
    """

    # Determine the path to the cell properties file
    cellprops_path = fs.join(output_path, prefix + "_ds_cellprops.dat")

    # Get the optical depth for which 90% of the cells have a smaller value
    optical_depth = None
    for line in reversed(open(cellprops_path).readlines()):
        if "of the cells have optical depth smaller than" in line:
            optical_depth = float(line.split("than: ")[1])
            break

    # Return the optical depth
    return optical_depth

# -----------------------------------------------------------------

def get_statistics(ski, simulation_path, input_path, prefix):

    """
    This function ...
    :param ski:
    :param simulation_path:
    :param input_path:
    :param prefix:
    :return:
    """

    # Local SKIRT execution environment
    skirt = SkirtExec()

    # Create output path
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
    simulation = skirt.run(arguments, show_progress=True)

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
