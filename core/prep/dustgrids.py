#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.dustgrids Contains the DustGridGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import the relevant PTS classes and modules
from ..simulation.grids import BinaryTreeDustGrid, OctTreeDustGrid, CartesianDustGrid
from ...core.tools.logging import log
from ...core.basics.range import zip_linear
from ..basics.configurable import Configurable
from ..basics.table import SmartTable
from ..basics.range import RealRange, QuantityRange, IntegerRange

# -----------------------------------------------------------------

class DustGridsTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(DustGridsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("Type", str, None, "grid type")
        self.add_column_info("Min x", float, "pc", "minimum x")
        self.add_column_info("Max x", float, "pc", "maximum x")
        self.add_column_info("Min y", float, "pc", "minimum y")
        self.add_column_info("Max y", float, "pc", "maximum y")
        self.add_column_info("Min z", float, "pc", "minimum z")
        self.add_column_info("Max z", float, "pc", "maximum z")
        self.add_column_info("Smallest scale", float, "pc", "Smallest scale")
        self.add_column_info("Min level", int, None, "Minimum level")
        self.add_column_info("Max mass fraction", float, None, "Maximum mass fraction")

    # -----------------------------------------------------------------

    def add_entry(self, grid_type, x_range, y_range, z_range, scale, min_level, max_mass_fraction):

        """
        This function ...
        :param grid_type:
        :param x_range:
        :param y_range:
        :param z_range:
        :param scale:
        :param min_level:
        :param max_mass_fraction:
        :return:
        """

        # Add a row to the table
        self.add_row([grid_type, x_range.min, x_range.max, y_range.min, y_range.max, z_range.min, z_range.max, scale, min_level, max_mass_fraction])

# -----------------------------------------------------------------

class DustGridGenerator(Configurable):
    
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
        super(DustGridGenerator, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Settings
        self.scale_range = None
        self.level_range = None
        self.mass_fraction_range = None
        self.ngrids = None

        self._grid_type = None
        self.x_radius = None
        self.y_radius = None
        self.z_radius = None

        # The dust grids
        self.grids = []

        # The dust grid property table
        self.table = None

    # -----------------------------------------------------------------

    @property
    def grid_type(self):

        """
        This function ...
        :return:
        """

        return self._grid_type

    # -----------------------------------------------------------------

    @grid_type.setter
    def grid_type(self, grid_type):

        """
        This function ...
        :return:
        """

        if not grid_type in ["cartesian", "bintree", "octtree"]: raise RuntimeError("Grid type '" + str(grid_type) + "' invalid. Must be either 'cartesian', 'bintree', or 'octtree'.")
        self._grid_type = grid_type

    # -----------------------------------------------------------------

    @property
    def x_min(self):

        """
        This function ...
        :return:
        """

        return - self.x_radius

    # -----------------------------------------------------------------

    @property
    def x_max(self):

        """
        This function ...
        :return:
        """

        return self.x_radius

    # -----------------------------------------------------------------

    @property
    def x_extent(self):

        """
        This function ...
        :return:
        """

        return self.x_max - self.x_min

    # -----------------------------------------------------------------

    @property
    def y_min(self):

        """
        This function ...
        :return:
        """

        return - self.y_radius

    # -----------------------------------------------------------------

    @property
    def y_max(self):

        """
        This function ...
        :return:
        """

        return self.y_radius

    # -----------------------------------------------------------------

    @property
    def y_extent(self):

        """
        This function ...
        :return:
        """

        return self.y_max - self.y_min

    # -----------------------------------------------------------------

    @property
    def z_min(self):

        """
        This function ...
        :return:
        """

        return - self.z_radius

    # -----------------------------------------------------------------

    @property
    def z_max(self):

        """
        This function ...
        :return:
        """

        return self.z_radius

    # -----------------------------------------------------------------

    @property
    def z_extent(self):

        """
        This function ...
        :return:
        """

        return self.z_max - self.z_min

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Generate the dust grids
        self.generate()

        # 3. Show
        if self.config.show: self.show()

        # 4. Write
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    @property
    def single_grid(self):

        """
        This function ...
        :return:
        """

        if len(self.grids) == 0: raise RuntimeError("No grid")
        elif len(self.grids) == 1: return self.grids[0]
        else: raise RuntimeError("More than one grid")

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DustGridGenerator, self).setup(**kwargs)

        # Get settings
        self.ngrids = kwargs.pop("ngrids")
        if self.ngrids == 1:
            self.scale_range = QuantityRange.infinitesimal(kwargs.pop("scale"))
            self.level_range = IntegerRange.infinitesimal(kwargs.pop("level"))
            self.mass_fraction_range = RealRange.infinitesimal(kwargs.pop("mass_fraction"))
        else:
            self.scale_range = kwargs.pop("scale_range")
            self.level_range = kwargs.pop("level_range")
            self.mass_fraction_range = kwargs.pop("mass_fraction_range")

        # Initialize the table
        self.table = DustGridsTable()

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the grids ...")

        # Loop over the different grid parameter values
        for scale, level, mass_fraction in zip_linear(self.scale_range, self.level_range, self.mass_fraction_range, npoints=self.ngrids):

            # Create the grid and add it to the list
            if self.grid_type == "cartesian": self.create_cartesian_dust_grid(scale)
            elif self.grid_type == "bintree": self.create_binary_tree_dust_grid(scale, level, mass_fraction)
            elif self.grid_type == "octtree": self.create_octtree_dust_grid(scale, level, mass_fraction)
            else: raise ValueError("Invalid grid type: " + self.grid_type)

    # -----------------------------------------------------------------

    def create_cartesian_dust_grid(self, scale):

        """
        This function ...
        :param scale:
        :return:
        """

        # Create the grid
        grid = create_one_cartesian_dust_grid(scale, self.x_extent, self.y_extent, self.z_extent, self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max)

        # Add the grid
        self.grids.append(grid)

        # Debugging
        log.debug("Created a cartesian dust grid with:")
        log.debug("")
        log.debug(" - x_min: " + str(self.x_min))
        log.debug(" - x_max: " + str(self.x_max))
        log.debug(" - y_min: " + str(self.y_min))
        log.debug(" - y_max: " + str(self.y_max))
        log.debug(" - z_min: " + str(self.z_min))
        log.debug(" - z_max: " + str(self.z_max))
        log.debug(" - scale: " + str(scale))
        log.debug(" - x_bins: " + str(grid.x_bins))
        log.debug(" - y_bins: " + str(grid.y_bins))
        log.debug(" - z_bins: " + str(grid.z_bins))
        log.debug(" - mesh type: " + str(grid.mesh_type))
        log.debug(" - ratio: " + str(grid.real))
        log.debug("")

        # Add a row to the table
        self.table.add_row([self.grid_type, self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max, scale, None, None])

    # -----------------------------------------------------------------

    def create_binary_tree_dust_grid(self, scale, min_level, max_mass_fraction):

        """
        This function ...
        :param scale:
        :param min_level:
        :param max_mass_fraction:
        :return:
        """

        # Create the grid
        grid = create_one_bintree_dust_grid(scale, self.x_extent, self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max, min_level, max_mass_fraction)

        # Add the grid
        self.grids.append(grid)

        # Debugging
        log.debug("Created a binary tree dust grid with:")
        log.debug("")
        log.debug(" - x_min: " + str(self.x_min))
        log.debug(" - x_max: " + str(self.x_max))
        log.debug(" - y_min: " + str(self.y_min))
        log.debug(" - y_max: " + str(self.y_max))
        log.debug(" - z_min: " + str(self.z_min))
        log.debug(" - z_max: " + str(self.z_max))
        log.debug(" - scale: " + str(scale))
        log.debug(" - min_level: " + str(grid.min_level))
        log.debug(" - max_level: " + str(grid.max_level))
        log.debug(" - search_method: " + str(grid.search_method))
        log.debug(" - sample_count: " + str(grid.sample_count))
        log.debug(" - max_optical_depth: " + str(grid.max_optical_depth))
        log.debug(" - max_mass_fraction: " + str(grid.max_mass_fraction))
        log.debug(" - max_dens_disp_fraction: " + str(grid.max_dens_disp_fraction))
        log.debug(" - direction_method: " + str(grid.direction_method))
        log.debug("")

        # Add a row to the table
        self.table.add_row([self.grid_type, self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max, scale, min_level, max_mass_fraction])

    # -----------------------------------------------------------------

    def create_octtree_dust_grid(self, scale, min_level, max_mass_fraction):

        """
        This function ...
        :param scale:
        :param min_level:
        :param max_mass_fraction:
        :return:
        """

        # Create the grid
        grid = create_one_octtree_dust_grid(scale, self.x_extent, self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max, min_level, max_mass_fraction)

        # Add the grid
        self.grids.append(grid)

        # Debugging
        log.debug("Created an octtree dust grid with:")
        log.debug("")
        log.debug(" - x_min: " + str(self.x_min))
        log.debug(" - x_max: " + str(self.x_max))
        log.debug(" - y_min: " + str(self.y_min))
        log.debug(" - y_max: " + str(self.y_max))
        log.debug(" - z_min: " + str(self.z_min))
        log.debug(" - z_max: " + str(self.z_max))
        log.debug(" - scale: " + str(scale))
        log.debug(" - min_level: " + str(grid.min_level))
        log.debug(" - max_level: " + str(grid.max_level))
        log.debug(" - search_method: " + str(grid.search_method))
        log.debug(" - sample_count: " + str(grid.sample_count))
        log.debug(" - max_optical_depth: " + str(grid.max_optical_depth))
        log.debug(" - max_mass_fraction: " + str(grid.max_mass_fraction))
        log.debug(" - max_dens_disp_fraction: " + str(grid.max_dens_disp_fraction))
        log.debug(" - barycentric: " + str(grid.barycentric))
        log.debug("")

        # Add entry to the table
        x_range = RealRange(self.x_min, self.x_max)
        y_range = RealRange(self.y_min, self.y_max)
        z_range = RealRange(self.z_min, self.z_max)
        self.table.add_entry(self, self.grid_type, x_range, y_range, z_range, scale, min_level, max_mass_fraction)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ..
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the grids
        self.write_grids()

        # Write table
        self.write_table()

    # -----------------------------------------------------------------

    def write_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing grids ...")

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing table ...")

# -----------------------------------------------------------------

def create_one_dust_grid(grid_type, scale, x_extent, y_extent, z_extent, x_min, x_max, y_min, y_max, z_min, z_max, min_level, max_mass_fraction):

    """
    This function ...
    :param grid_type:
    :param scale:
    :param x_extent:
    :param y_extent:
    :param z_extent:
    :param x_min:
    :param x_max:
    :param y_min:
    :param y_max:
    :param z_min:
    :param z_max:
    :param min_level:
    :param max_mass_fraction:
    :return:
    """

    # Create the specified type of grid
    if grid_type == "cartesian": return create_one_cartesian_dust_grid(scale, x_extent, y_extent, z_extent, x_min, x_max, y_min, y_max, z_min, z_max)
    elif grid_type == "bintree": return create_one_bintree_dust_grid(scale, x_extent, x_min, x_max, y_min, y_max, z_min, z_max, min_level, max_mass_fraction)
    elif grid_type == "octtree": return create_one_octtree_dust_grid(scale, x_extent, x_min, x_max, y_min, y_max, z_min, z_max, min_level, max_mass_fraction)
    else: raise ValueError("Unknown dust grid type: " + grid_type)

# -----------------------------------------------------------------

def create_one_cartesian_dust_grid(scale, x_extent, y_extent, z_extent, x_min, x_max, y_min, y_max, z_min, z_max):

    """
    This function ...
    :param scale:
    :param x_extent:
    :param y_extent:
    :param z_extent:
    :param x_min:
    :param x_max:
    :param y_min:
    :param y_max:
    :param z_min:
    :param z_max:
    :return:
    """

    # Inform the user
    log.info("Creating a cartesian dust grid with a smallest physical scale of " + str(scale) + " ...")

    # Calculate the number of bins in each direction
    x_bins = int(math.ceil(x_extent.to("pc").value / scale.to("pc").value))
    y_bins = int(math.ceil(y_extent.to("pc").value / scale.to("pc").value))
    z_bins = int(math.ceil(z_extent.to("pc").value / scale.to("pc").value))

    # Create the grid
    grid = CartesianDustGrid(x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, z_min=z_min, z_max=z_max,
                             x_bins=x_bins, y_bins=y_bins, z_bins=z_bins)

    # Return the grid
    return grid

# -----------------------------------------------------------------

def create_one_bintree_dust_grid(scale, x_extent, x_min, x_max, y_min, y_max, z_min, z_max, min_level, max_mass_fraction):

    """
    This function ...
    :param scale:
    :param x_extent:
    :param x_min:
    :param x_max:
    :param y_min:
    :param y_max:
    :param z_min:
    :param z_max:
    :param min_level:
    :param max_mass_fraction:
    :return:
    """

    # Inform the user
    log.info("Creating a binary tree dust grid with a smallest physical scale of " + str(scale) + ", with a minimum division level of " + str(min_level) + " and a maximum mass fraction of " + str(max_mass_fraction) + " ...")

    # Calculate the maximum division level that is necessary to resolve the smallest scale of the input maps
    extent_x = x_extent.to("pc").value
    smallest_scale = scale.to("pc").value
    max_level = min_level_for_smallest_scale_bintree(extent_x, smallest_scale)

    # Create the dust grid
    grid = BinaryTreeDustGrid(min_x=x_min, max_x=x_max, min_y=y_min, max_y=y_max, min_z=z_min, max_z=z_max,
                              min_level=min_level, max_level=max_level, max_mass_fraction=max_mass_fraction)

    # Return the grid
    return grid

# -----------------------------------------------------------------

def create_one_octtree_dust_grid(scale, x_extent, x_min, x_max, y_min, y_max, z_min, z_max, min_level, max_mass_fraction):

    """
    This function ...
    :param scale:
    :param x_extent:
    :param x_min:
    :param x_max:
    :param y_min:
    :param y_max:
    :param z_min:
    :param z_max:
    :param min_level:
    :param max_mass_fraction:
    :return:
    """

    # Inform the user
    log.info("Creating a octtree dust grid with a smallest physical scale of " + str(scale) + ", with a minimum division level of " + str(min_level) + " and a maximum mass fraction of " + str(max_mass_fraction) + " ...")

    # Calculate the minimum division level that is necessary to resolve the smallest scale of the input maps
    extent_x = x_extent.to("pc").value
    smallest_scale = scale.to("pc").value
    max_level = min_level_for_smallest_scale_octtree(extent_x, smallest_scale)

    # Create the dust grid
    grid = OctTreeDustGrid(min_x=x_min, max_x=x_max, min_y=y_min, max_y=y_max, min_z=z_min, max_z=z_max,
                           min_level=min_level, max_level=max_level, max_mass_fraction=max_mass_fraction)

    # Return the grid
    return grid

# -----------------------------------------------------------------

def min_level_for_smallest_scale_bintree(extent, smallest_scale):

    """
    This function ...
    :param extent:
    :param smallest_scale:
    :return:
    """

    ratio = extent / smallest_scale
    octtree_level = int(math.ceil(math.log(ratio, 2)))
    level = int(3 * octtree_level)
    return level

# -----------------------------------------------------------------

def min_level_for_smallest_scale_octtree(extent, smallest_scale):

    """
    This function ...
    :param extent:
    :param smallest_scale:
    :return:
    """

    ratio = extent / smallest_scale
    octtree_level = int(math.ceil(math.log(ratio, 2)))
    return octtree_level

# -----------------------------------------------------------------
