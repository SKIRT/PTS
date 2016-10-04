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

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS classes and modules
from ..basics.grids import BinaryTreeDustGrid, OctTreeDustGrid, CartesianDustGrid
from ...core.tools.logging import log
from ...core.basics.range import zip_linear, zip_log

# -----------------------------------------------------------------

class DustGridGenerator(object):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(DustGridGenerator, self).__init__()

        # -- Attributes --

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

    def run(self, scale_range, level_range, mass_fraction_range, ngrids):

        """
        This function ...
        :param scale_range:
        :param level_range:
        :param mass_fraction_range:
        :param ngrids:
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Generate the dust grids
        self.generate(scale_range, level_range, mass_fraction_range, ngrids)

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Initialize the table
        names = ["Type", "Min x", "Max x", "Min y", "Max y", "Min z", "Max z", "Smallest scale", "Min level", "Max mass fraction"]
        dtypes = ["S9", "f8", "f8", "f8", "f8", "f8", "f8", "f8", int, "f8"]
        self.table = Table(names=names, dtype=dtypes)

    # -----------------------------------------------------------------

    def generate(self, scale_range, level_range, mass_fraction_range, ngrids):

        """
        This function ...
        :param scale_range:
        :param level_range:
        :param mass_fraction_range:
        :param ngrids:
        :return:
        """

        # Inform the user
        log.info("Creating the grids ...")

        # Loop over the different grid parameter values
        for scale, level, mass_fraction in zip_linear(scale_range, level_range, mass_fraction_range, npoints=ngrids):

            # Create the grid and add it to the list
            if self.grid_type == "cartesian": self.create_cartesian_dust_grid(scale)
            elif self.grid_type == "bintree": self.create_binary_tree_dust_grid(scale, level, mass_fraction)
            elif self.grid_type == "octtree": self.create_octtree_dust_grid(scale, level, mass_fraction)

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

        # Add a row to the table
        self.table.add_row([self.grid_type, self.x_min, self.x_max, self.y_min, self.y_max, self.z_min, self.z_max, scale, min_level, max_mass_fraction])

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
    grid = CartesianDustGrid(x_min, x_max, y_min, y_max, z_min, z_max, x_bins, y_bins, z_bins)

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

    # Calculate the minimum division level that is necessary to resolve the smallest scale of the input maps
    extent_x = x_extent.to("pc").value
    smallest_scale = scale.to("pc").value
    max_level = min_level_for_smallest_scale_bintree(extent_x, smallest_scale)

    # Create the dust grid
    grid = BinaryTreeDustGrid(x_min, x_max, y_min, y_max, z_min, z_max, min_level=min_level, max_level=max_level, max_mass_fraction=max_mass_fraction)

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
    grid = OctTreeDustGrid(x_min, x_max, y_min, y_max, z_min, z_max, min_level=min_level, max_level=max_level, max_mass_fraction=max_mass_fraction)

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
