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
from astropy.units import Unit, dimensionless_angles

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
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

        assert self.grid_type in ["cartesian", "bintree", "octtree"]
        self._grid_type = grid_type

    # -----------------------------------------------------------------

    def run(self, scale_range, level_range, mass_fraction_range, ngrids, grid_type="bintree"):

        """
        This function ...
        :param scale_range:
        :param level_range:
        :param mass_fraction_range:
        :param ngrids:
        :param grid_type:
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Generate the dust grids
        self.generate(scale_range, level_range, mass_fraction_range, ngrids, grid_type)

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        pass

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

            # Add a row to the table


    # -----------------------------------------------------------------

    def create_cartesian_dust_grid(self, scale):

        """
        This function ...
        :param scale:
        :return:
        """

        # Inform the user
        log.info("Configuring the cartesian dust grid ...")

        # Calculate the boundaries of the dust grid
        min_x = - self.x_radius
        max_x = self.x_radius
        min_y = - self.y_radius
        max_y = self.y_radius
        min_z = - self.z_radius
        max_z = self.z_radius

        # Calculate the number of bins in each direction
        x_bins = int(math.ceil((max_x - min_x).to("pc").value / scale.to("pc").value))
        y_bins = int(math.ceil((max_y - min_y).to("pc").value / scale.to("pc").value))
        z_bins = int(math.ceil((max_z - min_z).to("pc").value / scale.to("pc").value))

        # Create the grid
        grid = CartesianDustGrid(min_x, max_x, min_y, max_y, min_z, max_z, x_bins, y_bins, z_bins)

        # Add the grid
        self.grids.append(grid)

    # -----------------------------------------------------------------

    def create_binary_tree_dust_grid(self, scale, min_level, max_mass_fraction):

        """
        This function ...
        :param scale:
        :param min_level:
        :param max_mass_fraction:
        :return:
        """

        # Inform the user
        log.info("Configuring the bintree dust grid ...")

        # Calculate the boundaries of the dust grid
        min_x = - self.x_radius
        max_x = self.x_radius
        min_y = - self.y_radius
        max_y = self.y_radius
        min_z = - self.z_radius
        max_z = self.z_radius

        # Calculate the minimum division level that is necessary to resolve the smallest scale of the input maps
        extent_x = (max_x - min_x).to("pc").value
        smallest_scale = scale.to("pc").value
        max_level = min_level_for_smallest_scale_bintree(extent_x, smallest_scale)

        # Create the dust grid
        grid = BinaryTreeDustGrid(min_x, max_x, min_y, max_y, min_z, max_z, min_level=min_level, max_level=max_level, max_mass_fraction=max_mass_fraction)

        # Add the grid
        self.grids.append(grid)

    # -----------------------------------------------------------------

    def create_octtree_dust_grid(self, smallest_cell_pixels, min_level, max_mass_fraction):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the octtree dust grid ...")

        # Calculate the boundaries of the dust grid
        min_x = - radius_physical
        max_x = radius_physical
        min_y = - radius_physical
        max_y = radius_physical
        min_z = -3. * Unit("kpc")
        max_z = 3. * Unit("kpc")

        # Calculate the minimum division level that is necessary to resolve the smallest scale of the input maps
        extent_x = (max_x - min_x).to("pc").value
        smallest_scale = smallest_scale.to("pc").value
        max_level = min_level_for_smallest_scale_octtree(extent_x, smallest_scale)

        # Create the dust grid
        grid = OctTreeDustGrid(min_x, max_x, min_y, max_y, min_z, max_z, min_level=min_level, max_level=max_level, max_mass_fraction=max_mass_fraction)

        # Add the grid
        self.grids.append(grid)

# -----------------------------------------------------------------
