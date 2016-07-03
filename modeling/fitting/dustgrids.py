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

# -----------------------------------------------------------------

class DustGridGenerator(object):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(DustGridGenerator, self).__init__()

        # -- Attributes --

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # ...

        # 17. Writing
        #self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def create_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the grids ...")

    # -----------------------------------------------------------------

    def create_low_res_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the low-resolution dust grid ...")

        # Set the dust grid
        if self.config.dust_grid == "cartesian": self.lowres_dust_grid = self.create_cartesian_dust_grid(10.)
        elif self.config.dust_grid == "bintree": self.lowres_dust_grid = self.create_binary_tree_dust_grid(10., 6, 1e-5)
        elif self.config.dust_grid == "octtree": self.lowres_dust_grid = self.create_octtree_dust_grid(10., 2, 1e-5)
        else: raise ValueError("Invalid option for dust grid type")

    # -----------------------------------------------------------------

    def create_high_res_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the high-resolution dust grid ...")

        # Set the dust grid
        if self.config.dust_grid == "cartesian": self.highres_dust_grid = self.create_cartesian_dust_grid(1.)
        elif self.config.dust_grid == "bintree": self.highres_dust_grid = self.create_binary_tree_dust_grid(0.5, 9, 0.5e-6)
        elif self.config.dust_grid == "octtree": self.highres_dust_grid = self.create_octtree_dust_grid(0.5, 3, 0.5e-6)
        else: raise ValueError("Invalid option for dust grid type")

    # -----------------------------------------------------------------

    def create_cartesian_dust_grid(self, smallest_cell_pixels):

        """
        This function ...
        :param smallest_cell_pixels:
        :return:
        """

        # Inform the user
        log.info("Configuring the cartesian dust grid ...")

        # Calculate the major radius of the truncation ellipse in physical coordinates (pc)
        major_angular = self.ellipse.major  # major axis length of the sky ellipse
        radius_physical = (major_angular * self.parameters.distance).to("pc", equivalencies=dimensionless_angles())

        # Calculate the boundaries of the dust grid
        min_x = - radius_physical
        max_x = radius_physical
        min_y = - radius_physical
        max_y = radius_physical
        min_z = -3. * Unit("kpc")
        max_z = 3. * Unit("kpc")

        # Get the pixelscale in physical units
        distance = self.parameters.distance
        pixelscale_angular = self.reference_wcs.xy_average_pixelscale * Unit("pix")  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        # Because we (currently) can't position the grid exactly as the 2D pixels,
        # take half of the pixel size to avoid too much interpolation
        # smallest_scale = 0.5 * pixelscale
        smallest_scale = smallest_cell_pixels * pixelscale  # limit the number of cells

        # Calculate the number of bins in each direction
        x_bins = int(math.ceil((max_x - min_x).to("pc").value / smallest_scale.to("pc").value))
        y_bins = int(math.ceil((max_y - min_y).to("pc").value / smallest_scale.to("pc").value))
        z_bins = int(math.ceil((max_z - min_z).to("pc").value / smallest_scale.to("pc").value))

        # Create and return the dust grid
        return CartesianDustGrid(min_x, max_x, min_y, max_y, min_z, max_z, x_bins, y_bins, z_bins)

    # -----------------------------------------------------------------

    def create_binary_tree_dust_grid(self, smallest_cell_pixels, min_level, max_mass_fraction):

        """
        This function ...
        :param smallest_cell_pixels:
        :param min_level:
        :param max_mass_fraction:
        :return:
        """

        # Inform the user
        log.info("Configuring the bintree dust grid ...")

        # Calculate the major radius of the truncation ellipse in physical coordinates (pc)
        major_angular = self.ellipse.major  # major axis length of the sky ellipse
        radius_physical = (major_angular * self.parameters.distance).to("pc", equivalencies=dimensionless_angles())

        # Calculate the boundaries of the dust grid
        min_x = - radius_physical
        max_x = radius_physical
        min_y = - radius_physical
        max_y = radius_physical
        min_z = -3. * Unit("kpc")
        max_z = 3. * Unit("kpc")

        # Get the pixelscale in physical units
        distance = self.parameters.distance
        pixelscale_angular = self.reference_wcs.xy_average_pixelscale * Unit("pix")  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        # Because we (currently) can't position the grid exactly as the 2D pixels (rotation etc.),
        # take half of the pixel size to avoid too much interpolation
        smallest_scale = smallest_cell_pixels * pixelscale

        # Calculate the minimum division level that is necessary to resolve the smallest scale of the input maps
        extent_x = (max_x - min_x).to("pc").value
        smallest_scale = smallest_scale.to("pc").value
        max_level = min_level_for_smallest_scale_bintree(extent_x, smallest_scale)

        # Create the dust grid
        return BinaryTreeDustGrid(min_x, max_x, min_y, max_y, min_z, max_z, min_level=min_level, max_level=max_level, max_mass_fraction=max_mass_fraction)

    # -----------------------------------------------------------------

    def create_octtree_dust_grid(self, smallest_cell_pixels, min_level, max_mass_fraction):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the octtree dust grid ...")

        # Calculate the major radius of the truncation ellipse in physical coordinates (pc)
        major_angular = self.ellipse.major  # major axis length of the sky ellipse
        radius_physical = (major_angular * self.parameters.distance).to("pc", equivalencies=dimensionless_angles())

        # Calculate the boundaries of the dust grid
        min_x = - radius_physical
        max_x = radius_physical
        min_y = - radius_physical
        max_y = radius_physical
        min_z = -3. * Unit("kpc")
        max_z = 3. * Unit("kpc")

        # Get the pixelscale in physical units
        distance = self.parameters.distance
        pixelscale_angular = self.reference_wcs.xy_average_pixelscale * Unit("pix")  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        # Because we (currently) can't position the grid exactly as the 2D pixels (rotation etc.),
        # take half of the pixel size to avoid too much interpolation
        smallest_scale = smallest_cell_pixels * pixelscale

        # Calculate the minimum division level that is necessary to resolve the smallest scale of the input maps
        extent_x = (max_x - min_x).to("pc").value
        smallest_scale = smallest_scale.to("pc").value
        max_level = min_level_for_smallest_scale_octtree(extent_x, smallest_scale)

        # Create the dust grid and return it
        return OctTreeDustGrid(min_x, max_x, min_y, max_y, min_z, max_z, min_level=min_level, max_level=max_level, max_mass_fraction=max_mass_fraction)

    # -----------------------------------------------------------------

    def write_dust_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust grids ...")

        # Write the low-resolution dust grid
        path = fs.join(self.fit_grid_lowres_path, "lowres.grid")
        self.lowres_dust_grid.save(path)

        # Write the high-resolution dust grid
        path = fs.join(self.fit_grid_highres_path, "highres.grid")
        self.highres_dust_grid.save(path)

# -----------------------------------------------------------------
