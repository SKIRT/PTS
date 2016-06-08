#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.basics.grids Contains the BinaryTreeGrid, OctTreeGrid and CartesianGrid classes.

# -----------------------------------------------------------------

# Import standard modules
import copy

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class BinaryTreeDustGrid(object):

    """
    This class ...
    """

    def __init__(self, min_x, max_x, min_y, max_y, min_z, max_z, write=True, min_level=6, max_level=30,
                 search_method="Neighbor", sample_count=100, max_optical_depth=0, max_mass_fraction=1e-6,
                 max_dens_disp_fraction=0, direction_method="Alternating"):

        """
        The constructor ...
        :param min_x:
        :param max_x:
        :param min_y:
        :param max_y:
        :param min_z:
        :param max_z:
        :param write:
        :param min_level:
        :param max_level:
        :param search_method:
        :param sample_count:
        :param max_optical_depth:
        :param max_mass_fraction:
        :param max_dens_disp_fraction:
        :param direction_method:
        """

        self.min_x = min_x
        self.max_x = max_x
        self.min_y = min_y
        self.max_y = max_y
        self.min_z = min_z
        self.max_z = max_z
        self.write = write
        self.min_level = min_level
        self.max_level = max_level
        self.search_method = search_method
        self.sample_count = sample_count
        self.max_optical_depth = max_optical_depth
        self.max_mass_fraction = max_mass_fraction
        self.max_dens_disp_fraction = max_dens_disp_fraction
        self.direction_method = direction_method

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """



    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Open the file and write the properties of this model
        with open(path, 'w') as grid_file:

            print("Min x:", str(self.min_x), file=grid_file)
            print("Max x:", str(self.max_x), file=grid_file)
            print("Min y:", str(self.min_y), file=grid_file)
            print("Max y:", str(self.max_x), file=grid_file)
            print("Min z:", str(self.min_z), file=grid_file)
            print("Max z:", str(self.max_z), file=grid_file)
            print("Write:", str(self.write), file=grid_file)
            print("")

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

# -----------------------------------------------------------------

class OctTreeDustGrid(object):

    """
    This class ...
    """

    def __init__(self, min_x, max_x, min_y, max_y, min_z, max_z, write=True, min_level=2, max_level=6,
                 search_method="Neighbor", sample_count=100, max_optical_depth=0, max_mass_fraction=1e-6,
                 max_dens_disp_fraction=0, barycentric=False):

        """
        The constructor ...
        :param min_x:
        :param max_x:
        :param min_y:
        :param max_y:
        :param min_z:
        :param max_z:
        :param write:
        :param min_level:
        :param max_level:
        :param search_method:
        :param sample_count:
        :param max_optical_depth:
        :param max_mass_fraction:
        :param max_dens_disp_fraction:
        :param barycentric:
        """

        self.min_x = min_x
        self.max_x = max_x
        self.min_y = min_y
        self.max_y = max_y
        self.min_z = min_z
        self.max_z = max_z
        self.write = write
        self.min_level = min_level
        self.max_level = max_level
        self.search_method = search_method
        self.sample_count = sample_count
        self.max_optical_depth = max_optical_depth
        self.max_mass_fraction = max_mass_fraction
        self.max_dens_disp_fraction = max_dens_disp_fraction
        self.barycentric = barycentric

# -----------------------------------------------------------------

class CartesianDustGrid(object):

    """
    This class ...
    """

    def __init__(self, min_x, max_x, min_y, max_y, min_z, max_z, x_bins, y_bins, z_bins, mesh_type="linear", ratio=1.,
                 write=True):

        """
        The constructor ...
        :param min_x:
        :param max_x:
        :param min_y:
        :param max_y:
        :param min_z:
        :param max_z:
        :param x_bins:
        :param y_bins:
        :param z_bins:
        :param mesh_type:
        :param ratio:
        :param write:
        """

        self.min_x = min_x
        self.max_x = max_x
        self.min_y = min_y
        self.max_y = max_y
        self.min_z = min_z
        self.max_z = max_z
        self.x_bins = x_bins
        self.y_bins = y_bins
        self.z_bins = z_bins
        self.mesh_type = mesh_type
        self.ratio = ratio
        self.write = write

# -----------------------------------------------------------------
