#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.basics.grids Contains the BinaryTreeGrid, OctTreeGrid and CartesianGrid classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy

# Import the relevant PTS classes and modules
from ...core.tools import parsing

# -----------------------------------------------------------------

def load_grid(path):

    """
    This function ...
    :param path:
    :return:
    """

    return None

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

        effective_radius = None
        index = None
        flattening = None
        tilt = None

        # Read the parameter file
        with open(path, 'r') as model_file:

            # Loop over all lines in the file
            for line in model_file:

                # Split the line
                splitted = line.split(": ")
                splitted[1] = splitted[1].split("\n")[0]

                first = splitted[0]
                second = splitted[1]

                if first == "Min x": min_x = parsing.quantity(second)
                elif first == "Max x": max_x = parsing.quantity(second)
                elif first == "Min y": min_y = parsing.quantity(second)
                elif first == "Max y": max_y = parsing.quantity(second)
                elif first == "Min z": min_z = parsing.quantity(second)
                elif first == "Max z": max_z = parsing.quantity(second)
                elif first == "Write": write = parsing.boolean(second)
                elif first == "Min level": min_level = int(second)
                elif first == "Max level": max_level = int(second)
                elif first == "Search method": search_method = second

        # Creaete the SersicModel and return it
        return cls(effective_radius, index, flattening, tilt)

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Open the file and write the properties of this grid
        with open(path, 'w') as grid_file:

            print("Min x:", str(self.min_x), file=grid_file)
            print("Max x:", str(self.max_x), file=grid_file)
            print("Min y:", str(self.min_y), file=grid_file)
            print("Max y:", str(self.max_x), file=grid_file)
            print("Min z:", str(self.min_z), file=grid_file)
            print("Max z:", str(self.max_z), file=grid_file)
            print("Write:", str(self.write), file=grid_file)
            print("Min level:", str(self.min_level), file=grid_file)
            print("Max level:", str(self.max_level), file=grid_file)
            print("Search method:", str(self.search_method), file=grid_file)
            print("Sample count:", str(self.sample_count), file=grid_file)
            print("Max optical depth:", str(self.max_optical_depth), file=grid_file)
            print("Max mass fraction:", str(self.max_mass_fraction), file=grid_file)
            print("Max dens disp fraction:", str(self.max_dens_disp_fraction), file=grid_file)
            print("Direction method:", str(self.direction_method), file=grid_file)

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

        # Open the file and write the properties of this grid
        with open(path, 'w') as grid_file:

            print("Min x:", str(self.min_x), file=grid_file)
            print("Max x:", str(self.max_x), file=grid_file)
            print("Min y:", str(self.min_y), file=grid_file)
            print("Max y:", str(self.max_x), file=grid_file)
            print("Min z:", str(self.min_z), file=grid_file)
            print("Max z:", str(self.max_z), file=grid_file)
            print("Write:", str(self.write), file=grid_file)
            print("Min level:", str(self.min_level), file=grid_file)
            print("Max level:", str(self.max_level), file=grid_file)
            print("Search method:", str(self.search_method), file=grid_file)
            print("Sample count:", str(self.sample_count), file=grid_file)
            print("Max optical depth:", str(self.max_optical_depth), file=grid_file)
            print("Max mass fraction:", str(self.max_mass_fraction), file=grid_file)
            print("Max dens disp fraction:", str(self.max_dens_disp_fraction), file=grid_file)
            print("Barycentric", str(self.barycentric), file=grid_file)

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

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

        # Open the file and write the properties of this grid
        with open(path, 'w') as grid_file:

            print("Min x:", str(self.min_x), file=grid_file)
            print("Max x:", str(self.max_x), file=grid_file)
            print("Min y:", str(self.min_y), file=grid_file)
            print("Max y:", str(self.max_x), file=grid_file)
            print("Min z:", str(self.min_z), file=grid_file)
            print("Max z:", str(self.max_z), file=grid_file)
            print("X bins:", str(self.x_bins), file=grid_file)
            print("Y bins:", str(self.y_bins), file=grid_file)
            print("Z bins:", str(self.z_bins), file=grid_file)
            print("Mesh type:", str(self.mesh_type), file=grid_file)
            print("Ratio:", str(self.ratio), file=grid_file)
            print("Write:", str(self.write), file=grid_file)

# -----------------------------------------------------------------
