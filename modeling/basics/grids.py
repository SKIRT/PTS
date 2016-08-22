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
from abc import ABCMeta

# Import the relevant PTS classes and modules
from ...core.basics.composite import SimplePropertyComposite

# -----------------------------------------------------------------

def load_grid(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Get the first line of the file
    with open(path, 'r') as f: first_line = f.readline()

    # Create and return the appropriate dust grid
    if "BinaryTreeDustGrid" in first_line: return BinaryTreeDustGrid.from_file(path)
    elif "OctTreeDustGrid" in first_line: return OctTreeDustGrid.from_file(path)
    elif "CartesianDustGrid" in first_line: return CartesianDustGrid.from_file(path)
    else: raise ValueError("Unrecognized dust grid file")

# -----------------------------------------------------------------

class DustGrid(SimplePropertyComposite):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

# -----------------------------------------------------------------

class BinaryTreeDustGrid(DustGrid):

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

class OctTreeDustGrid(DustGrid):

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

class CartesianDustGrid(DustGrid):

    """
    This class ...
    """

    def __init__(self, min_x, max_x, min_y, max_y, min_z, max_z, x_bins, y_bins, z_bins, mesh_type="linear", ratio=1., write=True):

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
