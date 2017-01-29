#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.grids Contains the BinaryTreeGrid, OctTreeGrid and CartesianGrid classes.

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

    def __init__(self):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(DustGrid, self).__init__()

        # Define properties
        self.add_property("min_x", "quantity", "minimum x")
        self.add_property("max_x", "quantity", "maximum x")
        self.add_property("min_y", "quantity", "minimum y")
        self.add_property("max_y", "quantity", "maximum y")
        self.add_property("min_z", "quantity", "minimum z")
        self.add_property("max_z", "quantity", "maximum z")
        self.add_property("write", "boolean", "write grid", True)

# -----------------------------------------------------------------

class BinaryTreeDustGrid(DustGrid):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(BinaryTreeDustGrid, self).__init__()

        # Define properties
        self.add_property("min_level", "positive_integer", "minimum level", 6)
        self.add_property("max_level", "positive_integer", "maximum level", 30)
        self.add_property("search_method", "string", "search method", "Neighbor", choices=["Neighbor", "TopDown", "Bookkeeping"])
        self.add_property("sample_count", "positive_integer", "number of samples", 100)
        self.add_property("max_optical_depth", "real", "maximum optical depth", 0)
        self.add_property("max_mass_fraction", "real", "maximum mass fraction", 1e-6)
        self.add_property("max_dens_disp_fraction", "real", "maximum density dispersion fraction", 0)
        self.add_property("direction_method", "string", "direction method", "Alternating", choices=["Alternating"])

        # Set properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------

class OctTreeDustGrid(DustGrid):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(OctTreeDustGrid, self).__init__()

        # Define properties
        self.add_property("min_level", "positive_integer", "minimum level", 2)
        self.add_property("max_level", "positive_integer", "maximum level", 6)
        self.add_property("search_method", "string", "search method", "Neighbor", choices=["Neighbor", "TopDown", "Bookkeeping"])
        self.add_property("sample_count", "positive_integer", "number of samples", 100)
        self.add_property("max_optical_depth", "real", "maximum optical depth", 0)
        self.add_property("max_mass_fraction", "real", "maximum mass fraction", 1e-6)
        self.add_property("max_dens_disp_fraction", "real", "maximum density dispersion fraction", 0)
        self.add_property("barycentric", "boolean", "barycentric", False)

        # Set properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------

class CartesianDustGrid(DustGrid):

    """
    This class ...
    """

    def __init__(self, *kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(CartesianDustGrid, self).__init__()

        # Define properties
        self.add_property("x_bins", "positive_integer", "number of x bins")
        self.add_property("y_bins", "positive_integer", "number of y bins")
        self.add_property("z_bins", "positive_integer", "number of z bins")
        self.add_property("mesh_type", "string", "mesh type", "linear")
        self.add_property("ratio", "real", "ratio", 1.)

        # Set properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------
