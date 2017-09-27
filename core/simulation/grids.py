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

def load_grid(path, remote=None):

    """
    This function ...
    :param path:
    :param remote:
    :return:
    """

    from ..tools import filesystem as fs

    # Get the first line of the file
    #with open(path, 'r') as f: first_line = f.readline()
    if remote is not None: first_line = remote.get_first_line(path)
    else: first_line = fs.get_first_line(path)

    # Create and return the appropriate dust grid
    if "BinaryTreeDustGrid" in first_line: return BinaryTreeDustGrid.from_file(path, remote=remote)
    elif "OctTreeDustGrid" in first_line: return OctTreeDustGrid.from_file(path, remote=remote)
    elif "CartesianDustGrid" in first_line: return CartesianDustGrid.from_file(path, remote=remote)
    elif "FileTreeDustGrid" in first_line: return FileTreeDustGrid.from_file(path, remote=remote)
    else: raise ValueError("Unrecognized dust grid file")

# -----------------------------------------------------------------

mesh_types = ["linear", "power", "symmetric_power", "logarithmic"]

# -----------------------------------------------------------------

search_methods = ["Neighbor", "TopDown", "Bookkeeping"]
default_search_method = "Neighbor"

# -----------------------------------------------------------------

cartesian = "cartesian"
bintree = "bintree"
octtree = "octtree"

# -----------------------------------------------------------------

class DustGrid(SimplePropertyComposite):

    """
    This function ...
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
        self.add_property("write", "boolean", "write grid", True)

# -----------------------------------------------------------------

class FileTreeDustGrid(DustGrid):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(FileTreeDustGrid, self).__init__()

        # Define properties
        self.add_property("filename", "string", "tree grid file")
        self.add_property("search_method", "string", "search method", default_search_method, choices=search_methods)

        # Set properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------

class DustGrid3D(DustGrid):

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
        super(DustGrid3D, self).__init__()

        # Define properties
        self.add_property("min_x", "quantity", "minimum x")
        self.add_property("max_x", "quantity", "maximum x")
        self.add_property("min_y", "quantity", "minimum y")
        self.add_property("max_y", "quantity", "maximum y")
        self.add_property("min_z", "quantity", "minimum z")
        self.add_property("max_z", "quantity", "maximum z")

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
    def y_extent(self):

        """
        This function ...
        :return:
        """

        return self.y_max - self.y_min

    # -----------------------------------------------------------------

    @property
    def z_extent(self):

        """
        This function ...
        :return:
        """

        return self.z_max - self.z_min

    # -----------------------------------------------------------------

    @property
    def smallest_extent(self):

        """
        This fucntion ...
        :return:
        """

        return min(self.x_extent, self.y_extent, self.z_extent)

# -----------------------------------------------------------------

class BinaryTreeDustGrid(DustGrid3D):

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
        self.add_property("search_method", "string", "search method", default_search_method, choices=search_methods)
        self.add_property("sample_count", "positive_integer", "number of samples", 100)
        self.add_property("max_optical_depth", "real", "maximum optical depth", 0)
        self.add_property("max_mass_fraction", "real", "maximum mass fraction", 1e-6)
        self.add_property("max_dens_disp_fraction", "real", "maximum density dispersion fraction", 0)
        self.add_property("direction_method", "string", "direction method", "Alternating", choices=["Alternating"])
        self.add_property("write_tree", "boolean", "write tree data", False)

        # Set properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------

class OctTreeDustGrid(DustGrid3D):

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
        self.add_property("search_method", "string", "search method", default_search_method, choices=search_methods)
        self.add_property("sample_count", "positive_integer", "number of samples", 100)
        self.add_property("max_optical_depth", "real", "maximum optical depth", 0)
        self.add_property("max_mass_fraction", "real", "maximum mass fraction", 1e-6)
        self.add_property("max_dens_disp_fraction", "real", "maximum density dispersion fraction", 0)
        self.add_property("barycentric", "boolean", "barycentric", False)
        self.add_property("write_tree", "boolean", "write tree data", False)

        # Set properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------

class CartesianDustGrid(DustGrid3D):

    """
    This class ...
    """

    def __init__(self, **kwargs):

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

class DustGrid2D(DustGrid):

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
        super(DustGrid2D, self).__init__()

# -----------------------------------------------------------------

class CylindricalGrid(DustGrid):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(CylindricalGrid, self).__init__()

        # Add properties
        self.add_property("max_r", "quantity", "maximum radius")
        self.add_property("min_z", "quantity", "minimum z")
        self.add_property("max_z", "quantity", "maximum z")
        self.add_property("nbins_r", "positive_integer", "number of bins in the radial direction")
        self.add_property("nbins_z", "positive_integer", "number of bins in the axial direction")

        # Specific
        self.add_property("type_r", "string", "type of mesh in radial direction", choices=mesh_types)
        self.add_property("type_z", "string", "type of mesh in axial direction", choices=mesh_types)

        self.add_property("central_bin_fraction_r", "real", "central bin fraction in radial direction (for log mesh)")
        self.add_property("central_bin_fraction_z", "real", "central bin fraction in axial direction (for log mesh)")

        self.add_property("ratio_r", "real", "ratio in radial direction (for symmetric_power mesh)")
        self.add_property("ratio_z", "real", "ratio in axial direction (for symmetric_power mesh")

        # Set properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------
