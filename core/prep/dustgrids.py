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
from collections import OrderedDict

# Import astronomical modules
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from ..simulation.grids import BinaryTreeDustGrid, OctTreeDustGrid, CartesianDustGrid
from ...core.basics.log import log
from ...core.basics.range import zip_linear
from ..basics.configurable import Configurable
from ..basics.table import SmartTable
from ..basics.range import RealRange, QuantityRange, IntegerRange
from ..tools import types

# -----------------------------------------------------------------

class DustGridsTable(SmartTable):

    """
    This class ...
    """

    # Add column info
    _column_info = OrderedDict()
    _column_info["Type"] = (str, None, "grid type")
    _column_info["Min x"] = (float, "pc", "minimum x")
    _column_info["Max x"] = (float, "pc", "maximum x")
    _column_info["Min y"] = (float, "pc", "minimum y")
    _column_info["Max y"] = (float, "pc", "maximum y")
    _column_info["Min z"] = (float, "pc", "minimum z")
    _column_info["Max z"] = (float, "pc", "maximum z")
    _column_info["Smallest scale"] = (float, "pc", "Smallest scale")
    _column_info["Min level"] = (int, None, "Minimum level")
    _column_info["Max mass fraction"] = (float, None, "Maximum mass fraction")

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(DustGridsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_all_column_info(self._column_info)

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
        return - self.x_radius

    # -----------------------------------------------------------------

    @property
    def x_max(self):
        return self.x_radius

    # -----------------------------------------------------------------

    @property
    def x_extent(self):
        return self.x_max - self.x_min

    # -----------------------------------------------------------------

    @property
    def y_min(self):
        return - self.y_radius

    # -----------------------------------------------------------------

    @property
    def y_max(self):
        return self.y_radius

    # -----------------------------------------------------------------

    @property
    def y_extent(self):
        return self.y_max - self.y_min

    # -----------------------------------------------------------------

    @property
    def z_min(self):
        return - self.z_radius

    # -----------------------------------------------------------------

    @property
    def z_max(self):
        return self.z_radius

    # -----------------------------------------------------------------

    @property
    def z_extent(self):
        return self.z_max - self.z_min

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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
        for scale, min_level, mass_fraction in zip_linear(self.scale_range, self.level_range, self.mass_fraction_range, npoints=self.ngrids):

            # Create the grid and add it to the list
            if self.grid_type == "cartesian": self.create_cartesian_dust_grid(scale)
            elif self.grid_type == "bintree": self.create_binary_tree_dust_grid(scale, min_level, mass_fraction)
            elif self.grid_type == "octtree": self.create_octtree_dust_grid(scale, min_level, mass_fraction)
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
        if log.is_debug:
            print("")
            print(grid)
            print("")

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
        if log.is_debug:
            print("")
            print(grid)
            print("")

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
        if log.is_debug:
            print("")
            print(grid)
            print("")

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

def create_one_dust_grid_for_galaxy_from_deprojection(grid_type, deprojection, distance, sky_ellipse, min_level,
                                                      max_mass_fraction, max_ndivisions_per_pixel=2, nscaleheights=10.):

    """
    This function ...
    :param grid_type: 
    :param deprojection:
    :param distance:
    :param sky_ellipse:
    :param min_level:
    :param max_mass_fraction:
    :param max_ndivisions_per_pixel:
    :param nscaleheights:
    :return: 
    """

    if sky_ellipse is not None:
        # Calculate the major radius of the truncation ellipse in physical coordinates (pc)
        semimajor_angular = sky_ellipse.semimajor  # semimajor axis length of the sky ellipse
        radius_physical = (semimajor_angular * distance).to("pc", equivalencies=dimensionless_angles())
    else:
        x_radius_physical = deprojection.x_range.radius
        y_radius_physical = deprojection.y_range.radius
        radius_physical = max(x_radius_physical, y_radius_physical)

    # Get properties
    average_pixelscale = deprojection.pixelscale
    scaleheight = deprojection.scale_height

    # Get the pixelscale in physical units
    if types.is_angle(average_pixelscale):
        pixelscale_angular = average_pixelscale.to("deg")
        # pixelscale_angular = self.reference_wcs.average_pixelscale.to("deg")  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())
    elif types.is_length_quantity(average_pixelscale): pixelscale = average_pixelscale.to("pc") # normally it should be this case (deprojections should have their pixelscale defined in physical units)
    else: raise ValueError("Pixelscale should be an angle or a length quantity")

    # Determine the minimum physical scale
    min_scale = pixelscale / float(max_ndivisions_per_pixel)

    # Create the dust grid
    return create_one_dust_grid_for_galaxy(grid_type, radius_physical, scaleheight, min_scale, min_level, max_mass_fraction, nscaleheights=nscaleheights)

# -----------------------------------------------------------------

def create_one_dust_grid_for_galaxy(grid_type, radius, scaleheight, min_scale, min_level, max_mass_fraction, nscaleheights=10.):

    """
    This function ...
    :param grid_type: 
    :param radius: IN PHYSICAL COORDINATES
    :param scaleheight: IN PHYSICAL COORDINATES
    :param min_scale:
    :param min_level:
    :param max_mass_fraction:
    :param nscaleheights: REAL NUMBER
    :return: 
    """

    # Determine x, y and z radius
    x_radius = radius
    y_radius = radius
    z_radius = scaleheight * nscaleheights

    # X
    x_min = - x_radius
    x_max = x_radius
    x_extent = x_max - x_min

    # Y
    y_min = - y_radius
    y_max = y_radius
    y_extent = y_max - y_min

    # Z
    z_min = - z_radius
    z_max = z_radius
    z_extent = z_max - z_min

    # Create the dust grid
    return create_one_dust_grid(grid_type, min_scale, x_extent, y_extent, z_extent, x_min, x_max, y_min, y_max, z_min, z_max, min_level, max_mass_fraction)

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
    log.info("Creating a cartesian dust grid with a physical scale of " + str(scale) + " ...")

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
    max_level = max_level_for_smallest_scale_bintree(extent_x, smallest_scale)

    # Check arguments
    if x_min is None: raise ValueError("'x_min' is undefined")
    if x_max is None: raise ValueError("'x_max' is undefined")
    if y_min is None: raise ValueError("'y_min' is undefined")
    if y_max is None: raise ValueError("'y_max' is undefined")
    if z_min is None: raise ValueError("'z_min' is undefined")
    if z_max is None: raise ValueError("'z_max' is undefined")
    if min_level is None: raise ValueError("'min_level' is undefined")
    if max_mass_fraction is None: raise ValueError("'max_mass_fraction' is undefined")

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
    max_level = max_level_for_smallest_scale_octtree(extent_x, smallest_scale)

    # Check arguments
    if x_min is None: raise ValueError("'x_min' is undefined")
    if x_max is None: raise ValueError("'x_max' is undefined")
    if y_min is None: raise ValueError("'y_min' is undefined")
    if y_max is None: raise ValueError("'y_max' is undefined")
    if z_min is None: raise ValueError("'z_min' is undefined")
    if z_max is None: raise ValueError("'z_max' is undefined")
    if min_level is None: raise ValueError("'min_level' is undefined")
    if max_mass_fraction is None: raise ValueError("'max_mass_fraction' is undefined")

    # Create the dust grid
    grid = OctTreeDustGrid(min_x=x_min, max_x=x_max, min_y=y_min, max_y=y_max, min_z=z_min, max_z=z_max,
                           min_level=min_level, max_level=max_level, max_mass_fraction=max_mass_fraction)

    # Return the grid
    return grid

# -----------------------------------------------------------------

def max_level_for_smallest_scale_bintree(extent, smallest_scale):

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

def smallest_scale_for_max_level_bintree(extent, max_level):

    """
    This function ...
    :param extent:
    :param max_level:
    :return:
    """

    octtree_level = max_level / 3
    max_ratio = 2**octtree_level
    min_scale = extent / max_ratio
    return min_scale

# -----------------------------------------------------------------

def max_level_for_smallest_scale_octtree(extent, smallest_scale):

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

def smallest_scale_for_max_level_octtree(extent, max_level):

    """
    This function ...
    :param extent:
    :param max_level:
    :return:
    """

    max_ratio = 2**max_level
    min_scale = extent / max_ratio
    return min_scale

# -----------------------------------------------------------------

def smallest_scale_for_dust_grid(grid):

    """
    This function ...
    :param grid:
    :return:
    """

    # Cartesian grid
    if isinstance(grid, CartesianDustGrid):

        min_x_scale = grid.x_extent / float(grid.x_bins)
        min_y_scale = grid.y_extent / float(grid.y_bins)
        min_z_scale = grid.z_extent / float(grid.z_bins)

        # Return the minimum scale
        return min(min_x_scale, min_y_scale, min_z_scale)

    # Octtree
    elif isinstance(grid, OctTreeDustGrid):

        extent = grid.smallest_extent
        max_level = grid.max_level
        return smallest_scale_for_max_level_octtree(extent, max_level)

    # Binary tree
    elif isinstance(grid, BinaryTreeDustGrid):

        extent = grid.smallest_extent
        max_level = grid.max_level
        return smallest_scale_for_max_level_bintree(extent, max_level)

    # Other
    else: raise NotImplementedError("Other dust grids not implemented")

# -----------------------------------------------------------------
