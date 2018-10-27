#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.projection.data Contains the DataProjections class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...magic.core.frame import Frame
from ...magic.core.image import Image
from ...magic.basics.vector import PixelShape
from ...core.tools import numbers
from ...core.basics.range import RealRange
from .general import is_faceon, is_edgeon, get_faceon_projection, get_edgeon_projection
from ...core.tools.utils import lazyproperty
from ...core.tools import sequences
from ...core.basics.log import log
from ...core.tools import types
from ...magic.basics.vector import Extent
from ...core.tools.stringify import tostr
from ...core.tools.progress import Bar
from ..basics.projection import FaceOnProjection, EdgeOnProjection
from ..core.data import Data3D

# -----------------------------------------------------------------

# Number of photon packages
default_npackages = 5e7

# Instruments/orientations
earth_name = "earth"
faceon_name = "faceon"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

def project_3d(name, x_coordinates, y_coordinates, z_coordinates, values, projection, length_unit, unit=None,
               return_stddev=False, return_ncells=False, weights=None, description=None, as_image=False):

    """
    This function projects 3D data into a 2D map
    :param name: name of the kind of data
    :param x_coordinates:
    :param y_coordinates:
    :param z_coordinates:
    :param values:
    :param projection:
    :param length_unit:
    :param unit:
    :param return_stddev:
    :param return_ncells:
    :param weights:
    :param description:
    :param as_image:
    :return:
    """

    # Create the data
    data = Data3D(name, x_coordinates, y_coordinates, z_coordinates, values, weights=weights, length_unit=length_unit,
                  unit=unit, description=description)

    # Create projection
    return project_data(name, data, projection, return_stddev=return_stddev, return_ncells=return_ncells,
                        description=description, as_image=as_image)

# -----------------------------------------------------------------

def project_data(name, data, projection, return_stddev=False, return_ncells=False, description=None, as_image=False,
                 height=None, width=None, cell_based=True, interpolate=False): # default for 'cell_based' was False before

    """
    This function ...
    :param name:
    :param data:
    :param projection:
    :param return_stddev:
    :param return_ncells:
    :param description:
    :param as_image:
    :param height:
    :param width:
    :param cell_based:
    :param interpolate:
    :return:
    """

    # Create faceon or edgeon projection
    faceon = is_faceon(projection)
    edgeon = is_edgeon(projection)
    if faceon:
        faceon_projection = get_faceon_projection(projection, strict=True)
        edgeon_projection = None
        if width is not None: raise ValueError("Width cannot be specified when projecting face-on")
    elif edgeon:
        faceon_projection = None
        edgeon_projection = get_edgeon_projection(projection, strict=True)
        if height is not None: raise ValueError("Height cannot be specified when projecting edge-on")
    else: raise ValueError("Projection must be face-on or edge-on")

    # Create the projections object
    projections = DataProjections(data, name=name, description=description, faceon=faceon, edgeon=edgeon,
                                  projection_faceon=faceon_projection, projection_edgeon=edgeon_projection,
                                  faceon_height=height, edgeon_width=width, cell_based=cell_based)

    # Get the projected map
    if faceon: frame = projections.faceon
    elif edgeon: frame = projections.edgeon
    else: raise RuntimeError("Something went wrong")
    if frame is None: raise RuntimeError("Something went wrong: projected frame not created")

    # Get the stddev
    #if return_stddev:
    if faceon: stddev_frame = projections.faceon_stddev
    elif edgeon: stddev_frame = projections.edgeon_stddev
    else: raise RuntimeError("Something went wrong")
    #else: stddev_frame = None

    # Get the ncells
    #if return_ncells:
    if faceon: ncells_frame = projections.faceon_ncells
    elif edgeon: ncells_frame = projections.edgeon_ncells
    else: raise RuntimeError("Something went wrong")
    #else: ncells_frame = None

    # Create interpolated version?
    original_frame = frame
    if interpolate: frame = create_interpolated_map(frame, ncells_frame)

    # Return as image
    if as_image:

        # Create the image
        image = Image.from_frame(frame, name=name)
        if interpolate: image.add_frame(original_frame, name="original")
        if return_stddev: image.add_frame(stddev_frame, name="stddev")
        if return_ncells: image.add_frame(ncells_frame, name="ncells")

        # Return the image
        return image

    # Return as individual frames
    else:

        if return_stddev:
            if return_ncells: return frame, stddev_frame, ncells_frame
            else: return frame, stddev_frame
        elif return_ncells: return frame, ncells_frame
        else: return frame

# -----------------------------------------------------------------

def create_interpolated_map(frame, ncells=None, min_ncells=10, replace_nans=True, interpolate_below=None, interpolate_above=None,
                            not_nans_dilation_radius=3., kernel_sigma=2.):

    """
    This function ...
    :param frame:
    :param ncells:
    :param min_ncells:
    :param replace_nans:
    :param interpolate_below:
    :param interpolate_above:
    :param not_nans_dilation_radius:
    :param kernel_sigma:
    :return:
    """

    # Copy
    interpolated = frame.copy()

    #from ...magic.tools import plotting

    # Get outside nans
    outside_nans = interpolated.nans.largest()
    not_nans = outside_nans.inverse()

    #plotting.plot_mask(not_nans, title="not nans")

    not_nans.disk_dilate(radius=not_nans_dilation_radius)

    #plotting.plot_mask(not_nans, title="not nans (2)")

    do_nans = not_nans.largest().inverse()

    #plotting.plot_mask(do_nans, title="do nans")

    # Get mask
    #
    #where = where * do_nans.inverse()  # don't interpolate outside (where ncells = 0)

    # Interpolate based on ncells
    if ncells is not None: where = ncells.where_smaller_than(min_ncells)

    # Interpolate NaNs
    else: where = interpolated.nans #+ interpolated.where_smaller_than(1e-12)

    # Mask below?
    if interpolate_below is not None: where = where + interpolated.where_smaller_than(interpolate_below)

    # Mask above?
    if interpolate_above is not None: where = where + interpolated.where_greater_than(interpolate_above)

    # Create mask of where to interpolate
    where = where * do_nans.inverse()

    #plotting.plot_mask(where, title="where")
    # plotting.plot_mask(where, title="where smaller than " + str(self.config.min_ncells))
    # plotting.plot_mask(self.map_interpolated.nans, title="nans")

    # Replace NaNs to zero that have to stay NaNs (don't interpolate)
    if replace_nans:
        interpolated[do_nans] = 0.0
        do_nans.disk_dilate(radius=not_nans_dilation_radius)
        #plotting.plot_mask(do_nans, title="do nans")

    # Put pixels to NaN
    interpolated.replace_by_nans(where)

    # plotting.plot_mask(self.map_interpolated.nans, title="nans")
    # exit()

    #plotting.plot_frame(interpolated)

    # Interpolate nans
    interpolated.interpolate_nans(sigma=kernel_sigma, error_on_max=replace_nans)
    interpolated.replace_by_nans(do_nans)

    # Return the interpolated frame
    return interpolated

# -----------------------------------------------------------------

def project_faceon(name, data, return_stddev=False, return_ncells=False, description=None, spacing="mean",
                   spacing_factor=1., height=None, as_image=False, cell_based=True):

    """
    This function ...
    :param name:
    :param data:
    :param return_stddev:
    :param return_ncells:
    :param description:
    :param spacing:
    :param spacing_factor:
    :param height:
    :param as_image:
    :param cell_based:
    :return:
    """

    # Create the projections object
    projections = DataProjections(data, name=name, description=description, faceon=True, edgeon=False,
                                  faceon_spacing=spacing, faceon_spacing_factor=spacing_factor, faceon_height=height,
                                  cell_based=cell_based)

    # Get results
    frame = projections.faceon
    stddev_frame = projections.faceon_stddev
    ncells_frame = projections.faceon_ncells
    if frame is None: raise RuntimeError("Something went wrong: projected frame not created")
    if stddev_frame is None: raise RuntimeError("Something went wrong: stddev frame not created")
    if ncells_frame is None: raise RuntimeError("Something went wrong: ncells frame not created")

    # Return as image
    if as_image:

        # Create the image
        image = Image.from_frame(frame, name=name)
        if return_stddev: image.add_frame(stddev_frame, name="stddev")
        if return_ncells: image.add_frame(ncells_frame, name="ncells")

        # Return the image
        return image

    # Return as seperate frames
    else:

        if return_stddev:
            if return_ncells: return frame, stddev_frame, ncells_frame
            else: return frame, stddev_frame
        elif return_ncells: return frame, ncells_frame
        else: return frame

# -----------------------------------------------------------------

def project_edgeon(name, data, return_stddev=False, return_ncells=False, description=None, spacing="mean",
                   spacing_factor=1., width=None, as_image=False, cell_based=True):

    """
    This function ...
    :param name:
    :param data:
    :param return_stddev:
    :param return_ncells:
    :param description:
    :param spacing:
    :param spacing_factor:
    :param width:
    :param as_image:
    :param cell_based:
    :return:
    """

    # Create the projections object
    projections = DataProjections(data, name=name, description=description, faceon=False, edgeon=True,
                                  edgeon_spacing=spacing, edgeon_spacing_factor=spacing_factor, edgeon_width=width,
                                  cell_based=cell_based)

    # Get results
    frame = projections.edgeon
    stddev_frame = projections.edgeon_stddev
    ncells_frame = projections.edgeon_ncells
    if frame is None: raise RuntimeError("Something went wrong: projected frame not created")
    if stddev_frame is None: raise RuntimeError("Something went wrong: projected stddev frame not created")
    if ncells_frame is None: raise RuntimeError("Something went wrong: projected ncells frame not created")

    # Return as image
    if as_image:

        # Create the image
        image = Image.from_frame(frame, name=name)
        if return_stddev: image.add_frame(stddev_frame, name="stddev")
        if return_ncells: image.add_frame(ncells_frame, name="ncells")

        # Return the image
        return image

    # Return as separate frames
    else:

        if return_stddev:
            if return_ncells: return frame, stddev_frame, ncells_frame
            else: return frame, stddev_frame
        elif return_ncells: return frame, ncells_frame
        else: return frame

# -----------------------------------------------------------------

class DataProjections(object):

    """
    This class ...
    """

    def __init__(self, data, name=None, projection_faceon=None, projection_edgeon=None,
                 faceon=True, edgeon=True, description=None, distance=None, faceon_height=None, edgeon_width=None,
                 faceon_spacing="mean", edgeon_spacing="mean", faceon_spacing_factor=1., edgeon_spacing_factor=1.,
                 logfreq=100, cell_based=True): # default of 'cell_based' was False before

        """
        The constructor ...
        :param name:
        :param data:
        :param projection_faceon:
        :param projection_edgeon:DataProjections
        :param path:
        :param faceon:
        :param edgeon:
        :param description:
        :param distance:
        :param faceon_height:
        :param edgeon_width:
        :param faceon_spacing:
        :param edgeon_spacing:
        :param faceon_spacing_factor:
        :param edgeon_spacing_factor:
        :param logfreq:
        :param cell_based:
        """

        # Set the data
        self.data = data

        # Name and description
        if name is not None: self.name = name
        self.description = description

        # Other
        self.distance = distance

        # Set the projections
        self.projection_faceon = projection_faceon
        self.projection_edgeon = projection_edgeon

        # Set the height and depth (to make cuts of midplane or vertical plane)
        if faceon_height is not None:
            if types.is_length_quantity(faceon_height): faceon_height = faceon_height.to(self.length_unit).value
            elif types.is_real_type(faceon_height): pass
            else: raise ValueError("Invalid type for 'faceon_height'")
            self.faceon_height = faceon_height
        else: self.faceon_height = None
        if edgeon_width is not None:
            if types.is_length_quantity(edgeon_width): edgeon_width = edgeon_width.to(self.length_unit).value
            elif types.is_real_type(edgeon_width): pass
            else: raise ValueError("Invalid type for 'edgeon_width'")
            self.edgeon_width = edgeon_width
        else: self.edgeon_width = None

        # Create projections?
        if faceon and not self.has_projection_faceon: self.create_projection_faceon(spacing=faceon_spacing, spacing_factor=faceon_spacing_factor)
        if edgeon and not self.has_projection_edgeon: self.create_projection_edgeon(spacing=edgeon_spacing, spacing_factor=edgeon_spacing_factor)

        # The maps
        self.faceon = None
        self.edgeon = None

        # The stddev maps
        self.faceon_stddev = None
        self.edgeon_stddev = None

        # The ncells maps
        self.faceon_ncells = None
        self.edgeon_ncells = None

        # Create?
        if faceon: self.project_faceon(logfreq=logfreq, cell_based=cell_based)
        if edgeon: self.project_edgeon(logfreq=logfreq, cell_based=cell_based)

    # -----------------------------------------------------------------

    @lazyproperty
    def name(self): # for when name is not a required argument in the constructor anymore..?
        return self.data.name

    # -----------------------------------------------------------------

    @lazyproperty
    def distance(self):
        if self.has_projection_faceon: return self.projection_faceon.distance
        elif self.has_projection_edgeon: return self.projection_edgeon.distance
        else: return None

    # -----------------------------------------------------------------
    # DATA
    # -----------------------------------------------------------------

    @property
    def has_weights(self):
        return self.data.has_weights

    # -----------------------------------------------------------------

    @lazyproperty
    def ncoordinates(self):
        return len(self.x)

    # -----------------------------------------------------------------

    @property
    def x(self):
        return self.data.valid_x

    # -----------------------------------------------------------------

    @lazyproperty
    def absolute_x(self):
        return abs(self.x)

    # -----------------------------------------------------------------

    @property
    def y(self):
        return self.data.valid_y

    # -----------------------------------------------------------------

    @lazyproperty
    def absolute_y(self):
        return abs(self.y)

    # -----------------------------------------------------------------

    @property
    def z(self):
        return self.data.valid_z

    # -----------------------------------------------------------------

    @lazyproperty
    def absolute_z(self):
        return abs(self.z)

    # -----------------------------------------------------------------

    @property
    def radii(self):
        return self.data.valid_radii

    # -----------------------------------------------------------------

    @property
    def values(self):
        return self.data.valid_values

    # -----------------------------------------------------------------

    @property
    def weights(self):
        return self.data.valid_weights

    # -----------------------------------------------------------------
    # DATA FOR FACE-ON MAP
    # -----------------------------------------------------------------

    @property
    def has_faceon_height(self):
        return self.faceon_height is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def within_faceon_height(self):
        return self.absolute_z <= self.faceon_height

    # -----------------------------------------------------------------

    @lazyproperty
    def nfaceon_coordinates(self):
        if self.has_faceon_height: return np.sum(self.within_faceon_height)
        else: return self.ncoordinates

    # -----------------------------------------------------------------

    @lazyproperty
    def x_faceon(self):
        if self.has_faceon_height: return self.x[self.within_faceon_height]
        else: return self.x

    # -----------------------------------------------------------------

    @lazyproperty
    def y_faceon(self):
        if self.has_faceon_height: return self.y[self.within_faceon_height]
        else: return self.y

    # -----------------------------------------------------------------

    @lazyproperty
    def z_faceon(self):
        if self.has_faceon_height: return self.z[self.within_faceon_height]
        else: return self.z

    # -----------------------------------------------------------------

    @lazyproperty
    def radii_faceon(self):
        if self.has_faceon_height: return self.radii[self.within_faceon_height]
        else: return self.radii

    # -----------------------------------------------------------------

    @lazyproperty
    def values_faceon(self):
        if self.has_faceon_height: return self.values[self.within_faceon_height]
        else: return self.values

    # -----------------------------------------------------------------

    @lazyproperty
    def weights_faceon(self):
        if self.has_faceon_height: return self.values[self.within_faceon_height]
        else: return self.weights

    # -----------------------------------------------------------------
    # DATA FOR EDGE-ON MAP
    # -----------------------------------------------------------------

    @property
    def has_edgeon_width(self):
        return self.edgeon_width is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def within_edgeon_width(self):
        return self.absolute_x <= self.edgeon_width

    # -----------------------------------------------------------------

    @lazyproperty
    def nedgeon_coordinates(self):
        if self.has_edgeon_width: return np.sum(self.within_edgeon_width)
        else: return self.ncoordinates

    # -----------------------------------------------------------------

    @lazyproperty
    def x_edgeon(self):
        if self.has_edgeon_width: return self.x[self.within_edgeon_width]
        else: return self.x

    # -----------------------------------------------------------------

    @lazyproperty
    def y_edgeon(self):
        if self.has_edgeon_width: return self.y[self.within_edgeon_width]
        else: return self.y

    # -----------------------------------------------------------------

    @lazyproperty
    def z_edgeon(self):
        if self.has_edgeon_width: return self.z[self.within_edgeon_width]
        else: return self.z

    # -----------------------------------------------------------------

    @lazyproperty
    def radii_edgeon(self):
        if self.has_edgeon_width: return self.radii[self.within_edgeon_width]
        else: return self.radii

    # -----------------------------------------------------------------

    @lazyproperty
    def values_edgeon(self):
        if self.has_edgeon_width: return self.values[self.within_edgeon_width]
        else: return self.values

    # -----------------------------------------------------------------

    @lazyproperty
    def weights_edgeon(self):
        if self.has_edgeon_width: return self.weights[self.within_edgeon_width]
        else: return self.weights

    # -----------------------------------------------------------------
    # UNITS
    # -----------------------------------------------------------------

    @property
    def unit(self):
        return self.data.unit

    # -----------------------------------------------------------------

    @property
    def length_unit(self):
        return self.data.length_unit

    # -----------------------------------------------------------------
    # SORTED UNIQUE COORDINATES
    # -----------------------------------------------------------------

    @lazyproperty
    def sorted_unique_x_coordinates(self):
        return np.sort(np.unique(self.x))

    # -----------------------------------------------------------------

    @lazyproperty
    def sorted_unique_y_coordinates(self):
        return np.sort(np.unique(self.y))

    # -----------------------------------------------------------------

    @lazyproperty
    def sorted_unique_z_coordinates(self):
        return np.sort(np.unique(self.z))

    # -----------------------------------------------------------------
    # MIN COORDINATES
    # -----------------------------------------------------------------

    @lazyproperty
    def min_x_coordinate(self):
        return self.sorted_unique_x_coordinates[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def min_y_coordinate(self):
        return self.sorted_unique_y_coordinates[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def min_z_coordinate(self):
        return self.sorted_unique_z_coordinates[0]

    # -----------------------------------------------------------------
    # MAX COORDINATES
    # -----------------------------------------------------------------

    @lazyproperty
    def max_x_coordinate(self):
        return self.sorted_unique_x_coordinates[-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def max_y_coordinate(self):
        return self.sorted_unique_y_coordinates[-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def max_z_coordinate(self):
        return self.sorted_unique_z_coordinates[-1]

    # -----------------------------------------------------------------
    # SPACINGS
    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_spacings(self):
        return np.diff(self.sorted_unique_x_coordinates)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_spacings(self):
        return np.diff(self.sorted_unique_y_coordinates)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_z_coordinates_spacings(self):
        return np.diff(self.sorted_unique_z_coordinates)

    # -----------------------------------------------------------------
    # AVERAGE SPACING
    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_average_spacing(self):
        return np.mean(self.unique_x_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_average_spacing(self):
        return np.mean(self.unique_y_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_z_coordinates_average_spacing(self):
        return np.mean(self.unique_z_coordinates_spacings)

    # -----------------------------------------------------------------
    # MEDIAN SPACING
    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_median_spacing(self):
        return np.median(self.unique_x_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_median_spacing(self):
        return np.median(self.unique_y_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_z_coordinates_median_spacing(self):
        return np.median(self.unique_z_coordinates_spacings)

    # -----------------------------------------------------------------
    # STD SPACING
    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_stddev_spacing(self):
        return np.std(self.unique_x_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_stddev_spacing(self):
        return np.std(self.unique_y_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_z_coordinates_stddev_spacing(self):
        return np.std(self.unique_z_coordinates_spacings)

    # -----------------------------------------------------------------
    # MIN SPACING
    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_min_spacing(self):
        return np.min(self.unique_x_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_min_spacing(self):
        return np.min(self.unique_y_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_z_coordinates_min_spacing(self):
        return np.min(self.unique_z_coordinates_spacings)

    # -----------------------------------------------------------------
    # MAX SPACING
    # -----------------------------------------------------------------

    @lazyproperty
    def unique_x_coordinates_max_spacing(self):
        return np.max(self.unique_x_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_y_coordinates_max_spacing(self):
        return np.max(self.unique_y_coordinates_spacings)

    # -----------------------------------------------------------------

    @lazyproperty
    def unique_z_coordinates_max_spacing(self):
        return np.max(self.unique_z_coordinates_spacings)

    # -----------------------------------------------------------------
    # RADIUS & HEIGHT
    # -----------------------------------------------------------------

    @lazyproperty
    def radius(self):
        return max(abs(self.min_x_coordinate), self.max_x_coordinate, abs(self.min_y_coordinate), self.max_y_coordinate)

    # -----------------------------------------------------------------

    @lazyproperty
    def height(self):
        return max(abs(self.min_z_coordinate), self.max_z_coordinate)

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def has_projection_faceon(self):
        return self.projection_faceon is not None

    # -----------------------------------------------------------------

    @property
    def has_projection_edgeon(self):
        return self.projection_edgeon is not None

    # -----------------------------------------------------------------

    def get_xy_spacing(self, measure, factor=1.):

        """
        This function ...
        :param measure:
        :param factor:
        :return:
        """

        # Determine
        if measure == "min": spacing = np.mean([self.unique_x_coordinates_min_spacing, self.unique_y_coordinates_min_spacing])
        elif measure == "max": spacing = np.mean([self.unique_x_coordinates_max_spacing, self.unique_y_coordinates_max_spacing])
        elif measure == "mean": spacing = np.mean([self.unique_x_coordinates_average_spacing, self.unique_y_coordinates_average_spacing])
        elif measure == "median": spacing = np.mean([self.unique_x_coordinates_median_spacing, self.unique_y_coordinates_median_spacing])
        else: raise ValueError("Invalid measure '" + measure + "'")

        # Return
        return spacing * factor * self.length_unit

    # -----------------------------------------------------------------

    def create_projection_faceon(self, spacing="mean", spacing_factor=1.):

        """
        This function ...
        :param spacing:
        :param spacing_factor:
        :return:
        """

        # Inform the user
        log.info("Creating the face-on projection ...")

        # Get scalar value of spacing
        if types.is_string_type(spacing): spacing = self.get_xy_spacing(spacing, factor=spacing_factor)
        elif types.is_length_quantity(spacing): spacing = spacing.to(self.length_unit)
        elif types.is_real_type(spacing): spacing = spacing * self.length_unit
        else: raise ValueError("Invalid value for 'spacing'")

        # Set field
        field_x = 2. * self.radius * self.length_unit
        field_y = 2. * self.radius * self.length_unit

        # Set center of galaxy in frame
        #center_x = self.radius * self.length_unit
        #center_y = self.radius * self.length_unit
        #center_x = 0.0 * self.length_unit
        #center_y = 0.0 * self.length_unit

        # Set number of pixels
        nx = numbers.round_up_to_odd_integer(field_x / spacing)
        ny = numbers.round_up_to_odd_integer(field_y / spacing)

        # Determine pixelscale
        pixelscale_x = field_x / nx
        pixelscale_y = field_y / ny

        # Determine center
        center_x = -0.5 * pixelscale_x
        center_y = -0.5 * pixelscale_y

        # Create the projection
        self.projection_faceon = FaceOnProjection(distance=self.distance, pixels_x=nx,
                                                  pixels_y=ny, center_x=center_x,
                                                  center_y=center_y, field_x=field_x,
                                                  field_y=field_y)

    # -----------------------------------------------------------------

    def get_yz_spacing(self, measure, factor=1.):

        """
        This function ...
        :param measure:
        :param factor:
        :return:
        """

        # Determine
        if measure == "min": spacing = np.mean([self.unique_y_coordinates_min_spacing, self.unique_z_coordinates_min_spacing])
        elif measure == "max": spacing = np.mean([self.unique_y_coordinates_max_spacing, self.unique_z_coordinates_max_spacing])
        elif measure == "mean": spacing = np.mean([self.unique_y_coordinates_average_spacing, self.unique_z_coordinates_average_spacing])
        elif measure == "median": spacing = np.mean([self.unique_y_coordinates_median_spacing, self.unique_z_coordinates_median_spacing])
        else: raise ValueError("Invalid measure '" + measure + "'")

        # Return
        return spacing * factor * self.length_unit

    # -----------------------------------------------------------------

    def create_projection_edgeon(self, spacing="mean", spacing_factor=1.):

        """
        This function ...
        :param spacing:
        :param spacing_factor:
        :return:
        """

        # Inform the user
        log.info("Creating the edge-on projection ...")

        # Get scalar value of spacing
        if types.is_string_type(spacing): spacing = self.get_yz_spacing(spacing, factor=spacing_factor)
        elif types.is_length_quantity(spacing): spacing = spacing.to(self.length_unit)
        elif types.is_real_type(spacing): spacing = spacing * self.length_unit
        else: raise ValueError("Invalid value for 'spacing'")

        # Set field
        field_y = 2. * self.radius * self.length_unit
        field_z = 2. * self.height * self.length_unit

        # Set center of galaxy in frame
        #center_y = self.radius * self.length_unit
        #center_z = self.height * self.length_unit
        #center_y = 0.0 * self.length_unit
        #center_z = 0.0 * self.length_unit

        # Set number of pixels
        ny = numbers.round_up_to_odd_integer(field_y / spacing)
        nz = numbers.round_up_to_odd_integer(field_z / spacing)

        # Determine pixelscale
        pixelscale_y = field_y / ny
        pixelscale_z = field_z / nz

        # Determine center
        center_y = -0.5 * pixelscale_y
        center_z = -0.5 * pixelscale_z

        # Create projection
        self.projection_edgeon = EdgeOnProjection(distance=self.distance, pixels_x=ny,
                                                  pixels_y=nz, center_x=center_y,
                                                  center_y=center_z, field_x=field_y,
                                                  field_y=field_z)

    # -----------------------------------------------------------------

    @property
    def faceon_nx(self):
        return self.projection_faceon.pixels_x

    # -----------------------------------------------------------------

    @property
    def faceon_ny(self):
        return self.projection_faceon.pixels_y

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_npixels(self):
        return self.faceon_nx * self.faceon_ny

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_shape(self):
        return PixelShape.from_xy(self.faceon_nx, self.faceon_ny)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_x_min(self):
        return self.projection_faceon.center_x - 0.5 * self.projection_faceon.field_x

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_x_max(self):
        return self.projection_faceon.center_x + 0.5 * self.projection_faceon.field_y

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_y_min(self):
        return self.projection_faceon.center_y - 0.5 * self.projection_faceon.field_y

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_y_max(self):
        return self.projection_faceon.center_y + 0.5 * self.projection_faceon.field_y

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_x_min_value(self):
        #print("center X", self.projection_faceon.center_x)
        #print("field X", self.projection_faceon.field_x)
        #print("xmin", self.faceon_x_min)
        return self.faceon_x_min.to(self.length_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_x_max_value(self):
        return self.faceon_x_max.to(self.length_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_y_min_value(self):
        return self.faceon_y_min.to(self.length_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_y_max_value(self):
        return self.faceon_y_max.to(self.length_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixelscale(self):
        return self.projection_faceon.physical_pixelscale

    # -----------------------------------------------------------------

    @property
    def faceon_pixelscale_x(self):
        return self.faceon_pixelscale.x

    # -----------------------------------------------------------------

    @property
    def faceon_pixelscale_y(self):
        return self.faceon_pixelscale.y

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixelscale_value(self):
        return Extent(self.faceon_pixelscale_x.to(self.length_unit).value, self.faceon_pixelscale_y.to(self.length_unit).value)

    # -----------------------------------------------------------------

    @property
    def faceon_pixelscale_value_x(self):
        return self.faceon_pixelscale_value.x

    # -----------------------------------------------------------------

    @property
    def faceon_pixelscale_value_y(self):
        return self.faceon_pixelscale_value.y

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_half_pixelscale(self):
        return self.faceon_pixelscale * 0.5

    # -----------------------------------------------------------------

    @property
    def faceon_half_pixelscale_x(self):
        return self.faceon_half_pixelscale.x

    # -----------------------------------------------------------------

    @property
    def faceon_half_pixelscale_y(self):
        return self.faceon_half_pixelscale.y

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_half_pixelscale_value(self):
        return Extent(self.faceon_half_pixelscale_x.to(self.length_unit).value, self.faceon_half_pixelscale_y.to(self.length_unit).value)

    # -----------------------------------------------------------------

    @property
    def faceon_half_pixelscale_value_x(self):
        return self.faceon_half_pixelscale_value.x

    # -----------------------------------------------------------------

    @property
    def faceon_half_pixelscale_value_y(self):
        return self.faceon_half_pixelscale_value.y

    # -----------------------------------------------------------------

    @property
    def faceon_distance(self):
        return self.projection_faceon.distance

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_metadata(self):

        """
        Thisf unction ...
        :return:
        """

        meta = OrderedDict()
        meta["x_min"] = repr(self.faceon_x_min_value) + " " + tostr(self.length_unit)
        meta["x_max"] = repr(self.faceon_x_max_value) + " " + tostr(self.length_unit)
        meta["y_min"] = repr(self.faceon_y_min_value) + " " + tostr(self.length_unit)
        meta["y_max"] = repr(self.faceon_y_max_value) + " " + tostr(self.length_unit)
        return meta

    # -----------------------------------------------------------------

    def get_faceon_coordinate_ranges(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        # Determine range of x and y
        x_range = RealRange(x - self.faceon_half_pixelscale_value_x, x + self.faceon_half_pixelscale_value_x)
        y_range = RealRange(y - self.faceon_half_pixelscale_value_y, y + self.faceon_half_pixelscale_value_y)

        # Return
        return x_range, y_range

    # -----------------------------------------------------------------

    def get_faceon_coordinate_mask(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        # Get range
        x_range, y_range = self.get_faceon_coordinate_ranges(x, y)

        # Create masks
        x_mask = (x_range.min < self.x_faceon) * (self.x_faceon <= x_range.max)
        y_mask = (y_range.min < self.y_faceon) * (self.y_faceon <= y_range.max)

        # Create and return the combined mask
        return x_mask * y_mask

    # -----------------------------------------------------------------

    def get_faceon_coordinate_indices(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        # Get mask
        mask = self.get_faceon_coordinate_mask(x, y)

        # Get and return
        return np.where(mask)[0]

    # -----------------------------------------------------------------
    # FACEON PIXEL X COORDINATES
    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixel_x_min_value(self):
        return self.faceon_x_min_value + 0.5 * self.faceon_pixelscale_value_x # coordinate for the center of the pixel

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixel_x_max_value(self):
        return self.faceon_pixel_x_min_value + (self.faceon_nx-1) * self.faceon_pixelscale_value_x

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixel_x_coordinates(self):
        return np.linspace(self.faceon_pixel_x_min_value, self.faceon_pixel_x_max_value, num=self.faceon_nx)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixel_x_edges(self):
        return np.linspace(self.faceon_x_min_value, self.faceon_x_max_value, num=self.faceon_nx+1)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixel_x_left_edges(self):
        return self.faceon_pixel_x_edges[:-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixel_x_right_edges(self):
        return self.faceon_pixel_x_edges[1:]

    # -----------------------------------------------------------------
    # FACEON PIXEL Y COORDINATES
    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixel_y_min_value(self):
        return self.faceon_y_min_value + 0.5 * self.faceon_pixelscale_value_y # coordinate for the center of the pixel

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixel_y_max_value(self):
        return self.faceon_pixel_y_min_value + (self.faceon_ny-1) * self.faceon_pixelscale_value_y

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixel_y_coordinates(self):
        return np.linspace(self.faceon_pixel_y_min_value, self.faceon_pixel_y_max_value, num=self.faceon_ny)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixel_y_edges(self):
        return np.linspace(self.faceon_y_min_value, self.faceon_y_max_value, num=self.faceon_ny+1)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixel_y_left_edges(self):
        return self.faceon_pixel_y_edges[:-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_pixel_y_right_edges(self):
        return self.faceon_pixel_y_edges[1:]

    # -----------------------------------------------------------------

    def project_faceon(self, logfreq=100, cell_based=True):

        """
        This function ...
        :param logfreq:
        :param cell_based:
        :return:
        """

        # Cell-based
        if cell_based: self.project_faceon_cells(logfreq=logfreq*100)

        # Cell
        else: self.project_faceon_pixels(logfreq=logfreq)

    # -----------------------------------------------------------------

    def project_faceon_pixels(self, logfreq=100):

        """
        This function ...
        :param logfreq:
        :return:
        """

        # Inform the user
        log.info("Creating the face-on map (pixel-based) ...")

        # Initialize maps
        self.faceon = Frame.initialize_nans(self.faceon_shape, unit=self.unit)
        self.faceon_stddev = Frame.initialize_nans(self.faceon_shape, unit=self.unit)
        self.faceon_ncells = Frame.initialize_nans(self.faceon_shape)

        # Set the pixelscale and the coordinate info
        self.faceon.pixelscale = self.faceon_pixelscale
        self.faceon.distance = self.faceon_distance
        self.faceon.metadata.update(self.faceon_metadata)

        # Show progress bar
        with Bar(label='', expected_size=self.faceon_npixels, every=1, add_datetime=True) as bar:

            # Loop over the pixels of the map
            index = 0
            for i, x, j, y in sequences.iterate_enumerated_combinations(self.faceon_pixel_x_coordinates, self.faceon_pixel_y_coordinates):

                # Debugging
                if index % logfreq == 0:
                    log.debug("Calculating projected value in pixel " + str(index) + " of " + str(self.faceon_npixels) + " (" + tostr(float(index) / self.faceon_npixels * 100, decimal_places=1, round=True) + "%) ...")
                    log.debug("Pixel position: x = " + tostr(x) + " " + tostr(self.length_unit) + ", y = " + tostr(y) + " " + tostr(self.length_unit))
                    log.debug("Pixel index: (" + str(i) + ", " + str(j) + ")")

                # Show progress
                #progress = int(float(index+1) / float(self.faceon_npixels))
                bar.show(float(index+1))

                # Get the indices
                indices = self.get_faceon_coordinate_indices(x, y)
                nindices = indices.shape[0]

                # Set number of cells
                self.faceon_ncells[j, i] = nindices

                # If any cells
                if nindices > 0:

                    # Calculate the heating fraction
                    vals = self.values_faceon[indices]

                    # With weights
                    if self.has_weights:

                        wghts = self.weights_faceon[indices]
                        fraction = numbers.weighed_arithmetic_mean_numpy(vals, weights=wghts)
                        fraction_stddev = numbers.weighed_standard_deviation_numpy(vals, weights=wghts, mean=fraction)

                    # Without weights
                    else:

                        # Calculate the mean heating fraction
                        fraction = numbers.arithmetic_mean_numpy(vals)
                        fraction_stddev = numbers.standard_deviation_numpy(vals)

                    # Set fraction
                    self.faceon[j, i] = fraction
                    self.faceon_stddev[j, i] = fraction_stddev

                # Increment the absolute pixel index
                index += 1

    # -----------------------------------------------------------------

    def get_faceon_pixel_index(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        # Check
        mask_x = (self.faceon_pixel_x_left_edges < x) * (x < self.faceon_pixel_x_right_edges)
        mask_y = (self.faceon_pixel_y_left_edges < y) * (y < self.faceon_pixel_y_right_edges)

        # Get indices
        indices_x = np.where(mask_x)[0]
        indices_y = np.where(mask_y)[0]
        nx = len(indices_x)
        ny = len(indices_y)

        # Should be single or None
        if nx == 0 or ny == 0: return None, None

        if nx > 1: raise RuntimeError("Something went wrong")
        if ny > 1: raise RuntimeError("Something went wrong")

        # Return
        return indices_x[0], indices_y[0]

    # -----------------------------------------------------------------

    def project_faceon_cells(self, logfreq=10000):

        """
        This function ...
        :param logfreq:
        :return:
        """

        # Inform the user
        log.info("Creating the face-on map (cell-based) ...")

        # Initialize maps
        self.faceon = Frame.initialize_zeroes(self.faceon_shape, unit=self.unit)
        self.faceon_stddev = Frame.initialize_nans(self.faceon_shape, unit=self.unit)
        self.faceon_ncells = Frame.initialize_zeroes(self.faceon_shape)

        # Set the pixelscale and the coordinate info
        self.faceon.pixelscale = self.faceon_pixelscale
        self.faceon.distance = self.faceon_distance
        self.faceon.metadata.update(self.faceon_metadata)

        # Show progress bar
        with Bar(label='', expected_size=self.nfaceon_coordinates, every=100, add_datetime=True) as bar:

            # Loop over the cells
            for index in range(self.nfaceon_coordinates):

                # Get coordinates
                x = self.x_faceon[index]
                y = self.y_faceon[index]
                value = self.values_faceon[index]
                #if value == 0: continue

                # Show progress
                bar.show(float(index + 1))

                # Get the pixel coordinates
                #for i, j in self.get_faceon_pixel_indices(x, y):
                i, j = self.get_faceon_pixel_index(x, y)
                if i is None: continue

                # Debugging
                if index % logfreq == 0:
                    log.debug("Calculating projected value from cell " + str(index) + " of " + str(self.nfaceon_coordinates) + " (" + tostr(float(index) / self.nfaceon_coordinates * 100, decimal_places=1, round=True) + "%) ...")
                    log.debug("Cell position: x = " + tostr(x) + " " + tostr(self.length_unit) + ", y = " + tostr(y) + " " + tostr(self.length_unit))
                    log.debug("Pixel index: (" + str(i) + ", " + str(j) + ")")

                #  Set fraction
                self.faceon[j, i] += value
                #self.faceon_stddev[j, i] = fraction_stddev
                self.faceon_ncells[j, i] += 1

        # Normalize
        self.faceon /= self.faceon_ncells # make average

    # -----------------------------------------------------------------

    @property
    def edgeon_ny(self):
        return self.projection_edgeon.pixels_x

    # -----------------------------------------------------------------

    @property
    def edgeon_nz(self):
        return self.projection_edgeon.pixels_y

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_npixels(self):
        return self.edgeon_ny * self.edgeon_nz

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_shape(self):
        return PixelShape.from_xy(self.edgeon_ny, self.edgeon_nz)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_y_min(self):
        return self.projection_edgeon.center_x - 0.5 * self.projection_edgeon.field_x

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_y_max(self):
        return self.projection_edgeon.center_x + 0.5 * self.projection_edgeon.field_x

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_z_min(self):
        return self.projection_edgeon.center_y - 0.5 * self.projection_edgeon.field_y

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_z_max(self):
        return self.projection_edgeon.center_y + 0.5 * self.projection_edgeon.field_y

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_y_min_value(self):
        return self.edgeon_y_min.to(self.length_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_y_max_value(self):
        return self.edgeon_y_max.to(self.length_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_z_min_value(self):
        return self.edgeon_z_min.to(self.length_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_z_max_value(self):
        return self.edgeon_z_max.to(self.length_unit).value

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixelscale(self):
        return self.projection_edgeon.physical_pixelscale

    # -----------------------------------------------------------------

    @property
    def edgeon_pixelscale_y(self):
        return self.edgeon_pixelscale.x

    # -----------------------------------------------------------------

    @property
    def edgeon_pixelscale_z(self):
        return self.edgeon_pixelscale.y

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixelscale_value(self):
        return Extent(self.edgeon_pixelscale_y.to(self.length_unit).value, self.edgeon_pixelscale_z.to(self.length_unit).value)

    # -----------------------------------------------------------------

    @property
    def edgeon_pixelscale_value_y(self):
        return self.edgeon_pixelscale_value.x

    # -----------------------------------------------------------------

    @property
    def edgeon_pixelscale_value_z(self):
        return self.edgeon_pixelscale_value.y

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_half_pixelscale(self):
        return self.edgeon_pixelscale * 0.5

    # -----------------------------------------------------------------

    @property
    def edgeon_half_pixelscale_y(self):
        return self.edgeon_half_pixelscale.x

    # -----------------------------------------------------------------

    @property
    def edgeon_half_pixelscale_z(self):
        return self.edgeon_half_pixelscale.y

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_half_pixelscale_value(self):
        return Extent(self.edgeon_half_pixelscale_y.to(self.length_unit).value, self.edgeon_half_pixelscale_z.to(self.length_unit).value)

    # -----------------------------------------------------------------

    @property
    def edgeon_half_pixelscale_value_y(self):
        return self.edgeon_half_pixelscale_value.x

    # -----------------------------------------------------------------

    @property
    def edgeon_half_pixelscale_value_z(self):
        return self.edgeon_half_pixelscale_value.y

    # -----------------------------------------------------------------

    @property
    def edgeon_distance(self):
        return self.projection_edgeon.distance

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_metadata(self):

        """
        Thisf unction ...
        :return:
        """

        meta = OrderedDict()
        meta["y_min"] = repr(self.edgeon_y_min_value) + " " + tostr(self.length_unit)
        meta["y_max"] = repr(self.edgeon_y_max_value) + " " + tostr(self.length_unit)
        meta["z_min"] = repr(self.edgeon_z_min_value) + " " + tostr(self.length_unit)
        meta["z_max"] = repr(self.edgeon_z_max_value) + " " + tostr(self.length_unit)
        return meta

    # -----------------------------------------------------------------

    def get_edgeon_coordinate_ranges(self, y, z):

        """
        This function ...
        :param y:
        :param z:
        :return:
        """

        # Determine range of y and z
        y_range = RealRange(y - self.edgeon_half_pixelscale_value_y, y + self.edgeon_half_pixelscale_value_y)
        z_range = RealRange(z - self.edgeon_half_pixelscale_value_z, z + self.edgeon_half_pixelscale_value_z)

        # Return
        return y_range, z_range

    # -----------------------------------------------------------------

    def get_edgeon_coordinate_mask(self, y, z):

        """
        This function ...
        :param y:
        :param z:
        :return:
        """

        # Get range
        y_range, z_range = self.get_edgeon_coordinate_ranges(y, z)

        # Create masks
        y_mask = (y_range.min < self.y_edgeon) * (self.y_edgeon <= y_range.max)
        z_mask = (z_range.min < self.z_edgeon) * (self.z_edgeon <= z_range.max)

        # Create combined mask
        return y_mask * z_mask

    # -----------------------------------------------------------------

    def get_edgeon_coordinate_indices(self, y, z):

        """
        This function ...
        :param y:
        :param z:
        :return:
        """

        # Get mask
        mask = self.get_edgeon_coordinate_mask(y, z)

        # Get and return
        return np.where(mask)[0]

    # -----------------------------------------------------------------
    # EDGEON PIXEL Y COORDINATES
    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixel_y_min_value(self):
        return self.edgeon_y_min_value + 0.5 * self.edgeon_pixelscale_value_y  # coordinate for the center of the pixel

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixel_y_max_value(self):
        return self.edgeon_pixel_y_min_value + (self.edgeon_ny - 1) * self.edgeon_pixelscale_value_y

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixel_y_coordinates(self):
        return np.linspace(self.edgeon_pixel_y_min_value, self.edgeon_pixel_y_max_value, num=self.edgeon_ny)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixel_y_edges(self):
        return np.linspace(self.edgeon_y_min_value, self.edgeon_y_max_value, num=self.edgeon_ny+1)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixel_y_left_edges(self):
        return self.edgeon_pixel_y_edges[:-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixel_y_right_edges(self):
        return self.edgeon_pixel_y_edges[1:]

    # -----------------------------------------------------------------
    # EDGEON PIXEL Z COORDINATES
    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixel_z_min_value(self):
        return self.edgeon_z_min_value + 0.5 * self.edgeon_pixelscale_value_z # coordinate for the center of the pixel

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixel_z_max_value(self):
        return self.edgeon_pixel_z_min_value + (self.edgeon_nz - 1) * self.edgeon_pixelscale_value_z

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixel_z_coordinates(self):
        return np.linspace(self.edgeon_pixel_z_min_value, self.edgeon_pixel_z_max_value, num=self.edgeon_nz)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixel_z_edges(self):
        return np.linspace(self.edgeon_z_min_value, self.edgeon_z_max_value, num=self.edgeon_nz+1)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixel_z_left_edges(self):
        return self.edgeon_pixel_z_edges[:-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_pixel_z_right_edges(self):
        return self.edgeon_pixel_z_edges[1:]

    # -----------------------------------------------------------------

    def project_edgeon(self, logfreq=100, cell_based=True):

        """
        This function ...
        :param logfreq:
        :param cell_based:
        :return:
        """

        # Cell-based
        if cell_based: self.project_edgeon_cells(logfreq=logfreq*100)

        # Cell
        else: self.project_edgeon_pixels(logfreq=logfreq)

    # -----------------------------------------------------------------

    def project_edgeon_pixels(self, logfreq=100):

        """
        This function ...
        :param logfreq:
        :return:
        """

        # Inform the user
        log.info("Creating the edge-on map (pixel-based) ...")

        # Initialize maps
        self.edgeon = Frame.initialize_nans(self.edgeon_shape, unit=self.unit)
        self.edgeon_stddev = Frame.initialize_nans(self.edgeon_shape, unit=self.unit)
        self.edgeon_ncells = Frame.initialize_nans(self.edgeon_shape)

        # Set the pixelscale and the coordinate info
        self.edgeon.pixelscale = self.edgeon_pixelscale
        self.edgeon.distance = self.edgeon_distance
        self.edgeon.metadata.update(self.edgeon_metadata)

        # Show progress bar
        with Bar(label='', expected_size=self.edgeon_npixels, every=1, add_datetime=True) as bar:

            # Loop over the pixels of the map
            index = 0
            for i, y, j, z in sequences.iterate_enumerated_combinations(self.edgeon_pixel_y_coordinates, self.edgeon_pixel_z_coordinates):

                # Debugging
                if index % logfreq == 0:
                    log.debug("Calculating projected value in pixel " + str(index) + " of " + str(self.edgeon_npixels) + " (" + tostr(float(index) / self.edgeon_npixels * 100, decimal_places=1, round=True) + "%) ...")
                    log.debug("Pixel position: y = " + tostr(y) + " " + tostr(self.length_unit) + ", z = " + tostr(z) + " " + tostr(self.length_unit))
                    log.debug("Pixel index: (" + str(i) + ", " + str(j) + ")")

                # Show progress
                bar.show(float(index+1))

                # Get the indices
                indices = self.get_edgeon_coordinate_indices(y, z)
                nindices = indices.shape[0]

                # Set number of cells
                self.edgeon_ncells[j, i] = nindices

                # If any cells
                if nindices > 0:

                    # Calculate the heating fraction
                    vals = self.values_edgeon[indices]

                    # With weights
                    if self.has_weights:

                        wghts = self.weights_edgeon[indices]
                        fraction = numbers.weighed_arithmetic_mean_numpy(vals, weights=wghts)
                        fraction_stddev = numbers.weighed_standard_deviation_numpy(vals, weights=wghts, mean=fraction)

                    # Without weights
                    else:

                        fraction = numbers.arithmetic_mean_numpy(vals)
                        fraction_stddev = numbers.standard_deviation_numpy(vals)

                    # Set fraction
                    self.edgeon[j, i] = fraction
                    self.edgeon_stddev[j, i] = fraction_stddev

                # Increment the absolute pixel index
                index += 1

    # -----------------------------------------------------------------

    def get_edgeon_pixel_index(self, y, z):

        """
        This function ...
        :param y:
        :param z:
        :return:
        """

        # Check
        mask_y = (self.edgeon_pixel_y_left_edges < y) * (y < self.edgeon_pixel_y_right_edges)
        mask_z = (self.edgeon_pixel_z_left_edges < z) * (z < self.edgeon_pixel_z_right_edges)

        # Get indices
        indices_y = np.where(mask_y)[0]
        indices_z = np.where(mask_z)[0]
        ny = len(indices_y)
        nz = len(indices_z)

        # Should be single or None
        if ny == 0 or nz == 0: return None, None

        if ny > 1: raise RuntimeError("Something went wrong")
        if nz > 1: raise RuntimeError("Something went wrong")

        # Return
        return indices_y[0], indices_z[0]

    # -----------------------------------------------------------------

    def project_edgeon_cells(self, logfreq=10000):

        """
        This function ...
        :param logfreq:
        :return:
        """

        # Inform the user
        log.info("Creating the edge-on map (cell-based) ...")

        # Initialize maps
        self.edgeon = Frame.initialize_zeroes(self.edgeon_shape, unit=self.unit)
        self.edgeon_stddev = Frame.initialize_nans(self.edgeon_shape, unit=self.unit)
        self.edgeon_ncells = Frame.initialize_zeroes(self.edgeon_shape)

        # Set the pixelscale and the coordinate info
        self.edgeon.pixelscale = self.edgeon_pixelscale
        self.edgeon.distance = self.edgeon_distance
        self.edgeon.metadata.update(self.edgeon_metadata)

        # Show progress bar
        with Bar(label='', expected_size=self.nedgeon_coordinates, every=100, add_datetime=True) as bar:

            # Loop over the cells
            for index in range(self.nedgeon_coordinates):

                # Get coordinates
                y = self.y_edgeon[index]
                z = self.z_edgeon[index]
                value = self.values_edgeon[index]
                # if value == 0: continue

                # Show progress
                bar.show(float(index + 1))

                # Get the pixel coordinates
                i, j = self.get_edgeon_pixel_index(y, z)
                if i is None: continue

                # Debugging
                if index % logfreq == 0:
                    log.debug("Calculating projected value from cell " + str(index) + " of " + str(self.nedgeon_coordinates) + " (" + tostr(float(index) / self.nedgeon_coordinates * 100, decimal_places=1, round=True) + "%) ...")
                    log.debug("Cell position: y = " + tostr(y) + " " + tostr(self.length_unit) + ", z = " + tostr(z) + " " + tostr(self.length_unit))
                    log.debug("Pixel index: (" + str(i) + ", " + str(j) + ")")

                #  Set fraction
                self.edgeon[j, i] += value
                # self.faceon_stddev[j, i] = fraction_stddev
                self.edgeon_ncells[j, i] += 1

        # Normalize
        self.edgeon /= self.edgeon_ncells  # make average

# -----------------------------------------------------------------
