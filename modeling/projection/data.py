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

# Import the relevant PTS classes and modules
from ...magic.core.frame import Frame
from ...magic.basics.vector import PixelShape
from ...core.tools import numbers
from ...core.basics.range import RealRange

# -----------------------------------------------------------------

# Number of photon packages
default_npackages = 5e7

# Instruments/orientations
earth_name = "earth"
faceon_name = "faceon"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

def project_3d(name, x_coordinates, y_coordinates, z_coordinates, values, projection, length_unit, unit=None,
               return_stddev=False, return_ncells=False, weights=None):

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
    :return:
    """

    # Get the projection
    #projection = get_projection(projection)
    #projection = get_faceon_projection(projection, strict=True)
    if is_faceon(projection): projection = get_faceon_projection(projection, strict=True)
    elif is_edgeon(projection): projection = get_edgeon_projection(projection, strict=True)
    else: raise ValueError("Projection must be face-on or edge-on")

    # Get the shape of the map
    nx, ny = projection.pixels_x, projection.pixels_y
    shape = PixelShape.from_xy(nx, ny)

    # Initialize maps
    frame = Frame.initialize_nans(shape, unit=unit)
    stddev_frame = Frame.initialize_nans(shape, unit=unit)
    ncells_frame = Frame.initialize_nans(shape)

    # Set boundaries
    x_min = projection.center_x - 0.5 * projection.field_x
    x_max = projection.center_x + 0.5 * projection.field_y
    y_min = projection.center_y - 0.5 * projection.field_y
    y_max = projection.center_y + 0.5 * projection.field_y
    x_min = x_min.to(length_unit)
    x_max = x_max.to(length_unit)
    y_min = y_min.to(length_unit)
    y_max = y_max.to(length_unit)

    # Get pixelscale
    pixelscale = projection.physical_pixelscale
    x_pixelscale = pixelscale.x.to(length_unit).value
    y_pixelscale = pixelscale.y.to(length_unit).value
    half_x_pixelscale = 0.5 * x_pixelscale
    half_y_pixelscale = 0.5 * y_pixelscale

    # Set the pixelscale and the coordinate info
    frame.pixelscale = pixelscale
    frame.distance = projection.distance
    frame.set_meta("x_min", repr(x_min.value) + " " + length_unit)
    frame.set_meta("x_max", repr(x_max.value) + " " + length_unit)
    frame.set_meta("y_min", repr(y_min.value) + " " + length_unit)
    frame.set_meta("y_max", repr(y_max.value) + " " + length_unit)

    # Mask invalid values
    values_nans = np.isnan(values)
    values_infs = np.isinf(values)
    values_mask = values_nans + values_infs #+ values_unphysical
    valid_x_coordinates = np.ma.MaskedArray(x_coordinates, mask=values_mask).compressed()
    valid_y_coordinates = np.ma.MaskedArray(y_coordinates, mask=values_mask).compressed()
    valid_values = np.ma.MaskedArray(values, mask=values_mask).compressed()
    valid_weights = np.ma.MaskedArray(weights, mask=values_mask).compressed() if weights is not None else None

    # Loop over the pixels
    x = x_min
    y = y_min
    for i in range(nx):
        for j in range(ny):

            # Determine range of x and y
            x_range = RealRange(x - half_x_pixelscale, x + half_x_pixelscale)
            y_range = RealRange(y - half_y_pixelscale, y + half_y_pixelscale)

            # Show
            #if index % 100 == 0: log.debug("Calculating heating fraction in the pixel " + str(index) + " of " + str(self.map_npixels) + " (" + tostr(float(index) / self.map_npixels * 100, decimal_places=1, round=True) + "%) ...")

            x_mask = (x_range.min < valid_x_coordinates) * (valid_x_coordinates <= x_range.max)
            y_mask = (y_range.min < valid_y_coordinates) * (valid_y_coordinates <= y_range.max)

            #x_mask = self.get_coordinate_mask_x_for_map(x_range)
            #y_mask = self.get_coordinate_mask_y_for_map(y_range)
            #return x_mask * y_mask
            mask = x_mask * y_mask

            # Get the indices
            #indices = self.get_coordinate_indices_in_column_for_map(x_range, y_range)
            #mask = get_coordinate_mask_for_map(x_range, y_range)
            indices = np.where(mask)[0]
            nindices = indices.shape[0]

            # Set number of cells
            ncells_frame[j, i] = nindices

            # If any cells
            if nindices > 0:

                # Calculate the heating fraction
                #fractions = self.valid_heating_fractions[indices]
                #weights = self.valid_cell_weights[indices]
                vals = valid_values[indices]
                wghts = valid_weights[indices] if valid_weights is not None else None

                # Calculate the mean heating fraction
                fraction = numbers.weighed_arithmetic_mean_numpy(vals, weights=wghts)
                fraction_stddev = numbers.weighed_standard_deviation_numpy(vals, weights=wghts, mean=fraction)

                # Set fraction
                frame[j, i] = fraction
                stddev_frame[j, i] = fraction_stddev

            # Increment the x and y coordinate
            x += x_pixelscale
            y += y_pixelscale

    # Return
    if return_stddev:
        if return_ncells: return frame, stddev_frame, ncells_frame
        else: return frame, stddev_frame
    elif return_ncells: return frame, ncells_frame
    else: return frame

# -----------------------------------------------------------------
