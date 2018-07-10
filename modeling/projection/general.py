#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.projection.genereal Contains general projection functions.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.simulation.simulation import createsimulations
from ...core.tools.utils import lazyproperty
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...magic.core.frame import Frame
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.simulation.execute import run_simulation
from ...core.prep.smile import get_oligochromatic_template
from ...core.tools import introspection
from ..basics.projection import GalaxyProjection, FaceOnProjection, EdgeOnProjection
from ..basics.instruments import FrameInstrument, FullInstrument
from ..build.representations.galaxy import create_projection_from_deprojection, create_faceon_projection_from_earth_projection, create_edgeon_projection_from_earth_projection
from ..basics.models import DeprojectionModel3D
from ...magic.basics.coordinatesystem import CoordinateSystem
from ..basics.instruments import Instrument
from ...core.tools import numbers

# -----------------------------------------------------------------

# Number of photon packages
default_npackages = 5e7

# Instruments/orientations
earth_name = "earth"
faceon_name = "faceon"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

default_scale_heights = 15.

# -----------------------------------------------------------------

def get_faceon_projection(argument, distance=None, center=None, radial_factor=1, strict=False):

    """
    This function ...
    :param argument:
    :param distance:
    :param center:
    :param radial_factor:
    :param strict:
    :return:
    """

    # Already a projection
    #if isinstance(argument, GalaxyProjection) and not isinstance(argument, FaceOnProjection): raise ValueError("Not a face-on projection")

    # Is already face-on?
    if is_faceon(argument):
        if is_projection(argument): return FaceOnProjection.from_projection(argument)
        elif is_instrument(argument): return FaceOnProjection.from_instrument(argument)
        else: raise RuntimeError("Something went wrong")
    elif strict: raise ValueError("Not face-on")

    # Create projection
    projection = get_projection(argument, distance=distance, center=center)

    # Convert into face-on
    return create_faceon_projection_from_earth_projection(projection, radial_factor=radial_factor)

# -----------------------------------------------------------------

def get_edgeon_projection(argument, distance=None, center=None, scaleheight=None, radial_factor=1,
                          scale_heights=default_scale_heights, strict=False):

    """
    This function ...
    :param argument:
    :param distance:
    :param center:
    :param scaleheight:
    :param radial_factor:
    :param scale_heights:
    :param strict:
    :return:
    """

    # Is already edge-on?
    if is_edgeon(argument):
        if is_projection(argument): return EdgeOnProjection.from_projection(argument)
        elif is_instrument(argument): return EdgeOnProjection.from_instrument(argument)
        else: raise RuntimeError("Something went wrong")
    elif strict: raise ValueError("Not edge-on")

    # Create projection
    projection = get_projection(argument, distance=distance, center=center)

    # Convert into edge-on
    if scaleheight is None: raise ValueError("Scaleheight must be passed for conversion into edge-on projection")
    z_extent = 2. * scaleheight * scale_heights
    return create_edgeon_projection_from_earth_projection(projection, z_extent, radial_factor=radial_factor)

# -----------------------------------------------------------------

def get_projection(argument, distance=None, center=None, inclination=None, position_angle=None, azimuth=0.0):

    """
    This function returns a galaxy projection object from various kinds of input
    :param argument:
    :param distance:
    :param center:
    :param inclination:
    :param position_angle:
    :param azimuth:
    :return:
    """

    # Already a projection
    if is_projection(argument): return argument

    # Deprojection
    elif is_deprojection(argument): return GalaxyProjection.from_deprojection(argument, distance=distance, azimuth=azimuth)

    # Coordinate system
    elif is_coordinate_system(argument): return GalaxyProjection.from_wcs(argument, center, distance, inclination, azimuth, position_angle)

    # Instrument
    elif is_instrument(argument): return GalaxyProjection.from_instrument(argument)

    # Invalid
    else: raise ValueError("Invalid argument of type '" + str(type(argument)) + "'")

# -----------------------------------------------------------------

def is_projection(argument):
    return isinstance(argument, GalaxyProjection)

# -----------------------------------------------------------------

def is_faceon_projection(argument):
    return isinstance(argument, FaceOnProjection)

# -----------------------------------------------------------------

def is_edgeon_projection(argument):
    return isinstance(argument, EdgeOnProjection)

# -----------------------------------------------------------------

def is_instrument(argument):
    return isinstance(argument, Instrument)

# -----------------------------------------------------------------

def is_faceon(argument):
    if is_projection(argument):
        if is_faceon_projection(argument): return True
        else: return has_faceon_angles(argument)
    elif is_instrument(argument): return has_faceon_angles(argument)
    else: return False

# -----------------------------------------------------------------

def is_edgeon(argument):
    if is_projection(argument):
        if is_edgeon_projection(argument): return True
        else: return has_edgeon_angles(argument)
    elif is_instrument(argument): return has_edgeon_angles(argument)
    else: return False

# -----------------------------------------------------------------

def is_coordinate_system(argument):
    return isinstance(argument, CoordinateSystem)

# -----------------------------------------------------------------

def is_deprojection(argument):
    return isinstance(argument, DeprojectionModel3D)

# -----------------------------------------------------------------

def has_faceon_angles(argument):
    return numbers.is_close_to_zero(argument.inclination.to("deg").value) and numbers.is_close_to_zero(argument.azimuth.to("deg").value) and numbers.is_close(argument.position_angle.to("deg").value, 90)

# -----------------------------------------------------------------

def has_edgeon_angles(argument):
    return numbers.is_close(argument.inclination.to("deg").value, 90) and numbers.is_close_to_zero(argument.azimuth.to("deg").value) and numbers.is_close_to_zero(argument.position_angle.to("deg").value)

# -----------------------------------------------------------------
