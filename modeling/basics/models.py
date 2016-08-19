#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.basics.models Contains the SersicModel, ExponentialDiskModel and DeprojectionModel classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
from abc import ABCMeta

# Import astronomical modules
from astropy.coordinates import Angle
from astropy.units import Unit
from astropy.modeling.models import Sersic2D

# Import the relevant PTS classes and modules
from ...core.basics.composite import SimplePropertyComposite

# -----------------------------------------------------------------

class Model(SimplePropertyComposite):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

# -----------------------------------------------------------------
#
# 3D MODELS
#
# -----------------------------------------------------------------

def load_3d_model(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Get the first line of the file
    with open(path, 'r') as f: first_line = f.readline()

    # Create the appropriate model
    if "SersicModel3D" in first_line: return SersicModel3D.from_file(path)
    elif "ExponentialDiskModel3D" in first_line: return ExponentialDiskModel3D.from_file(path)
    elif "DeprojectionModel3D" in first_line: return DeprojectionModel3D.from_file(path)
    else: raise ValueError("Unrecognized model file")

# -----------------------------------------------------------------

class SersicModel3D(Model):

    """
    This function ...
    """

    def __init__(self, effective_radius, index, flattening, tilt=0):

        """
        This function ...
        :return:
        """

        self.effective_radius = effective_radius
        self.index = index
        self.flattening = flattening
        self.tilt = tilt

    # -----------------------------------------------------------------

    @classmethod
    def from_2d(cls, sersic2d, inclination, position_angle):

        """
        :param sersic2d:
        :param inclination:
        :param position_angle:
        :return:
        """

        # Get effective radius and Sersic index
        effective_radius = sersic2d.effective_radius
        index = sersic2d.index

        # Calculate the intrinsic flattening
        flattening = intrinsic_flattening(sersic2d.axial_ratio, inclination)

        # Calculate the tilt angle of the bulge (tilt w.r.t. the x-axis)
        tilt = deproject_pa_to_tilt(sersic2d.position_angle - position_angle, inclination)
        tilt = Angle(90., "deg") - tilt

        # Create a new Sersic model and return it
        return cls(effective_radius, index, flattening, tilt)

# -----------------------------------------------------------------

class ExponentialDiskModel3D(Model):

    """
    This function ...
    """

    def __init__(self, radial_scale, axial_scale, radial_truncation=0, axial_truncation=0, inner_radius=0, tilt=0):

        """
        This function ...
        :return:
        """

        # Set the properties
        self.radial_scale = radial_scale
        self.axial_scale = axial_scale
        self.radial_truncation = radial_truncation
        self.axial_truncation = axial_truncation
        self.inner_radius = inner_radius
        self.tilt = tilt

    # -----------------------------------------------------------------

    @classmethod
    def from_2d(cls, exponentialdiskmodel2d, inclination, position_angle):

        """
        This function ...
        :param exponentialdiskmodel2d:
        :param inclination:
        :param position_angle:
        :return:
        """

        # Get the radial scale
        radial_scale = exponentialdiskmodel2d.scalelength

        # Calculate the intrinsic flattening
        flattening = intrinsic_flattening(exponentialdiskmodel2d.axial_ratio, inclination)

        # Calculate the axial scale
        axial_scale = flattening * radial_scale

        # Calculate the tilt angle of the disk (tilt w.r.t. the x-axis)
        #tilt = deproject_pa_to_tilt(parameters.PA - position_angle, inclination)
        tilt = deproject_pa_to_tilt(position_angle - exponentialdiskmodel2d.position_angle, inclination)
        tilt = Angle(90., "deg") - tilt

        # Create a new exponential disk model and return it
        return cls(radial_scale, axial_scale, tilt=tilt)

# -----------------------------------------------------------------

class DeprojectionModel3D(Model):

    """
    This class ...
    """

    def __init__(self, filename, pixelscale, position_angle, inclination, x_size, y_size, x_center, y_center, scale_height):

        """
        This function ...
        :return:
        """

        # position angle: -360 deg to 360 deg
        # inclination: 0 deg to 90 deg
        # center in pixel coordinates!

        self.filename = filename
        self.pixelscale = pixelscale
        self.position_angle = position_angle
        self.inclination = inclination
        self.x_size = x_size
        self.y_size = y_size
        self.x_center = x_center
        self.y_center = y_center
        self.scale_height = scale_height

# -----------------------------------------------------------------

def intrinsic_flattening(qprime, inclination):

    """
    This function ...
    :param qprime:
    :param inclination:
    """

    # Get the inclination angle in radians
    i = inclination.to("radian").value

    # Calculate the intrinsic flattening
    q = math.sqrt((qprime**2 - math.cos(i)**2)/math.sin(i)**2)

    # Return the intrinsic flattening
    return q

# -----------------------------------------------------------------

def project_azimuth_to_pa(azimuth, inclination):

    """
    This function ...
    :param azimuth:
    :param inclination:
    """

    # Get the azimuth angle and inclination in radians
    azimuth_radian = azimuth.to("radian").value
    i_radian = inclination.to("radian").value

    denominator = math.sqrt(math.cos(azimuth_radian)**2 * math.cos(i_radian)**2 + math.sin(azimuth_radian)**2)

    cos_pa = math.cos(azimuth_radian) * math.cos(i_radian) / denominator
    sin_pa = math.sin(azimuth_radian) / denominator

    pa_radian = math.atan2(sin_pa, cos_pa) * Unit("radian")

    return pa_radian.to("deg")

# -----------------------------------------------------------------

def deproject_pa_to_azimuth(pa, inclination):

    """
    This function ...
    :param pa:
    :param inclination:
    :return:
    """

    # Get the PA and inclination in radians
    pa_radian = pa.to("radian").value
    i_radian = inclination.to("radian").value

    denominator = math.sqrt(math.cos(pa_radian)**2 + math.sin(pa_radian)**2 * math.cos(i_radian)**2)

    cos_azimuth = math.cos(pa_radian) / denominator
    sin_azimuth = math.sin(pa_radian) * math.cos(inclination) / denominator

    azimuth_radian = math.atan2(sin_azimuth, cos_azimuth) * Unit("radian")

    return azimuth_radian.to("deg")

# -----------------------------------------------------------------

def project_tilt_to_pa(tilt, inclination):

    """
    This function ...
    :param tilt:
    :param inclination:
    :return:
    """

    # Get the tilt angle and inclination in radians
    tilt_radian = tilt.to("radian").value
    i_radian = inclination.to("radian").value

    denominator = math.sqrt(math.sin(tilt_radian)**2 * math.sin(i_radian)**2 + math.cos(tilt_radian)**2)

    cos_pa = math.sin(tilt_radian) * math.sin(i_radian) / denominator
    sin_pa = math.cos(tilt_radian) / denominator

    pa_radian = math.atan2(sin_pa, cos_pa) * Unit("radian")

    return pa_radian.to("deg")

# -----------------------------------------------------------------

def deproject_pa_to_tilt(pa, inclination):

    """
    This function ...
    :param pa:
    :param inclination:
    :return:
    """

    # Get the PA and inclination in radians
    pa_radian = pa.to("radian").value
    i_radian = inclination.to("radian").value

    denominator = math.sqrt(math.sin(pa_radian)**2 * math.sin(i_radian)**2 + math.cos(pa_radian)**2)

    cos_tilt = math.sin(pa_radian) * math.sin(i_radian) / denominator
    sin_tilt = math.cos(pa_radian) / denominator

    tilt_radian = math.atan2(sin_tilt, cos_tilt) * Unit("radian")

    return tilt_radian.to("deg")

# -----------------------------------------------------------------

# Test the deprojection functions:

#test_angles = []
#test_angles.append(Angle(33, "deg"))
#test_angles.append(Angle(45, "deg"))
#test_angles.append(Angle(0, "deg"))
#test_angles.append(Angle(90, "deg"))
#test_angles.append(Angle(189, "deg"))

#for test_angle in test_angles:

#    result = deproject_pa_to_tilt(test_angle, Angle(0.0, "deg"))
#    print("Should fail/be zero:", Angle(90., "deg") - result, test_angle)

#    result = deproject_pa_to_tilt(test_angle, Angle(90., "deg"))
#    print("Should be the same:", Angle(90., "deg") - result, test_angle)

    #result = project_tilt_to_pa(test_angle, Angle(0.0, "deg"))
    #print("Should fail/be zero:", Angle(90., "deg") - result, test_angle)

    #result = project_tilt_to_pa(test_angle, Angle(90., "deg"))
    #print("Should be the same:", Angle(90., "deg") - result, test_angle)

# -----------------------------------------------------------------
#
# 2D MODELS
#
# -----------------------------------------------------------------

def load_2d_model(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Get the first line of the file
    with open(path, 'r') as f: first_line = f.readline()

    # Create the appropriate model
    if "SersicModel2D" in first_line: return SersicModel2D.from_file(path)
    elif "ExponentialDiskModel2D" in first_line: return ExponentialDiskModel2D.from_file(path)
    else: raise ValueError("Unrecognized model file")

# -----------------------------------------------------------------

class SersicModel2D(Model):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        self.rel_contribution = kwargs.pop("rel_contribution", None)
        self.fluxdensity = kwargs.pop("fluxdensity", None)
        self.axial_ratio = kwargs.pop("axial_ratio", None)
        self.position_angle = kwargs.pop("position_angle", None)  # (degrees ccw from North)
        self.effective_radius = kwargs.pop("effective_radius", None)
        self.index = kwargs.pop("index", None)

# -----------------------------------------------------------------

class ExponentialDiskModel2D(Model):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        """

        self.rel_contribution = kwargs.pop("relative_contribution", None)
        self.fluxdensity = kwargs.pop("fluxdensity", None)
        self.axial_ratio = kwargs.pop("axial_ratio", None)
        self.position_angle = kwargs.pop("position_angle", None) # (degrees ccw from North)
        self.mu0 = kwargs.pop("mu0", None)
        self.scalelength = kwargs.pop("scalelength", None)

# -----------------------------------------------------------------
