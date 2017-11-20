#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.equivalent_lengths Show equivalent lenghts.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.tools import types
from pts.core.tools.stringify import tostr

# -----------------------------------------------------------------

environment = load_modeling_environment_cwd()

# -----------------------------------------------------------------

reference_names = ["disk_major", "disk_minor", "disk_semimajor", "disk_semiminor", "truncation_major", "truncation_minor", "truncation_semimajor", "truncation_semiminor"]

# -----------------------------------------------------------------

reference_lengths = dict()

reference_lengths["disk_major"] = environment.disk_ellipse.major
reference_lengths["disk_minor"] = environment.disk_ellipse.minor
reference_lengths["disk_semimajor"] = 0.5 * reference_lengths["disk_major"]
reference_lengths["disk_semiminor"] = 0.5 * reference_lengths["disk_minor"]

reference_lengths["truncation_major"] = environment.truncation_ellipse.major
reference_lengths["truncation_minor"] = environment.truncation_ellipse.minor
reference_lengths["truncation_semimajor"] = 0.5 * reference_lengths["truncation_major"]
reference_lengths["truncation_semiminor"] = 0.5 * reference_lengths["truncation_minor"]

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_required("length", "real_or_angle_or_length_quantity", "length (relative, pixels, angle or length quantity)")
definition.add_positional_optional("reference", "string", "reference length", choices=reference_names)
definition.add_optional("image", "string", "image name")
definition.add_optional("pixelscale", "pixelscale", "image pixelscale")
config = parse_arguments("equivalent_lengths", definition, "Show equivalent lengths")

# -----------------------------------------------------------------

def show_npixels(npixels, distance, pixelscale):

    """
    This function ...
    :param npixels:
    :param distance:
    :param pixelscale:
    :return:
    """

    print("")
    print("Input has been interpreted as a number of pixels: " + str(npixels))
    print("")

    # As angle
    angle = npixels * pixelscale.average
    print(" - angle: " + tostr(angle))
    #print("")

    # As quantity
    quantity = (angle * distance).to("kpc", equivalencies=dimensionless_angles())
    print(" - length quantity: " + tostr(quantity))

    print("")

# -----------------------------------------------------------------

def show_angle(angle, distance, pixelscale=None):

    """
    This function ...
    :param angle:
    :param distance:
    :param pixelscale:
    :return:
    """

    print("")
    print("Input has been interpreted as an angle: " + tostr(angle))
    print("")

    # As quantity
    quantity = (angle * distance).to("kpc", equivalencies=dimensionless_angles())
    print(" - length quantity: " + tostr(quantity))
    #print("")

    # As number of pixels?
    if pixelscale is not None:
        npixels = (angle / pixelscale.average).to("").value
        print(" - npixels: " + str(npixels))

    print("")

# -----------------------------------------------------------------

def show_quantity(quantity, distance, pixelscale=None):

    """
    This function ...
    :param quantity:
    :param distance:
    :param pixelscale:
    :return:
    """

    print("")
    print("Input has been interpreted as a length quantity: " + tostr(quantity))
    print("")

    # As angle
    angle = (quantity / distance).to("arcsec", equivalencies=dimensionless_angles())
    print(" - angle: " + tostr(angle))
    #print("")

    # As number of pixels?
    if pixelscale is not None:
        npixels = (angle / pixelscale.average).to("").value
        print(" - npixels: " + str(npixels))

    print("")

# -----------------------------------------------------------------

# Real: number of pixels or relative length
if types.is_real_type(config.length):

    # Relative to reference
    if config.reference is not None:

        # Get the quantity
        length = config.length * reference_lengths[config.reference]

        # Get the pixelscale
        if config.pixelscale is not None: pixelscale = config.pixelscale

        # Get pixelscale from image
        elif config.image is not None: pixelscale = environment.get_frame(config.image).pixelscale

        # No pixelscale
        else: pixelscale = None

        # Show
        show_angle(length, environment.galaxy_distance, pixelscale=pixelscale)

    # Number of pixels
    else:

        # Get the pixelscale
        if config.pixelscale is not None: pixelscale = config.pixelscale

        # Get pixelscale from image
        elif config.image is not None: pixelscale = environment.get_frame(config.image).pixelscale

        # No pixelscale?
        else: raise ValueError("Image name or pixelscale has to be specified")

        # Show
        show_npixels(config.length, environment.galaxy_distance, pixelscale)

# Angle
elif types.is_angle(config.length):

    # Get the pixelscale
    if config.pixelscale is not None: pixelscale = config.pixelscale

    # Get pixelscale from image
    elif config.image is not None: pixelscale = environment.get_frame(config.image).pixelscale

    # No pixcelscale
    else: pixelscale = None

    # Show
    show_angle(config.length, environment.galaxy_distance, pixelscale=pixelscale)

# Length quantity
elif types.is_length_quantity(config.length):

    # Get the pixelscale
    if config.pixelscale is not None: pixelscale = config.pixelscale

    # Get pixelscale from image
    elif config.image is not None: pixelscale = environment.get_frame(config.image).pixelscale

    # No pixelscale
    else: pixelscale = None

    # Show
    show_quantity(config.length, environment.galaxy_distance, pixelscale=pixelscale)

# Invalid
else: raise ValueError("Invalid length")

# -----------------------------------------------------------------
