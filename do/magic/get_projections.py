#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.log import log
from pts.core.tools import filesystem as fs
from pts.modeling.basics.models import load_2d_model
from pts.modeling.basics.properties import GalaxyProperties
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.modeling.basics.projection import GalaxyProjection, FaceOnProjection, EdgeOnProjection
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Required settings
definition.add_required("wcs", "file_path", "FITS file with the desired WCS")
definition.add_required("parameters", "directory_path", "path parameters directory")

definition.add_flag("show", "show the projections", True)

# The output directory
definition.add_optional("output", "directory_path", "output directory", letter="o")

# -----------------------------------------------------------------

config = parse_arguments("get_projections", definition, "Create projections for creating a simulation of a galaxy according to the galaxy parameters")

# -----------------------------------------------------------------

# The projections
projections = dict()

# -----------------------------------------------------------------

# Determine path to properties file
properties_path = fs.join(config.parameters, "properties.dat")

# Create the galaxy properties object
properties = GalaxyProperties.from_file(properties_path)

# Load the disk model
disk_path = fs.join(config.parameters, "disk.mod")
disk = load_2d_model(disk_path)

# Load the WCS
wcs = CoordinateSystem.from_file(config.wcs)

# -----------------------------------------------------------------

def create_projections():

    """
    This function ...
    :return:
    """

    # Inform the user
    log.info("Creating the projection systems ...")

    disk_pa = disk.position_angle

    # Create the 'earth' projection system
    azimuth = 0.0
    projections["earth"] = GalaxyProjection.from_wcs(wcs, properties.center, properties.distance, properties.inclination, azimuth, disk_pa)

    # Create the face-on projection system
    projections["faceon"] = FaceOnProjection.from_wcs(wcs, properties.center, properties.distance)

    # Create the edge-on projection system
    projections["edgeon"] = EdgeOnProjection.from_wcs(wcs, properties.center, properties.distance)

# -----------------------------------------------------------------

# def create_instruments():
#
#     """
#     This function ...
#     :return:
#     """
#
#     # Inform the user
#     log.info("Creating the instruments ...")
#
#     # Loop over the projection systems
#     for name in projections:
#
#         # Create the instrument from the projection system
#         instruments[name] = SimpleInstrument.from_projection(self.projections[name])

# -----------------------------------------------------------------

def show():

    """
    This function ...
    :return:
    """

    print(properties)
    #print(components)

# -----------------------------------------------------------------

def write_projections():

    """
    This function ...
    :return:
    """

    # Loop over the components
    for name in projections:

        # Determine the path
        path = fs.join(config.output, name + ".proj")

        # Save the projection
        projections[name].saveto(path)

# -----------------------------------------------------------------

# Create the projections
create_projections()

# Show
if config.show: show()

# Writing
if config.output is not None: write_projections()

# -----------------------------------------------------------------
