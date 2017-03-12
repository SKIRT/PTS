#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import inspect

# Import astronomical modules
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.unit import parse_unit as u
from pts.core.basics.configuration import Configuration
from pts.core.simulation.skifile import LabeledSkiFile
from pts.core.basics.range import QuantityRange, RealRange
from pts.core.basics.map import Map
from pts.do.commandline import Command
from pts.core.basics.quantity import parse_quantity
from pts.modeling.basics.instruments import FullInstrument

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "fitting a mock spiral galaxy"

# -----------------------------------------------------------------

# Determine the ski path
ski_path = fs.join(this_dir_path, "galaxy.ski")

# Get the initial dust mass of the exponential disk with spiral structure
#ski = LabeledSkiFile(ski_path)
#dust_mass = ski.get_labeled_value("exp_dustmass")

# -----------------------------------------------------------------

# Initialize list for the commands
commands = []

# -----------------------------------------------------------------
# SETUP FUNCTION
# -----------------------------------------------------------------

fraction = 0.5
count = 1000
radius = parse_quantity("10 pc")
cutoff = False
kernel_type = "uniform"

# -----------------------------------------------------------------

def setup(temp_path):

    """
    This function ...
    :param temp_path:
    """

    # Create clumpy ski file
    ski = LabeledSkiFile(ski_path)

    # Add clumpiness to all stellar components
    for component_id in ski.get_stellar_component_ids():
        ski.add_stellar_component_clumpiness(component_id, fraction, count, radius, cutoff, kernel_type)

    # Add clumpiness to all dust components
    for component_id in ski.get_dust_component_ids():
        ski.add_dust_component_clumpiness(component_id, fraction, count, radius, cutoff, kernel_type)

    # Add full instrument that writes out photon counts
    ski.remove_all_instruments()

    # Create a full instrument
    distance = parse_quantity("3.63 Mpc")
    inclination = Angle(59, "deg")
    azimuth = Angle(0, "deg")
    position_angle = Angle(67, "deg")
    # <SimpleInstrument azimuth="0.0 deg" centerX="129.164269644 pc" centerY="-181.633825059 pc" distance="3.69199991226 Mpc" fieldOfViewX="56063.3808772 pc" fieldOfViewY="56114.3939627 pc" inclination="60.0 deg" instrumentName="earth" pixelsX="1099" pixelsY="1100" positionAngle="66.3 deg"/>
    field_x =  parse_quantity("55000 pc")
    field_y = parse_quantity("5500 pc")
    pixels_x = 1000
    pixels_y = 1000
    center_x = parse_quantity("130 pc")
    center_y = parse_quantity("-181 pc")
    scattering_levels = 0
    counts = True # write photon counts
    instrument = FullInstrument(distance=distance, inclination=inclination, azimuth=azimuth, position_angle=position_angle,
                                field_x=field_x, field_y=field_y, pixels_x=pixels_x, pixels_y=pixels_y, center_x=center_x,
                                center_y=center_y, scattering_levels=scattering_levels, counts=counts)

    # Add the instrument
    ski.add_instrument("earth", instrument)

    # Save as new ski file
    new_path = fs.join(temp_path, "galaxy_clumpy.ski")
    ski.saveto(new_path)

# -----------------------------------------------------------------
# LAUNCH REFERENCE SIMULATION
# -----------------------------------------------------------------

# Determine the simulation output path
simulation_output_path = "./ref"

# Settings
settings_launch = dict()
settings_launch["ski"] = ski_path
settings_launch["output"] = simulation_output_path
settings_launch["create_output"] = True

# Input
input_launch = dict()

# Construct the command
launch = Command("launch_simulation", "launch the reference simulation", settings_launch, input_launch, cwd=".")

# Add the command
commands.append(launch)

# -----------------------------------------------------------------
# CREATE THE MOCK SED
# -----------------------------------------------------------------

# Settings
settings_sed = dict()
settings_sed["spectral_convolution"] = False

# Input
input_sed = dict()
input_sed["simulation_output_path"] = simulation_output_path
input_sed["output_path"] = "."

# Construct the command
create_sed = Command("observed_fluxes", "create the mock SED", settings_sed, input_sed, cwd=".")

# Add the command
commands.append(create_sed)

# Determine the path to the mock SED
mock_sed_path = "spiral_earth_fluxes.dat"

# -----------------------------------------------------------------
# SETUP THE MODELLING
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# PERFORM THE MODELLING
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# TEST FUNCTION
# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    """

    return

# -----------------------------------------------------------------
