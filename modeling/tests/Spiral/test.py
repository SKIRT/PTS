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

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.unit import parse_unit as u
from pts.core.basics.configuration import Configuration
from pts.core.simulation.skifile import LabeledSkiFile
from pts.core.basics.range import QuantityRange, RealRange
from pts.core.basics.map import Map
from pts.do.commandline import Command

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "determining parameters based on mock observations of a simple spiral galaxy model"

# -----------------------------------------------------------------

# Initialize list for the commands
commands = []

# -----------------------------------------------------------------
# SETUP FUNCTION
# -----------------------------------------------------------------

def setup(temp_path):

    """
    This function ...
    :param temp_path:
    """

    return

# -----------------------------------------------------------------
# LAUNCH REFERENCE SIMULATION
# -----------------------------------------------------------------

# Determine the ski path
ski_path = fs.join(this_dir_path, "spiral.ski")

# Determine the simulation output path
simulation_output_path = "ref"

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

# Create simulation
#prefix = name = "spiral"
output_path = simulation_output_path
#simulation = SkirtSimulation(prefix, outpath=output_path, ski_path=ski_path, name=name)

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

# Settings
settings_setup = dict()
settings_setup["type"] = "sed"
settings_setup["name"] = "Spiral"
settings_setup["fitting_host_ids"] = None

# Create object config
object_config = dict()
object_config["ski"] = ski_path

# Create input dict for setup
input_setup = dict()
input_setup["object_config"] = object_config
input_setup["sed"] = mock_sed_path

# Construct the command
stp = Command("setup", "setup the modeling", settings_setup, input_setup, cwd=".")

# Add the command
commands.append(stp)

# -----------------------------------------------------------------
# PEFORM THE MODELLING
# -----------------------------------------------------------------

# Settings
settings_model = dict()
settings_model["ngenerations"] = 4
settings_model["nsimulations"] = 20
settings_model["fitting_settings"] = {"spectral_convolution": False}

# Input

# Get free parameter names
ski = LabeledSkiFile(ski_path)
free_parameter_names = ski.labels

# Get fitting filter names
#filter_names = sed.filter_names()

# Set descriptions
descriptions = Map()
descriptions["luminosity"] = "total luminosity of the SN"
descriptions["dustmass"] = "total dust mass"
descriptions["grainsize"] = "dust grain size"
descriptions["fsil"] = "dust silicate fraction"

# Set types
types = Map()
types["luminosity"] = "luminosity"
types["dustmas"] = "mass"
types["grainsize"] = "grainsize"
types["fsil"] = "dimless"

# Set units
units = Map()
units["luminosity"] = u("Lsun")
units["dustmass"] = u("Msun")
units["grainsize"] = u("micron")
units["fsil"] = None

# Set ranges
luminosity_range = QuantityRange(100, 1000, "Lsun")
dustmass_range = QuantityRange(0.3, 5, "Msun")
grainsize_range = QuantityRange(0.1, 5, "micron")
fsil_range = RealRange(0.1, 100)

# Create input dict for model
input_model = dict()
input_model["parameters_config"] = Configuration(free_parameters=free_parameter_names)
input_model["descriptions_config"] = Configuration(descriptions=descriptions)
input_model["types_config"] = Configuration(types=types)
input_model["units_config"] = Configuration(units=units)
input_model["ranges_config"] = Configuration(luminosity_range=luminosity_range, dustmass_range=dustmass_range, grainsize_range=grainsize_range, fsil_range=fsil_range)
#input_model["filters_config"] = Configuration(filters=filter_names)

# Fitting initializer config
input_model["initialize_config"] = Configuration(npackages=1e4)

# Add dict of input for 'model' command to the list
#input_dicts.append(input_model)

# Construct the command
command = Command("model", "perform the modelling", settings_model, input_model, "./Spiral")

# Add the command
commands.append(command)

# -----------------------------------------------------------------
# TEST FUNCTION
# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    """

    return

# -----------------------------------------------------------------
