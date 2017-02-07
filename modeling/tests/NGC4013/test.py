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
from pts.core.tools import network
from pts.core.tools.logging import log
from pts.do.commandline import Command

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "Fitting the galaxy NGC4013 using a flattened Sersic profile for the central bulge and a double exponential for the stellar disk and dust disk"

# -----------------------------------------------------------------

# Initialize list for the commands
commands = []

# -----------------------------------------------------------------
# SETUP FUNCTION
# -----------------------------------------------------------------

input_url = "http://www.skirt.ugent.be/downloads/tutorial_NGC4013.tar.gz"
all_url = "http://www.skirt.ugent.be/downloads/tutorial_NGC4013_complete.tar.gz"

# -----------------------------------------------------------------

def setup(temp_path):

    """
    This function ...
    """

    # Inform the user
    log.info("Downloading the input ...")

    # Download the input
    network.download_and_decompress_file(input_url, temp_path, progress_bar=True)

# -----------------------------------------------------------------
# SETUP MODELLING
# -----------------------------------------------------------------

# Settings
settings_setup = dict()
settings_setup["type"] = "images"
settings_setup["name"] = "NGC4013"
settings_setup["fitting_host_ids"] = None

# -----------------------------------------------------------------

# Input

# Create object config
object_config = dict()
ski_path = fs.join(this_dir_path, "NGC4013.ski")
object_config["ski"] = ski_path

# Create input dict for setup
input_setup = dict()
input_setup["object_config"] = object_config
#input_setup["sed"] = sed

# Construct the command
setup_command = Command("setup", "setup the modelling", settings_setup, input_setup, ".")

# Add the command
commands.append(setup_command)

# -----------------------------------------------------------------
# PERFORM THE MODELLING
# -----------------------------------------------------------------

# Settings for 'model_sed'
settings_model = dict()
settings_model["ngenerations"] = "4"

# -----------------------------------------------------------------

# Create input dict for model
input_model = dict()
#input_model["parameters_config"] = Configuration(free_parameters=free_parameter_names)
#input_model["descriptions_config"] = Configuration(descriptions=descriptions)
#input_model["types_config"] = Configuration(types=types)
#input_model["units_config"] = Configuration(units=units)
#input_model["ranges_config"] = Configuration(luminosity_range=luminosity_range, dustmass_range=dustmass_range, grainsize_range=grainsize_range, fsil_range=fsil_range)
#input_model["filters_config"] = Configuration(filters=filter_names)

# Fitting initializer config
#input_model["initialize_config"] = Configuration(npackages=1e4)

# Construct the command
model_command = Command("model", "perform the modelling", settings_model, input_model, cwd="./NGC4013")

# Add the command
commands.append(model_command)

# -----------------------------------------------------------------
# TEST FUNCTION
# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    :param temp_path:
    """

    return

# -----------------------------------------------------------------
