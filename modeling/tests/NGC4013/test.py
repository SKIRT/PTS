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

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "Fitting the galaxy NGC4013 using a flattened Sersic profile for the central bulge and a double exponential for the stellar disk and dust disk"

# -----------------------------------------------------------------

# Initialize lists
commands = []
input_dicts = []
settings = []
cwds = []

# -----------------------------------------------------------------
# COMMANDS
# -----------------------------------------------------------------

commands.append("setup")
commands.append("model")

# -----------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------

# Settings for 'setup'
settings_setup = dict()
settings_setup["type"] = "images"
settings_setup["name"] = "NGC4013"
settings_setup["fitting_host_ids"] = None
settings.append(settings_setup)

# Settings for 'model_sed'
settings_model = dict()
settings_model["ngenerations"] = "4"
settings.append(settings_model)

# -----------------------------------------------------------------
# INPUT DICTS
# -----------------------------------------------------------------

# Create object config
object_config = dict()
ski_path = fs.join(this_dir_path, "NGC4013.ski")
object_config["ski"] = ski_path

# Create input dict for setup
input_setup = dict()
input_setup["object_config"] = object_config
#input_setup["sed"] = sed
input_dicts.append(input_setup)

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

input_dicts.append(input_model)

# -----------------------------------------------------------------
# WORKING DIRECTORIES
# -----------------------------------------------------------------

cwds.append(".")
cwds.append("./NGC4013")

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
# TEST FUNCTION
# -----------------------------------------------------------------

def test():

    """
    This function ...
    """

    return

# -----------------------------------------------------------------
