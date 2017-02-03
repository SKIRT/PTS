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
from pts.core.data.sed import ObservedSED
from pts.core.basics.unit import parse_unit as u

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "modeling of the Supernova 1987A with genetic algorithms and 4 free parameters"

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
settings_setup["type"] = "sed"
settings_setup["name"] = "SN1987A"
settings_setup["fitting_host_ids"] = None
settings.append(settings_setup)

# Settings for 'model_sed'
settings_model = dict()
settings_model["ngenerations"] = "4"
settings.append(settings_model)

# -----------------------------------------------------------------
# INPUT DICTS
# -----------------------------------------------------------------

# Construct the observed SED
sed = ObservedSED(photometry_unit="Jy")
sed.add_point("Pacs 70", 0.0455 * u("Jy"), 0.0034 * u("Jy"))
sed.add_point("Pacs 100", 0.0824 * u("Jy"), 0.0045 * u("Jy"))
sed.add_point("Pacs 160", 0.1530 * u("Jy"), 0.0090 * u("Jy"))
sed.add_point("SPIRE 250", 0.1107 * u("Jy"), 0.0252 * u("Jy"))
sed.add_point("SPIRE 350", 0.0693 * u("Jy"), 0.0228 * u("Jy"))
sed.add_point("ALMA 440mu", 0.0500 * u("Jy"), 0.0150 * u("Jy"))
sed.add_point("ALMA 870mu", 0.0050 * u("Jy"), 0.0010 * u("Jy"))

# Create object config
object_config = dict()
object_config["ski"] = fs.join(this_dir_path, "SN1987A.ski")

# Create input dict for setup
input_setup = dict()
input_setup["object_config"] = object_config
input_setup["sed"] = sed
input_dicts.append(input_setup)

# Create input dict for model
input_model = dict()
input_dicts.append(input_model)

# -----------------------------------------------------------------
# WORKING DIRECTORIES
# -----------------------------------------------------------------

cwds.append(".")
cwds.append("./SN1987A")

# -----------------------------------------------------------------
# SETUP FUNCTION
# -----------------------------------------------------------------

def setup():

    """
    This function ...
    """

    return

# -----------------------------------------------------------------
# TEST FUNCTION
# -----------------------------------------------------------------

def test():

    """
    This function ...
    """

    return

# -----------------------------------------------------------------
