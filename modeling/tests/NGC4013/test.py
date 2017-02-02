#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

description = "fitting the galaxy NGC4013 using a flattened Sersic profile for the central bulge and a double exponential for the stellar disk and dust disk"

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
commands.append("model_images")

# -----------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------

# Settings for 'setup'
settings_setup = dict()
settings_setup["type"] = "other"
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

input_dicts.append(input_setup)
input_dicts.append(input_model)

# -----------------------------------------------------------------
# WORKING DIRECTORIES
# -----------------------------------------------------------------

cwds.append(".")
cwds.append("./NGC4013")

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
