#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------

description = "SED modeling of the Supernova 1987A with genetic algorithms and 4 free parameters"

# -----------------------------------------------------------------

# Initialize lists
commands = []
settings = []
cwds = []

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

commands.append("setup")
commands.append("model_sed")

# -----------------------------------------------------------------

cwds.append(".")
cwds.append("./SN1987A")

# -----------------------------------------------------------------

def setup():

    """
    This function ...
    """

    return

# -----------------------------------------------------------------

def test():

    """
    This function ...
    """

    return

# -----------------------------------------------------------------
