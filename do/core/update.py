#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.update Update SKIRT and/or PTS locally and/or remotely.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.prep.update import SkirtUpdater, PTSUpdater
from pts.core.basics.configuration import ConfigurationDefinition, ConfigurationReader
from pts.core.basics.errors import ConfigurationError

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("skirt_or_pts", str, "choose to update SKIRT or PTS", choices=["skirt", "pts"])

# Add optional
definition.add_optional("remote", str, "update SKIRT on a remote system")

# Get the configuration
reader = ConfigurationReader("update")
config = reader.read(definition)

# -----------------------------------------------------------------

# SKIRT
if config.skirt_or_pts == "skirt":

    # Create a SkirtUpdater instance
    updater = SkirtUpdater(config)

# PTS
elif config.skirt_or_pts == "pts":

    # Create a PTSUpdater instance
    updater = PTSUpdater(config)

# Invalid
else: raise ConfigurationError("Invalid option")

# Run the updater
updater.run()

# -----------------------------------------------------------------
