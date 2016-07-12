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
from pts.core.basics.configuration import Configuration
from pts.core.basics.errors import ConfigurationError

# -----------------------------------------------------------------

# Create configuration instance
config = Configuration("update")

# Add required
config.add_required("skirt_or_pts", str, "choose to update SKIRT or PTS", to_instance=False, choices=["skirt", "pts"])

# Add optional
config.add_optional("remote", str, "update SKIRT on a remote system")

# Read
config.read()

# -----------------------------------------------------------------

# SKIRT
if config.arguments.skirt_or_pts == "skirt":

    # Create a SkirtUpdater instance
    updater = SkirtUpdater(config)

# PTS
elif config.arguments.skirt_or_pts == "pts":

    # Create a PTSUpdater instance
    updater = PTSUpdater(config)

# Invalid
else: raise ConfigurationError("Invalid option")

# Run the updater
updater.run()

# -----------------------------------------------------------------
