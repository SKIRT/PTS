#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.dustpedia.get_poisson_errors Calculate poisson error maps for DustPedia UV and optical images.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

import sys

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import Configuration
from pts.dustpedia.core.dataprocessing import DustPediaDataProcessing
from pts.core.basics.remote import Remote

# -----------------------------------------------------------------

# Configuration
config = Configuration("get_poisson_errors", "Calculate poisson error maps for DustPedia UV and optical images")

# Galaxy name
config.add_required("galaxy_name", str, "the name of the galaxy")
config.add_required("band", str, "the band (GALEX or SDSS u/g/r/i/z)")

# Optional
config.add_optional("remote", str, "the remote host name", None)

# Read the command line options
config.read()

# Get settings
settings = config.get_settings()

# Get command line arguments
arguments = config.get_arguments()

# -----------------------------------------------------------------

# If remote
if settings.remote is not None:

    # Create the remote
    remote = Remote()

    # Setup the remote
    remote.setup(settings.remote)

    # Remove specification of the remote
    arguments_no_remote = arguments[:-2]

    # Send the appropriate command
    remote.launch_pts_command("get_poisson_errors", arguments_no_remote)

else:

    # Create the DustPedia data processing instance
    dpdp = DustPediaDataProcessing()

    # GALEX
    if "GALEX" in settings.band: pass

    # SDSS
    elif "SDSS" in settings.band: pass

    # Invalid option
    else: raise ValueError("Invalid option for 'band'")

# -----------------------------------------------------------------
