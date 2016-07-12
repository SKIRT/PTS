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

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import Configuration
from pts.dustpedia.core.dataprocessing import DustPediaDataProcessing
from pts.core.basics.remote import Remote
from pts.core.tools import filesystem as fs
from pts.core.tools import logging
from pts.core.tools import time

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

# Determine the log file path
logfile_path = fs.join(fs.cwd(), time.unique_name("log") + ".txt") if config.arguments.report else None

# Determine the log level
level = "DEBUG" if config.arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting get_poisson_errors ...")

# -----------------------------------------------------------------

temp_name = time.unique_name(settings.band.replace(" ", ""))

# -----------------------------------------------------------------

# Remotely
if settings.remote is not None:

    # Create the remote
    remote = Remote()

    # Setup the remote
    remote.setup(settings.remote)

    # Remove specification of the remote
    arguments_no_remote = arguments[:-2]

    # Make a new temporary directory remotely
    remote_path = fs.join(remote.home_directory, temp_name)
    remote.create_directory(remote_path)
    remote.change_cwd(remote_path)

    # Send the appropriate command
    remote.launch_pts_command("get_poisson_errors", arguments_no_remote)

    # Retrieve the remote directory
    remote.download(remote_path, fs.cwd(), show_output=True)
    local_path = fs.join(fs.cwd(), fs.name(remote_path))

    # Remove the temporary remote directory
    remote.remove_directory(remote_path)

# Locally
else:

    # Make a local directory
    local_path = fs.join(fs.cwd(), temp_name)
    fs.create_directory(local_path)

    # Create the DustPedia data processing instance
    dpdp = DustPediaDataProcessing()

    # GALEX
    if "GALEX" in settings.band:

        dpdp.make_galex_mosaic_and_poisson_frame(settings.galaxy_name, local_path)

    # SDSS
    elif "SDSS" in settings.band:

        # Make ...
        dpdp.make_sdss_mosaic_and_poisson_frame(settings.galaxy_name, settings.band, local_path)

    # Invalid option
    else: raise ValueError("Invalid option for 'band'")

# Inform the user
log.info("Results are placed in '" + local_path)

# -----------------------------------------------------------------
