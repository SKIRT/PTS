#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.core.poisson Contains the PoissonErrorCalculator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from .dataprocessing import DustPediaDataProcessing

# -----------------------------------------------------------------

temp_name = time.unique_name(config.band.replace(" ", ""))

# -----------------------------------------------------------------

# Remotely
if config.remote is not None:

    # Create the remote
    remote = Remote()

    # Setup the remote
    remote.setup(config.remote)

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
    if "GALEX" in config.band:

        dpdp.make_galex_mosaic_and_poisson_frame(config.galaxy_name, local_path)

    # SDSS
    elif "SDSS" in config.band:

        band = config.band.split(" ")[1]

        # Make ...
        dpdp.make_sdss_mosaic_and_poisson_frame(config.galaxy_name, band, local_path)

    # Invalid option
    else: raise ValueError("Invalid option for 'band'")

# Inform the user
log.info("Results are placed in '" + local_path)

# -----------------------------------------------------------------

class PoissonErrorCalculator(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        """

        super(PoissonErrorCalculator, self).__init__(config)

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        self.setup()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        super(PoissonErrorCalculator, self).setup()

# -----------------------------------------------------------------
