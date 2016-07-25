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
from ...core.tools import filesystem as fs
from ...core.tools import time

# -----------------------------------------------------------------

class PoissonErrorCalculator(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(PoissonErrorCalculator, self).__init__(config)

        # The DustPedia data processing instance
        self.dpdp = None

        # The path to the temporary directory
        self.path = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # Call the setup function
        self.setup()

        # Create the temporary directory
        self.create_directory()

        # GALEX
        if "GALEX" in self.config.band: self.get_galex()

        # SDSS
        elif "SDSS" in self.config.band: self.get_sdss()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(PoissonErrorCalculator, self).setup()

        # Create the DustPedia data processing instance
        self.dpdp = DustPediaDataProcessing()

    # -----------------------------------------------------------------

    def create_directory(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the output directory ...")

        temp_name = time.unique_name(self.config.band.replace(" ", ""))

        # Make a local directory
        self.path = fs.join(fs.cwd(), temp_name)
        fs.create_directory(self.path)

    # -----------------------------------------------------------------

    def get_galex(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing GALEX data ...")

        self.dpdp.make_galex_mosaic_and_poisson_frame(self.config.galaxy_name, self.path)

    # -----------------------------------------------------------------

    def get_sdss(self):

        """
        This function ...
        :return:
        """

        # Determine which SDSS band
        band = self.config.band.split(" ")[1]

        # Inform the user
        log.info("Processing SDSS " + band + " data ...")

        # Make ...
        self.dpdp.make_sdss_mosaic_and_poisson_frame(self.config.galaxy_name, band, self.path)

# -----------------------------------------------------------------
