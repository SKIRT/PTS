#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.core.sdss Contains the GALEXMosaicMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from .dataprocessing import DustPediaDataProcessing

# -----------------------------------------------------------------

class GALEXMosaicMaker(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(GALEXMosaicMaker, self).__init__(config)

        # The DustPedia data processing instance
        self.dpdp = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Do the mosaicing
        self.mosaic()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(GALEXMosaicMaker, self).setup()

        # Create the DustPedia data processing instance
        self.dpdp = DustPediaDataProcessing()

    # -----------------------------------------------------------------

    def mosaic(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the GALEX mosaics ...")

        # Determine the output path
        output_path = self.output_path

        # Perform the mosaicing
        self.dpdp.make_galex_mosaic_and_poisson_frame(self.config.galaxy_name, output_path)

# -----------------------------------------------------------------
