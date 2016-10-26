#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.prepare.batch Contains the BatchImagePreparer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.basics.filter import Filter
from ...core.basics.configurable import Configurable
from ..misc.calibration import CalibrationError
from ..misc.extinction import GalacticExtinction

# -----------------------------------------------------------------

class BatchImagePreparer(Configurable):
    
    """
    This class ...
    """
    
    def __init__(self, config=None):
        
        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(BatchImagePreparer, self).__init__(config)

        # The input dataset
        self.dataset = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Get the galactic attenuation
        self.get_extinction()

        # Getting calibration error
        #CalibrationError.from_filter(image.filter)

        # 3. Prepare the images
        self.prepare()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BatchImagePreparer, self).setup()

        # Get the dataset
        self.dataset = kwargs.pop("dataset")

    # -----------------------------------------------------------------

    def get_extinction(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the galactic extinction ...")

        # From modeling preparation:

        # Create the galactic extinction calculator
        extinction = GalacticExtinction(self.center_coordinate)

        # Loop over all image paths
        for image_path in self.paths:

            # Get the filter name
            filter_name = fs.name(fs.directory_of(image_path))

            # Create a filter instance
            fltr = Filter.from_string(filter_name)

            # Get the exintinction
            self.attenuations[filter_name] = extinction.extinction_for_filter(fltr)

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the images ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
