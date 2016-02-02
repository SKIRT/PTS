#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.photometry Contains the PhotoMeter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import time

# -----------------------------------------------------------------

class PhotoMeter(Configurable):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(PhotoMeter, self).__init__(config, "modeling")

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new PhotoMeter instance
        photometer = cls(arguments.config)

        # Logging
        if arguments.debug:

            photometer.config.logging.level = "DEBUG"
            photometer.config.logging.cascade = True

        # Set the input and output path
        photometer.config.path = arguments.path
        photometer.config.input_path = os.path.join(arguments.path, "data")
        photometer.config.output_path = os.path.join(arguments.path, "prep")

        # A single image can be specified so the photometry is only calculated for that image
        photometer.config.single_image = arguments.image

        # Set logging path
        if arguments.report: photometer.config.logging.path = os.path.join(photometer.config.output_path, time.unique_name("log") + ".txt")

        # Return the new instance
        return photometer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :param image:
        :return:
        """

        # 1. Call the setup function
        self.setup()

# -----------------------------------------------------------------
