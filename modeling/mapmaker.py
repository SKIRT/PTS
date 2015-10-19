#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import os.path
import inspect

# Import astronomical modules
import astropy.units as u
from astropy import log
import astropy.logger

# Import Astromagic modules
from astromagic import Image
from astromagic.magic.starextraction import StarExtractor
from astromagic.magic.galaxyextraction import GalaxyExtractor
from astromagic.magic.skyextraction import SkyExtractor
from astromagic.tools import configuration

# *****************************************************************

class MapMaker(object):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "imagepreparation.cfg")

        # Open the default configuration if no configuration file is specified, otherwise adjust the default
        # settings according to the user defined configuration file
        if config is None: self.config = configuration.open(default_config)
        else: self.config = configuration.open(config, default_config)

        ### SET-UP LOGGING SYSTEM

        # Set the log level
        log.setLevel(self.config.logging.level)

        # Set log file path
        if self.config.logging.path is not None: astropy.logger.conf.log_file_path = self.config.logging.path.decode('unicode--escape')

        ### TEMPORARY

        self.config.extract_sky = False
        self.config.correct_for_extinction = False
        self.config.convert_unit = False
        self.config.convolve = False

        ###

        # Set extractors to None
        self.galaxyex = None
        self.starex = None
        self.skyex = None

        # Set the image reference to None initially
        self.image = None

    # *****************************************************************

    def run(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        # Cache a reference to the image
        self.image = image

        # Extract galaxies
        self.extract_galaxies()

        # If requested, extract the stars
        if self.config.extract_stars: self.extract_stars()

        # If requested, extract the sky
        if self.config.extract_sky: self.extract_sky()

        # If requested, correct for galactic extinction
        if self.config.correct_for_extinction: self.correct_for_extinction()

        # If requested, convert the unit
        if self.config.convert_unit: self.convert_unit()

        # If requested, convolve
        if self.config.convolve: self.convolve()

        # If requested, rebin
        if self.config.rebin: self.rebin()

        # If requested, crop
        if self.config.crop: self.crop()

    # *****************************************************************


# *****************************************************************
