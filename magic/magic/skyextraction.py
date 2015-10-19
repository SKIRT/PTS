#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import os.path
import inspect
import numpy as np

# Import Astromagic modules
from ..core import masks
from ..core.masks import Mask
from ..tools import statistics
from ..tools import interpolation
from ..tools import configuration
from ..core.frames import Frame

# Import astronomical modules
from astropy import log
import astropy.logger

# *****************************************************************

class SkyExtractor(object):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "skyextractor.cfg")

        # Open the default configuration if no configuration file is specified, otherwise adjust the default
        # settings according to the user defined configuration file
        if config is None: self.config = configuration.open(default_config)
        else: self.config.configuration.open(config, default_config)

        ### SET-UP LOGGING SYSTEM

        # Set the log level
        log.setLevel(self.config.logging.level)

        # Set log file path
        if self.config.logging.path is not None: astropy.logger.conf.log_file_path = self.config.logging.path.decode('unicode--escape')

        ###

        # Set the galaxy and star mask to None initialy
        self.galaxy_mask = None
        self.star_mask = None

        # Set the frame to None initially
        self.frame = None

        # Set the sky frame to None intially
        self.sky = None

    # *****************************************************************

    def run(self, frame, galaxyextractor, starextractor=None):

        """
        This function ...
        :return:
        """

        # Make a local reference to the frame
        self.frame = frame

        # Set the galaxy mask
        self.galaxy_mask = galaxyextractor.mask()

        # Set the star mask
        if starextractor is not None: self.star_mask = starextractor.mask()
        else: self.star_mask = Mask(np.zeros_like(self.frame))

        # Sigma-clipping
        #if self.config.sigma_clip: self.mask = statistics.sigma_clip_mask(frame, self.config.sigma_level, self.mask)

        # TODO: allow different estimation methods

        # Estimate the sky
        #data = interpolation.low_res_interpolation(frame, self.config.downsample_factor, self.mask)

        # Create sky map
        #self.sky = Frame(data, frame.wcs, frame.pixelscale, frame.description, frame.selected, frame.unit)

        #self.filtered_sky = self.sky

    # *****************************************************************

    @property
    def mask(self):

        """
        This function ...
        :return:
        """

        return masks.union(self.galaxy_mask, self.star_mask)

    # *****************************************************************

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the masks to None
        self.galaxy_mask = None
        self.star_mask = None

        # Set the frames to None
        self.frame = None

# *****************************************************************
