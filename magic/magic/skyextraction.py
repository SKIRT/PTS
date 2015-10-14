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
from config import Config

# Import Astromagic modules
from ..core import masks
from ..tools import statistics
from ..tools import interpolation
from ..tools import configuration
from ..core.frames import Frame

# *****************************************************************

class SkyExtractor(object):

    """
    This class ...
    """

    def __init__(self, config_file=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Load the configuration
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_path = os.path.join(directory, "config", "skyextractor.cfg")

        # Open the default configuration if no configuration file is specified, otherwise adjust the default
        # settings according to the user defined configuration file
        if config_file is None: self.config = configuration.open(default_path)
        else: self.config = configuration.open(config_file, default=default_path)

        # Set the mask to None initialy
        self.mask = None

        # Set the sky frame to None intially
        self.sky = None

    # *****************************************************************

    def run(self, frame, galaxyextractor, starextractor):

        """
        This function ...
        :return:
        """

        # Create a mask that covers the galaxies and stars (including saturation)
        self.mask = masks.union(galaxyextractor.mask(frame), starextractor.mask(frame))

        # Sigma-clipping
        if self.config.sigma_clip: self.mask = statistics.sigma_clip_mask(frame, self.config.sigma_level, self.mask)

        # TODO: allow different estimation methods

        # Estimate the sky
        data = interpolation.low_res_interpolation(frame, self.config.downsample_factor, self.mask)

        # Create sky map
        self.sky = Frame(data, frame.wcs, frame.pixelscale, frame.description, frame.selected, frame.unit)

        #self.filtered_sky = self.sky

    # *****************************************************************

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the mask to None
        self.mask = None

# *****************************************************************
