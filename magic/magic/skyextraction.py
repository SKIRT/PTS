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

        if config is None:

            directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))

            # Load the default configurations for the star remover
            config_path = os.path.join(directory, "config", "skyextractor.cfg")
            self.config = Config(file(config_path))

        else: self.config = config

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
        self.mask = masks.union(galaxyextractor.create_mask(frame), starextractor.create_mask(frame))

        # Sigma-clipping
        if self.config.sigma_clip: self.mask = statistics.sigma_clip_mask(frame, self.config.sigma_level, self.mask)

        # TODO: allow different estimation methods

        # Estimate the sky
        self.sky = interpolation.low_res_interpolation(frame, self.config.downsample_factor, self.mask)

        self.filtered_sky = self.sky

    # *****************************************************************

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the mask to None
        self.mask = None

# *****************************************************************
