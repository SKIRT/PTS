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
import matplotlib.pyplot as plt
import copy

# Import Astromagic modules
from ..core import masks
from ..core.masks import Mask
from ..tools import statistics
from ..tools import interpolation
from ..tools import configuration
from ..core.frames import Frame
from ..tools import plotting

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
        else: self.config = configuration.open(config, default_config)

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
        self.galaxy_mask = galaxyextractor.mask

        # Set the star mask
        if starextractor is not None: self.star_mask = starextractor.mask
        else: self.star_mask = Mask(np.zeros_like(self.frame))

        # Sigma-clipping
        #if self.config.sigma_clip: self.mask = statistics.sigma_clip_mask(frame, self.config.sigma_level, self.mask)

        # TODO: allow different estimation methods

        # Estimate the sky
        #data = interpolation.low_res_interpolation(frame, self.config.downsample_factor, self.mask)

        # Create sky map
        #self.sky = Frame(data, frame.wcs, frame.pixelscale, frame.description, frame.selected, frame.unit)

        #self.filtered_sky = self.sky

        # If requested, save the frame where the galaxies are masked
        if self.config.save_masked_frame: self.save_masked_frame()

    # *****************************************************************

    def plot(self):

        """
        This function ...
        :return:
        """

        # Plot the frame with the sky subtracted
        #plotting.plot_box(np.ma.masked_array(self, mask=self.mask))
        plotting.plot_difference(self.frame, self.sky)

    # *****************************************************************

    def histogram(self):

        """
        This function ...
        :return:
        """

        # Create a masked array
        masked = np.ma.masked_array(self.frame, mask=self.mask)

        # Calculate minimum and maximum value
        min = np.ma.min(masked)
        max = np.ma.max(masked)

        # Plot the histogram
        #b: blue, g: green, r: red, c: cyan, m: magenta, y: yellow, k: black, w: white
        plt.hist(masked.flatten(), 200, range=(min,max), fc='k', ec='k')
        plt.show()

    # *****************************************************************

    @property
    def mask(self):

        """
        This function ...
        :return:
        """

        return masks.union(self.galaxy_mask, self.star_mask)

    # *****************************************************************

    def save_masked_frame(self):

        """
        This function ...
        :return:
        """

        # Create a frame where the objects are masked
        frame = copy.deepcopy(self.frame)
        frame[self.mask] = 0.0

        # Save the masked frame
        frame.save(self.config.saving.masked_frame_path)

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
