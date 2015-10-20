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

import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")

import matplotlib.pyplot as plt
import copy

# Import Astromagic modules
from ..core import masks
from ..core.masks import Mask
from ..core.regions import Region
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

        # Set the extra mask to None initially
        self.extra_mask = None

        # Set the clipped mask to None initially
        self.clipped_mask = None

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

        # Set the extra mask
        if self.config.extra_region is not None: self.set_extra_mask()

        # Sigma-clipping
        if self.config.sigma_clip_mask: self.sigma_clip_mask()

        # Save histogram
        if self.config.save_histogram: self.save_histogram()

        # If requested, save the frame where the galaxies and stars are masked
        if self.config.save_masked_frame: self.save_masked_frame()

        # If requested, save the frame where pixels covered by the sigma-clipped mask are zero
        if self.config.save_clipped_masked_frame: self.save_clipped_masked_frame()

        # If requested, estimate the sky
        if self.config.estimate: self.estimate()

        # If requested, subtract the sky
        if self.config.subtract: self.subtract()

    # *****************************************************************

    def sigma_clip_mask(self):

        """
        This function ...
        :return:
        """

        # Create the sigma-clipped mask
        self.clipped_mask = statistics.sigma_clip_mask(self.frame, self.config.sigma_clipping.sigma_level, self.mask)

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

    def save_histogram(self):

        """
        This function ...
        :return:
        """

        # Create a masked array
        masked = np.ma.masked_array(self.frame, mask=self.mask)
        masked_clipped = np.ma.masked_array(self.frame, mask=self.clipped_mask)

        # Create a figure
        plt.figure(1)

        # Plot the histograms
        #b: blue, g: green, r: red, c: cyan, m: magenta, y: yellow, k: black, w: white
        plt.subplot(211)
        plt.hist(masked.compressed(), 200, range=(-10,20), alpha=0.5, normed=1, facecolor='g', histtype='stepfilled', label='not clipped')
        if self.config.histogram.log_scale: plt.semilogy()

        plt.subplot(212)
        plt.hist(masked_clipped.compressed(), 200, range=(-10,20), alpha=0.5, normed=1, facecolor='g', histtype='stepfilled', label='clipped')
        if self.config.histogram.log_scale: plt.semilogy()

        # Save the figure
        plt.savefig(self.config.saving.histogram_path, bbox_inches='tight', pad_inches=0.25)

    # *****************************************************************

    @property
    def mask(self):

        """
        This function ...
        :return:
        """

        # Create a total mask consisting only of the galaxy mask initially
        mask = self.galaxy_mask

        # Check whether a star mask is present, if so add it to the total mask
        if self.star_mask is not None: mask += self.star_mask

        # Check whether an extra mask is present, if so add it to the total mask
        if self.extra_mask is not None: mask += self.extra_mask

        # Return the total mask
        return mask

    # *****************************************************************

    @property
    def mean(self):

        """
        This function ...
        :return:
        """

        # Return the sigma-clipped mean
        return np.ma.mean(np.ma.masked_array(self.frame, mask=self.clipped_mask))

    # *****************************************************************

    @property
    def median(self):

        """
        This function ...
        :return:
        """

        # Return the sigma-clipped median
        return np.ma.median(np.ma.masked_array(self.frame, mask=self.clipped_mask))

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

    def save_clipped_masked_frame(self):

        """
        This function ...
        :return:
        """

        # Create a frame with masked pixels
        frame = copy.deepcopy(self.frame)
        frame[self.clipped_mask] = 0.0

        # Save the masked frame
        frame.save(self.config.saving.clipped_masked_frame_path)

    # *****************************************************************

    def set_extra_mask(self):

        """
        This function ...
        :return:
        """

        # Load the region and create a mask from it
        extra_region = Region.from_file(self.config.extra_region, self.frame.wcs)
        self.extra_mask = Mask.from_region(extra_region, self.frame.shape)

    # *****************************************************************

    def estimate(self):

        """
        This function ...
        :return:
        """

        # TODO: allow different estimation methods

        # If the mean sky level should be used
        if self.config.estimation.method == "mean":

            # Create a frame filled with the mean value
            self.sky = Frame(np.full(self.frame.shape, self.mean), self.frame.wcs, self.frame.pixelscale, "estimated sky", False, self.frame.unit)

        # If the median sky level should be used
        elif self.config.estimation.method == "median":

            # Create a frame filled with the median value
            self.sky = Frame(np.full(self.frame.shape, self.median), self.frame.wcs, self.frame.pixelscale, "estimated sky", False, self.frame.unit)

        # If the sky should be estimated by using low-resolution interpolation
        elif self.config.estimation.method == "low-res-interpolation":

            # Estimate the sky
            data = interpolation.low_res_interpolation(self.frame, self.config.estimation.downsample_factor, self.mask)

            # Create sky map
            self.sky = Frame(data, self.frame.wcs, self.frame.pixelscale, "estimated sky", False, self.frame.unit)

        # if the sky should be approximated by a polynomial function
        elif self.config.estimatation.method == "polynomial":

            pass

    # *****************************************************************

    def subtract(self):

        """
        This function ...
        :return:
        """

        pass

    # *****************************************************************

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the masks to None
        self.galaxy_mask = None
        self.star_mask = None
        self.clipped_mask = None

        # Set the frames to None
        self.frame = None

# *****************************************************************
