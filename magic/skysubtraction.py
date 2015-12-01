#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for Astronomers        **
# *****************************************************************

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import copy

# Import astronomical modules
from astropy import log

# Import the relevant AstroMagic classes and modules
from .basics import Mask, Region
from .core import Frame
from .tools import statistics, interpolation, plotting

# Import the relevant PTS classes and modules
from ..core.basics import Configurable

# -----------------------------------------------------------------

class SkySubtractor(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(SkySubtractor, self).__init__(config)

        ## Attributes

        # Set the galaxy and star mask to None initialy
        self.galaxy_mask = None
        self.star_mask = None

        # Set the extra mask to None initially
        self.extra_mask = None

        # Set the clipped mask to None initially
        self.clipped_mask = None

        # Set the frame to None initially
        self.frame = None

        # Set the sky frame to None initially
        self.sky = None

        # Set the logger to None initially
        self.log = None

    # -----------------------------------------------------------------

    def run(self, frame, galaxyextractor, starextractor=None):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, galaxyextractor, starextractor)

        # Set the extra mask
        if self.config.extra_region is not None: self.set_extra_mask()

        # Perform sigma-clipping
        if self.config.sigma_clip_mask: self.sigma_clip_mask()

        # If requested, estimate the sky
        if self.config.estimate: self.estimate()

        # If requested, subtract the sky
        if self.config.subtract: self.subtract()

        # Write out the results
        self.write()

    # -----------------------------------------------------------------

    def setup(self, frame, galaxyextractor, starextractor=None):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SkySubtractor, self).setup()

        # Make a local reference to the frame
        self.frame = frame

        # Set the galaxy mask
        self.galaxy_mask = galaxyextractor.mask

        # Set the star mask
        if starextractor is not None: self.star_mask = starextractor.mask

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Write out a histogram of the sky pixels
        if self.config.write_histogram: self.write_histogram()

        # If requested, write out the frame where the galaxies and stars are masked
        if self.config.write_masked_frame: self.write_masked_frame()

        # If requested, write out the frame where pixels covered by the sigma-clipped mask are zero
        if self.config.write_clipped_masked_frame: self.write_clipped_masked_frame()

    # -----------------------------------------------------------------

    def sigma_clip_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating a sigma-clipped mask for the sky")

        # Create the sigma-clipped mask
        self.clipped_mask = statistics.sigma_clip_mask(self.frame, self.config.sigma_clipping.sigma_level, self.mask)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Plot the frame with the sky subtracted
        #plotting.plot_box(np.ma.masked_array(self, mask=self.mask))
        plotting.plot_difference(self.frame, self.sky)

    # -----------------------------------------------------------------

    def write_histogram(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing sky histogram to " + self.config.writing.histogram_path)

        # Create a masked array
        masked = np.ma.masked_array(self.frame, mask=self.mask)
        masked_clipped = np.ma.masked_array(self.frame, mask=self.clipped_mask)

        # Create the PDF figure
        with PdfPages(self.config.writing.histogram_path) as pdf:

            # Create a figure
            fig = plt.figure()

            min = self.mean - 4.0*self.stddev
            max = self.mean + 4.0*self.stddev

            # Plot the histograms
            #b: blue, g: green, r: red, c: cyan, m: magenta, y: yellow, k: black, w: white
            plt.subplot(211)
            plt.hist(masked.compressed(), 200, range=(min,max), alpha=0.5, normed=1, facecolor='g', histtype='stepfilled', label='not clipped')
            if self.config.histogram.log_scale: plt.semilogy()

            plt.subplot(212)
            plt.hist(masked_clipped.compressed(), 200, range=(min,max), alpha=0.5, normed=1, facecolor='g', histtype='stepfilled', label='clipped')
            if self.config.histogram.log_scale: plt.semilogy()

            # Save the figure
            #plt.savefig(self.config.writing.histogram_path, bbox_inches='tight', pad_inches=0.25)

            pdf.savefig(fig)

            # Clear and close
            #plt.close()

    # -----------------------------------------------------------------

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

    # -----------------------------------------------------------------

    @property
    def mean(self):

        """
        This function ...
        :return:
        """

        # Return the sigma-clipped mean
        return np.ma.mean(np.ma.masked_array(self.frame, mask=self.clipped_mask))

    # -----------------------------------------------------------------

    @property
    def median(self):

        """
        This function ...
        :return:
        """

        # Return the sigma-clipped median
        return np.median(np.ma.masked_array(self.frame, mask=self.clipped_mask).compressed())

    # -----------------------------------------------------------------

    @property
    def stddev(self):

        """
        This function ...
        :return:
        """

        # Return the standard deviation of the sigma-clipped frame
        return np.ma.masked_array(self.frame, mask=self.clipped_mask).std()

    # -----------------------------------------------------------------

    def write_masked_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the masked frame to " + self.config.writing.masked_frame_path)

        # Create a frame where the objects are masked
        frame = copy.deepcopy(self.frame)
        frame[self.mask] = float(self.config.writing.mask_value)

        # Save the masked frame
        frame.save(self.config.writing.masked_frame_path)

    # -----------------------------------------------------------------

    def write_clipped_masked_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Saving sigma-clipped masked frame to " + self.config.writing.clipped_masked_frame_path)

        # Create a frame with masked pixels
        frame = copy.deepcopy(self.frame)
        frame[self.clipped_mask] = float(self.config.writing.mask_value)

        # Save the masked frame
        frame.save(self.config.writing.clipped_masked_frame_path)

    # -----------------------------------------------------------------

    def set_extra_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating an extra mask from " + self.config.extra_region)

        # Load the region and create a mask from it
        extra_region = Region.from_file(self.config.extra_region, self.frame.wcs)
        self.extra_mask = Mask.from_region(extra_region, self.frame.shape)

    # -----------------------------------------------------------------

    def estimate(self):

        """
        This function ...
        :return:
        """

        # TODO: allow different estimation methods

        # Inform the user
        log.info("Estimating the sky by using " + self.config.estimation.method)

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

        # Unkown estimation method
        else: raise ValueError("Unkown sky estimation method")

    # -----------------------------------------------------------------

    def subtract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the sky from the frame")

        # Check whether the median sky level exceeds the standard deviation
        if self.median > self.stddev:

            # Inform the user
            log.info("The median sky level exceeds the standard deviation")

            # Subtract the estimated sky from the image frame
            self.frame -= self.sky

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the sky extractor")

        # Set the masks to None
        self.galaxy_mask = None
        self.star_mask = None
        self.clipped_mask = None

        # Set the frames to None
        self.frame = None

# -----------------------------------------------------------------
