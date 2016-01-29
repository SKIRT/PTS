#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.subtract.skysubtraction Contains the SkySubtractor class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
import matplotlib.pyplot as plt
import copy

# Import the relevant AstroMagic classes and modules
from ..core import Frame
from ..basics import Region, Mask
from ..tools import interpolation, plotting, masks, statistics, regions

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable

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
        super(SkySubtractor, self).__init__(config, "magic")

        # -- Attributes --

        # The image frame
        self.frame = None

        # The input mask (galaxies, stars and other objects)
        self.input_mask = None

        # The mask of bad pixels
        self.bad_mask = None

        # The galaxy and saturation region
        self.galaxy_region = None
        self.saturation_region = None

        # The output mask (combined input + bad mask + sigma-clipping mask)
        self.mask = None

        # The estimated sky frame
        self.sky = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new SkySubtractor instance
        if arguments.config is not None: subtractor = cls(arguments.config)
        elif arguments.settings is not None: subtractor = cls(arguments.settings)
        else: subtractor = cls()

        # Debug mode
        if arguments.debug:

            subtractor.config.logging.level = "DEBUG"
            subtractor.config.logging.cascade = True

        # Report file
        if arguments.report: subtractor.config.logging.path = "log.txt"

        # Set the input and output path
        if arguments.input_path is not None: subtractor.config.input_path = arguments.input_path
        if arguments.output_path is not None: subtractor.config.output_path = arguments.output_path

        # Return the new instance
        return subtractor

    # -----------------------------------------------------------------

    def run(self, frame, input_mask, galaxy_region_path, saturation_region_path=None, bad_mask=None):

        """
        This function ...
        :param frame:
        :param input_mask:
        :param bad_mask:
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, input_mask, galaxy_region_path, saturation_region_path, bad_mask)

        # 2. Create mask
        self.create_mask()

        # 2. Perform sigma-clipping
        if self.config.sigma_clip_mask: self.sigma_clip()

        # 3. Estimate the sky
        #self.estimate()

        # 4. If requested, subtract the sky
        #if self.config.subtract: self.subtract()

        # 5. Write out the results
        self.write()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Clearing the sky extractor ...")

        # Set all attributes to None
        self.frame = None
        self.input_mask = None
        self.bad_mask = None
        self.galaxy_region = None
        self.mask = None
        self.sky = None

    # -----------------------------------------------------------------

    def setup(self, frame, input_mask, galaxy_region_path, saturation_region_path=None, bad_mask=None):

        """
        This function ...
        :param frame:
        :param input_mask:
        :param bad_mask:
        :return:
        """

        # Call the setup function of the base class
        super(SkySubtractor, self).setup()

        # Set the paths to the resulting frame and the total mask
        self.config.write_result = True
        self.config.write_masked_frame = True
        self.config.writing.result_path = "subtracted.fits"
        self.config.writing.masked_frame_path = "masked_sky_frame.fits"

        # Make a local reference to the frame
        self.frame = frame

        # The input mask and bad pixel mask
        self.input_mask = input_mask
        self.bad_mask = bad_mask

        # Set the galaxy region and saturation region
        self.galaxy_region = Region.from_file(galaxy_region_path)
        if saturation_region_path is not None: self.saturation_region = Region.from_file(saturation_region_path)

    # -----------------------------------------------------------------

    def create_mask(self):

        """
        This function ...
        :return:
        """

        # Get the ellipse describing the principal galaxy
        principal_ellipse = None
        for shape in self.galaxy_region:

            if not shape.name == "ellipse" and not shape.name == "circle": continue

            # Get the center and radius of the shape (can be a circle or an ellipse)
            ellipse = regions.ellipse(shape)
            if principal_ellipse is None or ellipse.major > principal_ellipse.major: principal_ellipse = ellipse

        # Create a mask from the ellipse
        annulus_outer_factor = 1.5
        annulus_inner_factor = 1.0
        annulus_mask = Mask.from_ellipse(self.frame.xsize, self.frame.ysize, principal_ellipse * annulus_outer_factor).inverse() + \
                       Mask.from_ellipse(self.frame.xsize, self.frame.ysize, principal_ellipse * annulus_inner_factor)

        # Set the mask, make a copy of the input mask initially
        self.mask = self.input_mask + annulus_mask

        if self.bad_mask is not None: self.mask += self.bad_mask

        if self.saturation_region is not None:

            for shape in self.saturation_region:

                # Expand
                shape.coord_list[2] *= 1.5
                shape.coord_list[3] *= 1.5

            saturation_mask = Mask.from_region(self.saturation_region, self.frame.shape)

            self.mask += saturation_mask

    # -----------------------------------------------------------------

    def sigma_clip(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Performing sigma-clipping on the pixel values ...")

        ### TEMPORARY: WRITE OUT MASK BEFORE CLIPPING

        # Create a frame where the objects are masked
        frame = copy.deepcopy(self.frame)
        frame[self.mask] = float(self.config.writing.mask_value)

        # Save the masked frame
        frame.save("masked_sky_frame_notclipped.fits")

        ###

        # Create the sigma-clipped mask
        self.mask = statistics.sigma_clip_mask(self.frame, self.config.sigma_clipping.sigma_level, self.mask)

    # -----------------------------------------------------------------

    def estimate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Estimating the sky by using " + self.config.estimation.method + " ...")

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
        self.log.info("Subtracting the sky from the frame ...")

        # Check whether the median sky level exceeds the standard deviation
        if self.median > self.stddev:

            # Inform the user
            self.log.info("The median sky level exceeds the standard deviation")

            # Subtract the estimated sky from the image frame
            self.frame -= self.sky

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Writing ...")

        # Write out the result
        #if self.config.write_result: self.write_result()

        # Write out a histogram of the sky pixels
        #if self.config.write_histogram: self.write_histogram()

        # If requested, write out the masked frame
        if self.config.write_masked_frame: self.write_masked_frame()

        # If requested, write out the sky statistics
        if self.config.write_statistics: self.write_statistics()

    # -----------------------------------------------------------------

    def write_result(self, header=None):

        """
        This function ...
        :return:
        """

        # Determine the full path to the result file
        path = self.full_output_path(self.config.writing.result_path)

        # Inform the user
        self.log.info("Writing resulting frame to " + path + " ...")

        # Write out the resulting frame
        self.frame.save(path, header)

    # -----------------------------------------------------------------

    def write_masked_frame(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the masked frame file
        path = self.full_output_path(self.config.writing.masked_frame_path)

        # Inform the user
        self.log.info("Writing masked frame to " + path + " ...")

        # Create a frame where the objects are masked
        frame = self.frame.copy()
        frame[self.mask] = float(self.config.writing.mask_value)

        # Write out the masked frame
        frame.save(path)

    # -----------------------------------------------------------------

    def write_statistics(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the statistics file
        path = self.full_output_path(self.config.writing.statistics_path)

        # Inform the user
        self.log.info("Writing statistics to " + path + " ...")

        

    # -----------------------------------------------------------------

    def write_histogram(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Writing sky histogram to " + self.config.writing.histogram_path +  " ...")

        # Create a masked array
        masked = np.ma.masked_array(self.frame, mask=self.mask)
        masked_clipped = np.ma.masked_array(self.frame, mask=self.clipped_mask)

        # Create a figure
        fig = plt.figure()

        min = self.mean - 4.0 * self.stddev
        max = self.mean + 4.0 * self.stddev

        # Plot the histograms
        #b: blue, g: green, r: red, c: cyan, m: magenta, y: yellow, k: black, w: white
        plt.subplot(211)
        plt.hist(masked.compressed(), 200, range=(min,max), alpha=0.5, normed=1, facecolor='g', histtype='stepfilled', label='not clipped')
        if self.config.histogram.log_scale: plt.semilogy()

        plt.subplot(212)
        plt.hist(masked_clipped.compressed(), 200, range=(min,max), alpha=0.5, normed=1, facecolor='g', histtype='stepfilled', label='clipped')
        if self.config.histogram.log_scale: plt.semilogy()

        # Save the figure
        plt.savefig(self.config.writing.histogram_path, bbox_inches='tight', pad_inches=0.25)
        plt.close()

    # -----------------------------------------------------------------

    @property
    def mean(self):

        """
        This function ...
        :return:
        """

        # Return the sigma-clipped mean
        return np.ma.mean(np.ma.masked_array(self.frame, mask=self.mask))

    # -----------------------------------------------------------------

    @property
    def median(self):

        """
        This function ...
        :return:
        """

        # Return the sigma-clipped median
        return np.median(np.ma.masked_array(self.frame, mask=self.mask).compressed())

    # -----------------------------------------------------------------

    @property
    def stddev(self):

        """
        This function ...
        :return:
        """

        # Return the standard deviation of the sigma-clipped frame
        return np.ma.masked_array(self.frame, mask=self.mask).std()

# -----------------------------------------------------------------
