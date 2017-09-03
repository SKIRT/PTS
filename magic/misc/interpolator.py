#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.interpolate Interpolate an image within regions defined by the user.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from pts.magic.misc.imageimporter import ImageImporter
from pts.magic.tools import interpolation
from pts.magic.region.list import load_as_pixel_region_list
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.basics.log import log
from pts.magic.core.cutout import interpolation_methods
from pts.magic.core.detection import Detection
from pts.magic.tools import plotting

# -----------------------------------------------------------------

class Interpolator(Configurable):

    """
    This class
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(Interpolator, self).__init__(*args, **kwargs)

        # The frame
        self.frame = None

        # The regions
        self.regions = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Load the sources
        self.load_sources()

        # Interpolate
        self.interpolate()

        # Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Interpolator, self).setup(**kwargs)

        # Set input path
        if config.input is not None:
            input_path = fs.absolute_or_in_cwd(config.input)
        else:
            input_path = fs.cwd()

        # Set output path
        if config.output is not None:
            output_path = fs.absolute_or_in_cwd(config.output)
        else:
            output_path = fs.cwd()

        # Load the frame
        if "frame" in kwargs: self.frame = kwargs.pop("frame")
        else: self.load_frame()

        # Load the regions
        if "region" in kwargs: self.regions = [kwargs.pop("region")]
        elif "regions" in kwargs: self.regions = kwargs.pop("regions")
        else: self.load_regions()

    # -----------------------------------------------------------------

    def load_frame(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the frame ...")

        # Determine the full path to the image
        #image_path = fs.absolute_path(config.image)

        # Import the image
        importer = ImageImporter()
        importer.run(image_path)

        # Get the primary image frame
        frame = importer.image.primary

        # Get the original header
        header = importer.image.original_header

        # Create a mask of the pixels that are NaNs
        nans = frame.nans()

        # Set the NaN pixels to zero in the frame
        frame[nans] = 0.0

    # -----------------------------------------------------------------

    def load_regions(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the regions ...")

        # Load the region
        region_path = fs.join(input_path, config.regions)
        self.regions = load_as_pixel_region_list(region_path, frame.wcs, only=config.shapes, color=config.color,
                                            ignore_color=config.ignore_color)

    # -----------------------------------------------------------------

    def load_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the sources ...")

        sources = []

        # Create sources
        for shape in regions:
            # Create a source
            source = Detection.from_shape(frame, shape, config.source_outer_factor)
            sources.append(source)

    # -----------------------------------------------------------------

    def interpolate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the frame within the masked pixels ...")

        # Inform the user
        log.info("Interpolating ...")

        # Loop over the sources
        for source in sources:

            # Estimate the background
            source.estimate_background(config.interpolation_method, sigma_clip=config.sigma_clip)

            # Replace the pixels by the background
            source.background.replace(frame, where=source.mask)

        # Set the original NaN pixels back to NaN
        frame[nans] = float("nan")

        # Inform the user
        #log.info("Creating a mask from the region ...")

        # Create a mask from the region list
        mask = regions.to_mask(frame.xsize, frame.ysize)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the frmae
        self.write_frame()

        # Write the mask
        self.write_mask()

    # -----------------------------------------------------------------

    def write_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Saving the result ...")

        # Determine the path
        path = fs.join(output_path, config.image)
        if fs.is_file(path):
            if config.replace:
                if config.backup:
                    fs.backup_file(path, suffix=config.backup_suffix)
                    fs.remove_file(path)
                else:
                    fs.remove_file(path)
            else:
                raise ValueError("The image already exists")

        # Save
        frame.saveto(path, header=header)

    # -----------------------------------------------------------------

    def write_mask(self):

        """
        Thisfunction ...
        :return:
        """

        path = fs.join(output_path, "mask.fits")
        mask.saveto(path, header=header)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Plot
        if self.config.plot: plotting.plot_box(frame)

# -----------------------------------------------------------------
