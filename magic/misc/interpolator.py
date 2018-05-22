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
from pts.magic.region.list import load_as_pixel_region_list
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
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

        # The original header
        self.header = None

        # The nans mask
        self.nans = None

        # The regions
        self.regions = None

        # The sources
        self.sources = []

        # The interpolation mask
        self.mask = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the sources
        self.load_sources()

        # 3. Interpolate
        self.interpolate()

        # 4. Create the mask
        self.create_mask()

        # 5. Write
        self.write()

        # 6. SHow
        self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Interpolator, self).setup(**kwargs)

        # Load the frame
        if "frame" in kwargs: self.frame = kwargs.pop("frame")
        else: self.load_frame()

        # Load the regions
        if "region" in kwargs: self.regions = [kwargs.pop("region")]
        elif "regions" in kwargs: self.regions = kwargs.pop("regions")
        else: self.load_regions()

    # -----------------------------------------------------------------

    @property
    def wcs(self):

        """
        This function ...
        :return:
        """

        return self.frame.wcs

    # -----------------------------------------------------------------

    def load_frame(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the frame ...")

        # Import the image
        importer = ImageImporter()
        importer.run(self.config.image)

        # Get the primary image frame
        frame = importer.image.primary

        # Get the original header
        self.header = importer.image.original_header

        # Get the original nan pixels
        self.nans = frame.nans

        # Set the NaN pixels to zero in the frame
        frame[self.nans] = 0.0

        # Set the frame
        self.frame = frame

    # -----------------------------------------------------------------

    def load_regions(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the regions ...")

        # Load the region
        self.regions = load_as_pixel_region_list(self.config.regions, self.wcs, only=self.config.shapes, color=self.config.color, ignore_color=self.config.ignore_color)

    # -----------------------------------------------------------------

    def load_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the sources ...")

        # Create sources
        for shape in self.regions:

            # Create a source
            source = Detection.from_shape(self.frame, shape, self.config.source_outer_factor)

            # Add the source
            self.sources.append(source)

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
        for source in self.sources:

            # Estimate the background
            source.estimate_background(self.config.interpolation_method, sigma_clip=self.config.sigma_clip)

            # Replace the pixels by the background
            source.background.replace(self.frame, where=source.mask)

        # Set the original NaN pixels back to NaN
        if self.nans is not None: self.frame[self.nans] = float("nan")

    # -----------------------------------------------------------------

    def create_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the mask ...")

        # Create a mask from the region list
        self.mask = self.regions.to_mask(self.frame.xsize, self.frame.ysize)

        # Set WCS
        self.mask.wcs = self.wcs

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

    @property
    def image_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.config.image)

    # -----------------------------------------------------------------

    def write_frame(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Saving the result ...")

        # Determine the path
        path = self.output_path_file(self.image_name)

        # Check if already existing
        if fs.is_file(path):

            if self.config.replace:
                if self.config.backup: fs.backup_file(path, suffix=self.config.backup_suffix)
                fs.remove_file(path)
            else:
                #raise ValueError("The image already exists")
                path = self.output_path_file("result.fits")

        # Save
        self.frame.saveto(path, header=self.header)

    # -----------------------------------------------------------------

    def write_mask(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Writing the mask ...")

        # Determine the path
        path = self.output_path_file("mask.fits")

        # Write the mask
        self.mask.saveto(path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the result ...")

        # Plot
        if self.config.plot: plotting.plot_box(self.frame)

# -----------------------------------------------------------------
