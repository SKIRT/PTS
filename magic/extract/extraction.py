#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.extraction Contains the Extractor class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant AstroMagic classes and modules
from .core import Frame
from .tools import masks
from .galaxyextraction import GalaxyExtractor
from .starextraction import StarExtractor
from .trainedextractor import TrainedExtractor

# Import the relevant PTS classes and modules
from ..core.tools import filesystem
from ..core.basics.configurable import Configurable

# -----------------------------------------------------------------

class Extractor(Configurable):

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
        super(Extractor, self).__init__(config, "magic")

        # -- Attributes --

        # The image frame on which to perform the extraction
        self.frame = None

        # The mask covering pixels that should be ignored throughout the entire extraction procedure
        self.input_mask = None

        # The galaxy and star extractors
        self.galaxy_extractor = None
        self.star_extractor = None
        self.trained_extractor = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        """

        # Create a new Extractor instance
        extractor = cls(arguments.config)

        # Set the input and output path
        extractor.config.input_path = arguments.input_path
        extractor.config.output_path = arguments.output_path

        # Set options for writing out regions or masks
        extractor.config.write_regions = arguments.regions
        extractor.config.write_masked_frames = arguments.masks

        # Return the new instance
        return extractor

    # -----------------------------------------------------------------

    def run(self, frame, input_mask):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(frame, input_mask)

        # 2. Create the directory that will contain the output
        self.create_output_path()

        # 3. Run the galaxy extraction
        self.extract_galaxies()
        
        # 4. Run the star extraction
        self.extract_stars()

        # 5. Look for other sources
        self.find_other_sources()

        # 5. Writing phase
        self.write()

    # -----------------------------------------------------------------

    def setup(self, frame, input_mask):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(Extractor, self).setup()

        # Inform the user
        self.log.info("Setting up the extractor...")

        # Make a local reference to the frame and mask
        self.frame = frame
        self.input_mask = input_mask

        # Set the paths to the resulting frame and the total mask
        self.config.write_result = False
        self.config.write_mask = True
        self.config.writing.result_path = "result.fits"
        self.config.writing.mask_path = "mask.fits"

        # Initialize a galaxy extractor according to the settings defined in the provided configuration file
        self.galaxy_extractor = GalaxyExtractor(self.config.galaxies)

        # Initialize a star extractor according to the settings defined in the provided configuration file
        self.star_extractor = StarExtractor(self.config.stars)

        # Initialize a trained extractor according to the settings defined in the provided configuration file
        self.trained_extractor = TrainedExtractor(self.config.other_sources)

        # Set the input and output path for the galaxy and star extractor
        self.galaxy_extractor.config.input_path = self.config.input_path
        self.galaxy_extractor.config.output_path = self.config.output_path
        self.star_extractor.config.input_path = self.config.input_path
        self.star_extractor.config.output_path = self.config.output_path
        self.trained_extractor.config.input_path = self.config.input_path
        self.trained_extractor.config.output_path = self.config.output_path

        # Set the appropriate configuration settings for writing out the region files
        if self.config.write_regions:

            # Galaxy extractor
            self.galaxy_extractor.config.write_region = True
            self.galaxy_extractor.config.writing.region_path = "galaxies.reg"

            # Star extractor
            self.star_extractor.config.write_region = True
            self.star_extractor.config.writing.region_path = "stars.reg"

        # Set the appropriate configuration settings for writing out the masked frames
        if self.config.save_masked_frames:

            # Galaxy extractor
            self.galaxy_extractor.config.write_masked_frame = True
            self.galaxy_extractor.config.writing.masked_frame_path = "masked_galaxies.fits"

            # Star extractor
            self.star_extractor.config.write_masked_frame = True
            self.star_extractor.config.writing.masked_frame_path = "masked_stars.fits"

        # Options for logging
        self.galaxy_extractor.config.logging.level = "WARNING"
        self.galaxy_extractor.config.logging.path = self.config.logging.path
        self.star_extractor.config.logging.level = "WARNING"
        self.star_extractor.config.logging.path = self.config.logging.path
        self.trained_extractor.config.logging.level = "WARNING"
        self.trained_extractor.config.logging.path = self.config.logging.path

        # Set the log level and path for the different children of this extractor, if cascading is enabled
        if self.config.logging.cascade:

            # Galaxy extractor
            self.galaxy_extractor.config.loggin.cascade = True
            self.galaxy_extractor.config.logging.level = self.config.logging.level

            # Star extractor
            self.star_extractor.config.logging.cascade = True
            self.star_extractor.config.logging.level = self.config.logging.level

            # Trained extractor
            self.trained_extractor.config.logging.cascade = True
            self.trained_extractor.config.logging.level = self.config.logging.level

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Clearing the extractor")

        # Clear the frame and input mask
        self.frame = None
        self.input_mask = None

        # Clear the extractors
        self.galaxy_extractor.clear()
        self.star_extractor.clear()
        self.trained_extractor.clear()

    # -----------------------------------------------------------------

    def create_output_path(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Creating the output directory if necessary")

        # Create the directory if necessary
        filesystem.create_directory(self.config.output_path)

    # -----------------------------------------------------------------

    def extract_galaxies(self):
        
        """
        This function ...
        """

        # Inform the user
        self.log.info("Extracting the galaxies...")

        # Run the galaxy extractor
        self.galaxy_extractor.run(self.frame, self.input_mask)
    
    # -----------------------------------------------------------------
    
    def extract_stars(self):
        
        """
        This function ...
        """

        # Inform the user
        self.log.info("Extracting the stars...")

        # Run the star extractor
        self.star_extractor.run(self.frame, self.input_mask, self.galaxy_extractor)

    # -----------------------------------------------------------------

    def find_other_sources(self):

        """
        This function ...
        :return:
        """

        # Total mask
        mask = self.input_mask + self.star_extractor.mask + self.galaxy_extractor.mask

        # Run the trained extractor just to find sources
        self.trained_extractor.run(self.frame, mask, self.galaxy_extractor, self.star_extractor)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        self.log.info("Writing...")

        # If requested, write out the result
        if self.config.write_result: self.write_result()

        # If requested, write out the total mask
        if self.config.write_mask: self.write_mask()

    # -----------------------------------------------------------------

    def write_result(self, header=None):

        """
        This function ...
        :return:
        """

        # Determine the full path to the result file
        path = self.full_output_path(self.config.writing.result_path)

        # Inform the user
        self.log.info("Writing resulting frame to " + path)

        # Write out the resulting frame
        self.frame.save(path, header)

    # -----------------------------------------------------------------

    def write_mask(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the mask file
        path = self.full_output_path(self.config.writing.mask_path)

        # Inform the user
        self.log.info("Writing the total mask to " + path)

        # Set the galaxy mask
        galaxy_mask = self.galaxy_extractor.mask

        # Set the star mask
        star_mask = self.star_extractor.mask

        # Create a frame for the total mask
        frame = Frame(masks.union(galaxy_mask, star_mask).astype(int))

        # Write out the total mask
        frame.save(path)

# -----------------------------------------------------------------
