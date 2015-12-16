#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       AstroMagic -- the image editor for astronomers        **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant AstroMagic classes and modules
from . import GalaxyExtractor, StarExtractor

# Import the relevant PTS classes and modules
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

        ## Attributes

        # Set the frame to None initially
        self.frame = None

        # Create the galaxy and star extractors
        self.galaxy_extractor = GalaxyExtractor()
        self.star_extractor = StarExtractor()

    # -----------------------------------------------------------------

    def run(self, frame):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(frame)

        # 2. Create a temporary directory for the input and output of SExtractor
        self.create_output_path()

        # 3. Run the galaxy extraction
        self.extract_galaxies()
        
        # 4. Run the star extraction
        self.extract_stars()

        # 5. Write the mask
        self.write()

    # -----------------------------------------------------------------
    
    @classmethod
    def from_arguments(cls, arguments):
        
        """
        This function ...
        :param arguments:
        """

        # Create a new Extractor instance
        extractor = cls(arguments.config)

        # Set the output path
        extractor.config.output_path = arguments.output_path

        # Set options for saving regions or masks
        extractor.config.save_regions = arguments.regions
        extractor.config.save_masks = arguments.masks

        # Return the new instance
        return extractor

    # -----------------------------------------------------------------

    def setup(self, frame):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(Extractor, self).setup()

        # Make a local reference to the frame
        self.frame = frame

        # Set the appropriate configuration settings for saving the region files
        if self.config.save_regions:

            self.galaxy_extractor.config.save_region = True
            self.galaxy_extractor.config.saving.region_path = os.path.join(self.config.output_path, "galaxies.reg")
            self.star_extractor.config.save_region = True
            self.star_extractor.config.saving.region_path = os.path.join(self.config.output_path, "stars.reg")

        # Set the appropriate configuration settings for saving the masked frames
        if self.config.save_masks:

            self.galaxy_extractor.config.save_masked_frame = True
            self.galaxy_extractor.config.saving.masked_frame_path = os.path.join(self.config.output_path, "masked_galaxies.fits")
            self.star_extractor.config.save_masked_frame = True
            self.star_extractor.config.saving.masked_frame_path = os.path.join(self.config.output_path, "masked_stars.fits")
        
    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def create_output_path(self):

        """
        This function ...
        :return:
        """

        # Create the directory
        os.mkdir(self.config.output_path)

    # -----------------------------------------------------------------

    def extract_galaxies(self):
        
        """
        This function ...
        """

        # Run the galaxy extractor
        self.galaxy_extractor.run(self.frame)
    
    # -----------------------------------------------------------------
    
    def extract_stars(self):
        
        """
        This function ...
        """

        # Run the star extractor
        self.star_extractor.run(self.frame, self.galaxy_extractor)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Set the galaxy mask
        galaxy_mask = self.galaxy_extractor.mask

        # Set the star mask
        star_mask = self.star_extractor.mask

# -----------------------------------------------------------------
