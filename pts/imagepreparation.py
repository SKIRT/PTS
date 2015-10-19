#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import os.path
import inspect

# Import astronomical modules
import astropy.units as u
from astropy import log
import astropy.logger

# Import Astromagic modules
from astromagic import Image
from astromagic.magic.starextraction import StarExtractor
from astromagic.magic.galaxyextraction import GalaxyExtractor
from astromagic.magic.skyextraction import SkyExtractor
from astromagic.tools import configuration

# *****************************************************************

class ImagePreparation(object):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "imagepreparation.cfg")

        # Open the default configuration if no configuration file is specified, otherwise adjust the default
        # settings according to the user defined configuration file
        if config is None: self.config = configuration.open(default_config)
        else: self.config = configuration.open(config, default_config)

        ### SET-UP LOGGING SYSTEM

        # Set the log level
        log.setLevel(self.config.logging.level)

        # Set log file path
        if self.config.logging.path is not None: astropy.logger.conf.log_file_path = self.config.logging.path.decode('unicode--escape')

        ### TEMPORARY

        self.config.extract_sky = False
        self.config.correct_for_extinction = False
        self.config.convert_unit = False
        self.config.convolve = False

        ###

        # Set extractors to None
        self.galaxyex = None
        self.starex = None
        self.skyex = None

        # Set the image reference to None initially
        self.image = None

    # *****************************************************************

    def run(self, image):

        """
        This function ...
        :param image:
        :return:
        """

        # Cache a reference to the image
        self.image = image

        # Extract galaxies
        self.extract_galaxies()

        # If requested, extract the stars
        if self.config.extract_stars: self.extract_stars()

        # If requested, extract the sky
        if self.config.extract_sky: self.extract_sky()

        # If requested, correct for galactic extinction
        if self.config.correct_for_extinction: self.correct_for_extinction()

        # If requested, convert the unit
        if self.config.convert_unit: self.convert_unit()

        # If requested, convolve
        if self.config.convolve: self.convolve()

        # If requested, rebin
        if self.config.rebin: self.rebin()

        # If requested, crop
        if self.config.crop: self.crop()

    # *****************************************************************

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set extractors to None
        self.galaxyex = None
        self.starex = None
        self.skyex = None

        # Set the image reference to None
        self.image = None

    # *****************************************************************

    def extract_galaxies(self):

        """
        This function ...
        :return:
        """

        # Create a galaxy extractor
        #if "galaxyextraction" in self.config:
        self.galaxyex = GalaxyExtractor(self.config.galaxy_extraction)
        #else: self.galaxyex = GalaxyExtractor()

        # Run the galaxy extractor
        self.galaxyex.run(self.image.frames[self.config.primary])

        # Print table
        print(self.galaxyex.table)

    # *****************************************************************

    def extract_stars(self):

        """
        This function ...
        :return:
        """

        # Create a star extractor
        #if "starextraction" in self.config:
        self.starex = StarExtractor(self.config.star_extraction)
        #else: self.starex = StarExtractor()

        # Special region for stars
        #special_region_path = join("special", basename(path)[:-5] + "_special.reg")
        #if isfile(special_region_path): starex.config.special_region = special_region_path

        # Ignore region for stars
        #ignore_region_path = join("ignore", basename(path)[:-5] + "_ignore.reg")
        #if isfile(ignore_region_path): starex.config.ignore_region = ignore_region_path

        # Run the star extractor
        self.starex.run(self.image.frames[self.config.primary], self.galaxyex)

        #image.deselect_all()
        #image.frames.primary.select()

        #mask = starex.aperture_mask(image.frames.primary, expansion_factor=1.1)
        #image.masks["apertures_expanded"] = mask
        #image.masks.apertures_expanded.select()
        #image.apply_masks(0.0)

        #image.frames.primary = image.frames.primary.interpolate(image.masks.apertures_expanded)

        #image.frames.primary.select()

        # Export the image
        #image.save(join("out", "apertures_interpolated_"+basename(path)))

    # *****************************************************************

    def extract_sky(self):

        """
        This function ...
        :return:
        """

        # Create a sky extractor
        #if "skyextraction" in self.config:
        self.skyex = SkyExtractor(self.config.sky_extraction)
        #else: self.skyex = SkyExtractor()

        # Run the sky extraction
        self.skyex.run(self.image.frames[self.config.primary], self.galaxyex, self.starex)

        # Create a sky frame
        #image.add_mask(skyex.mask, "sky")
        #image.masks.sky.select()
        #image.apply_masks(0.0)

        # Export the sky frame
        #image.save(join("out", "sky_" + basename(path)))

    # *****************************************************************

    def correct_for_extinction(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Correct the primary frame for galactic extinction
        self.image.frames[self.config.primary] *= 10**(0.4*self.config.attenuation)

    # *****************************************************************

    def convert_unit(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Create a unit object
        unit = u.Unit(self.config.unit_conversion.to_unit)

        # Convert the image to different units (primary and errors frame)
        self.image.convert_to(unit)

    # *****************************************************************

    def convolve(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Open the kernel image
        kernel_path = "Kernel_HiRes_" + self.config.convolution.aniano_name + "_to_" + self.config.convolution.convolve_to + ".fits"
        kernel = Image(kernel_path)

        # Convolve the image (the primary and errors frame)
        self.image.convolve(kernel.frames.primary)

    # *****************************************************************

    def rebin(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Open the reference image
        reference = Image(self.config.rebinning.rebin_to)

        # Rebin the image (the primary and errors frame)
        self.image.rebin(reference.frames.primary)

    # *****************************************************************

    def crop(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Get the cropping limits
        x_min = self.config.cropping.limits[0]
        x_max = self.config.cropping.limits[1]
        y_min = self.config.cropping.limits[2]
        y_max = self.config.cropping.limits[3]

        # Crop the image (the primary and errors frame)
        self.image.crop(x_min, x_max, y_min, y_max)

# *****************************************************************
