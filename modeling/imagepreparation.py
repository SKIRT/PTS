#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import numpy as np
import os.path
import inspect

# Import astronomical modules
import astropy.units as u
from astropy import log
import astropy.logger

# Import Astromagic modules
from astromagic.core.frames import Frame
from astromagic.magic.starextraction import StarExtractor
from astromagic.magic.galaxyextraction import GalaxyExtractor
from astromagic.magic.skyextraction import SkyExtractor
from astromagic.tools import configuration
from astromagic.core import regions
from astromagic.tools import cropping

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

        # If requested, set the uncertainties
        if self.config.set_uncertainties: self.set_uncertainties()

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
        self.galaxyex = GalaxyExtractor(self.config.galaxy_extraction)

        # Run the galaxy extractor
        self.galaxyex.run(self.image.frames[self.config.primary])

    # *****************************************************************

    def extract_stars(self):

        """
        This function ...
        :return:
        """

        # Create a star extractor
        self.starex = StarExtractor(self.config.star_extraction)

        # Run the star extractor
        self.starex.run(self.image.frames[self.config.primary], self.galaxyex)

    # *****************************************************************

    def extract_sky(self):

        """
        This function ...
        :return:
        """

        # Create a sky extractor
        self.skyex = SkyExtractor(self.config.sky_extraction)

        # Run the sky extraction
        self.skyex.run(self.image.frames[self.config.primary], self.galaxyex, self.starex)

        # Print the statistics of the sky frame
        log.info("Mean sky level = " + str(self.skyex.mean))
        log.info("Median sky level = " + str(self.skyex.median))
        log.info("Standard deviation of sky = " + str(self.skyex.stddev))

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

        #### TODO: make this function automatic

        # Create a unit object
        unit = u.Unit(self.config.unit_conversion.to_unit)

        # Inform the user
        log.info("Converting image to " + str(unit) + ", but not yet automatic")

        # Convert the image to different units (primary and errors frame)
        #self.image.convert_to(unit)

        ####

        # THE GALEX FUV IMAGE
        if self.image.name == "GALEXFUV":

            # Get the pixelscale of the image
            pixelscale = self.image.frames.primary.pixelscale

            # Get the wavelength of the image
            wavelength = self.image.frames.primary.filter.centerwavelength()
            wavelength = (wavelength * u.micron).to(u.AA)

            # Speed of light in Angstrom per seconds
            c = (299792458.0 * u.m / u.s).to(u.AA / u.s)
            spectralfactor = wavelength.value**2 / c.value
            pixelfactor = (206264.806247 / pixelscale)**2
            factor = 1.4e-15 * spectralfactor * 1e17 * pixelfactor

            # Multiply the image (primary and errors frame) by the conversion factor
            self.image *= factor

        # THE 2MASS H IMAGE
        elif self.image.name == "2MASSH":

            # Conversion factor to magnitudes
            m0 = 20.6473999

            # Conversion factor back to fluxes
            F0 = 1024.0

            # Get the pixelscale of the image
            pixelscale = self.image.frames.primary.pixelscale

            # Calculate the conversion factor
            pixelfactor = (206264.806247 / pixelscale)**2
            factor = 10e-6 * pixelfactor
            factor *= F0 * np.power(10.0, -m0/2.5)

            # Multiply the image (primary and errors frame) by the conversion factor
            self.image *= factor

        # THE Halpha IMAGE
        elif self.image.name == "Ha":

            # Get the pixelscale of the image
            pixelscale = self.image.frames.primary.pixelscale

            pixelfactor = (206264.806247 / pixelscale)**2

            # The wavelength (in meter)
            wavelength = 0.657894736 * 1e-6

            # The speed of light (in m / s)
            c = 299792458

            # The frequency (in Hz)
            frequency = c / wavelength

            # Calculate the conversion factor
            factor = 1e23 * 1e-6 * pixelfactor / frequency

            # Multiply the image (primary and errors frame) by the conversion factor
            self.image *= factor

        # THE IRAC I1 IMAGE IS ALREADY IN MJY/SR
        elif self.image.name == "IRACI1": pass

        # THE MIPS 24 IMAGE IS ALREADY IN MJY/SR
        elif self.image.name == "MIPS24": pass

        # THE PACS 70 IMAGE
        elif self.image.name == "PACS70":

            # Get the pixelscale of the image
            pixelscale = self.image.frames.primary.pixelscale

            # Calculate the conversion factor
            pixelfactor = (206264.806247 / pixelscale)**2
            factor = 1e-6 * pixelfactor

            # Multiply the image (primary and errors frame) by the conversion factor
            self.image *= factor

        # THE PACS 160 IMAGE
        elif self.image.name == "PACS160":

            # Get the pixelscale of the image
            pixelscale = self.image.frames.primary.pixelscale

            # Calculate the conversion factor
            pixelfactor = (206264.806247 / pixelscale)**2
            factor = 1e-6 * pixelfactor

            # Multiply the image (primary and errors frame) by the conversion factor
            self.image *= factor

        # UNKOWN IMAGE NAME
        else: raise ValueError("Unkown image: " + self.image.name)

    # *****************************************************************

    def convolve(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Open the kernel frame
        kernel_path = "Kernel_HiRes_" + self.config.convolution.aniano_name + "_to_" + self.config.convolution.convolve_to + ".fits"
        kernel = Frame.from_file(kernel_path)

        # Convolve the image (the primary and errors frame)
        self.image.convolve(kernel)

    # *****************************************************************

    def rebin(self):

        """
        This function ...
        :param frame:
        :return:
        """

        # Open the reference frame
        reference = Frame.from_file(self.config.rebinning.rebin_to)

        # Rebin the image (the primary and errors frame)
        self.image.rebin(reference)

    # *****************************************************************

    def set_uncertainties(self):

        """
        This function ...
        :return:
        """

        # Create noise region
        region = regions.Region.from_file(self.config.uncertainties.noise_path, self.image.frames[self.config.primary].wcs)

        # Initialize lists for the mean value and standard deviation in the different shapes
        means = []
        stddevs = []

        # Loop over all shapes
        for shape in region:

            x_center, y_center, x_radius, y_radius = regions.ellipse_parameters(shape)

            box, x_min, x_max, y_min, y_max = cropping.crop(self.image.frames[self.config.primary], x_center, y_center, x_radius, y_radius)

            # Calculate the mean and standard deviation in the box
            mean = np.mean(box)
            stddev = np.std(box)

            means.append(mean)
            stddevs.append(stddev)

        means = np.array(means)
        stddevs = np.array(stddevs)

        # Calculate the global uncertainty
        uncertainty = np.sqrt(np.std(means)**2 + np.median(stddevs))

        # If there is no errors frame
        if self.image.frames[self.config.errors] is None:

            wcs = self.image.frames[self.config.primary].wcs
            pixelscale = self.image.frames[self.config.primary].pixelscale
            unit = self.image.frames[self.config.primary].unit
            filter = self.image.frames[self.config.primary].filter
            selected = True

            # Create the errors frame (select it)
            frame = Frame(np.full(self.image.frames[self.config.primary].shape, uncertainty), wcs, pixelscale, None, selected, unit, self.config.errors, filter)

            # Add the errors frame to the image
            self.image.add_frame(frame, self.config.errors)

        # If there is an errors frame, add the uncertainty to the errors quadratically
        else: self.image.frames[self.config.errors] = np.sqrt(np.power(self.image.frames[self.config.errors], 2) + uncertainty**2)

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
