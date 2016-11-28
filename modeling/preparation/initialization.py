#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.preparation.initialization Contains the PreparationInitializer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from .component import PreparationComponent
from ...magic.sources.finder import SourceFinder
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.misc.imageimporter import ImageImporter
from ...magic.core.image import Image
from ...magic.core.frame import Frame
from ...core.basics.animation import Animation
from ...core.tools import time
from ...core.tools import parsing
from ...magic.core.dataset import DataSet

# -----------------------------------------------------------------

# Reference: Common-Resolution Convolution Kernels for Space- and Ground-Based Telescopes (G. Aniano et. al)
fwhms = {"GALEX FUV": 4.48 * Unit("arcsec"),
         "GALEX NUV": 5.05 * Unit("arcsec"),
         "Mosaic Halpha": 2.0 * Unit("arcsec"),
         "IRAC I1": 1.90 * Unit("arcsec"),
         "IRAC I2": 1.81 * Unit("arcsec"),
         "IRAC I3": 2.11 * Unit("arcsec"),
         "IRAC I4": 2.82 * Unit("arcsec"),
         "WISE W1": 5.79 * Unit("arcsec"),
         "WISE W2": 6.37 * Unit("arcsec"),
         "WISE W3": 6.60 * Unit("arcsec"),
         "WISE W4": 11.89 * Unit("arcsec"),
         "MIPS 24mu": 6.43 * Unit("arcsec"),
         "MIPS 70mu": 18.74 * Unit("arcsec"),
         "MIPS 160mu": 38.78 * Unit("arcsec"),
         "Pacs blue": 5.67 * Unit("arcsec"),
         "Pacs green": 7.04 * Unit("arcsec"),
         "Pacs red": 11.18 * Unit("arcsec"),
         "SPIRE PSW": 18.15 * Unit("arcsec"),
         "SPIRE PMW": 24.88 * Unit("arcsec"),
         "SPIRE PLW": 36.09 * Unit("arcsec")}

# -----------------------------------------------------------------

# For M81:
#fwhms_from_finding = {"2MASS H": 4.640929858306589 * Unit("arcsec"),
#                      "2MASS J": 4.580828087551186 * Unit("arcsec"),
#                      "2MASS Ks": 4.662813601376219 * Unit("arcsec"),
#                      "SDSS g": 2.015917936060279 * Unit("arcsec"),
#                      "SDSS i": 1.85631074608032 * Unit("arcsec"),
#                      "SDSS r": 2.026862297071852 * Unit("arcsec"),
#                      "SDSS u": 2.327165667182196 * Unit("arcsec"),
#                      "SDSS z": 1.841443699129355 * Unit("arcsec")}

# -----------------------------------------------------------------

class PreparationInitializer(PreparationComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(PreparationInitializer, self).__init__(config)

        # -- Attributes --

        # The frame paths
        self.paths = dict()
        self.error_paths = dict()

        # The source finder
        self.finder = None

        # The initial dataset
        self.set = DataSet()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Get the image paths
        self.get_paths()

        # 3. Process the images (identify sources, create error frames)
        self.process_images()

        # 4. Create the dataset
        self.create_dataset()

        # 5. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(PreparationInitializer, self).setup()

        # Create the source finder
        self.finder = SourceFinder(self.config.sources)

    # -----------------------------------------------------------------

    def get_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Looking for images and error frames ...")

        # Loop over the different image origins
        for path, origin in fs.directories_in_path(self.data_images_path, returns=["path", "name"]):

            # Ignore the Planck data (for now)
            if origin == "Planck": continue

            # Loop over the FITS files in the current directory
            for image_path, image_name in fs.files_in_path(path, extension="fits", not_contains="poisson", returns=["path", "name"]):

                # Open the image frame
                frame = Frame.from_file(image_path)

                # Determine the preparation name
                if frame.filter is not None: prep_name = str(frame.filter)
                else: prep_name = image_name

                # Add the image path
                self.paths[prep_name] = image_path

                # Determine path to poisson error map
                poisson_path = fs.join(path, image_name + "_poisson.fits")

                # Set the path to the poisson error map
                if fs.is_file(poisson_path): self.error_paths[prep_name] = poisson_path

    # -----------------------------------------------------------------

    def process_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing the images ...")

        # Loop over all image paths
        for prep_name in self.paths:

            # Get the image path
            image_path = self.paths[prep_name]

            # Determine the output path for this image
            output_path = self.get_prep_path(prep_name)

            # Check whether this image already has an initialized image
            initialized_path = fs.join(output_path, "initialized.fits")
            if fs.is_file(initialized_path): continue

            # Debugging
            log.debug("Processing image '" + image_path + "' ...")

            # Set the path to the region of bad pixels
            bad_region_path = fs.join(self.data_path, "bad", prep_name + ".reg")
            if not fs.is_file(bad_region_path): bad_region_path = None

            # Set the FWHM if the instrument has a fixed PSF
            if prep_name in fwhms: fwhm = fwhms[prep_name]
            else: fwhm = None

            # Debugging
            log.debug("Loading image " + image_path + " as " + prep_name + " ...")

            # Import the image
            importer = ImageImporter()
            importer.run(image_path, bad_region_path, fwhm=fwhm, find_error_frame=False) # don't look for error frames

            # Get the imported image
            image = importer.image

            # Set the image name
            image.name = prep_name

            # -----------------------------------------------------------------

            # Remove all frames except for the primary frame
            image.remove_frames_except("primary")

            # -----------------------------------------------------------------

            # Determine the path to the "sources" directory within the output path for this image
            sources_output_path = fs.join(output_path, "sources")

            # If the source finding step has already been performed on this image, don't do it again
            if fs.is_directory(sources_output_path):

                # Debugging
                log.debug("Source finder output has been found for this image, skipping source finding step")

                # Set the FWHM of the image from the output of the SourceFinder if it is undefined
                if image.fwhm is None: self.set_fwhm_from_source_finder(image, sources_output_path)

            # The source finding step has yet to be performed
            else:

                # Debugging
                log.debug("Source finder will be run for this image")

                # Create the directory if necessary
                fs.create_directory(sources_output_path)

                # Find sources for the current image
                self.find_sources_for_image(image, sources_output_path)

            # -----------------------------------------------------------------

            # Save the image
            image.save(initialized_path)

    # -----------------------------------------------------------------

    def create_dataset(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the initial dataset ...")

        # Loop over the image paths
        for prep_name in self.paths:

            # Add entry to the dataset
            self.set.add_path(prep_name, self.paths[prep_name])

            # Set the path to the poisson error map
            if prep_name in self.error_paths: self.set.add_error_path(prep_name, self.error_paths[prep_name])

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the dataset
        self.write_dataset()

    # -----------------------------------------------------------------

    def write_dataset(self):

        """
        This function ...
        :return:
        """

        # Save the dataset
        self.set.saveto(self.initial_dataset_path)

    # -----------------------------------------------------------------

    def set_fwhm_from_source_finder(self, image, sources_output_path):

        """
        This function ...
        :param image:
        :param sources_output_path:
        :return:
        """

        # Determine the path to the sources/statistics file
        statistics_path = fs.join(sources_output_path, "statistics.dat")

        # Check whether the file exists
        if not fs.is_file(statistics_path): raise RuntimeError("The statistics file could not be found")

        # Get the FWHM from the statistics file
        fwhm = None
        with open(statistics_path) as statistics_file:
            for line in statistics_file:
                if "FWHM" in line: fwhm = parsing.quantity(line.split("FWHM: ")[1].replace("\n", ""))

        # Check whether the FWHM is valid
        if fwhm is None: raise RuntimeError("The FWHM could not be found")

        # Set the FWHM of the image
        image.fwhm = fwhm

    # -----------------------------------------------------------------

    def find_sources_for_image(self, image, sources_output_path):

        """
        This function ...
        :param image:
        :param sources_output_path:
        :return:
        """

        # Get the mask of bad pixels
        bad_mask = image.masks.bad if "bad" in image.masks else None

        # Don't look for stars in the Halpha image
        if "Halpha" in image.name: self.finder.config.find_stars = False
        else: self.finder.config.find_stars = True  # still up to the SourceFinder to decide whether stars should be found (based on the filter)

        # Fix: don't look for other sources in the IRAC images
        if "IRAC" in image.name: self.finder.config.find_other_sources = False
        else: self.finder.config.find_other_sources = True

        # Create an animation for the source finder
        if self.config.visualise: animation = Animation()
        else: animation = None

        # Run the source finder on this image
        self.finder.run(image.frames.primary, self.galactic_catalog, self.stellar_catalog, bad_mask=bad_mask, animation=animation)

        # Write the animation
        if self.config.visualise:

            # Determine the path to the animation
            path = fs.join(self.visualisation_path, time.unique_name(self.image.name + "_sourcefinding") + ".gif")

            # Debugging
            log.debug("Writing animation of the source finding to '" + path + "' ...")

            # Save the animation
            animation.save(path)

        # Save the galaxy region
        galaxy_region = self.finder.galaxy_region
        galaxy_region_path = fs.join(sources_output_path, "galaxies.reg")
        galaxy_region.save(galaxy_region_path)

        # Save the star region
        star_region = self.finder.star_region
        star_region_path = fs.join(sources_output_path, "stars.reg")
        if star_region is not None: star_region.save(star_region_path)

        # Save the saturation region
        saturation_region = self.finder.saturation_region
        saturation_region_path = fs.join(sources_output_path, "saturation.reg")
        if saturation_region is not None: saturation_region.save(saturation_region_path)

        # Save the region of other sources
        other_region = self.finder.other_region
        path = fs.join(sources_output_path, "other_sources.reg")
        if other_region is not None: other_region.save(path)

        # -----------------------------------------------------------------

        # Create an image with the segmentation maps
        segments = Image("segments")

        # Add the segmentation map of the galaxies
        segments.add_frame(self.finder.galaxy_segments, "galaxies")

        # Add the segmentation map of the saturated stars
        if self.finder.star_segments is not None: segments.add_frame(self.finder.star_segments, "stars")

        # Add the segmentation map of the other sources
        if self.finder.other_segments is not None: segments.add_frame(self.finder.other_segments, "other_sources")

        # Save the FITS file with the segmentation maps
        path = fs.join(sources_output_path, "segments.fits")
        segments.save(path)

        # -----------------------------------------------------------------

        # Get the FWHM
        fwhm = self.finder.fwhm

        # Debugging
        log.debug("The FWHM as determined by the source finder is " + str(fwhm) + " ...")

        # Set the FWHM of the image
        if image.fwhm is None: image.fwhm = fwhm

        # -----------------------------------------------------------------

        # Write statistics file
        statistics_path = fs.join(sources_output_path, "statistics.dat")
        self.finder.write_statistics(statistics_path)

        # -----------------------------------------------------------------

        # Clear the source finder
        self.finder.clear()

# -----------------------------------------------------------------
