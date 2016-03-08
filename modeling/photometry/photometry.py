#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.photometry.photometry Contains the PhotoMeter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant AstroMagic classes and modules
from ...magic.core import Image
from ...magic.tools import headers
from ...magic.basics import Mask
from ...magic.basics.skyregion import SkyRegion

# Import the relevant PTS classes and modules
from .component import PhotometryComponent
from .sedfetching import SEDFetcher
from ...core.tools import filesystem
from ...core.tools.logging import log
from ..core.sed import ObservedSED
from ...core.basics.errorbar import ErrorBar
from ...core.tools import tables
from ..plotting.sed import SEDPlotter

# -----------------------------------------------------------------

class PhotoMeter(PhotometryComponent):
    
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
        super(PhotoMeter, self).__init__(config)

        # The list of images
        self.images = []

        # The disk ellipse
        self.disk_ellipse = None

        # The SED
        self.sed = None

        # The SEDFetcher
        self.sed_fetcher = None

        # The differences
        self.differences = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new PhotoMeter instance
        photometer = cls(arguments.config)

        # Set the input and output path
        photometer.config.path = arguments.path
        photometer.config.input_path = filesystem.join(arguments.path, "prep")
        photometer.config.output_path = filesystem.join(arguments.path, "phot")

        # A single image can be specified so the photometry is only calculated for that image
        photometer.config.single_image = arguments.image

        # Return the new instance
        return photometer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the prepared images
        self.load_images()

        # 3. Get the disk ellipse
        self.get_disk_ellipse()

        # 4. Do the photometry
        self.do_photometry()

        # 5. Get the photometric flux points from the literature for comparison
        self.get_references()

        # 6. Calculate the differences between the calculated photometry and the reference SEDs
        self.calculate_differences()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(PhotoMeter, self).setup()

        # Create an observed SED
        self.sed = ObservedSED()

        # Create an SEDFetcher instance
        self.sed_fetcher = SEDFetcher()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images ...")

        # Loop over all directories in the preparation directory
        for directory_path, directory_name in filesystem.directories_in_path(self.prep_path, returns="both"):

            # If only a single image has to be processeds, skip the other images
            if self.config.single_image is not None and directory_name != self.config.single_image: continue

            # Determine the filter
            filter_name = directory_name
            if filter_name == "Ha": filter_name = "H alpha"
            filter = headers.get_filter(filter_name)

            # Look for a file called 'result.fits'
            image_path = filesystem.join(directory_path, "result.fits")
            if not filesystem.is_file(image_path):
                log.warning("Prepared image could not be found for " + directory_name)
                continue

            # Open the prepared image
            image = Image.from_file(image_path)

            # Set the image name
            image.name = directory_name

            # Set the filter
            image.filter = filter

            # Add the image to the list
            self.images.append(image)

    # -----------------------------------------------------------------

    def get_disk_ellipse(self):

        """
        This function ...
        :return:
        """

        # Get the path to the disk region
        path = filesystem.join(self.components_path, "disk.reg")

        # Open the region
        region = SkyRegion.from_file(path)

        # Get ellipse in sky coordinates
        self.disk_ellipse = region[0]

    # -----------------------------------------------------------------

    def do_photometry(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing the photometry calculation ...")

        # Create directory to contain the images in Jansky with truncation
        truncated_path = self.full_output_path("truncated")
        filesystem.create_directory(truncated_path)

        # Loop over all the images
        for image in self.images:

            # Inform the user
            log.debug("Performing photometry for " + image.name + " image ...")

            # Inform the user
            log.debug("Creating mask of truncated pixels ...")

            # Create ellipse in image coordinates from ellipse in sky coordinates for this image
            ellipse = self.disk_ellipse.to_ellipse(image.wcs)

            # Create mask from ellipse
            inverted_mask = Mask.from_shape(ellipse, image.xsize, image.ysize, invert=True)

            # Get mask of padded (and bad) pixels, for example for WISE, this mask even covers pixel within the elliptical region
            mask = inverted_mask + image.masks.bad
            if "padded" in image.masks: mask += image.masks.padded

            # Inform the user
            log.debug("Converting image to Jansky ...")

            # Convert from MJy/sr to Jy/sr
            conversion_factor = 1.0
            conversion_factor *= 1e6

            # Conversion from Jy / sr to Jy / pixel
            pixelscale = image.frames.primary.xy_average_pixelscale
            pixel_factor = (1.0/pixelscale**2).to("pix2/sr").value
            conversion_factor /= pixel_factor

            # Frame in Janskys
            jansky_frame = image.frames.primary * conversion_factor
            jansky_frame[mask] = 0.0

            # Debugging
            log.debug("Calculating the total flux and flux error ...")

            # Calculate the total flux in Jansky
            flux = np.sum(jansky_frame)

            jansky_errors_frame = image.frames.errors * conversion_factor
            jansky_errors_frame[mask] = 0.0
            flux_error = np.sum(jansky_errors_frame)

            # Create new image with primary and errors frame in Jansky and save it
            new_image = Image()
            new_image.add_frame(jansky_frame, "primary")
            new_image.add_frame(jansky_errors_frame, "errors")
            new_image_path = filesystem.join(truncated_path, image.name + ".fits")
            new_image.save(new_image_path)

            # Create errorbar
            errorbar = ErrorBar(float(flux_error))

            # Add this entry to the SED
            self.sed.add_entry(image.filter, flux, errorbar)

    # -----------------------------------------------------------------

    def get_references(self):

        """
        This function ...
        :return:
        """

        # Fetch the reference SEDs
        self.sed_fetcher.run(self.galaxy_name)

    # -----------------------------------------------------------------

    def calculate_differences(self):

        """
        This function ...
        :return:
        """

        # Get list of instruments, bands and fluxes of the calculated SED
        instruments = self.sed.instruments()
        bands = self.sed.bands()
        fluxes = self.sed.fluxes(unit="Jy")

        # The number of data points
        number_of_points = len(instruments)

        # Initialize data and names
        reference_labels = self.sed_fetcher.seds.keys()
        data = [[] for _ in range(len(reference_labels)+3)]
        names = ["Instrument", "Band", "Flux"]
        for label in reference_labels:
            names.append(label)

        # Loop over the different points in the calculated SED
        for i in range(number_of_points):

            # Add instrument, band and flux
            data[0].append(instruments[i])
            data[1].append(bands[i])
            data[2].append(fluxes[i])

            column_index = 3

            # Loop over the different reference SEDs
            for label in reference_labels:

                relative_difference = None

                # Loop over the data points in the reference SED
                for j in range(len(self.sed_fetcher.seds[label].table["Wavelength"])):

                    if self.sed_fetcher.seds[label].table["Instrument"][j] == instruments[i] and self.sed_fetcher.seds[label].table["Band"][j] == bands[i]:

                        difference = fluxes[i] - self.sed_fetcher.seds[label].table["Flux"][j]
                        relative_difference = difference / fluxes[i] * 100.

                        # Break because a match has been found within this reference SED
                        break

                # Add percentage to the table (or None if no match was found in this reference SED)
                data[column_index].append(relative_difference)

                column_index += 1

        # Create table of differences
        self.differences = tables.new(data, names=names)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write SED table
        self.write_sed()

        # Write the differences
        self.write_differences()

        # Plot the SED
        self.plot_sed()

        # Plot the SED with references
        self.plot_sed_with_references()

    # -----------------------------------------------------------------

    def write_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing SED to a data file ...")

        # Determine the full path to the output file
        path = self.full_output_path("fluxes.dat")

        # Save the SED
        self.sed.save(path)

    # -----------------------------------------------------------------

    def write_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the percentual differences with reference fluxes to a data file ...")

        # Determine the full path to the output file
        path = self.full_output_path("differences.dat")

        # Save the differences table
        tables.write(self.differences, path)

    # -----------------------------------------------------------------

    def plot_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SED ...")

        # Create a new SEDPlotter instance
        plotter = SEDPlotter("M81")

        # Add the SED
        plotter.add_observed_sed(self.sed, "PTS")

        # Determine the full path to the plot file
        path = self.full_output_path("sed.pdf")
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_sed_with_references(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SED with reference fluxes ...")

        # Create a new SEDPlotter instance
        plotter = SEDPlotter("M81")

        # Add the SED
        plotter.add_observed_sed(self.sed, "PTS")

        # Add the reference SEDs
        for label in self.sed_fetcher.seds: plotter.add_observed_sed(self.sed_fetcher.seds[label], label)

        # Determine the full path to the plot file
        path = self.full_output_path("sed_with_references.pdf")
        plotter.run(path)

# -----------------------------------------------------------------
