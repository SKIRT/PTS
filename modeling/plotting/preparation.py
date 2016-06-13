#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.preparation Contains the PreparationPlotter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.core.frame import Frame
from ...magic.basics.mask import Mask
from ...magic.basics.region import Region
from ...magic.plot.imagegrid import StandardImageGridPlotter
from ...core.plot.distribution import DistributionGridPlotter
from ...core.basics.distribution import Distribution

# -----------------------------------------------------------------

class PreparationPlotter(PlottingComponent):
    
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
        super(PlottingComponent, self).__init__(config)

        # -- Attributes --

        # The paths to the resulting FITS files
        self.result_paths = dict()

        # The paths to the sky directories
        self.sky_paths = dict()

        # The dictionary of prepared image frames
        self.images = dict()

        # The dictionary of sources masks
        self.sources_masks = dict()

        # The dictionary of sky masks
        self.sky_masks = dict()

        # The dictionary of sky values
        self.sky_values = dict()

        # The dictionary of sky annuli
        self.annuli = dict()

        # The dictionary of sky apertures
        self.apertures = dict()

        # The dictionary of error frames
        self.errors = dict()

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

        # 3. Load the source and sky masks
        self.load_masks()

        # 4. Load the sky values
        self.load_sky()

        # 5. Load the galaxy and sky annuli
        self.load_annuli()

        # 6. Load the sky apertures
        self.load_apertures()

        # 7. Load the error frames
        self.load_errors()

        # 8. Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(PreparationPlotter, self).setup()

        # Loop over all directories in the preparation directory
        for directory_path, directory_name in fs.directories_in_path(self.prep_path, returns=["path", "name"]):

            # Look for a file called 'result.fits'
            image_path = fs.join(directory_path, "result.fits")
            if not fs.is_file(image_path):
                log.warning("Prepared image could not be found for " + directory_name)
                continue

            # Add the image path to the dictionary
            self.result_paths[directory_name] = image_path

            # Look for the 'sky' directory
            sky_path = fs.join(directory_path, "sky")
            if not fs.is_directory(sky_path):
                log.warning("Sky directory is not present for " + directory_name)
                continue

            # Add the sky directory path to the dictionary
            self.sky_paths[directory_name] = sky_path

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the prepared images ...")

        # Loop over the image paths
        for label in self.result_paths:

            # Open the prepared image frame
            frame = Frame.from_file(self.result_paths[label])

            # Set the image name
            frame.name = label

            # Add the image to the dictionary
            self.images[label] = frame

    # -----------------------------------------------------------------

    def load_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the masks ...")

        # Load sources masks
        self.load_sources_masks()

        # Load sky masks
        self.load_sky_masks()

    # -----------------------------------------------------------------

    def load_sources_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the sources masks ...")

        # Loop over the image paths
        for label in self.result_paths:

            # Open the sources mask
            mask = Mask.from_file(self.result_paths[label], plane="sources")

            # Add the mask to the dictionary
            self.sources_masks[label] = mask

    # -----------------------------------------------------------------

    def load_sky_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the sky masks ...")

        # Loop over the image paths
        for label in self.result_paths:

            # Open the sky mask
            mask = Mask.from_file(self.result_paths[label], plane="sky")

            # Add the sky mask to the dictionary
            self.sky_masks[label] = mask

    # -----------------------------------------------------------------

    def load_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the sky values ...")

        # Loop over the image paths
        for label in self.result_paths:

            # Open the sky frame
            sky = Frame.from_file(self.result_paths[label], plane="sky")

            # Get the sky value (assuming the sky frame is constant)
            value = sky[0,0]

            # Add the sky value to the dictionary
            self.sky_values[label] = value

    # -----------------------------------------------------------------

    def load_annuli(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the galaxy and sky annuli ...")

        # Loop over the sky paths
        for label in self.sky_paths:

            # Look for the annulus region file
            region_path = fs.join(self.sky_paths[label], "annulus.reg")
            if not fs.is_file(region_path):
                log.warning("The annulus region could not be found for " + label)
                continue

            # Open the annulus region
            region = Region.from_file(region_path)

            # Add the region to the dictionary
            self.annuli[label] = region

    # -----------------------------------------------------------------

    def load_apertures(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the sky apertures ...")

        # Loop over the sky paths
        for label in self.sky_paths:

            # Look for the apertures FITS file
            apertures_path = fs.join(self.sky_paths[label], "apertures.fits")
            if not fs.is_file(apertures_path):
                log.warning("The apertures image could not be found for " + label)
                continue

            # Open the apertures image
            apertures = Frame.from_file(apertures_path)

            # Add the apertures image to the dictionary
            self.apertures[label] = apertures

    # -----------------------------------------------------------------

    def load_errors(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Loading the error frames ...")

        # Loop over the image paths
        for label in self.result_paths:

            # Open the prepared image frame
            frame = Frame.from_file(self.result_paths[label], plane="errors")

            # Set the image name
            frame.name = label

            # Add the error frame to the dictionary
            self.errors[label] = frame

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot a grid of the prepared images
        self.plot_images()

        # Plot the grid of images with the sources masks and sky annuli overlayed
        self.plot_masks_and_annuli()

        # Plot the sky values
        self.plot_sky()

        # Plot the distributions of the relative errors
        self.plot_errors()

    # -----------------------------------------------------------------

    def plot_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the images ...")

        # Create the image plotter
        plotter = StandardImageGridPlotter()

        # Add the images
        for label in self.sorted_labels: plotter.add_image(self.images[label], label)

        # Determine the path to the plot file
        path = fs.join(self.plot_path, "preparation.pdf")

        # plotter.colormap = "hot"

        plotter.vmin = 0.0

        plotter.set_title("Prepared images")

        # Make the plot
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the sky values ...")

        # Create the distribution grid plotter
        plotter = DistributionGridPlotter()

        # Loop over the different images
        for label in self.sorted_labels:

            # Create the distribution from the image pixel values
            distribution = Distribution.from_values(self.images[label].flatten() + self.sky_values[label])

            # Create an array of all the pixels used for estimating the sky
            notnan = np.logical_not(np.isnan(self.apertures))
            sky_values = self.apertures[notnan]

            # Create the distribution of pixel values used for the sky estimation
            sky_distribution = Distribution.from_values(sky_values)

            # Add the distributions
            plotter.add_distribution(distribution, label)
            plotter.add_distribution(sky_distribution, label)

        # Determine the path to the plot file
        path = fs.join(self.plot_preparation_path, "sky_distribution.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_errors(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the errors ...")

        # Plot histograms of the absolute error values
        self.plot_error_histograms_absolute()

        # Plot histograms of the relative error values
        self.plot_error_histograms_relative()

        # Plot the relative errors of each pixel
        self.plot_errors_pixels()

    # -----------------------------------------------------------------

    def plot_error_histograms_absolute(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absolute error values in a histogram for each prepared image compared to the histogram of the actual image values ...")

        # Create the distribution grid plotter
        plotter = DistributionGridPlotter()

        # Loop over the different images
        for label in self.sorted_labels:

            # Create the distribution from the image pixel values
            distribution = Distribution.from_values(self.images[label].flatten())

            # Create the distribution from the error values
            error_distribution = Distribution.from_values(self.errors[label].flatten())

            # Add an entry to the distribution grid plotter
            plotter.add_distribution(distribution, label)
            plotter.add_distribution(error_distribution, label)

        # Determine the path to the plot file
        path = fs.join(self.plot_preparation_path, "absolute_errors.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_error_histograms_relative(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the relative error values in a histogram for each prepared image ...")

        # Create the distribution grid plotter
        plotter = DistributionGridPlotter()

        # Loop over the different images
        for label in self.sorted_labels:

            # Calculate the relative errors
            rel_errors = self.errors[label] / self.images[label]

            # Create a distribution from the relative errors
            rel_error_distribution = Distribution.from_values(rel_errors)

            # Add the distribution to the plotter
            plotter.add_distribution(rel_error_distribution, label)

        # Determine the path to the plot file
        path = fs.join(self.plot_preparation_path, "relative_errors.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_errors_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the relative error for each pixel in each prepared image ...")



    # -----------------------------------------------------------------

    def plot_masks_and_annuli(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the images with the sources masks and sky annuli overlayed ...")

        # Create the image plotter
        plotter = StandardImageGridPlotter()

        # Add the images
        for label in self.sorted_labels: plotter.add_image(self.images[label], label, mask=self.masks[label])

        # Determine the path to the plot file
        path = fs.join(self.plot_preparation_path, "preparation_masks_annuli.pdf")

        plotter.vmin = 0.0

        plotter.set_title("Prepared images with source masks and sky annuli")

        # Make the plot
        plotter.run(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def sorted_labels(self):

        """
        This function ...
        :return:
        """

        sorted_labels = sorted(self.images.keys(), key=lambda key: self.images[key].filter.pivotwavelength())
        return sorted_labels

# -----------------------------------------------------------------
