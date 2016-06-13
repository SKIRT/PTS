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

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.core.frame import Frame
from ...magic.basics.mask import Mask
from ...magic.basics.region import Region
from ...magic.plot.imagegrid import StandardImageGridPlotter

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

        # The dictionary of prepared image frames
        self.images = dict()

        # The dictionary of sources masks
        self.masks = dict()

        # The dictionary of sky annuli
        self.annuli = dict()

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

        # 3. Load the source masks
        self.load_masks()

        # 4. Load the galaxy and sky annuli
        self.load_annuli()

        # 5. Load the error frames
        self.load_errors()

        # 6. Plot
        self.plot()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the prepared images ...")

        # Loop over all directories in the preparation directory
        for directory_path, directory_name in fs.directories_in_path(self.prep_path, returns=["path", "name"]):

            # Look for a file called 'result.fits'
            image_path = fs.join(directory_path, "result.fits")
            if not fs.is_file(image_path):
                log.warning("Prepared image could not be found for " + directory_name)
                continue

            # Open the prepared image frame
            frame = Frame.from_file(image_path)

            # Set the image name
            frame.name = directory_name

            # Add the image to the dictionary
            self.images[directory_name] = frame

    # -----------------------------------------------------------------

    def load_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the source masks ...")

        # Loop over all directories in the preparation directory
        for directory_path, directory_name in fs.directories_in_path(self.prep_path, returns=["path", "name"]):

            # Look for the 'sources' directory
            #sources_path = fs.join(directory_path, "sources")
            #if not fs.is_directory(sources_path):
            #    log.warning("Sources directory is not present for " + directory_name)
            #    continue

            # Look for a file called 'result.fits'
            image_path = fs.join(directory_path, "result.fits")
            if not fs.is_file(image_path):
                log.warning("Prepared image could not be found for " + directory_name)
                continue

            # Open the sources mask
            mask = Mask.from_file(image_path, plane="sources")

            # Add the mask to the dictionary
            self.masks[directory_name] = mask

    # -----------------------------------------------------------------

    def load_annuli(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the galaxy and sky annuli ...")

        # Loop over all directories in the preparation directory
        for directory_path, directory_name in fs.directories_in_path(self.prep_path, returns=["path", "name"]):

            # Look for the 'sky' directory
            sky_path = fs.join(directory_path, "sky")
            if not fs.is_directory(sky_path):
                log.warning("Sky directory is not present for " + directory_name)
                continue

            # Look for the annulus region file
            region_path = fs.join(sky_path, "annulus.reg")
            if not fs.is_file(region_path):
                log.warning("The annulus region could not be found for " + directory_name)
                continue

            # Open the annulus region
            region = Region.from_file(region_path)

            # Add the region to the dictionary
            self.annuli[directory_name] = region

    # -----------------------------------------------------------------

    def load_errors(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Loading the error frames ...")

        # Loop over all directories in the preparation directory
        for directory_path, directory_name in fs.directories_in_path(self.prep_path, returns=["path", "name"]):

            # Look for a file called 'result.fits'
            image_path = fs.join(directory_path, "result.fits")
            if not fs.is_file(image_path):
                log.warning("Prepared image could not be found for " + directory_name)
                continue

            # Open the prepared image frame
            frame = Frame.from_file(image_path, plane="errors")

            # Set the image name
            frame.name = directory_name

            # Add the error frame to the dictionary
            self.errors[directory_name] = frame

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

        # Sort the image labels based on wavelength
        sorted_labels = sorted(self.images.keys(), key=lambda key: self.images[key].filter.pivotwavelength())

        # Add the images
        for label in sorted_labels: plotter.add_image(self.images[label], label)

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

        # Plot histogram of sky pixels

    # -----------------------------------------------------------------

    def plot_errors(self):

        """
        This function ...
        :return:
        """

        self.plot_error_histograms()

        self.plot_errors_pixels()

    # -----------------------------------------------------------------

    def plot_error_histogram(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_errors_pixels(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_masks_and_annuli(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the images with the sources masks and sky annuli overlayed ...")

# -----------------------------------------------------------------
