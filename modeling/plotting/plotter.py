#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.plotter Contains the Plotter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...magic.core.image import Image
from ...magic.core.frame import Frame
from ...magic.plot.imagegrid import StandardImageGridPlotter

# -----------------------------------------------------------------

class Plotter(PlottingComponent):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(Plotter, self).__init__(config)

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new Plotter instance
        plotter = cls(arguments.config)

        # Set the modeling path
        plotter.config.path = arguments.path

        # Set the modeling step for which to make the plots
        plotter.config.step = arguments.step

        # Return the new instance
        return plotter

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Make a plot of the data
        if self.config.step == "data": self.plot_data()

        # 3. Make a plot of the preparation step
        elif self.config.step == "preparation": self.plot_preparation()

        # 4. Make a plot of the decomposition step
        elif self.config.step == "decomposition": self.plot_decomposition()

        # 5. Make a plot of the truncation step
        elif self.config.step == "truncation": self.plot_truncation()

        # 6. Make a plot of the photometry step
        elif self.config.step == "photometry": self.plot_photometry()

        # 7. Make a plot of the map making step
        elif self.config.step == "maps": self.plot_maps()

        # 8. Make a plot of the fitting step
        elif self.config.step == "fit": self.plot_fit()

        # Invalid modelling step
        else: raise ValueError("Invalid modelling step")

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(Plotter, self).setup()

    # -----------------------------------------------------------------

    def plot_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images ...")

        # The dictionary of image frames
        images = dict()

        # Loop over all subdirectories of the preparation directory
        for directory_path, directory_name in fs.directories_in_path(self.prep_path, returns=["path", "name"]):

            # Debugging
            log.debug("Opening " + directory_name + " image ...")

            # Look if an initialized image file is present
            image_path = fs.join(directory_path, "initialized.fits")
            if not fs.is_file(image_path):

                log.warning("Initialized image could not be found for " + directory_name)
                continue

            # Open the prepared image frame
            frame = Frame.from_file(image_path)

            # Set the image name
            frame.name = directory_name

            # Add the image to the dictionary
            images[directory_name] = frame

        # Inform the user
        log.info("Plotting the images ...")

        # Create the image plotter
        plotter = StandardImageGridPlotter()

        # Sort the image labels based on wavelength
        sorted_labels = sorted(images.keys(), key=lambda key: images[key].filter.pivotwavelength())

        # Add the images
        for label in sorted_labels: plotter.add_image(images[label], label)

        # Determine the path to the plot file
        path = fs.join(self.plot_path, "data.pdf")

        # Set the plot title
        plotter.set_title("Images")

        # Make the plot
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_preparation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the prepared images ...")

        # The dictionary of image frames
        images = dict()

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
            images[directory_name] = frame

        # Inform the user
        log.info("Plotting the images ...")

        # Create the image plotter
        plotter = StandardImageGridPlotter()

        # Sort the image labels based on wavelength
        sorted_labels = sorted(images.keys(), key=lambda key: images[key].filter.pivotwavelength())

        # Add the images
        for label in sorted_labels: plotter.add_image(images[label], label)

        # Determine the path to the plot file
        path = fs.join(self.plot_path, "preparation.pdf")

        # plotter.colormap = "hot"

        plotter.vmin = 0.0

        plotter.set_title("Prepared images")

        # Make the plot
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_decomposition(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

    # -----------------------------------------------------------------

    def plot_truncation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the truncated images ...")

        # The dictionary of image frames
        images = dict()

        # Loop over all files found in the truncation directory
        for path, name in fs.files_in_path(self.truncation_path, extension="fits", returns=["path", "name"]):

            # Skip the H alpha image
            if "Halpha" in name: continue

            # Debugging
            log.debug("Loading the " + name + " image ...")

            # Open the truncated image
            image = Image.from_file(path)

            # Check that the image has a primary and and errors frame
            if "primary" not in image.frames:
                log.warning("The " + name + " image does not contain a primary frame: skipping")
                continue
            if "errors" not in image.frames:
                log.warning("The " + name + " image does not contain an errors frame: skipping")
                continue

            # Add the image to the dictionary
            images[name] = image



    # -----------------------------------------------------------------

    def plot_photometry(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

    # -----------------------------------------------------------------

    def plot_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

    # -----------------------------------------------------------------

    def plot_fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
