#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.data Contains the DataPlotter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.core.frame import Frame
from ...magic.plot.imagegrid import StandardImageGridPlotter

# -----------------------------------------------------------------

class DataPlotter(PlottingComponent):
    
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
        super(DataPlotter, self).__init__(config)

        # -- Attributes --

        # The dictionary of image frames
        self.images = dict()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the images
        self.load_images()

        # 3. Plot
        self.plot()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images ...")

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
            self.images[directory_name] = frame

    # -----------------------------------------------------------------

    def plot(self):

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
        path = fs.join(self.plot_data_path, "data.pdf")

        # Set the plot title
        plotter.set_title("Images")

        # Make the plot
        plotter.run(path)

# -----------------------------------------------------------------
