#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.truncation Contains the TruncationPlotter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.core.frame import Frame
from ...magic.plot.imagegrid import StandardImageGridPlotter

# -----------------------------------------------------------------

class TruncationPlotter(ModelingComponent):
    
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
        super(TruncationPlotter, self).__init__(config)

        # -- Attributes --

        # The image frames
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
        log.info("Loading the truncated images ...")

        # Loop over all files found in the truncation directory
        for path, name in fs.files_in_path(self.truncation_path, extension="fits", returns=["path", "name"]):

            # Skip the H alpha image
            if "Halpha" in name: continue

            # Debugging
            log.debug("Loading the " + name + " image ...")

            # Open the truncated image
            image = Frame.from_file(path)

            # Add the image to the dictionary
            self.images[name] = image

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
        sorted_labels = sorted(self.images.keys(), key=(lambda key: self.images[key].filter.pivotwavelength() if self.images[key].filter is not None else None))

        # Add the images
        for label in sorted_labels: plotter.add_image(self.images[label], label)

        # Determine the path to the plot file
        path = fs.join(self.plot_path, "truncation.pdf")

        plotter.colormap = "hot"
        plotter.vmin = 0.0

        plotter.set_title("Truncated images")

        # Make the plot
        plotter.run(path)

# -----------------------------------------------------------------
