#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.component Contains the PlottingComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log

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

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        self.setup()

        self.plot()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(TruncationPlotter, self).setup()

    # -----------------------------------------------------------------

    def plot(self):

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
            image = Frame.from_file(path)

            # Add the image to the dictionary
            images[name] = image

        # Inform the user
        log.info("Plotting the images ...")

        # Create the image plotter
        plotter = StandardImageGridPlotter()

        # Sort the image labels based on wavelength
        sorted_labels = sorted(images.keys(), key=(
        lambda key: images[key].filter.pivotwavelength() if images[key].filter is not None else None))

        # Add the images
        for label in sorted_labels: plotter.add_image(images[label], label)

        # Determine the path to the plot file
        path = fs.join(self.plot_path, "truncation.pdf")

        plotter.colormap = "hot"
        plotter.vmin = 0.0

        plotter.set_title("Truncated images")

        # Make the plot
        plotter.run(path)

# -----------------------------------------------------------------
