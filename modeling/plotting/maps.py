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
from .component import PlottingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log

# -----------------------------------------------------------------

class MapsPlotter(PlottingComponent):
    
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
        super(MapsPlotter, self).__init__(config)

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
        super(MapsPlotter, self).setup()

        # Set the output path
        self.config.output_path = self.plot_path

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the maps ...")

        # The dictionary of image frames
        images = dict()

        for name in ["dust", "ionizing_stars", "old_stars", "young_stars"]:
            # Determine the path to the image
            path = fs.join(self.maps_path, name + ".fits")

            # Debugging
            log.debug("Loading the " + name + " image ...")

            # Open the map
            image = Frame.from_file(path)

            # Add the image to the dictionary
            images[name] = image

        # Inform the user
        log.info("Plotting ...")

        # Create the image plotter
        plotter = StandardImageGridPlotter()

        # Add the images
        for label in images: plotter.add_image(images[label], label)

        # Determine the path to the plot file
        path = fs.join(self.plot_path, "maps.pdf")

        plotter.colormap = "hot"
        plotter.vmin = 0.0

        plotter.set_title("Input maps")

        # Make the plot
        plotter.run(path)

# -----------------------------------------------------------------
