#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.decomposition Contains the DecompositionPlotter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.core.frame import Frame
from ...magic.plot.imagegrid import ResidualImageGridPlotter

# -----------------------------------------------------------------

class DecompositionPlotter(PlottingComponent):
    
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
        super(DecompositionPlotter, self).__init__(config)

        # -- Attributes --

        self.frame = None
        self.bulge = None
        self.disk = None
        self.model = None

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
        log.info("Loading the IRAC 3.6 micron image ...")

        # Determine the path to the truncated 3.6 micron image
        path = fs.join(self.truncation_path, "IRAC I1.fits")
        self.frame = Frame.from_file(path)

        # Convert the frame to Jy/pix
        conversion_factor = 1.0
        conversion_factor *= 1e6

        # Convert the 3.6 micron image from Jy / sr to Jy / pixel
        pixelscale = self.frame.xy_average_pixelscale
        pixel_factor = (1.0 / pixelscale ** 2).to("pix2/sr").value
        conversion_factor /= pixel_factor
        self.frame *= conversion_factor
        self.frame.unit = "Jy"

        # frame.save(fs.join(self.truncation_path, "i1_jy.fits"))

        # Inform the user
        log.info("Loading the bulge image ...")

        # Determine the path to the truncated bulge image
        bulge_path = fs.join(self.truncation_path, "bulge.fits")
        self.bulge = Frame.from_file(bulge_path)

        # Inform the user
        log.info("Loading the disk image ...")

        # Determine the path to the truncated disk image
        disk_path = fs.join(self.truncation_path, "disk.fits")
        self.disk = Frame.from_file(disk_path)

        # Inform the user
        log.info("Loading the model image ...")

        # Determine the path to the truncated model image
        model_path = fs.join(self.truncation_path, "model.fits")
        self.model = Frame.from_file(model_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Create the image plotter
        plotter = ResidualImageGridPlotter()

        # Add the rows
        plotter.add_row(self.frame, self.bulge, "Bulge")
        plotter.add_row(self.frame, self.disk, "Disk")
        plotter.add_row(self.frame, self.model, "Bulge+disk")

        # Set the title
        plotter.set_title("Decomposition")

        # Determine the path to the plot file
        path = fs.join(self.plot_path, "decomposition.pdf")

        plotter.absolute = True
        plotter.colormap = "hot"

        # Make the plot
        plotter.run(path)

# -----------------------------------------------------------------
