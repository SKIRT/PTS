#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.plotting.decomposition.plotter Contains the DecompositionPlotter class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import PlottingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log

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
        super(DecompositionPlotter, self).setup()

    # -----------------------------------------------------------------

    def plot(self):

        # Inform the user
        log.info("Plotting ...")

        # Inform the user
        log.info("Loading the IRAC 3.6 micron image ...")

        # Determine the path to the truncated 3.6 micron image
        path = fs.join(self.truncation_path, "IRAC I1.fits")
        frame = Frame.from_file(path)

        # Convert the frame to Jy/pix
        conversion_factor = 1.0
        conversion_factor *= 1e6

        # Convert the 3.6 micron image from Jy / sr to Jy / pixel
        pixelscale = frame.xy_average_pixelscale
        pixel_factor = (1.0 / pixelscale ** 2).to("pix2/sr").value
        conversion_factor /= pixel_factor
        frame *= conversion_factor
        frame.unit = "Jy"

        # frame.save(fs.join(self.truncation_path, "i1_jy.fits"))

        # Inform the user
        log.info("Loading the bulge image ...")

        # Determine the path to the truncated bulge image
        bulge_path = fs.join(self.truncation_path, "bulge.fits")
        bulge = Frame.from_file(bulge_path)

        # Inform the user
        log.info("Loading the disk image ...")

        # Determine the path to the truncated disk image
        disk_path = fs.join(self.truncation_path, "disk.fits")
        disk = Frame.from_file(disk_path)

        # Inform the user
        log.info("Loading the model image ...")

        # Determine the path to the truncated model image
        model_path = fs.join(self.truncation_path, "model.fits")
        model = Frame.from_file(model_path)

        # -----------------------------------------------------------------

        # Calculate the bulge residual frame
        # bulge_residual = frame - bulge
        # bulge_residual_path = fs.join(residuals_path, "bulge_residual.fits")
        # bulge_residual.save(bulge_residual_path)

        # Calculate the disk residual frame
        # disk_residual = frame - disk
        # disk_residual_path = fs.join(residuals_path, "disk_residual.fits")
        # disk_residual.save(disk_residual_path)

        # Calculate the model residual frame
        # model_residual = frame - model
        # model_residual = frame - (bulge*1.3)
        # model_residual = model_residual - disk
        # model_residual_path = fs.join(residuals_path, "model_residual.fits")
        # model_residual.save(model_residual_path)

        # Inform the user
        log.info("Plotting ...")

        # Create the image plotter
        plotter = ResidualImageGridPlotter()

        plotter.add_row(frame, bulge, "Bulge")
        plotter.add_row(frame, disk, "Disk")
        plotter.add_row(frame, model, "Bulge+disk")

        plotter.set_title("Decomposition")

        # Determine the path to the plot file
        path = fs.join(self.plot_path, "decomposition.pdf")

        plotter.absolute = True
        plotter.colormap = "hot"

        # Make the plot
        plotter.run(path)

# -----------------------------------------------------------------
