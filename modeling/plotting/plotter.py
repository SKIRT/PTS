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
from ...magic.plot.imagegrid import StandardImageGridPlotter, ResidualImageGridPlotter

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

        #frame.save(fs.join(self.truncation_path, "i1_jy.fits"))

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
        #bulge_residual = frame - bulge
        #bulge_residual_path = fs.join(residuals_path, "bulge_residual.fits")
        #bulge_residual.save(bulge_residual_path)

        # Calculate the disk residual frame
        #disk_residual = frame - disk
        #disk_residual_path = fs.join(residuals_path, "disk_residual.fits")
        #disk_residual.save(disk_residual_path)

        # Calculate the model residual frame
        #model_residual = frame - model
        # model_residual = frame - (bulge*1.3)
        #model_residual = model_residual - disk
        #model_residual_path = fs.join(residuals_path, "model_residual.fits")
        #model_residual.save(model_residual_path)

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
            image = Frame.from_file(path)

            # Add the image to the dictionary
            images[name] = image

        # Inform the user
        log.info("Plotting the images ...")

        # Create the image plotter
        plotter = StandardImageGridPlotter()

        # Sort the image labels based on wavelength
        sorted_labels = sorted(images.keys(), key=(lambda key: images[key].filter.pivotwavelength() if images[key].filter is not None else None))

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

    def plot_fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
