#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.colours.colours Contains the ColourAnalyser class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import ColourAnalysisComponent
from ....core.tools.logging import log
from ....core.tools import filesystem as fs
from ....magic.core.frame import Frame
from ....magic.plot.imagegrid import ResidualImageGridPlotter
from ....core.basics.distribution import Distribution
from ....core.basics.filter import Filter

# -----------------------------------------------------------------

# Names that identify the interesting wavelengths/filters
colour_names = ["70/100", "100/160", "160/250", "250/350", "350/500"]
keys = ["70", "100", "160", "250", "350", "500"]
ids = {"Pacs blue": "70", "Pacs green": "100", "Pacs red": "160", "SPIRE PSW": "250", "SPIRE PMW": "350", "SPIRE PLW": "500"}
filter_names = {"70": "Pacs blue", "100": "Pacs green", "160": "Pacs red", "250": "SPIRE PSW", "350": "SPIRE PMW", "500": "SPIRE PLW"}

# -----------------------------------------------------------------

class ColourAnalyser(ColourAnalysisComponent):
    
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
        super(ColourAnalyser, self).__init__(config)

        # -- Attributes --

        # The dictionaries that contain the observed and simulated far-infrared images
        self.observed = dict()
        self.simulated = dict()

        # The dictionaries that contain the observed and simulated colour maps
        self.observed_colours = dict()
        self.simulated_colours = dict()

        # The residual colour maps
        self.residual_colours = dict()

        # The residual pixel value distributions
        self.residual_distributions = dict()

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

        # 2. Rebin the images to the same pixel grid
        self.rebin()

        # 3. Calculate the colour maps
        self.calculate_colours()

        # 4. Calculate the residual colour maps
        self.calculate_residuals()

        # 5. Writing
        self.write()

        # 6. Plotting
        self.plot()

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images ...")

        # 2. Load the observed images
        self.load_observed_images()

        # 3. Load the simulated images
        self.load_simulated_images()

    # -----------------------------------------------------------------

    def load_observed_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed images ...")

        # Loop over the appropriate observed images
        for key in keys:

            # Get the corresponding filter
            fltr = Filter(filter_names[key])

            # Get the observed image for this filter
            frame = self.dataset.get_frame_for_filter(fltr)

            # Add the frame to the dictionary
            self.observed[key] = frame

    # -----------------------------------------------------------------

    def load_simulated_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated images ...")

        # Loop over the appropriate simulated images
        for key in keys:

            # Determine the path to the image
            filter_name = filter_names[key]
            image_name = self.galaxy_name + "_earth_total__" + filter_name
            path = fs.join(self.analysis_run.misc_path, image_name + ".fits")

            # Check whether the image is present
            if not fs.is_file(path):
                log.warning("Simulated " + key + " micron image is not present")
                continue

            # Load the image frame
            frame = Frame.from_file(path)

            # Add the frame to the dictionary
            self.simulated[key] = frame

    # -----------------------------------------------------------------

    def rebin(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning all images to the pixel grid of the observed image with the lowest resolution ...")

        # Set the target coordinate system to None initially
        target_wcs = None

        # Get the WCS of the observed image with the largest pixelscale
        for key in self.observed:
            if target_wcs is None or self.observed[key].wcs.average_pixelscale > target_wcs.average_pixelscale: target_wcs = self.observed[key].wcs

        # Debugging
        log.debug("The target coordinate system has a pixelscale of " + str(target_wcs.average_pixelscale) + " and a shape of (nx=" + str(target_wcs.naxis1) + ", ny=" + str(target_wcs.naxis2) + ")")

        # Rebin all observed images to the target coordinate system
        for key in self.observed:

            # Debugging
            log.debug("Checking observed images ...")

            # Check whether the coordinate systems match
            if self.observed[key].wcs == target_wcs:

                # Debugging
                log.debug("The coordinate system of the observed " + key + " micron image matches the target coordinate system")

            else:

                # Debugging
                log.debug("The coordinate system of the observed " + key + " micron image does not match the target coordinate system: rebinning ...")

                # Rebin the image to the specific coordinate system
                self.observed[key] = self.observed[key].rebinned(target_wcs)

        # Rebin all simulated images to the target coordinate system
        for key in self.simulated:

            # Check whether the coordinate systems match
            if self.simulated[key].wcs == target_wcs:

                # Debugging
                log.debug("The coordinate system of the simulated " + key + " micron image matches the target coordinate system")

            else:

                # Debugging
                log.debug("The coordinate system of the simulated " + key + " micron image does not match the target coordinate system: rebinning ...")

                # Rebin the image to the standardized coordinate system
                self.simulated[key] = self.simulated[key].rebinned(target_wcs)

    # -----------------------------------------------------------------

    def calculate_colours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the colour maps")

        # Loop over the different colours
        for colour_name in colour_names:

            # Get the wavelengths to create the colour map
            wavelength1, wavelength2 = colour_name.split("/")

            # Check whether the corresponding observed and simulated images exist
            if not wavelength1 in self.observed:
                log.warning("Observed " + filter_names[wavelength1] + " not present, " + colour_name + " can not be calculated")
                continue
            if not wavelength1 in self.simulated:
                log.warning("Simulated " + filter_names[wavelength1] + " not present, " + colour_name + " can not be calculated")
                continue
            if not wavelength2 in self.observed:
                log.warning("Observed " + filter_names[wavelength2] + " not present, " + colour_name + " can not be calculated")
                continue
            if not wavelength2 in self.simulated:
                log.warning("Simulated " + filter_names[wavelength2] + " not present, " + colour_name + " can not be calculated")
                continue

            # Debugging
            log.debug("Calculating the observed and simulated " + colour_name + " colour maps ...")

            # Calculate the colour maps
            observed_colour = np.log10(self.observed[wavelength1]/self.observed[wavelength2])
            simulated_colour = np.log10(self.simulated[wavelength1]/self.simulated[wavelength2])

            # Add the colour maps to the appropriate dictionary
            self.observed_colours[colour_name] = observed_colour
            self.simulated_colours[colour_name] = simulated_colour

    # -----------------------------------------------------------------

    def calculate_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating residuals between the observed and simulated colour maps ...")

        # Loop over the different colours
        for colour_name in self.observed_colours:

            # Debugging
            log.debug("Calculating residuals between the observed and simulated " + colour_name + " colour maps ...")

            # Calculate the residual
            residual = (self.observed_colours[colour_name] - self.simulated_colours[colour_name])/self.observed_colours[colour_name]

            # Add the residual map to the dictionary
            self.residual_colours[colour_name] = residual

    # -----------------------------------------------------------------

    def calculate_residual_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating distributions of residual pixel values ...")

        # Loop over the different colours
        for colour_name in self.observed_colours:

            # Debugging
            log.debug("Calculating the distribution for the pixels of the " + colour_name + " residual map ...")

            # Get an 1D array of the valid pixel values
            pixel_values = None

            # Create the distribution
            distribution = Distribution.from_values(pixel_values)

            # Debugging
            #log.debug("Median " + colour_name + " residual: " + str(np.nanmedian(np.abs(residual))))
            #log.debug("Standard deviation of " + colour_name + " residual: " + str(np.nanstd(residual)))

            # Add the distribution to the dictionary
            self.residual_distributions[colour_name] = distribution

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the observed colour maps
        self.write_observed_colours()

        # Write the simulated colour maps
        self.write_simulated_colours()

        # Write the residual colour maps
        self.write_residual_colours()

        # Write the distributions of the residual pixel values
        self.write_residual_distributions()

    # -----------------------------------------------------------------

    def write_observed_colours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the observed colour maps ...")

        # Loop over the observed colour maps
        for colour_name in self.observed_colours:

            # Determine the path
            path = fs.join(self.colours_observed_path, colour_name.replace("/", "-") + ".fits")

            # Save the image
            self.observed_colours[colour_name].saveto(path)

    # -----------------------------------------------------------------

    def write_simulated_colours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the simulated colour maps ...")

        # Loop over the simulated colour maps
        for colour_name in self.simulated_colours:

            # Determine the path
            path = fs.join(self.colours_simulated_path, colour_name.replace("/", "-") + ".fits")

            # Save the image
            self.simulated_colours[colour_name].saveto(path)

    # -----------------------------------------------------------------

    def write_residual_colours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the residual colour maps ...")

        # Loop over the residual colour maps
        for colour_name in self.residual_colours:

            # Determine the path
            path = fs.join(self.colours_residuals_path, colour_name.replace("/", "-") + ".fits")

            # Save the image
            self.residual_colours[colour_name].saveto(path)

    # -----------------------------------------------------------------

    def write_residual_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the residual colour pixel distributions ...")

        # Loop over the residual pixel distributions
        for colour_name in self.residual_distributions:

            # Determine the path
            path = fs.join()

            # Save the distribution data
            self.residual_distributions[colour_name].saveto(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot a grid with the observed, simulated and residual colour maps
        self.plot_image_grid()

        # Plot the distributions of the pixel values of the residual maps
        self.plot_residual_distributions()

    # -----------------------------------------------------------------

    def plot_image_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a grid with the observed, simulated and residual colour maps ...")

        # Create the image grid plotter
        plotter = ResidualImageGridPlotter(title="Colours")

        # Add the rows
        for colour_name in self.observed_colours:

            observed_colour = self.observed_colours[colour_name]
            simulated_colour = self.simulated_colours[colour_name]

            plotter.add_row(observed_colour, simulated_colour, colour_name)

        # Set the bounding box for the plotter
        plotter.set_bounding_box(self.truncation_box)

        # Determine the path to the plot file
        path = fs.join(self.analysis_run.colours_path, "colours.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_residual_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the pixel value distributions for the residual maps ...")

        # Loop over the residual pixel distributions
        for colour_name in self.residual_distributions:

            # Determine the path
            path = fs.join(self.colours_residuals_path, colour_name + ".pdf")

            # Save the distribution data
            self.residual_distributions[colour_name].plot(title="Distribution of residual pixel values for " + colour_name + " colour", path=path)

# -----------------------------------------------------------------
