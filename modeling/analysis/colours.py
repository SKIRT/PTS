#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.colours Contains the ColourAnalyser class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...magic.core.frame import Frame
from ...magic.plot.imagegrid import ImageGridPlotter

# -----------------------------------------------------------------

class ColourAnalyser(AnalysisComponent):
    
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

        # Names that identify the interesting wavelengths/filters
        self.colour_names = ["70/100", "100/160", "160/250", "250/350", "350/500"]
        self.keys = ["70", "100", "160", "250", "350", "500"]
        self.ids = {"Pacs blue": "70", "Pacs green": "100", "Pacs red": "160", "SPIRE PSW": "250", "SPIRE PMW": "350", "SPIRE PLW": "500"}
        self.filter_names = {"70": "Pacs blue", "100": "Pacs green", "160": "Pacs red", "250": "SPIRE PSW", "350": "SPIRE PMW", "500": "SPIRE PLW"}

        # The dictionaries that contain the observed and simulated far-infrared images
        self.observed = dict()
        self.simulated = dict()

        # The dictionaries that contain the observed and simulated colour maps
        self.observed_colours = dict()
        self.simulated_colours = dict()

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new ColourAnalyser instance
        analyser = cls()

        # Set the modeling path
        analyser.config.path = arguments.path

        # Return the new instance
        return analyser

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the observed images
        self.load_observed_images()

        # 3. Load the simulated images
        self.load_simulated_images()

        # 4. Rebin the images to the same pixel grid
        self.rebin()

        # 4. Calculate the colour maps
        self.calculate_colours()

        # 5. Writing
        self.write()

        # 6. Plotting
        self.plot()

    # -----------------------------------------------------------------

    def load_observed_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed images ...")

        # Loop over the appropriate observed images
        for key in self.keys:

            # Determine the path to the image
            filter_name = self.filter_names[key]
            path = fs.join(self.truncation_path, filter_name + ".fits")

            # Check whether the image is present
            if not fs.is_file(path):
                log.warning("Observed " + filter_name + " image is not present")
                continue

            # Load the image frame
            frame = Frame.from_file(path)

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
        for key in self.keys:

            # Determine the path to the image
            filter_name = self.filter_names[key]
            image_name = self.galaxy_name + "_earth_total__" + filter_name
            path = fs.join(self.analysis_misc_path, image_name + ".fits")

            # Check whether the image is present
            if not fs.is_file(path):
                log.warning("Simulated " + filter_name + " image is not present")
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
            if target_wcs is None or self.observed[key].wcs.xy_average_pixelscale > target_wcs.xy_average_pixelscale: target_wcs = self.observed[key].wcs

        # Debugging
        log.debug("The target coordinate system has a pixelscale of " + str(target_wcs.xy_average_pixelscale) + " and a shape of (nx=" + str(target_wcs.naxis1) + ", ny=" + str(target_wcs.naxis2) + ")")

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

        # Loop over the different colours
        for colour_name in self.colour_names:

            # Get the wavelengths to create the colour map
            wavelength1, wavelength2 = colour_name.split("/")

            # Check whether the corresponding observed and simulated images exist
            if not wavelength1 in self.observed:
                log.warning("Observed " + self.filter_names[wavelength1] + " not present, " + colour_name + " can not be calculated")
                continue
            if not wavelength1 in self.simulated:
                log.warning("Simulated " + self.filter_names[wavelength1] + " not present, " + colour_name + " can not be calculated")
                continue
            if not wavelength2 in self.observed:
                log.warning("Observed " + self.filter_names[wavelength2] + " not present, " + colour_name + " can not be calculated")
                continue
            if not wavelength2 in self.simulated:
                log.warning("Simulated " + self.filter_names[wavelength2] + " not present, " + colour_name + " can not be calculated")
                continue

            # Calculate the colour maps
            observed_colour = np.log10(self.observed[wavelength1]/self.observed[wavelength2])
            simulated_colour = np.log10(self.simulated[wavelength1]/self.simulated[wavelength2])

            # Add the colour maps to the appropriate dictionary
            self.observed_colours[colour_name] = observed_colour
            self.simulated_colours[colour_name] = simulated_colour

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Directories
        observed_path = fs.join(self.analysis_colours_path, "observed")
        simulated_path = fs.join(self.analysis_colours_path, "simulated")
        fs.create_directory(observed_path)
        fs.create_directory(simulated_path)

        # Loop over the observed colour maps
        for colour_name in self.observed_colours:

            # Determine the path
            path = fs.join(observed_path, colour_name.replace("/", "-") + ".fits")

            # Save the image
            self.observed_colours[colour_name].save(path)

        # Loop over the simulated colour maps
        for colour_name in self.simulated_colours:

            # Determine the path
            path = fs.join(simulated_path, colour_name.replace("/", "-") + ".fits")

            # Save the image
            self.simulated_colours[colour_name].save(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Create an ImageGridPlotter instance
        plotter = ImageGridPlotter(title="Colours")

        # Add the rows
        for colour_name in self.observed_colours:

            observed_colour = self.observed_colours[colour_name]
            simulated_colour = self.simulated_colours[colour_name]
            plotter.add_row(observed_colour, simulated_colour, colour_name)

        # Determine the path to the plot file
        path = fs.join(self.analysis_colours_path, "colours.pdf")

        # Run the plotter
        plotter.run(path)

# -----------------------------------------------------------------
