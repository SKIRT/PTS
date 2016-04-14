#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.residuals Contains the ResidualAnalyser class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools import filesystem
from ...core.tools.logging import log
from ...magic.core.frame import Frame

# -----------------------------------------------------------------

class ResidualAnalyser(AnalysisComponent):
    
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
        super(ResidualAnalyser, self).__init__(config)

        # The simulated images
        self.simulated = dict()

        # The observed images
        self.observed = dict()

        # The residual images
        self.residuals = dict()

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new ResidualAnalyser instance
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

        # 2. Load the simulated images
        self.load_simulated_images()

        # 3. Load the observed images
        self.load_observations()

        # 4. Calculate the residual images
        self.calculate_residuals()

        # 5. Writing
        self.write()

    # -----------------------------------------------------------------

    def load_simulated_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated images ...")

        # Loop over all FITS files found in the analysis/misc directory
        for path, name in filesystem.files_in_path(self.analysis_misc_path, extension="fits", returns=["path", "name"], contains="__"):

            # Debugging
            log.debug("Loading the '" + name + "' image ...")

            # Get the filter name
            filter_name = name.split("__")[1]

            # Open the image
            frame = Frame.from_file(path)

            # Add the image frame to the dictionary
            self.simulated[filter_name] = frame

    # -----------------------------------------------------------------

    def load_observed_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed images ...")

        # Loop over all FITS files found in the 'truncated' directory
        for path, name in filesystem.files_in_path(self.truncation_path, extension="fits", returns=["path", "name"]):

            # Ignore the bulge, disk and model images
            if name == "bulge" or name == "disk" or name == "model": continue

            # Ignore the H alpha image
            if "Halpha" in name: continue

            # Debugging
            log.debug("Loading the '" + name + "' image ...")

            # The filter name is the image name
            filter_name = name

            # Open the image
            frame = Frame.from_file(path)

            # Add the image frame to the dictionary
            self.observed[filter_name] = frame

    # -----------------------------------------------------------------

    def calculate_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the residual images ...")

        # Get the filter names which appear in both the simulated and observed images
        filter_names = list(set(self.simulated.keys() + self.observed.keys()))

        # Loop over the filter names
        for filter_name in filter_names:

            # Calculate the residual image
            residual = (self.simulated[filter_name] - self.observed[filter_name]) / self.observed[filter_name]
            #residual.replace_infs(0.0)

            # Add the residual image to the dictionary
            self.residuals[filter_name] = residual

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Write the residual frames
        self.write_residuals()

    # -----------------------------------------------------------------

    def write_residuals(self):

        """
        This function ...
        :return:
        """

        # Loop over the residual frames
        for filter_name in self.residuals:

            # Determine the path for this residual image
            path = filesystem.join(self.analysis_residuals_path, filter_name + ".fits")

            # Write the image
            self.residuals[filter_name].save(path)

# -----------------------------------------------------------------
