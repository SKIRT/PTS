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
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.core.frame import Frame
from ...magic.plot.imagegrid import ImageGridPlotter
from ...magic.basics.skyregion import SkyRegion

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

        # -- Attributes --

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
        self.load_observed_images()

        # 4. Rebin the images to the same pixel grid
        self.rebin()

        # Load the truncation ellipse
        self.load_truncation_ellipse()

        # 5. Calculate the residual images
        #self.calculate_residuals()

        # 6. Writing
        #self.write()

        # 7. Plotting
        self.plot()

    # -----------------------------------------------------------------

    def load_simulated_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated images ...")

        # Loop over all FITS files found in the analysis/misc directory
        for path, name in fs.files_in_path(self.analysis_misc_path, extension="fits", returns=["path", "name"], contains="__"):

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
        for path, name in fs.files_in_path(self.truncation_path, extension="fits", returns=["path", "name"]):

            # Ignore the bulge, disk and model images
            if name == "bulge" or name == "disk" or name == "model": continue

            # Ignore the H alpha image
            if "Halpha" in name: continue

            # Check whether a simulated image exists for this band
            if name not in self.simulated:
                log.warning("The simulated version of the " + name + " image could not be found, skipping " + name + " data ...")
                continue

            # Debugging
            log.debug("Loading the '" + name + "' image ...")

            # The filter name is the image name
            filter_name = name

            # Open the image
            frame = Frame.from_file(path)

            # Add the image frame to the dictionary
            self.observed[filter_name] = frame

    # -----------------------------------------------------------------

    def rebin(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning the observed and simulated images to the same resolution ...")

        # Loop over the filter names
        for filter_name in self.filter_names:

            simulated = self.simulated[filter_name]
            observed = self.observed[filter_name]

            # Check whether the coordinate systems of the observed and simulated image match
            if simulated.wcs == observed.wcs:
                # Debugging
                log.debug("The coordinate system of the simulated and observed image for the " + filter_name + " filter matches")
                continue

            # Debugging
            log.debug("The coordinate system of the simulated and observed image for the " + filter_name + " does not match: rebinning the simulated image ...")

            # Rebin the simulated image to the coordinate system of the observed image
            simulated_rebinned = simulated.rebinned(observed.wcs)

            # Replace the simulated frame by the rebinned frame
            self.simulated[filter_name] = simulated_rebinned

    # -----------------------------------------------------------------

    def load_truncation_ellipse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ellipse region used for truncating the observed images ...")

        # Determine the path
        path = fs.join(self.truncation_path, "ellipse.reg")

        # Get the ellipse
        region = SkyRegion.from_file(path)
        self.ellipse = region[0]

    # -----------------------------------------------------------------

    def calculate_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the residual images ...")

        # Loop over the filter names
        for filter_name in self.filter_names:

            simulated = self.simulated[filter_name]
            observed = self.observed[filter_name]

            # Calculate the residual image
            residual = (simulated - observed) / observed

            #residual.replace_infs(0.0)

            # Add the residual image to the dictionary
            self.residuals[filter_name] = residual

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the residual frames
        self.write_residuals()

    # -----------------------------------------------------------------

    def write_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the residual frames ...")

        # Loop over the residual frames
        for filter_name in self.residuals:

            # Determine the path for this residual image
            path = fs.join(self.analysis_residuals_path, filter_name + ".fits")

            # Debugging
            log.debug("Writing the residual frame for the " + filter_name + " band to '" + path + "' ...")

            # Write the image
            self.residuals[filter_name].save(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot a grid with the observed, simulated and residual images
        self.plot_image_grid()

    # -----------------------------------------------------------------

    def plot_image_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a grid with the observed, simulated and residual images ...")

        # Create an ImageGridPlotter instance
        plotter = ImageGridPlotter(title="Image residuals")

        # Loop over the filter names, add a row to the image grid plotter for each filter
        for filter_name in self.filter_names_sorted:

            observed = self.observed[filter_name]
            simulated = self.simulated[filter_name]

            plotter.add_row(observed, simulated, filter_name)

        # Set the bounding box for the plotter
        plotter.set_bounding_box(self.ellipse.bounding_box)

        # Determine the path to the plot file
        path = fs.join(self.analysis_residuals_path, "residuals.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    @property
    def filter_names(self):

        """
        This function ...
        :return:
        """

        # Get the filter names which appear in both the simulated and observed images
        filter_names = list(set(self.simulated.keys() + self.observed.keys()))
        return filter_names

    # -----------------------------------------------------------------

    @property
    def filter_names_sorted(self):

        """
        This function returns a list of the filter names, sorted on wavelength
        :return:
        """

        return sorted(self.filter_names, key=lambda key: self.observed[key].filter.pivotwavelength())

# -----------------------------------------------------------------
