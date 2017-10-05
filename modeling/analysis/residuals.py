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
from ...core.basics.log import log
from ...magic.core.frame import Frame
from ...magic.plot.imagegrid import ResidualImageGridPlotter
from ...core.tools import sequences
from ...core.basics.distribution import Distribution

# -----------------------------------------------------------------

class ResidualAnalyser(AnalysisComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(ResidualAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The analysis run
        self.analysis_run = None

        # The simulated images
        self.simulated = dict()

        # The observed images
        self.observed = dict()

        # The error maps
        self.errors = dict()

        # The residual images
        self.residuals = dict()

        # The weighed residual images
        self.weighed = dict()

        # The distributions of residual values
        self.residual_distributions = dict()
        self.weighed_distributions = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the observed and simulated images
        self.load_images()

        # 3. Calculate the residual images
        self.calculate_residuals()

        # 4. Calculate the weighed residual images
        self.calculate_weighed()

        # 5. Create distributions of the residual values
        self.create_distributions()

        # 6. Writing
        self.write()

        # 7. Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ResidualAnalyser, self).setup(**kwargs)

        # Load the analysis run
        self.load_run()

    # -----------------------------------------------------------------

    def load_run(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the analysis run " + self.config.run + " ...")

        # Get the run
        self.analysis_run = self.get_run(self.config.run)

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images ...")

        # Observed images
        self.load_observed_images()

        # Simulated images
        self.load_simulated_images()

        # Error maps
        self.load_errors()

    # -----------------------------------------------------------------

    def load_observed_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed images ...")

        # Load all data
        for name in self.dataset.names:

            # Debugging
            log.debug("Loading the observed " + name + " image ...")

            # Load the frame, not truncated
            frame = self.dataset.get_frame(name, masked=False)

            # Define the name of the filter as the image name
            image_name = str(frame.filter)

            # Add the image frame to the dictionary
            self.observed[image_name] = frame

    # -----------------------------------------------------------------

    def load_simulated_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated images ...")

        # Loop over all FITS files found in the analysis run misc directory
        for path, name in fs.files_in_path(self.analysis_run.misc_path, extension="fits", returns=["path", "name"], contains="__"):

            # Debugging
            log.debug("Loading the simulated " + name + " image ...")

            # Open the image
            frame = Frame.from_file(path)

            # Get the filter name
            image_name = str(frame.filter)

            # Check whether a simulated image exists for this band
            if image_name not in self.observed:
                log.warning("The observed " + image_name + " image could not be found, skipping simulated " + image_name + " image ...")
                continue

            # Get the corresponding observed image
            observed = self.observed[image_name]

            # Check whether the coordinate systems of the observed and simulated image match
            if frame.wcs == observed.wcs: log.debug("The coordinate system of the simulated and observed image for the " + image_name + " filter matches")

            # The observed image has a smaller pixelscale as the simulated image -> rebin the observed image
            elif observed.average_pixelscale < frame.average_pixelscale:

                # Debugging
                log.debug("The observed image has a better resolution as the simulated image: rebinning the observed image ...")

                # Rebin
                observed.rebin(frame.wcs)

            # The simulated image has a smaller pixelscale as the observed image
            elif frame.average_pixelscale < observed.average_pixelscale:

                # Debugging
                log.debug("The simulated image has a better resolution as the observed image: rebinning the simulated image ...")

                # Rebin the simulated image to the coordinate system of the observed image
                frame.rebin(observed.wcs)

            # Error
            else: raise RuntimeError("Something unexpected happened")

            # Add the simulated image frame to the dictionary
            self.simulated[image_name] = frame

    # -----------------------------------------------------------------

    def load_errors(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the error maps ...")

        # Load all data
        for name in self.dataset.names:

            # Debugging
            log.debug("Loading the " + name + " error map ...")

            # Load the error map, not truncated
            errors = self.dataset.get_errormap(name, masked=False)

            # Define the name of the filter as the image name
            image_name = str(errors.filter)

            # Get the corresponding observed image
            observed = self.observed[image_name]

            # CHECK THE WCS:
            if errors.wcs != observed.wcs:

                # Debugging
                log.debug("The coordinate system of the error map is not identical to the coordinate system of the observed image: rebinning the error map ...")

                # Rebin
                errors.rebin(observed.wcs)

            # Add the error map to the dictionary
            self.errors[image_name] = errors

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

            # Get the observed and simulated image
            simulated = self.simulated[filter_name]
            observed = self.observed[filter_name]
            errors = self.errors[filter_name]

            # Calculate the residual image
            residual = (simulated - observed) / observed

            # Replace infs
            residual.replace_infs(0.0)

            # Get the truncation mask
            truncation_mask = self.get_truncation_mask(observed.wcs)

            # Get the significance mask
            significance_mask = self.get_significance_mask(observed, errors)

            # MASK
            residual[truncation_mask] = 0.0
            residual[significance_mask] = 0.0

            # Add the residual image to the dictionary
            self.residuals[filter_name] = residual

    # -----------------------------------------------------------------

    def calculate_weighed(self):

        """
        This ufnction ...
        :return:
        """

        # Inform the user
        log.info("Calculating the weighed residual images ...")

        # Loop over the filter names
        for filter_name in self.filter_names:

            # Get the observed image, simulated image and error map
            simulated = self.simulated[filter_name]
            observed = self.observed[filter_name]
            errors = self.errors[filter_name]

            # Calculate the weighed residual image
            residual = (simulated - observed) / errors

            # Replace infs
            residual.replace_infs(0.0)

            # Get the truncation mask
            truncation_mask = self.get_truncation_mask(observed.wcs)

            # Get the significance mask
            significance_mask = self.get_significance_mask(observed, errors)

            # MASK
            residual[truncation_mask] = 0.0
            residual[significance_mask] = 0.0

            # Add the weighed residual image to the dictionary
            self.weighed[filter_name] = residual

    # -----------------------------------------------------------------

    def create_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the distributions ...")

        # Residuals
        self.create_residuals_distributions()

        # Weighed residuals
        self.create_weighed_distributions()

    # -----------------------------------------------------------------

    def create_residuals_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the residuals distributions ...")

        # Loop over the filter names
        for filter_name in self.filter_names:

            # Get the values within the truncation ellipse
            values = self.residuals[filter_name].values_in(self.truncation_ellipse)

            # Create distribution
            distribution = Distribution.from_values(values, bins=self.config.nbins)

            # Add the distribution
            self.residual_distributions[filter_name] = distribution

    # -----------------------------------------------------------------

    def create_weighed_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the weighed residuals distributions ...")

        # Loop over the filter names
        for filter_name in self.filter_names:

            # Get the values within the truncation ellipse
            values = self.weighed[filter_name].values_in(self.truncation_ellipse)

            # Create distribution
            distribution = Distribution.from_values(values, bins=self.config.nbins)

            # Add the distribution
            self.weighed_distributions[filter_name] = distribution

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

        # Write the weighed residual frames
        self.write_weighed()

        # Write distributions
        self.write_distributions()

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
            path = fs.join(self.analysis_run.residuals_path, filter_name + ".fits")

            # Debugging
            log.debug("Writing the residual frame for the '" + filter_name + "' band to '" + path + "' ...")

            # Write the image
            self.residuals[filter_name].saveto(path)

    # -----------------------------------------------------------------

    def write_weighed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the weighed residual frames ...")

        # Loop over the residual frames
        for filter_name in self.weighed:

            # Determine the path for this residual image
            path = fs.join(self.analysis_run.weighed_residuals_path, filter_name + ".fits")

            # Debugging
            log.debug("Writing the weighed residual frame for the '" + filter_name + "' band to '" + path + "' ...")

            # Write
            self.weighed[filter_name].saveto(path)

    # -----------------------------------------------------------------

    def write_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the distributions ...")

        # Residuals
        self.write_residuals_distributions()

        # Weighed residuals
        self.write_weighed_distributions()

    # -----------------------------------------------------------------

    def write_residuals_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing residuals distributions ...")

        # Loop over the filter names
        for filter_name in self.filter_names:

            # Determine plot path
            path = fs.join(self.analysis_run.residuals_path, filter_name + ".dat")

            # Debugging
            log.debug("Writing the residuals distributions for the '" + filter_name + "' band to '" + path + "' ...")

            # Write
            self.residual_distributions[filter_name].saveto(path)

    # -----------------------------------------------------------------

    def write_weighed_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing weighed residuals distributions ...")

        # Loop over the filter names
        for filter_name in self.filter_names:

            # Determine plot path
            path = fs.join(self.analysis_run.weighed_residuals_path, filter_name + ".dat")

            # Debugging
            log.debug("Writing the weighed residuals distributions for the '" + filter_name + "' band to '" + path + "' ...")

            # Write
            self.weighed_distributions[filter_name].saveto(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the distributions
        self.plot_distributions()

        # Plot a grid with the observed, simulated and residual images
        self.plot_image_grid()

        # Plot a grid with the observed, simulated and weighed residual images
        self.plot_image_grid_weighed()

    # -----------------------------------------------------------------

    def plot_image_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a grid with the observed, simulated and residual images ...")

        # Create the image grid plotter
        plotter = ResidualImageGridPlotter(title="Image residuals")

        # Loop over the filter names, add a row to the image grid plotter for each filter
        for filter_name in self.filter_names_sorted:

            # Get input
            observed = self.observed[filter_name]
            simulated = self.simulated[filter_name]

            # Add row
            plotter.add_row(observed, simulated, filter_name)

        # Set the bounding box for the plotter
        plotter.set_bounding_box(self.truncation_box)

        # Determine the path to the plot file
        path = fs.join(self.analysis_run.residuals_path, "residuals.pdf")

        # Run the plotter
        plotter.run(path)

    # -----------------------------------------------------------------

    def plot_image_grid_weighed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting a grid with the observed, simulated and weighed residual images ...")

        # Create the image grid plotter
        plotter = ResidualImageGridPlotter(title="Weighed image residuals", weighed=True)

        # Loop over the filter names, add a row to the image grid plotter for each filter
        for filter_name in self.filter_names_sorted:

            # Get input
            observed = self.observed[filter_name]
            simulated = self.simulated[filter_name]
            errors = self.errors[filter_name]

            # Add row
            plotter.add_row(observed, simulated, filter_name, errors=errors)

        # Set the bounding box for the plotter
        plotter.set_bounding_box(self.truncation_box)

        # Determine the path to the plot file
        path = fs.join(self.analysis_run.weighed_residuals_path, "weighed_residuals.pdf")

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
        return sequences.intersection(self.simulated.keys(), self.observed.keys())

    # -----------------------------------------------------------------

    @property
    def filter_names_sorted(self):

        """
        This function returns a list of the filter names, sorted on wavelength
        :return:
        """

        return sorted(self.filter_names, key=lambda key: self.observed[key].filter.pivotwavelength())

# -----------------------------------------------------------------
