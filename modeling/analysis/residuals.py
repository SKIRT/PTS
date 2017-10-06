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

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...magic.core.frame import Frame
from ...magic.plot.imagegrid import ResidualImageGridPlotter
from ...core.tools import sequences
from ...core.basics.distribution import Distribution
from ...core.tools.utils import lazyproperty
from ...core.filter.filter import parse_filter

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

        # 1. Observed images
        self.load_observed_images()

        # 2. Simulated images
        self.load_simulated_images()

        # 3. Error maps
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

            # Get filter name
            filter_name = self.dataset.get_filter_name(name) # should be the same?

            # Check whether residual frames already created
            if self.has_residuals(filter_name) and self.has_weighed_residuals(filter_name):
                log.success("Residual and weighed residual map for the '" + filter_name + "' filter already created")
                continue

            # Debugging
            log.debug("Loading the observed " + name + " image ...")

            # Load the frame, not truncated
            frame = self.dataset.get_frame(name, masked=False)

            # Add the image frame to the dictionary
            self.observed[filter_name] = frame

    # -----------------------------------------------------------------

    def load_simulated_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated images ...")

        # Loop over all FITS files found in the analysis run misc directory
        for path, image_name in fs.files_in_path(self.analysis_run.misc_path, extension="fits", returns=["path", "name"], contains="__"):

            # Get filter name
            filter_name = image_name.split("__")[1]

            # Check whether residual frames already creeated
            if self.has_residuals(filter_name) and self.has_weighed_residuals(filter_name):
                log.success("Residual and weighed residual map for the '" + filter_name + "' filter already created")
                continue

            # Debugging
            log.debug("Loading the simulated " + filter_name + " image ...")

            # Open the image
            frame = Frame.from_file(path)

            # Check whether a simulated image exists for this band
            if filter_name not in self.observed:
                log.warning("The observed " + filter_name + " image could not be found, skipping simulated " + filter_name + " image ...")
                continue

            # Get the corresponding observed image
            observed = self.observed[filter_name]

            # Check whether the coordinate systems of the observed and simulated image match
            if frame.wcs == observed.wcs: log.debug("The coordinate system of the simulated and observed image for the " + filter_name + " filter matches")

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
            self.simulated[filter_name] = frame

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

            # Get filter name
            filter_name = self.dataset.get_filter_name(name)  # should be the same?

            # Check whether residual frames already created
            if self.has_residuals(filter_name) and self.has_weighed_residuals(filter_name):
                log.success("Residual and weighed residual map for the '" + filter_name + "' filter already created")
                continue

            # Debugging
            log.debug("Loading the " + name + " error map ...")

            # Load the error map, not truncated
            errors = self.dataset.get_errormap(name, masked=False)

            # Get the corresponding observed image
            observed = self.observed[filter_name]

            # CHECK THE WCS:
            if errors.wcs != observed.wcs:

                # Debugging
                log.debug("The coordinate system of the error map is not identical to the coordinate system of the observed image: rebinning the error map ...")

                # Rebin
                errors.rebin(observed.wcs)

            # Add the error map to the dictionary
            self.errors[filter_name] = errors

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

            # Already have residuals
            if self.has_residuals(filter_name):
                log.success("Residual frame was already created for the '" + filter_name + " filter': loading from file ...")
                self.residuals[filter_name] = self.load_residuals(filter_name)
                continue

            # Debugging
            log.debug("Creating the residual frame for the '" + filter_name + "' filter ...")

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
            significance_mask = self.get_significance_mask(observed, errors, min_npixels=self.config.min_npixels, connectivity=self.config.connectivity)

            # MASK
            residual[truncation_mask] = 0.0
            residual[significance_mask] = 0.0

            # Add the residual image to the dictionary
            self.residuals[filter_name] = residual

            # WRITE THE SIGNIFICANCE MASK
            mask_path = fs.join(self.analysis_run.residuals_path, filter_name + "_significance.fits")
            significance_mask.saveto(mask_path)

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

            # Already have residuals
            if self.has_weighed_residuals(filter_name):
                log.success("Weighed residual frame was already created for the '" + filter_name + "' filter: loading from file ...")
                self.weighed[filter_name] = self.load_weighed_residuals(filter_name)
                continue

            # Debugging
            log.debug("Creating the weighed residual frame for the '" + filter_name + "' filter ...")

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
            significance_mask = self.get_significance_mask(observed, errors, min_npixels=self.config.min_npixels, connectivity=self.config.connectivity)

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

            # Already have distribution
            if self.has_residuals_distribution(filter_name):
                log.success("Residuals distribution was already created for the '" + filter_name + "' filter: loading from file ...")
                self.residual_distributions[filter_name] = self.load_residuals_distribution(filter_name)
                continue

            # Debugging
            log.debug("Creating the residuals distribution for the '" + filter_name + "' filter ...")

            # Get the values within the truncation ellipse
            values = self.residuals[filter_name].values_in(self.truncation_ellipse)

            # REMOVE EXACT ZEROES
            indices = np.argwhere(values == 0)
            values = np.delete(values, indices)

            # Check
            if len(values) == 0 or sequences.all_equal(values):
                log.error("Cannot create distribution for the '" + filter_name + "' filter")
                continue

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

            # Already have distribution
            if self.has_weighed_residuals_distribution(filter_name):
                log.success("Weighed residuals distribution was already created for the '" + filter_name + "' filter: loading from file ...")
                self.weighed_distributions[filter_name] = self.load_weighed_residuals_distribution(filter_name)
                continue

            # Debugging
            log.debug("Creating the weighed residuals distribution for the '" + filter_name + "' filter ...")

            # Get the values within the truncation ellipse
            values = self.weighed[filter_name].values_in(self.truncation_ellipse)

            # REMOVE EXACT ZEROES
            indices = np.argwhere(values == 0)
            values = np.delete(values, indices)

            # Check
            if len(values) == 0 or sequences.all_equal(values):
                log.error("Cannot create distribution for the '" + filter_name + "' filter")
                continue

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

    def has_residuals(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # File present?
        if fs.is_file(self.get_residuals_path(filter_name)):

            # REMAKE?
            if self.config.remake_residuals:

                # Remove files
                fs.remove_file(self.get_residuals_path(filter_name))
                if self.has_residuals_distribution(filter_name): fs.remove_file(self.get_residuals_distribution_path(filter_name))
                if self.has_residuals_distribution_plot(filter_name): fs.remove_file(self.get_residuals_distribution_plot_path(filter_name))
                if self.has_residuals_plot(filter_name): fs.remove_file(self.get_residuals_plot_path(filter_name))
                return False

            # Don't remake
            else: return True

        # No file
        else:

            # Check whether derived things are present, if so, they need to be removed
            if self.has_residuals_distribution(filter_name): fs.remove_file(self.get_residuals_distribution_path(filter_name))
            if self.has_residuals_distribution_plot(filter_name): fs.remove_file(self.get_residuals_distribution_plot_path(filter_name))
            if self.has_residuals_plot(filter_name): fs.remove_file(self.get_residuals_plot_path(filter_name))
            return False

    # -----------------------------------------------------------------

    def load_residuals(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        return Frame.from_file(self.get_residuals_path(filter_name))

    # -----------------------------------------------------------------

    def get_residuals_path(self, filter_name):

        """
        This function ...
        :return:
        """

        # Determine the path for this residual image
        path = fs.join(self.analysis_run.residuals_path, filter_name + ".fits")
        return path

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

            # Get path
            path = self.get_residuals_path(filter_name)
            if fs.is_file(path): continue

            # Debugging
            log.debug("Writing the residual frame for the '" + filter_name + "' band to '" + path + "' ...")

            # Write the image
            self.residuals[filter_name].saveto(path)

    # -----------------------------------------------------------------

    def has_weighed_residuals(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # File present?
        if fs.is_file(self.get_weighed_residuals_path(filter_name)):

            # REMAKE?
            if self.config.remake_weighed:

                # Remove files
                fs.remove_file(self.get_weighed_residuals_path(filter_name))
                if self.has_weighed_residuals_distribution(filter_name): fs.remove_file(self.get_weighed_residuals_distribution_path(filter_name))
                if self.has_weighed_residuals_distribution_plot(filter_name): fs.remove_file(self.get_weighed_residuals_distribution_plot_path(filter_name))
                if self.has_weighed_residuals_plot(filter_name): fs.remove_file(self.get_weighed_residuals_plot_path(filter_name))
                return False

            # Don't remake
            else: return True

        # No file
        else:

            # Check whether derived things are present, if so, remove them
            if self.has_weighed_residuals_distribution(filter_name): fs.remove_file(self.get_weighed_residuals_distribution_path(filter_name))
            if self.has_weighed_residuals_distribution_plot(filter_name): fs.remove_file(self.get_weighed_residuals_distribution_plot_path(filter_name))
            if self.has_weighed_residuals_plot(filter_name): fs.remove_file(self.get_weighed_residuals_plot_path(filter_name))
            return False

    # -----------------------------------------------------------------

    def load_weighed_residuals(self, filter_name):

        """
        Thisf unction ...
        :param filter_name:
        :return:
        """

        return Frame.from_file(self.get_weighed_residuals_path(filter_name))

    # -----------------------------------------------------------------

    def get_weighed_residuals_path(self, filter_name):

        """
        Thisf unction ...
        :param filter_name:
        :return:
        """

        # Determine the path for this residual image
        path = fs.join(self.analysis_run.weighed_residuals_path, filter_name + ".fits")
        return path

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

            # Get path
            path = self.get_weighed_residuals_path(filter_name)
            if fs.is_file(path): continue

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

    def has_residuals_distribution(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # File present?
        if fs.is_file(self.get_residuals_distribution_path(filter_name)):

            # Remake?
            if self.config.remake_distributions:

                # Remove files
                fs.remove_file(self.get_residuals_distribution_path(filter_name))
                if self.has_residuals_distribution_plot(filter_name): fs.remove_file(self.get_residuals_distribution_plot_path(filter_name))
                return False

            # Don't remake
            else: return True

        # No file
        else:

            # Check if derived things are present, if so, remove
            if self.has_residuals_distribution_plot(filter_name): fs.remove_file(self.get_residuals_distribution_plot_path(filter_name))
            return False

    # -----------------------------------------------------------------

    def load_residuals_distribution(self, filter_name):

        """
        Thisf ucntion ...
        :param filter_name:
        :return:
        """

        return Distribution.from_file(self.get_residuals_distribution_path(filter_name))

    # -----------------------------------------------------------------

    def get_residuals_distribution_path(self, filter_name):

        """
        Thisnfunction ...
        :param filter_name:
        :return:
        """

        # Determine data file path
        path = fs.join(self.analysis_run.residuals_path, filter_name + ".dat")
        return path

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

            # Get path
            path = self.get_residuals_distribution_path(filter_name)
            if fs.is_file(path): continue

            # Debugging
            log.debug("Writing the residuals distribution for the '" + filter_name + "' band to '" + path + "' ...")

            # Write
            self.residual_distributions[filter_name].saveto(path)

    # -----------------------------------------------------------------

    def has_weighed_residuals_distribution(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # File present?
        if fs.is_file(self.get_weighed_residuals_distribution_path(filter_name)):

            # Remake?
            if self.config.remake_weighed_distributions:

                # Remove files
                fs.remove_file(self.get_weighed_residuals_distribution_path(filter_name))
                if self.has_weighed_residuals_distribution_plot(filter_name): fs.remove_file(self.get_weighed_residuals_distribution_plot_path(filter_name))
                return False

            # Don't remake
            else: return True

        # No file
        else:

            # Check if derived things are present, if so, remove them
            if self.has_weighed_residuals_distribution_plot(filter_name): fs.remove_file(self.get_weighed_residuals_distribution_plot_path(filter_name))
            return False

    # -----------------------------------------------------------------

    def load_weighed_residuals_distribution(self, filter_name):

        """
        This function ....
        :param filter_name:
        :return:
        """

        return Distribution.from_file(self.get_weighed_residuals_distribution_path(filter_name))

    # -----------------------------------------------------------------

    def get_weighed_residuals_distribution_path(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # Determine data file path
        path = fs.join(self.analysis_run.weighed_residuals_path, filter_name + ".dat")
        return path

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

            # Get the path
            path = self.get_weighed_residuals_distribution_path(filter_name)
            if fs.is_file(path): continue

            # Debugging
            log.debug("Writing the weighed residuals distribution for the '" + filter_name + "' band to '" + path + "' ...")

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

        # Plot the residuals
        self.plot_residuals()

        # Plot the weighed rsiduals
        self.plot_weighed_residuals()

        # Plot a grid with the observed, simulated and residual images
        #self.plot_image_grid()

        # Plot a grid with the observed, simulated and weighed residual images
        #self.plot_image_grid_weighed()

    # -----------------------------------------------------------------

    def plot_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the distributions ...")

        # Residuals
        self.plot_residuals_distributions()

        # Weighed residuals
        self.plot_weighed_distributions()

    # -----------------------------------------------------------------

    def has_residuals_distribution_plot(self, filter_name):

        """
        Thisn function ...
        :param filter_name:
        :return:
        """

        # File present?
        if fs.is_file(self.get_residuals_distribution_plot_path(filter_name)):

            # Replot?
            if self.config.replot_distributions:

                # Remove file
                fs.remove_file(self.get_residuals_distribution_plot_path(filter_name))
                return False

            # Don't replot
            else: return True

        # No file
        else: return False

    # -----------------------------------------------------------------

    def get_residuals_distribution_plot_path(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # Determine plot path
        path = fs.join(self.analysis_run.residuals_path, filter_name + ".pdf")
        return path

    # -----------------------------------------------------------------

    def plot_residuals_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting residuals distributions ...")

        # Loop over the filter names
        for filter_name in self.filter_names:

            # Determine the path
            path = self.get_residuals_distribution_plot_path(filter_name)
            if fs.is_file(path): continue

            # Debugging
            log.debug("Plotting the residuals distribution for the '" + filter_name + "' band to '" + path + "' ...")

            # Plot
            self.residual_distributions[filter_name].plot(path=path)

    # -----------------------------------------------------------------

    def has_weighed_residuals_distribution_plot(self, filter_name):

        """
        Thisn function ...
        :param filter_name:
        :return:
        """

        # File present?
        if fs.is_file(self.get_weighed_residuals_distribution_plot_path(filter_name)):

            # Replot?
            if self.config.replot_weighed_distributions:

                # Remove file
                fs.remove_file(self.get_weighed_residuals_distribution_plot_path(filter_name))
                return False

            # Don't replot
            else: return True

        # No file
        else: return False

    # -----------------------------------------------------------------

    def get_weighed_residuals_distribution_plot_path(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # Determine plot path
        path = fs.join(self.analysis_run.weighed_residuals_path, filter_name + ".pdf")
        return path

    # -----------------------------------------------------------------

    def plot_weighed_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting weighed residuals distributions ...")

        # Loop over the filter names
        for filter_name in self.filter_names:

            # Get the path
            path = self.get_weighed_residuals_distribution_plot_path(filter_name)
            if fs.is_file(path): continue

            # Debugging
            log.debug("Plotting the weighed residuals distribution for the '" + filter_name + "' band to '" + path + "' ...")

            # Plot
            self.weighed_distributions[filter_name].plot(path=path)

    # -----------------------------------------------------------------

    def has_residuals_plot(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # File present?
        if fs.is_file(self.get_residuals_plot_path(filter_name)):

            # Replot
            if self.config.replot_residuals:

                # Remove file
                fs.remove_file(self.get_residuals_plot_path(filter_name))
                return False

            # Don't replot
            else: return True

        # No file
        else: return False

    # -----------------------------------------------------------------

    def get_residuals_plot_path(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        return fs.join(self.analysis_run.residuals_path, filter_name + ".png")

    # -----------------------------------------------------------------

    def plot_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the residuals ...")

        # Loop over the filter names
        for filter_name in self.filter_names:

            # Get the path
            path = self.get_residuals_plot_path(filter_name)
            if fs.is_file(path): continue

            # Save as PNG
            residuals = self.residuals[filter_name]
            vmin, vmax = residuals.saveto_png(path, colours=self.config.colours,
                                          interval=self.config.interval,
                                          scale=self.config.scale, alpha=self.config.alpha_method,
                                          peak_alpha=self.config.peak_alpha)

    # -----------------------------------------------------------------

    def has_weighed_residuals_plot(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # File present?
        if fs.is_file(self.get_weighed_residuals_plot_path(filter_name)):

            # Replot?
            if self.config.replot_weighed:

                # Remove file
                fs.remove_file(self.get_weighed_residuals_plot_path(filter_name))
                return False

            # Don't replot
            else: return True

        # No file
        else: return False

    # -----------------------------------------------------------------

    def get_weighed_residuals_plot_path(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        return fs.join(self.analysis_run.weighed_residuals_path, filter_name + ".png")

    # -----------------------------------------------------------------

    def plot_weighed_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the weighed residuals ...")

        # Loop over the filter names
        for filter_name in self.filter_names:

            # Get the path
            path = self.get_weighed_residuals_plot_path(filter_name)
            if fs.is_file(path): continue

            # Save as PNG
            residuals = self.weighed[filter_name]
            vmin, vmax = residuals.saveto_png(path, colours=self.config.colours,
                                              interval=self.config.interval,
                                              scale=self.config.scale, alpha=self.config.alpha_method,
                                              peak_alpha=self.config.peak_alpha)

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

    # @property
    # def filter_names(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Get the filter names which appear in both the simulated and observed images
    #     return sequences.intersection(self.simulated.keys(), self.observed.keys())

    # -----------------------------------------------------------------

    @lazyproperty
    def filter_names(self):

        """
        This function ...
        :return:
        """

        image_names = fs.files_in_path(self.analysis_run.misc_path, extension="fits", returns="name", contains="__")
        return [name.split("__")[1] for name in image_names]

    # -----------------------------------------------------------------

    @property
    def filter_names_sorted(self):

        """
        This function returns a list of the filter names, sorted on wavelength
        :return:
        """

        #return sorted(self.filter_names, key=lambda key: self.observed[key].filter.pivotwavelength())
        return sorted(self.filter_names, key=lambda key: parse_filter(key).pivotwavelength())

# -----------------------------------------------------------------
