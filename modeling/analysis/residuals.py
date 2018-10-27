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
from .component import AnalysisRunComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...magic.core.frame import Frame
from ...core.tools import sequences
from ...core.basics.distribution import Distribution
from ...core.tools.utils import lazyproperty
from ...core.filter.filter import parse_filter
from ...magic.core.list import rebin_to_highest_pixelscale
from ...magic.tools import plotting
from ...core.plot.distribution import plot_distribution
from ...magic.core.mask import Mask

# -----------------------------------------------------------------

class ResidualAnalyser(AnalysisRunComponent):
    
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

        # The simulated images
        self.simulated = dict()

        # The observed images
        self.observed = dict()

        # The masks
        self.masks = dict()

        # The residual images
        self.residuals = dict()

        # The distributions of residual values
        self.distributions = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Get the observed and simulated images
        self.load_images()

        # Rebin the images
        self.rebin_images()

        # Crop the images
        self.crop_images()

        # Write the images
        self.write_images()

        # Make masks
        self.create_masks()

        # Load present residuals
        self.load_residuals()

        # Calculate residuals
        self.calculate_residuals()

        # Create distributions of the residual values
        self.create_distributions()

        # Writing
        self.write()

        # Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ResidualAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------
    # OBSERVED IMAGES
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_names_for_filters(self):
        return self.static_photometry_dataset.get_names_for_filters(self.filters, as_dict=True)

    # -----------------------------------------------------------------
    # SIMULATED IMAGES
    # -----------------------------------------------------------------

    @property
    def simulated_images_path(self):
        return self.analysis_run.images_path

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_image_names(self):
        return fs.files_in_path(self.simulated_images_path, returns="name", extension="fits")

    # -----------------------------------------------------------------

    @property
    def nsimulated_images(self):
        return len(self.simulated_image_names)

    # -----------------------------------------------------------------

    @property
    def has_simulated_images(self):
        return self.nsimulated_images > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_image_filters(self):
        return [parse_filter(name) for name in self.simulated_image_names]

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):
        return sequences.sorted_by_attribute(self.simulated_image_filters, "wavelength")

    # -----------------------------------------------------------------

    @lazyproperty
    def filter_names(self):
        return [str(fltr) for fltr in self.filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def to_calculate_filter_names(self):

        filter_names = []
        for filter_name in self.filter_names:

            # Already have residuals?
            has_residual = self.has_residuals(filter_name)
            if has_residual: log.success("Residual frame was already created for the '" + filter_name + "' filter: loading from file ...")

            # To calculate?
            if has_residual: continue
            else: filter_names.append(filter_name)

        return filter_names

    # -----------------------------------------------------------------

    @lazyproperty
    def to_calculate_filters(self):
        return [parse_filter(name) for name in self.to_calculate_filter_names]

    # -----------------------------------------------------------------

    @lazyproperty
    def present_residuals_filter_names(self):
        return [filter_name for filter_name in self.filter_names if self.has_residuals(filter_name)]

    # -----------------------------------------------------------------

    @lazyproperty
    def present_residuals_filters(self):
        return [parse_filter(filter_name) for filter_name in self.present_residuals_filter_names]

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

    # -----------------------------------------------------------------

    def load_observed_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed images ...")

        # Loop over the filters
        for fltr in self.to_calculate_filters:

            # Get filter name
            filter_name = str(fltr)

            # Get observed image name
            name = self.observed_names_for_filters[fltr]

            # Debugging
            log.debug("Loading the observed " + name + " image ...")

            # Load
            if self.has_observed_frame(filter_name): frame = self.load_observed_frame(filter_name)
            else: frame = self.static_photometry_dataset.get_frame(name, masked=False)

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

        # Loop over the filters
        for fltr in self.to_calculate_filters:

            # Get filter name
            filter_name = str(fltr)

            # Debugging
            log.debug("Loading the simulated " + filter_name + " image ...")

            # Load
            if self.has_simulated_frame(filter_name): frame = self.load_simulated_frame(filter_name)
            else: frame = Frame.from_file(fs.join(self.simulated_images_path, filter_name + ".fits"))

            # Add the simulated image frame to the dictionary
            self.simulated[filter_name] = frame

    # -----------------------------------------------------------------

    def rebin_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning the images to the same pixelscale ...")

        # Loop over the images
        for filter_name in self.to_calculate_filter_names:

            # Debugging
            log.debug("Rebinning the '" + filter_name + "' images ...")

            # Get the images
            simulated = self.simulated[filter_name]
            observed = self.observed[filter_name]

            # Rebin
            rebin_to_highest_pixelscale(simulated, observed, names=["simulated", "observed"], in_place=True)

    # -----------------------------------------------------------------

    def crop_images(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Cropping the images to the bounding box of the model ...")

        # Loop over the images
        for filter_name in self.to_calculate_filter_names:

            # Debugging
            log.debug("Cropping the '" + filter_name + "' images ...")

            # Get the images
            simulated = self.simulated[filter_name]
            observed = self.observed[filter_name]

            # Crop each of them
            if not simulated.get_meta("prepared", False): simulated.crop_to(self.truncation_box, factor=self.config.cropping_factor)
            if not observed.get_meta("prepared", False): observed.crop_to(self.truncation_box, factor=self.config.cropping_factor)

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the images ...")

        # Observed
        self.write_observed()

        # Simulated
        self.write_simulated()

    # -----------------------------------------------------------------

    def write_observed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the prepared observed frames ...")

        # Loop over the observed image frames
        for filter_name in self.observed: self.write_observed_frame(filter_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_frames_path(self):
        return fs.create_directory_in(self.analysis_run.residuals_path, "observed")

    # -----------------------------------------------------------------

    def get_observed_path(self, filter_name):

        """
        This function ...
        :return:
        """

        # Determine the path for this residual image
        path = fs.join(self.observed_frames_path, filter_name + ".fits")
        return path

    # -----------------------------------------------------------------

    def has_observed_frame(self, filter_name):
        return fs.is_file(self.get_observed_path(filter_name))

    # -----------------------------------------------------------------

    def load_observed_frame(self, filter_name):
        frame = Frame.from_file(self.get_observed_path(filter_name))
        frame.set_meta("prepared", True)
        return frame

    # -----------------------------------------------------------------

    def write_observed_frame(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # Get path
        path = self.get_observed_path(filter_name)
        if fs.is_file(path): return

        # Debugging
        log.debug("Writing the observed frame for the '" + filter_name + "' band to '" + path + "' ...")

        # Write the image
        self.observed[filter_name].saveto(path)

    # -----------------------------------------------------------------

    def write_simulated(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the prepared simulated frames ...")

        # Loop over the simulated image frames
        for filter_name in self.simulated: self.write_simulated_frame(filter_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_frames_path(self):
        return fs.create_directory_in(self.analysis_run.residuals_path, "simulated")

    # -----------------------------------------------------------------

    def get_simulated_path(self, filter_name):

        """
        This function ...
        :return:
        """

        # Determine the path for this frame
        path = fs.join(self.simulated_frames_path, filter_name + ".fits")
        return path

    # -----------------------------------------------------------------

    def has_simulated_frame(self, filter_name):
        return fs.is_file(self.get_simulated_path(filter_name))

    # -----------------------------------------------------------------

    def load_simulated_frame(self, filter_name):
        frame = Frame.from_file(self.get_simulated_path(filter_name))
        frame.set_meta("prepared", True)
        return frame

    # -----------------------------------------------------------------

    def write_simulated_frame(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # Get path
        path = self.get_simulated_path(filter_name)
        if fs.is_file(path): return

        # Debugging
        log.debug("Writing the simulated frame for the '" + filter_name + "' band to '" + path + "' ...")

        # Write the image
        self.simulated[filter_name].saveto(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def masks_path(self):
        return fs.create_directory_in(self.analysis_run.residuals_path, "masks")

    # -----------------------------------------------------------------

    def get_mask_path(self, filter_name):
        return fs.join(self.masks_path, filter_name + ".fits")

    # -----------------------------------------------------------------

    def has_mask(self, filter_name):
        return fs.is_file(self.get_mask_path(filter_name))

    # -----------------------------------------------------------------

    def load_mask(self, filter_name):
        return Mask.from_file(self.get_mask_path(filter_name))

    # -----------------------------------------------------------------

    def create_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating masks for the images ...")

        # Loop over the filters
        for filter_name in self.filter_names:

            # Check if present
            if self.has_mask(filter_name):
                self.masks[filter_name] = self.load_mask(filter_name)
                continue

            # Get the images
            simulated = self.simulated[filter_name]
            observed = self.observed[filter_name]

            # Get pixel distribution
            #distribution = Distribution.from_data("mock pixel values", simulated, logarithmic=True)
            #plot_distribution(distribution, title=filter_name, logscale=True)

            # Debugging
            log.debug("Creating mask for the '" + filter_name + "' image ...")

            # Get relative to maxium
            relative_simulated = simulated / simulated.max

            # Get mask of smallest pixel values
            simulated_mask = relative_simulated.where_smaller_than(1e-3)

            # Get the truncation mask
            truncation_mask = self.get_truncation_mask(observed.wcs)

            # Get the observation mask
            observation_mask = observed.nans + observed.zeroes

            # Combine the masks
            mask = simulated_mask + truncation_mask + observation_mask

            # Set the mask
            self.masks[filter_name] = mask

            # Save the mask
            mask.saveto(self.get_mask_path(filter_name))

    # -----------------------------------------------------------------

    def load_residuals(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Loading residual maps ...")

        # Loop over the filter names
        for filter_name in self.present_residuals_filter_names:

            # Debugging
            log.debug("Loading the residual frame for the '" + filter_name + "' filter ...")

            # Load
            self.residuals[filter_name] = Frame.from_file(self.get_residuals_path(filter_name))

    # -----------------------------------------------------------------

    def calculate_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the residual images ...")

        # Loop over the filter names
        for filter_name in self.to_calculate_filter_names:

            # Debugging
            log.debug("Creating the residual frame for the '" + filter_name + "' filter ...")

            # Get the observed and simulated image
            simulated = self.simulated[filter_name]
            observed = self.observed[filter_name]

            # Calculate the residual image
            residual = (simulated - observed) / observed

            # Replace infs by NaNs
            residual.replace_infs_by_nans()

            # MASK
            residual.replace_by_nans(self.masks[filter_name])

            # Add the residual image to the dictionary
            self.residuals[filter_name] = residual

            # Write
            self.write_residual_map(filter_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def significance_masks_path(self):
        return fs.create_directory_in(self.analysis_run.residuals_path, "significance")

    # -----------------------------------------------------------------

    def get_significance_mask_path(self, filter_name):
        return fs.join(self.significance_masks_path, filter_name + ".fits")

    # -----------------------------------------------------------------

    def has_significance_mask(self, filter_name):
        return fs.is_file(self.get_significance_mask_path(filter_name))

    # -----------------------------------------------------------------

    def create_distributions(self):

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
                self.distributions[filter_name] = self.load_residuals_distribution(filter_name)
                continue

            # Debugging
            log.debug("Creating the residuals distribution for the '" + filter_name + "' filter ...")

            # Get the values within the truncation ellipse
            values = self.residuals[filter_name].values_in(self.truncation_ellipse)

            # REMOVE EXACT ZEROES
            indices = np.argwhere(values == 0)
            values = np.delete(values, indices)

            # REMOVE NANS
            notnan = np.logical_not(np.isnan(values))
            values = values[notnan]

            # REMOVE TOO LOW OR TOO HIGH (PROBABLY NOISE)
            indices = np.argwhere(values < self.distribution_x_min)
            values = np.delete(values, indices)
            indices = np.argwhere(values > self.distribution_x_max)
            values = np.delete(values, indices)

            # Check
            if len(values) == 0 or sequences.all_equal(values):
                log.error("Cannot create distribution for the '" + filter_name + "' filter")
                continue

            # Create distribution
            distribution = Distribution.from_values("Residual", values, nbins=self.config.nbins)

            # Add the distribution
            self.distributions[filter_name] = distribution

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

    @lazyproperty
    def residual_maps_path(self):
        return fs.create_directory_in(self.analysis_run.residuals_path, "maps")

    # -----------------------------------------------------------------

    def get_residuals_path(self, filter_name):

        """
        This function ...
        :return:
        """

        # Determine the path for this residual image
        path = fs.join(self.residual_maps_path, filter_name + ".fits")
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
        for filter_name in self.residuals: self.write_residual_map(filter_name)

    # -----------------------------------------------------------------

    def write_residual_map(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # Get path
        path = self.get_residuals_path(filter_name)
        if fs.is_file(path): return

        # Debugging
        log.debug("Writing the residual frame for the '" + filter_name + "' band to '" + path + "' ...")

        # Write the image
        self.residuals[filter_name].saveto(path)

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

    @lazyproperty
    def distributions_path(self):
        return fs.create_directory_in(self.analysis_run.residuals_path, "distributions")

    # -----------------------------------------------------------------

    def get_residuals_distribution_path(self, filter_name):

        """
        Thisnfunction ...
        :param filter_name:
        :return:
        """

        # Determine data file path
        path = fs.join(self.distributions_path, filter_name + ".dat")
        return path

    # -----------------------------------------------------------------

    def write_distributions(self):

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
            self.distributions[filter_name].saveto(path)

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

        # Plot the absolute residual maps
        self.plot_absolute_residuals()

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

    @lazyproperty
    def distributions_plot_path(self):
        return fs.create_directory_in(self.analysis_run.residuals_path, "distribution_plots")

    # -----------------------------------------------------------------

    def get_residuals_distribution_plot_path(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # Determine plot path
        path = fs.join(self.distributions_plot_path, filter_name + ".pdf")
        return path

    # -----------------------------------------------------------------

    @property
    def distribution_x_min(self):
        return -4

    # -----------------------------------------------------------------

    @property
    def distribution_x_max(self):
        return 4

    # -----------------------------------------------------------------

    @property
    def distribution_x_limits(self):
        return (self.distribution_x_min, self.distribution_x_max,)  # because should be less than an order of magnitude (10)

    # -----------------------------------------------------------------

    def plot_distributions(self):

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
            plot_distribution(self.distributions[filter_name], path=path, x_limits=self.distribution_x_limits)

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

    @lazyproperty
    def residuals_plot_path(self):
        return fs.create_directory_in(self.analysis_run.residuals_path, "plots")

    # -----------------------------------------------------------------

    def get_residuals_plot_path(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        return fs.join(self.residuals_plot_path, filter_name + ".png")

    # -----------------------------------------------------------------

    @property
    def residuals_plot_interval(self):
        return (-100., 100.,)

    # -----------------------------------------------------------------

    @property
    def residuals_plot_cmap(self):
        return "RdBu_r"

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

            # Debugging
            log.debug("Plotting the residual map for the '" + filter_name + "' band to '" + path + "' ...")

            # Get residuals in percentage
            residuals = self.residuals[filter_name]
            residuals_percentage = residuals * 100

            # Save as PNG
            # vmin, vmax = residuals.saveto_png(path, colours=self.config.colours,
            #                               interval=self.config.interval,
            #                               scale=self.config.scale, alpha=self.config.alpha_method,
            #                               peak_alpha=self.config.peak_alpha)

            # Plot
            plotting.plot_box(residuals_percentage, interval=self.residuals_plot_interval, path=path, colorbar=True,
                              around_zero=True, scale="linear", cmap=self.residuals_plot_cmap, check_around_zero=False)

    # -----------------------------------------------------------------

    def has_absolute_residuals_plot(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # File present?
        if fs.is_file(self.get_absolute_residuals_plot_path(filter_name)):

            # Replot
            if self.config.replot_absolute_residuals:

                # Remove file
                fs.remove_file(self.get_absolute_residuals_plot_path(filter_name))
                return False

            # Don't replot
            else: return True

        # No file
        else: return False

    # -----------------------------------------------------------------

    @lazyproperty
    def absolute_residuals_plot_path(self):
        return fs.create_directory_in(self.analysis_run.residuals_path, "absolute_plots")

    # -----------------------------------------------------------------

    def get_absolute_residuals_plot_path(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        return fs.join(self.absolute_residuals_plot_path, filter_name + "_absolute.png")

    # -----------------------------------------------------------------

    @property
    def residuals_absolute_plot_interval(self):
        return (0, 100.)

    # -----------------------------------------------------------------

    @property
    def residuals_absolute_plot_cmap(self):
        return "viridis"

    # -----------------------------------------------------------------

    def plot_absolute_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absolute values of the residuals ...")

        # Loop over the filter names
        for filter_name in self.filter_names:

            # Get the path
            path = self.get_absolute_residuals_plot_path(filter_name)
            if fs.is_file(path): continue

            # Debugging
            log.debug("Plotting the absolute residual map for the '" + filter_name + "' band to '" + path + "' ...")

            # Get the absolute residual map in percentage
            residuals = self.residuals[filter_name].absolute
            residuals_percentage = residuals * 100

            # Plot
            plotting.plot_box(residuals_percentage, interval=self.residuals_absolute_plot_interval, path=path,
                              colorbar=True, scale="linear", cmap=self.residuals_absolute_plot_cmap)

# -----------------------------------------------------------------
