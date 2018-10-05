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

        # Load residuals
        self.load()

        # Calculate residuals
        self.calculate()

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

    @property
    def do_weighed(self):
        return self.config.errors

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
            has_weighed = self.has_weighed_residuals(filter_name)
            if has_residual: log.success("Residual frame was already created for the '" + filter_name + "' filter: loading from file ...")
            if has_weighed: log.success("Weighed residual map was already created for the '" + filter_name + "' filter: loading from file ...")

            # To calculate?
            if has_residual and (not self.do_weighed or has_weighed): continue
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

    @lazyproperty
    def present_weighed_filter_names(self):
        return [filter_name for filter_name in self.filter_names if self.has_weighed_residuals(filter_name)]

    # -----------------------------------------------------------------

    @lazyproperty
    def present_weiged_filters(self):
        return [parse_filter(filter_name) for filter_name in self.present_weighed_filter_names]

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
        if self.config.errors: self.load_errors()

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

            # Check whether residual frames already created
            #if self.has_residuals(filter_name) and self.has_weighed_residuals(filter_name):
            #    log.success("Residual and weighed residual map for the '" + filter_name + "' filter already created")
            #    continue

            # Get observed image name
            name = self.observed_names_for_filters[fltr]

            # Debugging
            log.debug("Loading the observed " + name + " image ...")

            # Load the frame, not truncated
            frame = self.static_photometry_dataset.get_frame(name, masked=False)

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

        # Loop over all FITS files found in the images directory
        for filepath, filter_name in fs.files_in_path(self.simulated_images_path, extension="fits", returns=["path", "name"]):

            # Check whether residual frames already creeated
            #if self.has_residuals(filter_name) and self.has_weighed_residuals(filter_name):
            #    log.success("Residual and weighed residual map for the '" + filter_name + "' filter already created")
            #    continue

            # Debugging
            log.debug("Loading the simulated " + filter_name + " image ...")

            # Check whether a simulated image exists for this band
            if filter_name not in self.observed:
                log.warning("The observed " + filter_name + " image could not be found, skipping simulated " + filter_name + " image ...")
                continue

            # Open the image
            frame = Frame.from_file(filepath)

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

            # With errors
            if self.config.errors:

                # Get error map
                errors = self.errors[filter_name]

                # Rebin in-place
                names = ["simulated", "observed", "errors"]
                rebin_to_highest_pixelscale(simulated, observed, errors, names=names, in_place=True)

            # No errors
            else: rebin_to_highest_pixelscale(simulated, observed, names=["simulated", "observed"], in_place=True)

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

            # With errors
            if self.config.errors:

                # Get error map
                errors = self.errors[filter_name]

                # Crop
                errors.crop_to(self.truncation_box, factor=self.config.cropping_factor)

            # Crop each of them
            simulated.crop_to(self.truncation_box, factor=self.config.cropping_factor)
            observed.crop_to(self.truncation_box, factor=self.config.cropping_factor)

    # -----------------------------------------------------------------

    def load(self):

        """
        This function ...
        :return:
        """

        # Residuals
        self.load_residuals()

        # Weighed residuals
        if self.do_weighed: self.load_weighed()

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

    def load_weighed(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Loading weighed residual maps ...")

        # Loop over the filter names
        for filter_name in self.present_weighed_filter_names:

            # Debugging
            log.debug("Loading the weighed residual frame for the '" + filter_name + "' filter ...")

            # Load
            self.weighed[filter_name] = Frame.from_file(self.get_weighed_residuals_path(filter_name))

    # -----------------------------------------------------------------

    def calculate(self):

        """
        This function ...
        :return:
        """

        # Calculate the residual images
        self.calculate_residuals()

        # Calculate the weighed residual images
        if self.do_weighed: self.calculate_weighed()

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

            # Replace infs
            #residual.replace_infs(0.0)
            residual.replace_infs_by_nans()

            # Get the truncation mask
            truncation_mask = self.get_truncation_mask(observed.wcs)

            # MASK
            #residual[truncation_mask] = nan_value
            residual.replace_by_nans(truncation_mask)

            # Get the significance mask
            if self.config.errors:

                errors = self.errors[filter_name]
                significance_mask = self.get_significance_mask(observed, errors, min_npixels=self.config.min_npixels, connectivity=self.config.connectivity)
                #residual[significance_mask] = nan_value
                residual.replace_by_nans(significance_mask)

                # WRITE THE SIGNIFICANCE MASK
                mask_path = self.get_significance_mask_path(filter_name)
                significance_mask.saveto(mask_path)

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
            #residual.replace_infs(0.0)
            residual.replace_infs_by_nans()

            # Get the truncation mask
            truncation_mask = self.get_truncation_mask(observed.wcs)

            # Get the significance mask
            significance_mask = self.get_significance_mask(observed, errors, min_npixels=self.config.min_npixels, connectivity=self.config.connectivity)

            # MASK
            #residual[truncation_mask] = 0.0
            residual.replace_by_nans(truncation_mask)
            #residual[significance_mask] = 0.0
            residual.replace_by_nans(significance_mask)

            # Add the weighed residual image to the dictionary
            self.weighed[filter_name] = residual

            # Write
            self.write_weighed_map(filter_name)

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
        if self.config.errors: self.create_weighed_distributions()

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

            # REMOVE NANS
            notnan = np.logical_not(np.isnan(values))
            values = values[notnan]

            # Check
            if len(values) == 0 or sequences.all_equal(values):
                log.error("Cannot create distribution for the '" + filter_name + "' filter")
                continue

            # Create distribution
            distribution = Distribution.from_values("Weighed residual", values, nbins=self.config.nbins)

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
        if self.do_weighed: self.write_weighed()

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

    @lazyproperty
    def weighed_residual_maps_path(self):
        return fs.create_directory_in(self.analysis_run.residuals_path, "weighed_maps")

    # -----------------------------------------------------------------

    def get_weighed_residuals_path(self, filter_name):

        """
        Thisf unction ...
        :param filter_name:
        :return:
        """

        # Determine the path for this residual image
        path = fs.join(self.weighed_residual_maps_path, filter_name + ".fits")
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
        for filter_name in self.weighed: self.write_weighed_map(filter_name)

    # -----------------------------------------------------------------

    def write_weighed_map(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # Get path
        path = self.get_weighed_residuals_path(filter_name)
        if fs.is_file(path): return

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
        if self.do_weighed: self.write_weighed_distributions()

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

    @lazyproperty
    def weighed_distributions_path(self):
        return fs.create_directory_in(self.analysis_run.residuals_path, "weighed_distributions")

    # -----------------------------------------------------------------

    def get_weighed_residuals_distribution_path(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # Determine data file path
        path = fs.join(self.weighed_distributions_path, filter_name + ".dat")
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

        # Plot the absolute residual maps
        self.plot_absolute_residuals()

        # Plot the weighed rsiduals
        if self.config.errors: self.plot_weighed_residuals()

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
        if self.config.errors: self.plot_weighed_distributions()

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
            plot_distribution(self.residual_distributions[filter_name], path=path, x_limits=self.distribution_x_limits)

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

    @lazyproperty
    def weighed_distributions_plot_path(self):
        return fs.create_directory_in(self.analysis_run.residuals_path, "weighed_distribution_plots")

    # -----------------------------------------------------------------

    def get_weighed_residuals_distribution_plot_path(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        # Determine plot path
        path = fs.join(self.weighed_distributions_plot_path, filter_name + ".pdf")
        return path

    # -----------------------------------------------------------------

    @property
    def weighed_distribution_x_limits(self):
        return None

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
            plot_distribution(self.weighed_distributions[filter_name], path=path, x_limits=self.weighed_distribution_x_limits)

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

    @lazyproperty
    def weighed_residuals_plot_path(self):
        return fs.create_directory_in(self.analysis_run.residuals_path, "weighed_plots")

    # -----------------------------------------------------------------

    def get_weighed_residuals_plot_path(self, filter_name):

        """
        This function ...
        :param filter_name:
        :return:
        """

        return fs.join(self.weighed_residuals_plot_path, filter_name + ".png")

    # -----------------------------------------------------------------

    @property
    def weighed_residuals_plot_cmap(self):
        return "RdBu_r"

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

            # Debugging
            log.debug("Plotting the weighed residual map for the '" + filter_name + "' band to '" + path + "' ...")

            # Get the map
            residuals = self.weighed[filter_name]

            # Save as PNG
            # vmin, vmax = residuals.saveto_png(path, colours=self.config.colours,
            #                                   interval=self.config.interval,
            #                                   scale=self.config.scale, alpha=self.config.alpha_method,
            #                                   peak_alpha=self.config.peak_alpha)

            # Plot
            plotting.plot_box(residuals, path=path, colorbar=True, around_zero=True, scale="linear", symmetric=True,
                              cmap=self.weighed_residuals_plot_cmap, check_around_zero=False)

    # -----------------------------------------------------------------

    # def plot_image_grid(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Plotting a grid with the observed, simulated and residual images ...")
    #
    #     # Create the image grid plotter
    #     plotter = ResidualImageGridPlotter(title="Image residuals")
    #
    #     # Loop over the filter names, add a row to the image grid plotter for each filter
    #     for filter_name in self.filter_names_sorted:
    #
    #         # Get input
    #         observed = self.observed[filter_name]
    #         simulated = self.simulated[filter_name]
    #
    #         # Add row
    #         plotter.add_row(observed, simulated, filter_name)
    #
    #     # Set the bounding box for the plotter
    #     plotter.set_bounding_box(self.truncation_box)
    #
    #     # Determine the path to the plot file
    #     path = fs.join(self.analysis_run.residuals_path, "residuals.pdf")
    #
    #     # Run the plotter
    #     plotter.run(path)
    #
    # # -----------------------------------------------------------------
    #
    # def plot_image_grid_weighed(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Plotting a grid with the observed, simulated and weighed residual images ...")
    #
    #     # Create the image grid plotter
    #     plotter = ResidualImageGridPlotter(title="Weighed image residuals", weighed=True)
    #
    #     # Loop over the filter names, add a row to the image grid plotter for each filter
    #     for filter_name in self.filter_names_sorted:
    #
    #         # Get input
    #         observed = self.observed[filter_name]
    #         simulated = self.simulated[filter_name]
    #         errors = self.errors[filter_name]
    #
    #         # Add row
    #         plotter.add_row(observed, simulated, filter_name, errors=errors)
    #
    #     # Set the bounding box for the plotter
    #     plotter.set_bounding_box(self.truncation_box)
    #
    #     # Determine the path to the plot file
    #     path = fs.join(self.analysis_run.weighed_residuals_path, "weighed_residuals.pdf")
    #
    #     # Run the plotter
    #     plotter.run(path)

    # -----------------------------------------------------------------

    # @lazyproperty
    # def filter_names(self):
    #     image_names = fs.files_in_path(self.analysis_run.misc_path, extension="fits", returns="name", contains="__")
    #     return [name.split("__")[1] for name in image_names]
    #
    # # -----------------------------------------------------------------
    #
    # @lazyproperty
    # def filter_names_sorted(self):
    #     return list(sorted(self.filter_names, key=lambda key: parse_filter(key).pivotwavelength()))

# -----------------------------------------------------------------
