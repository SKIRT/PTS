#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.evaluation Contains the AnalysisModelEvaluator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ..fitting.tables import WeightsTable
from ..fitting.modelanalyser import FluxDifferencesTable
from ...magic.tools import wavelengths
from ...core.tools.utils import lazyproperty
from ...magic.core.list import FrameList
from ...core.tools import sequences
from ...core.basics.distribution import Distribution
from ...core.basics.containers import FilterBasedList
from ...core.data.sed import ObservedSED

# -----------------------------------------------------------------

class AnalysisModelEvaluator(AnalysisComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        Thisn function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(AnalysisModelEvaluator, self).__init__(*args, **kwargs)

        # The analysis run
        self.analysis_run = None

        # Create the table to contain the weights
        self.weights = WeightsTable()

        # Initialize the differences table
        self.differences = FluxDifferencesTable()

        # The chi squared value
        self.chi_squared = None

        # The mock observed images and their errors
        self.images = FrameList()
        self.errors = FrameList()

        # Fluxes calculated based on images
        self.images_fluxes = ObservedSED() # based on simulated images
        self.images_sed = ObservedSED() # based on observed images
        self.images_differences = FluxDifferencesTable() # differences

        # The residual and weighed residual frames
        self.residuals = FrameList()
        self.weighed = FrameList()

        # The distributions
        self.residuals_distributions = FilterBasedList()
        self.weighed_distributions = FilterBasedList()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Calculate weight for each band
        self.calculate_weights()

        # 3. Calculate flux differences
        self.calculate_differences()

        # 4. Calculate chi squared
        self.calculate_chi_squared()

        # 6. Make images
        self.make_images()

        # 7. Calculate fluxes from the images
        self.calculate_image_fluxes()

        # 8. Calculate SED
        self.calculate_image_sed()

        # 9. Calcualte differences between fluxes calculated from the images and the observed fluxes
        self.calculate_image_differences()

        # 10. Calculate the residual images
        self.calculate_residuals()

        # 11. Calculate the weighed residual images
        self.calculate_weighed()

        # 12. Create distributions of the residual values
        self.create_residuals_distributions()

        # 13. Create distributions of the weighed residual values
        self.create_weighed_distributions()

        # 14. Writing
        self.write()

        # 15. Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisModelEvaluator, self).setup(**kwargs)

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

    @property
    def simulated_sed(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.simulated_sed

    # -----------------------------------------------------------------

    @property
    def simulated_fluxes(self):

        """
        This function ..
        :return:
        """

        return self.analysis_run.simulated_fluxes

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_filters(self):

        """
        This function ...
        :return:
        """

        return self.simulated_fluxes.filters()

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_filters_no_iras_planck(self):

        """
        This function ...
        :return:
        """

        return [fltr for fltr in self.simulated_filters if fltr not in self.iras_and_planck_filters]

    # -----------------------------------------------------------------

    @property
    def simulation_prefix(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.simulation_prefix

    # -----------------------------------------------------------------

    @property
    def datacube(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.simulated_datacube

    # -----------------------------------------------------------------

    def get_frame_for_filter(self, fltr, convolve=False):

        """
        This function ...
        :param fltr:
        :param convolve:
        :return:
        """

        return self.datacube.frame_for_filter(fltr, convolve=convolve)

    # -----------------------------------------------------------------

    def calculate_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the weight to give to each band ...")

        # Initialize lists to contain the filters of the different wavelength ranges
        uv_bands = []
        optical_bands = []
        nir_bands = []
        mir_bands = []
        fir_bands = []
        submm_bands = []

        # Loop over the observed SED filters
        for fltr in self.simulated_filters_no_iras_planck:

            # Get the central wavelength
            wavelength = fltr.center

            # Get a string identifying which portion of the wavelength spectrum this wavelength belongs to
            spectrum = wavelengths.name_in_spectrum(wavelength)

            # Determine to which group
            if spectrum[0] == "UV": uv_bands.append(fltr)
            elif spectrum[0] == "Optical": optical_bands.append(fltr)
            elif spectrum[0] == "Optical/IR": optical_bands.append(fltr)
            elif spectrum[0] == "IR":
                if spectrum[1] == "NIR": nir_bands.append(fltr)
                elif spectrum[1] == "MIR": mir_bands.append(fltr)
                elif spectrum[1] == "FIR": fir_bands.append(fltr)
                else: raise RuntimeError("Unknown IR range")
            elif spectrum[0] == "Submm": submm_bands.append(fltr)
            else: raise RuntimeError("Unknown wavelength range")

        # Set the number of groups
        number_of_groups = 0

        # Check which groups are present
        has_uv = len(uv_bands) > 0
        has_optical = len(optical_bands) > 0
        has_nir = len(nir_bands) > 0
        has_mir = len(mir_bands) > 0
        has_fir = len(fir_bands) > 0
        has_submm = len(submm_bands) > 0

        if has_uv: number_of_groups += 1
        if has_optical: number_of_groups += 1
        if has_nir: number_of_groups += 1
        if has_mir: number_of_groups += 1
        if has_fir: number_of_groups += 1
        if has_submm: number_of_groups += 1

        # Detemrine total number of data points
        number_of_data_points = len(self.simulated_filters_no_iras_planck)

        # Determine the weight for each group of filters
        uv_weight = 1. / (len(uv_bands) * number_of_groups) * number_of_data_points if has_uv else 0.0
        optical_weight = 1. / (len(optical_bands) * number_of_groups) * number_of_data_points if has_optical else 0.0
        nir_weight = 1. / (len(nir_bands) * number_of_groups) * number_of_data_points if has_nir else 0.0
        mir_weight = 1. / (len(mir_bands) * number_of_groups) * number_of_data_points if has_mir else 0.0
        fir_weight = 1. / (len(fir_bands) * number_of_groups) * number_of_data_points if has_fir else 0.0
        submm_weight = 1. / (len(submm_bands) * number_of_groups) * number_of_data_points if has_submm else 0.0

        # Debugging
        if has_uv: log.debug("UV: number of bands = " + str(len(uv_bands)) + ", weight = " + str(uv_weight))
        if has_optical: log.debug("Optical: number of bands = " + str(len(optical_bands)) + ", weight = " + str(optical_weight))
        if has_nir: log.debug("NIR: number of bands = " + str(len(nir_bands)) + ", weight = " + str(nir_weight))
        if has_mir: log.debug("MIR: number of bands = " + str(len(mir_bands)) + ", weight = " + str(mir_weight))
        if has_fir: log.debug("FIR: number of bands = " + str(len(fir_bands)) + ", weight = " + str(fir_weight))
        if has_submm: log.debug("Submm: number of bands = " + str(len(submm_bands)) + ", weight = " + str(submm_weight))

        # Loop over the bands in each group and set the weight in the weights table
        for fltr in uv_bands: self.weights.add_point(fltr, uv_weight)
        for fltr in optical_bands: self.weights.add_point(fltr, optical_weight)
        for fltr in nir_bands: self.weights.add_point(fltr, nir_weight)
        for fltr in mir_bands: self.weights.add_point(fltr, mir_weight)
        for fltr in fir_bands: self.weights.add_point(fltr, fir_weight)
        for fltr in submm_bands: self.weights.add_point(fltr, submm_weight)

    # -----------------------------------------------------------------

    def calculate_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the differences between the observed and simulated SED ...")

        # Loop over the entries in the fluxdensity table (SED) derived from the simulation
        for i in range(len(self.simulated_fluxes)):

            # Get the filter
            fltr = self.simulated_fluxes.get_filter(i)

            # Get the flux density
            fluxdensity = self.simulated_fluxes.get_photometry(i, unit="Jy", add_unit=False)

            # Find the corresponding flux in the SED derived from observation
            #observed_fluxdensity = self.observed_sed.photometry_for_band(instrument, band, unit="Jy", add_unit=False)
            observed_fluxdensity = self.observed_sed.photometry_for_filter(fltr, unit="Jy", add_unit=False)

            # Find the corresponding flux error in the SED derived from observation
            #observed_fluxdensity_error = self.observed_sed.error_for_band(instrument, band, unit="Jy").average.to("Jy").value
            observed_fluxdensity_error = self.observed_sed.error_for_filter(fltr, unit="Jy", add_unit=False)

            # If no match with (instrument, band) is found in the observed SED
            if observed_fluxdensity is None:
                log.warning("The observed flux density could not be found for the " + str(fltr) + " filter")
                continue

            # Calculate the difference
            difference = fluxdensity - observed_fluxdensity
            relative_difference = difference / observed_fluxdensity

            # Find the index of the current band in the weights table
            index = self.weights.index_for_filter(fltr, return_none=True)
            if index is None:
                log.warning("A weight is not calculated for the " + str(fltr) + " filter")
                continue # Skip this band if a weight is not found

            # Get the weight
            weight = self.weights["Weight"][index]

            # Calculate the chi squared term
            chi_squared_term = weight * difference ** 2 / observed_fluxdensity_error ** 2

            # Add an entry to the differences table
            self.differences.add_entry(fltr.instrument, fltr.band, difference, relative_difference, chi_squared_term)

    # -----------------------------------------------------------------

    def calculate_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the chi squared value for this model ...")

        # Calculate the degrees of freedom
        dof = len(self.differences) - 3. - 1.  # number of data points - number of fitted parameters - 1

        # The (reduced) chi squared value is the sum of all the terms (for each band),
        # divided by the number of degrees of freedom
        self.chi_squared = np.sum(self.differences["Chi squared term"]) / dof

        # Debugging
        log.debug("The chi squared value is " + str(self.chi_squared))

    # -----------------------------------------------------------------

    def make_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the mock observed images ...")

        # Loop over the filters for which we want to create images
        for fltr in self.simulated_filters_no_iras_planck:

            # Debugging
            log.debug("Making the mock observed '" + str(fltr) + "' image ...")

            # Get the appropriate simulated frame
            frame = self.get_frame_for_filter(fltr, convolve=False)

            # Add the frame
            self.images.append(frame)

            # Debugging
            log.debug("Making an approximate error map for the '" + str(fltr) + "' image ...")

            # Create an approximate error frame
            errors = np.sqrt(frame.data, wcs=frame.wcs, filter=frame.filter)

            # Add the errors frame
            self.errors.append(errors)

    # -----------------------------------------------------------------

    def calculate_image_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the fluxes based on the mock observed images ...")

        # Loop over the filters
        for fltr in self.simulated_filters_no_iras_planck:

            # Calculate the total flux
            flux = self.images[fltr].sum(add_unit=True)

            # Add to the SED
            self.images_fluxes.add_point(fltr, flux)

    # -----------------------------------------------------------------

    def calculate_image_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the fluxes based on the observed images ...")

        # Loop over the observed images
        for frame in self.frame_list:

            # Calculate the total flux
            flux = frame.sum(add_unit=True)

            # Add to the SED
            self.images_sed.add_point(frame.fltr, flux)

    # -----------------------------------------------------------------

    def calculate_image_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the flux differences based on the images ...")

        # Loop over the fluxes
        for i in range(len(self.images_fluxes)):

            # Get the filter
            fltr = self.images_fluxes.get_filter(i)

            # Get the flux density
            fluxdensity = self.images_fluxes.get_photometry(i, unit="Jy", add_unit=False)

            # Find the corresponding flux in the SED derived from observation
            # observed_fluxdensity = self.observed_sed.photometry_for_band(instrument, band, unit="Jy", add_unit=False)
            #observed_fluxdensity = self.observed_sed.photometry_for_filter(fltr, unit="Jy", add_unit=False)
            observed_fluxdensity = self.images_sed.photometry_for_filter(fltr, unit="Jy", add_unit=False)

            # If no match with (instrument, band) is found in the observed SED
            if observed_fluxdensity is None:
                log.warning("The observed flux density could not be found for the " + str(fltr) + " filter")
                continue

            # Calculate the difference
            difference = fluxdensity - observed_fluxdensity
            relative_difference = difference / observed_fluxdensity

            # Add an entry to the differences table
            self.images_differences.add_entry(fltr.instrument, fltr.band, difference, relative_difference)

    # -----------------------------------------------------------------

    def calculate_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the residual images ...")

        # Loop over the filters
        for fltr in self.simulated_filters_no_iras_planck:

            # # Already have residuals
            # if self.has_residuals(filter_name):
            #     log.success("Residual frame was already created for the '" + filter_name + " filter': loading from file ...")
            #     self.residuals[filter_name] = self.load_residuals(filter_name)
            #     continue

            # Debugging
            log.debug("Creating the residual frame for the '" + str(fltr) + "' filter ...")

            # Get the observed and simulated image
            simulated = self.images[fltr]
            observed = self.frame_list[fltr]
            errors = self.errormap_list[fltr]

            # Calculate the residual image
            residual = (simulated - observed) / observed

            # Set the filter
            residual.filter = fltr

            # Replace infs
            residual.replace_infs(0.0)

            # Get the truncation mask
            truncation_mask = self.get_truncation_mask(observed.wcs)

            # Get the significance mask
            significance_mask = self.get_significance_mask(observed, errors, min_npixels=self.config.min_npixels, connectivity=self.config.connectivity)

            # MASK
            residual[truncation_mask] = 0.0
            residual[significance_mask] = 0.0

            # Add the residual image
            self.residuals.append(residual)

    # -----------------------------------------------------------------

    def calculate_weighed(self):

        """
        This ufnction ...
        :return:
        """

        # Inform the user
        log.info("Calculating the weighed residual images ...")

        # Loop over the filters
        for fltr in self.simulated_filters_no_iras_planck:

            # # Already have residuals
            # if self.has_weighed_residuals(filter_name):
            #     log.success("Weighed residual frame was already created for the '" + filter_name + "' filter: loading from file ...")
            #     self.weighed[filter_name] = self.load_weighed_residuals(filter_name)
            #     continue

            # Debugging
            log.debug("Creating the weighed residual frame for the '" + str(fltr) + "' filter ...")

            # Get the observed and simulated image
            simulated = self.images[fltr]
            observed = self.frame_list[fltr]
            errors = self.errormap_list[fltr]

            # Calculate the weighed residual image
            residual = (simulated - observed) / errors

            # Set the filter
            residual.filter = fltr

            # Replace infs
            residual.replace_infs(0.0)

            # Get the truncation mask
            truncation_mask = self.get_truncation_mask(observed.wcs)

            # Get the significance mask
            significance_mask = self.get_significance_mask(observed, errors, min_npixels=self.config.min_npixels, connectivity=self.config.connectivity)

            # MASK
            residual[truncation_mask] = 0.0
            residual[significance_mask] = 0.0

            # Add the weighed residual image
            self.weighed.append(residual)

    # -----------------------------------------------------------------

    def create_residuals_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the residuals distributions ...")

        # Loop over the filters
        for fltr in self.simulated_filters_no_iras_planck:

            # # Already have distribution
            # if self.has_residuals_distribution(filter_name):
            #     log.success("Residuals distribution was already created for the '" + filter_name + "' filter: loading from file ...")
            #     self.residual_distributions[filter_name] = self.load_residuals_distribution(filter_name)
            #     continue

            # Debugging
            log.debug("Creating the residuals distribution for the '" + str(fltr) + "' filter ...")

            # Get the values within the truncation ellipse
            values = self.residuals[fltr].values_in(self.truncation_ellipse)

            # REMOVE EXACT ZEROES
            indices = np.argwhere(values == 0)
            values = np.delete(values, indices)

            # Check
            if len(values) == 0 or sequences.all_equal(values):
                log.error("Cannot create distribution for the '" + str(fltr) + "' filter")
                continue

            # Create distribution
            distribution = Distribution.from_values(values, bins=self.config.nbins)

            # Add the distribution
            self.residuals_distributions.append(fltr, distribution)

    # -----------------------------------------------------------------

    def create_weighed_distributions(self):

        """
        This function ...
        :return:
        """

        # Loop over the filters
        for fltr in self.simulated_filters_no_iras_planck:

            # Debugging
            log.debug("Creating the residuals distribution for the '" + str(fltr) + "' filter ...")

            # Get the values within the truncation ellipse
            values = self.weighed[fltr].values_in(self.truncation_ellipse)

            # REMOVE EXACT ZEROES
            indices = np.argwhere(values == 0)
            values = np.delete(values, indices)

            # Check
            if len(values) == 0 or sequences.all_equal(values):
                log.error("Cannot create distribution for the '" + str(fltr) + "' filter")
                continue

            # Create distribution
            distribution = Distribution.from_values(values, bins=self.config.nbins)

            # Add the distribution
            self.weighed_distributions.append(fltr, distribution)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Infomr the user
        log.info("Writing ...")

        # Write the weights
        self.write_weights()

        # Write the flux differences
        self.write_differences()

        # Write the image fluxes
        self.write_image_fluxes()

        # Write the image SED
        self.write_image_sed()

        # Write the image flux differences
        self.write_image_differences()

        # Write the images
        self.write_images()

        # Write the residual frames
        self.write_residuals()

        # Write the weighed residual frames
        self.write_weighed()

        # Write the residual distributions
        self.write_residuals_distributions()

        # Write the weighed residual distributions
        self.write_weighed_distributions()

    # -----------------------------------------------------------------

    @lazyproperty
    def fluxes_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.evaluation_path, "fluxes")

    # -----------------------------------------------------------------

    def write_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with weights ...")

        # Determine the path
        path = fs.join(self.fluxes_path, "weights.dat")

        # Write the table with weights
        self.weights.saveto(path)

    # -----------------------------------------------------------------

    def write_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with the flux density differences ...")

        # Determine the path
        path = fs.join(self.fluxes_path, "differences.dat")

        # Save the differences table
        self.differences.saveto(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def images_fluxes_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.evaluation_path, "fluxes_images")

    # -----------------------------------------------------------------

    def write_image_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing image fluxes ...")

        # Detemrine the path
        path = fs.join(self.images_fluxes_path, "fluxes.dat")

        # Save
        self.images_fluxes.saveto(path)

    # -----------------------------------------------------------------

    def write_image_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing image SED ...")

        # Determine the path
        path = fs.join(self.images_fluxes_path, "observed_sed.dat")

        # Save
        self.images_sed.saveto(path)

    # -----------------------------------------------------------------

    def write_image_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing image flux differences ...")

        # Determine the path
        path = fs.join(self.images_fluxes_path, "differences.dat")

        # Save
        self.images_differences.saveto(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def images_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.evaluation_path, "images")

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Infomr the user
        log.info("Writing the images ...")

        # Loop over the images
        for frame in self.images:

            # Determine the path
            path = fs.join(self.images_path, frame.filter_name + ".fits")

            # Save the frame
            frame.saveto(path)

    # -----------------------------------------------------------------

    def write_errors(self):

        """
        Thisj function ...
        :return:
        """

        # Inform the user
        log.info("Writing the error maps ...")

        # Loop over the error maps
        for errors in self.errors:

            # Determine the path
            path = fs.join(self.images_path, errors.filter_name + "_errors.fits")

            # Save the error map
            errors.saveto(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def residuals_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.evaluation_path, "residuals")

    # -----------------------------------------------------------------

    @lazyproperty
    def weighed_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.evaluation_path, "weighed")

    # -----------------------------------------------------------------

    def write_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the residuals ...")

        # Loop over the residual maps
        for residuals in self.residuals:

            # Determine the path
            path = fs.join(self.residuals_path, residuals.filter_name + ".fits")

            # Save
            residuals.saveto(path)

    # -----------------------------------------------------------------

    def write_weighed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the weighed residuals ...")

        # Loop over the weighed residuals maps
        for residuals in self.weighed:

            # Determine the path
            path = fs.join(self.weighed_path, residuals.filter_name + ".fits")

            # Save
            residuals.saveto(path)

    # -----------------------------------------------------------------

    def write_residuals_distributions(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing the residuals distributions ...")

        # Loop over the distributions
        for fltr in self.simulated_filters_no_iras_planck:

            # Get the distribution
            distribution = self.residuals_distributions[fltr]

            # Determine the path
            path = fs.join(self.residuals_path, str(fltr) + "_distribution.dat")

            # Save
            distribution.saveto(path)

    # -----------------------------------------------------------------

    def write_weighed_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the weighed distributions ...")

        # Loop over the distributions
        for fltr in self.simulated_filters_no_iras_planck:

            # Get the distribution
            distribution = self.weighed_distributions[fltr]

            # Determine the path
            path = fs.join(self.weighed_path, str(fltr) + "_distribution.dat")

            # Save
            distribution.saveto(path)

# -----------------------------------------------------------------
