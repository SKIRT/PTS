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
from ...core.data.sed import SED, ObservedSED
from ...magic.core.frame import Frame
from ...core.plot.sed import SEDPlotter
from ...magic.core.list import convert_to_same_unit
from ...magic.tools import plotting
from ...magic.core.list import rebin_to_highest_pixelscale

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

        # The simulated SED derived from the datacube
        self.simulated_datacube_sed = None

        # The mock observed images and their errors
        self.images = FrameList()
        self.errors = FrameList()

        # Rebinned?
        self.rebinned_observed = []
        self.rebinned_simulated = []

        # The proper mock observed images
        self.proper_images = FrameList()

        # The observed images and their errors
        self.observed_images = FrameList()
        self.observed_errors = FrameList()

        # Fluxes calculated based on images
        self.images_fluxes = ObservedSED(photometry_unit="Jy") # based on simulated images
        self.images_sed = ObservedSED(photometry_unit="Jy") # based on observed images
        self.images_differences = FluxDifferencesTable() # differences

        # Fluxes calculated based on properly made mock observed SEDs
        self.proper_images_fluxes = ObservedSED(photometry_unit="Jy")

        # Differences between either direct total fluxes or fluxes from images
        self.fluxes_differences = FluxDifferencesTable()
        self.sed_differences = FluxDifferencesTable()

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
        if not self.has_weights: self.calculate_weights()
        else: self.load_weights()

        # 3. Calculate flux differences
        if not self.has_differences: self.calculate_differences()
        else: self.load_differences()

        # 4. Calculate chi squared
        self.calculate_chi_squared()

        # 5. Create simulated SED derived from the datacubes
        if not self.has_simulated_datacube_sed: self.create_simulated_datacube_sed()
        else: self.load_simulated_datacube_sed()

        # 6. Make images
        self.make_images()

        # 7. Load images created with convolution etc.
        self.load_proper_images()

        # 8. Load observed images
        self.load_observed_images()

        # 9. Rebin the images to the same pixelscale
        self.rebin_images()

        # 10. Calculate fluxes from the images
        if not self.has_image_fluxes: self.calculate_image_fluxes()
        else: self.load_image_fluxes()

        # 11. Calculate fluxes from the proper images
        if not self.has_proper_image_fluxes: self.calculate_proper_image_fluxes()
        else: self.load_proper_image_fluxes()

        # 12. Calculate SED
        if not self.has_images_sed: self.calculate_image_sed()
        else: self.load_images_sed()

        # 13. Calcualte differences between fluxes calculated from the images and the observed fluxes
        if not self.has_images_differences: self.calculate_image_differences()
        else: self.load_image_differences()

        # 14. Calculate the differences between the simulated fluxes and the fluxes from the observed images
        if not self.has_fluxes_differences: self.calculate_fluxes_differences()
        else: self.load_fluxes_differences()

        # 15. Calculate the differences between the observed fluxes and the observed fluxes from the images
        if not self.has_sed_differences: self.calculate_sed_differences()
        else: self.load_sed_differences()

        # 16. Calculate the residual images
        self.calculate_residuals()

        # 17. Calculate the weighed residual images
        self.calculate_weighed()

        # 18. Create distributions of the residual values
        self.create_residuals_distributions()

        # 19. Create distributions of the weighed residual values
        self.create_weighed_distributions()

        # 20. Writing
        self.write()

        # 21. Plotting
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

        # Get the filters
        filters = [fltr for fltr in self.simulated_filters if fltr not in self.iras_and_planck_filters]

        # Sort them
        return list(sorted(filters, key=lambda fltr: fltr.wavelength.to("micron").value))

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

    def load_weights(self):

        """
        This function ...
        :return:
        """

        # Give message
        log.success("Loading the weights ...")

        # Load the weights
        self.weights = WeightsTable.from_file(self.weights_filepath)

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
            observed_fluxdensity_error = self.observed_sed.error_for_filter(fltr, unit="Jy", add_unit=False).average

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

    def load_differences(self):

        """
        This function ...
        :return:
        """

        # Give message
        log.success("Loading the differences ...")

        # Load
        self.differences = FluxDifferencesTable.from_file(self.differences_filepath)

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

    def create_simulated_datacube_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating simulated SED derived from the datacube ...")

        # Get datacube, but not in surface brigthness
        datacube = self.datacube.converted_to_corresponding_non_brightness_unit()

        # Create global SED from the datacube
        self.simulated_datacube_sed = datacube.global_sed()

    # -----------------------------------------------------------------

    def load_simulated_datacube_sed(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.success("Simulated SED derived from the datacube is already present: loading from file ...")

        # Load
        self.simulated_datacube_sed = SED.from_file(self.simulated_datacube_sed_filepath)

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

            # Present?
            if self.has_image_for_filter(fltr):

                # Needed?
                if self.has_residuals_for_filter(fltr) and self.has_weighed_for_filter(fltr) and self.has_image_fluxes and self.has_all_image_plots and self.has_all_relative_error_plots: continue

                # Success
                log.success("The '" + str(fltr) + "' image is already present: loading ...")

                # Load the frame
                frame = Frame.from_file(self.get_image_filepath_for_filter(fltr))

                # Add the frame
                self.images.append(frame)

            # Not yet present
            else:

                # Debugging
                log.debug("Making the mock observed '" + str(fltr) + "' image ...")

                # Get the appropriate simulated frame
                frame = self.get_frame_for_filter(fltr, convolve=False)

                # Print wavelengths to check
                #print(fltr.wavelength, fltr.center, fltr.pivot, frame._wavelength)

                # Give warning
                rel_difference = (frame._wavelength - fltr.pivot) / fltr.pivot
                if rel_difference > 0.05: log.warning("Wavelength of filter: " + str(fltr.pivot) + ", wavelength of frame: " + str(frame._wavelength))

                # Print unit BEFORE
                #print("UNIT BEFORE:", frame.unit)

                # Convert to non-brightness
                conversion_factor = frame.convert_to_corresponding_non_brightness_unit()

                #print("PIXELSCALE:", frame.average_pixelscale)
                #print("CONVERSION FACTOR", conversion_factor)

                # Print unit AFTER
                #print("UNIT AFTER:", frame.unit)

                # Add the frame
                self.images.append(frame)

            # Present?
            if self.has_errors_for_filter(fltr):

                # Needed?
                if self.has_residuals_for_filter(fltr) and self.has_weighed_for_filter(fltr) and self.has_all_relative_error_plots: continue

                # Success
                log.success("The '" + str(fltr) + "' image is already present: loading ...")

                # Load the error frame
                errors = Frame.from_file(self.get_errors_filepath_for_filter(fltr))

                # Add the errors frame
                self.errors.append(errors)

            # Not yet present
            else:

                # Debugging
                log.debug("Making an approximate error map for the '" + str(fltr) + "' image ...")

                # Create an approximate error frame
                errors = Frame(np.sqrt(frame.data), wcs=frame.wcs, filter=frame.filter, unit=frame.unit)

                # Add the errors frame
                self.errors.append(errors)

    # -----------------------------------------------------------------

    def load_proper_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the proper mock observed images ...")

        # Loop over the filters for which we want to create images
        for fltr in self.simulated_filters_no_iras_planck:
            
            # Load the image
            frame = self.analysis_run.get_simulated_frame_for_filter(fltr)

            # Add the image
            self.proper_images.append(frame)

    # -----------------------------------------------------------------

    def load_observed_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed images ...")

        # Loop over the appropriate observed images
        image_names_for_filters = self.dataset.get_names_for_filters(self.simulated_filters_no_iras_planck)
        for fltr, image_name in zip(self.simulated_filters_no_iras_planck, image_names_for_filters):

            # Needs to be loaded?
            if self.has_residuals_for_filter(fltr) and self.has_weighed_for_filter(fltr) and self.has_images_sed: continue

            # # Check if not None (in database)
            # if image_name is None:
            #     log.warning("No observed " + key + " image was found")
            #     continue

            # # Check whether residual frames already created
            # if self.has_residuals(filter_name) and self.has_weighed_residuals(filter_name):
            #     log.success("Residual and weighed residual map for the '" + filter_name + "' filter already created")
            #     continue

            # Debugging
            log.debug("Loading the observed '" + str(fltr) + "' image ...")

            # Load the frame, not truncated
            frame = self.dataset.get_frame(image_name, masked=False)

            # Add the frame
            self.observed_images.append(frame)

            # Debugging
            log.debug("Loading the observed '" + str(fltr) + "' error map ...")

            # Get the error map
            errors = self.dataset.get_errormap(image_name)

            # Add the error map
            self.observed_errors.append(errors)

    # -----------------------------------------------------------------

    def rebin_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning the images to the same pixelscale ...")

        # Loop over the filters
        for fltr in self.simulated_filters_no_iras_planck:

            # All already present
            if self.has_residuals_for_filter(fltr) and self.has_weighed_for_filter(fltr) and self.has_images_sed and self.has_image_fluxes: continue

            # Debugging
            log.debug("Rebinning the '" + str(fltr) + "' images if necessary ...")

            # Get the frames
            simulated = self.images[fltr]
            simulated_errors = self.errors[fltr]
            observed = self.observed_images[fltr]
            observed_errors = self.observed_errors[fltr]

            # Rebin in-place
            names = ["simulated", "simulated_errors", "observed", "observed_errors"]
            rebin_to_highest_pixelscale(simulated, simulated_errors, observed, observed_errors, names=names, in_place=True)

    # -----------------------------------------------------------------

    # def rebin_images_old(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info('Rebinning the images to the same pixelscale ...')
    #
    #     # Loop over the filters
    #     for fltr in self.simulated_filters_no_iras_planck:
    #
    #         # Debugging
    #         log.debug("Rebinning the '" + str(fltr) + "' image if necessary ...")
    #
    #         # All already present
    #         if self.has_residuals_for_filter(fltr) and self.has_weighed_for_filter(fltr) and self.has_images_sed and self.has_image_fluxes: continue
    #
    #         # Get the frames
    #         simulated = self.images[fltr]
    #         simulated_errors = self.errors[fltr]
    #         observed = self.observed_images[fltr]
    #         observed_errors = self.observed_errors[fltr]
    #
    #         # Check coordinate systems
    #         if simulated.wcs != simulated_errors.wcs: raise ValueError("The coordinate system of the simulated frame and error map are not the same")
    #         if observed.wcs != observed_errors.wcs: raise ValueError("The coordinate system of the observed frame and error map are not the same")
    #
    #         # Check whether the coordinate systems of the observed and simulated image match
    #         if simulated.wcs == observed.wcs:
    #             log.debug("The coordinate system of the simulated and observed image for the " + str(fltr) + " filter matches")
    #
    #         # The observed image has a smaller pixelscale as the simulated image -> rebin the observed image
    #         elif observed.average_pixelscale < simulated.average_pixelscale:
    #
    #             # Debugging
    #             log.debug("The observed '" + str(fltr) + "' image has a better resolution as the simulated image: rebinning the observed image ...")
    #
    #             # Rebin the observed images
    #             observed.rebin(simulated.wcs)
    #             observed_errors.rebin(simulated.wcs)
    #
    #             # Add
    #             self.rebinned_observed.append(fltr)
    #
    #         # The simulated image has a smaller pixelscale as the observed image
    #         elif simulated.average_pixelscale < observed.average_pixelscale:
    #
    #             # Debugging
    #             log.debug("The simulated '" + str(fltr) + "' image has a better resolution as the observed image: rebinning the simulated image ...")
    #
    #             # Rebin the simulated images
    #             simulated.rebin(observed.wcs)
    #             simulated_errors.rebin(observed.wcs)
    #
    #             # Add
    #             self.rebinned_simulated.append(fltr)
    #
    #         # Error
    #         else: raise RuntimeError("Something unexpected happened")
    #
    #     # DEBUGGING
    #     log.debug("")
    #     log.debug("The following simulated images (and error maps) have been rebinned:")
    #     log.debug("")
    #     for fltr in self.rebinned_simulated: log.debug(" - " + str(fltr))
    #     log.debug("")
    #
    #     # DEBUGGING
    #     log.debug("The following observed images (and error maps) have been rebinned:")
    #     log.debug("")
    #     for fltr in self.rebinned_observed: log.debug(" - " + str(fltr))
    #     log.debug("")

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

            # Get the image
            image = self.images[fltr]

            # Convert to non - per angular or intrinsic area
            image.convert_to_corresponding_non_angular_or_intrinsic_area_unit()

            # Calculate the total flux
            flux = image.sum_in(self.truncation_ellipse, add_unit=True)

            # Add to the SED
            self.images_fluxes.add_point(fltr, flux)

    # -----------------------------------------------------------------

    def load_image_fluxes(self):

        """
        This function ...
        :return:
        """

        # Success
        log.success("Image fluxes are already present: loading from file ...")

        # Load
        self.images_fluxes = ObservedSED.from_file(self.images_fluxes_filepath)

    # -----------------------------------------------------------------

    def calculate_proper_image_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the fluxes based on the proper mock observed images ...")

        # Loop over the filters
        for fltr in self.simulated_filters_no_iras_planck:

            # Get the image
            image = self.proper_images[fltr]

            # Convert to non - per angular or intrinsic area
            image.convert_to_corresponding_non_angular_or_intrinsic_area_unit()

            # Calculate the total flux
            flux = image.sum_in(self.truncation_ellipse, add_unit=True)

            # Add to the SED
            self.proper_images_fluxes.add_point(fltr, flux)

    # -----------------------------------------------------------------

    def load_proper_image_fluxes(self):

        """
        This function ...
        :return:
        """

        # Succes
        log.success("Proper image fluxes are already presnet: loading from file ...")

        # Load
        self.proper_images_fluxes = ObservedSED.from_file(self.proper_images_fluxes_filepath)

    # -----------------------------------------------------------------

    def calculate_image_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the fluxes based on the observed images ...")

        # Loop over the filters
        for fltr in self.simulated_filters_no_iras_planck:

            # Get the frame
            frame = self.observed_images[fltr]

            # Convert to non - per angular or intrinsic area unit
            frame.convert_to_corresponding_non_angular_or_intrinsic_area_unit()

            # Calculate the total flux
            flux = frame.sum_in(self.truncation_ellipse, add_unit=True)

            # Add to the SED
            self.images_sed.add_point(frame.filter, flux)

    # -----------------------------------------------------------------

    def load_images_sed(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.success("Images SED is already present: loading from file ...")

        # Load
        self.images_sed = ObservedSED.from_file(self.images_sed_filepath)

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

    def load_image_differences(self):

        """
        This function ...
        :return:
        """

        # MEssage
        log.success("Image differences are already present: loading from file ...")

        # Load
        self.images_differences = FluxDifferencesTable.from_file(self.images_differences_filepath)

    # -----------------------------------------------------------------

    def calculate_fluxes_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculate the differences between the simulated fluxes and the fluxes from the simulated images ...")

        # Loop over the fluxes from mock observed images
        for i in range(len(self.images_fluxes)):

            # Get the filter
            fltr = self.images_fluxes.get_filter(i)

            # Get the flux density
            flux = self.images_fluxes.get_photometry(i, unit="Jy", add_unit=False)

            # Get the calculated mock observed flux
            observed_flux = self.simulated_fluxes.get_photometry(i, unit="Jy", add_unit=False)

            # Add to table
            self.fluxes_differences.add_from_filter_and_fluxes(fltr, flux, observed_flux)

    # -----------------------------------------------------------------

    def load_fluxes_differences(self):

        """
        Thisfunction ...
        :return:
        """

        # Message
        log.success("Fluxes differences are already present: loading from file ...")

        # Load
        self.fluxes_differences = FluxDifferencesTable.from_file(self.fluxes_differences_filepath)

    # -----------------------------------------------------------------

    def calculate_sed_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the differences between the observed fluxes and the fluxes from the observed images ...")

        # Loop over the fluxes
        for i in range(len(self.images_sed)):

            # Get the filter
            fltr = self.images_sed.get_filter(i)

            # Get the flux density
            flux = self.images_sed.get_photometry(i, unit="Jy", add_unit=False)

            # Get the observed flux
            observed_flux = self.observed_sed.photometry_for_filter(fltr, unit="Jy", add_unit=False)

            # Add to table
            self.sed_differences.add_from_filter_and_fluxes(fltr, flux, observed_flux)

    # -----------------------------------------------------------------

    def load_sed_differences(self):

        """
        This function ...
        :return:
        """

        # Message
        log.success("SED differences are already present: loading from file ...")

        # Load
        self.sed_differences = FluxDifferencesTable.from_file(self.sed_differences_filepath)

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

            if self.has_residuals_for_filter(fltr): residual = Frame.from_file(self.get_residuals_filepath_for_filter(fltr))
            else:

                # Debugging
                log.debug("Creating the residual frame for the '" + str(fltr) + "' filter ...")

                # Get the images in the same units
                simulated, observed, errors = convert_to_same_unit(self.images[fltr], self.observed_images[fltr], self.observed_errors[fltr])

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

            if self.has_weighed_for_filter(fltr): residual = Frame.from_file(self.get_weighed_filepath_for_filter(fltr))
            else:

                # Debugging
                log.debug("Creating the weighed residual frame for the '" + str(fltr) + "' filter ...")

                # Get the observed and simulated image
                #simulated = self.images[fltr]
                #observed = self.observed_images[fltr]
                #errors = self.observed_errors[fltr]

                # Get the images in the same units
                simulated, observed, errors = convert_to_same_unit(self.images[fltr], self.observed_images[fltr], self.observed_errors[fltr])

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

    @property
    def distribution_x_min(self):

        """
        This function ...
        :return:
        """

        return -4

    # -----------------------------------------------------------------

    @property
    def distribution_x_max(self):

        """
        This function ...
        :return:
        """

        return 4

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

            if self.has_residuals_distribution_for_filter(fltr): distribution = Distribution.from_file(self.get_residuals_distribution_filepath_for_filter(fltr))
            else:

                # Debugging
                log.debug("Creating the residuals distribution for the '" + str(fltr) + "' filter ...")

                # Get the values within the truncation ellipse
                values = self.residuals[fltr].values_in(self.truncation_ellipse)

                # REMOVE EXACT ZEROES
                indices = np.argwhere(values == 0)
                values = np.delete(values, indices)

                # REMOVE TOO LOW OR TOO HIGH (PROBABLY NOISE)
                indices = np.argwhere(values < self.distribution_x_min)
                values = np.delete(values, indices)
                indices = np.argwhere(values > self.distribution_x_max)
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

            if self.has_weighed_distribution_for_filter(fltr): distribution = Distribution.from_file(self.get_weighed_distribution_filepath_for_filter(fltr))
            else:

                # Debugging
                log.debug("Creating the residuals distribution for the '" + str(fltr) + "' filter ...")

                # Get the values within the truncation ellipse
                values = self.weighed[fltr].values_in(self.truncation_ellipse)

                # REMOVE EXACT ZEROES
                indices = np.argwhere(values == 0)
                values = np.delete(values, indices)

                # REMOVE TOO LOW OR TOO HIGH (PROBABLY NOISE)
                indices = np.argwhere(values < self.distribution_x_min)
                values = np.delete(values, indices)
                indices = np.argwhere(values > self.distribution_x_max)
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

        # 1. Write the weights
        if not self.has_weights: self.write_weights()

        # 2. Write the flux differences
        if not self.has_differences: self.write_differences()

        # Write the simulated datacube SED
        if not self.has_simulated_datacube_sed: self.write_simulated_datacube_sed()

        # 3. Write the image fluxes
        if not self.has_image_fluxes: self.write_image_fluxes()

        # Write proper image fluxes
        if not self.has_proper_image_fluxes: self.write_proper_image_fluxes()

        # 4. Write the image SED
        if not self.has_images_sed: self.write_image_sed()

        # 5. Write the image flux differences
        if not self.has_images_differences: self.write_image_differences()

        # 6. Write fluxes differences
        if not self.has_fluxes_differences: self.write_fluxes_differences()

        # 7. Write SED differences
        if not self.has_sed_differences: self.write_sed_differences()

        # 8. Write the images
        self.write_images()

        # 9. Write the residual frames
        self.write_residuals()

        # 10. Write the weighed residual frames
        self.write_weighed()

        # 11. Write the residual distributions
        self.write_residuals_distributions()

        # 12. Write the weighed residual distributions
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

    @property
    def weights_filepath(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.fluxes_path, "weights.dat")

    # -----------------------------------------------------------------

    @property
    def has_weights(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.weights_filepath)

    # -----------------------------------------------------------------

    def write_weights(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with weights ...")

        # Write the table with weights
        self.weights.saveto(self.weights_filepath)

    # -----------------------------------------------------------------

    @property
    def differences_filepath(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.fluxes_path, "differences.dat")

    # -----------------------------------------------------------------

    @property
    def has_differences(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.differences_filepath)

    # -----------------------------------------------------------------

    def write_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with the flux density differences ...")

        # Save the differences table
        self.differences.saveto(self.differences_filepath)

    # -----------------------------------------------------------------

    @property
    def simulated_datacube_sed_filepath(self):

        """
        THis function ...
        :return:
        """

        return fs.join(self.fluxes_path, "simulated_datacube_sed.dat")

    # -----------------------------------------------------------------

    @property
    def has_simulated_datacube_sed(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.simulated_datacube_sed_filepath)

    # -----------------------------------------------------------------

    def write_simulated_datacube_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the simulated datacube SED ...")

        # Write
        self.simulated_datacube_sed.saveto(self.simulated_datacube_sed_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def images_fluxes_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.evaluation_path, "fluxes_images")

    # -----------------------------------------------------------------

    @lazyproperty
    def images_fluxes_simulated_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.images_fluxes_path, "simulated")

    # -----------------------------------------------------------------

    @lazyproperty
    def images_fluxes_observed_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.images_fluxes_path, "observed")

    # -----------------------------------------------------------------

    @property
    def images_fluxes_filepath(self):

        """
        Thisn function ...
        :return:
        """

        #return fs.join(self.images_fluxes_path, "fluxes.dat")
        return fs.join(self.images_fluxes_simulated_path, "fluxes.dat")

    # -----------------------------------------------------------------

    @property
    def has_image_fluxes(self):

        """
        Thisn function ...
        :return:
        """

        return fs.is_file(self.images_fluxes_filepath)

    # -----------------------------------------------------------------

    def write_image_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing image fluxes ...")

        # Save
        self.images_fluxes.saveto(self.images_fluxes_filepath)

    # -----------------------------------------------------------------

    @property
    def proper_images_fluxes_filepath(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.images_fluxes_simulated_path, "fluxes_proper.dat")

    # -----------------------------------------------------------------

    @property
    def has_proper_image_fluxes(self):

        """
        Thisfunction ...
        :return:
        """

        return fs.is_file(self.proper_images_fluxes_filepath)

    # -----------------------------------------------------------------

    def write_proper_image_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing proper image fluxes ...")

        # Save
        self.proper_images_fluxes.saveto(self.proper_images_fluxes_filepath)

    # -----------------------------------------------------------------

    @property
    def images_sed_filepath(self):

        """
        This function ...
        :return:
        """

        #return fs.join(self.images_fluxes_path, "observed_sed.dat")
        return fs.join(self.images_fluxes_observed_path, "observed_sed.dat")

    # -----------------------------------------------------------------

    @property
    def has_images_sed(self):

        """
        Thisn function ...
        :return:
        """

        return fs.is_file(self.images_sed_filepath)

    # -----------------------------------------------------------------

    def write_image_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing image SED ...")

        # Save
        self.images_sed.saveto(self.images_sed_filepath)

    # -----------------------------------------------------------------

    @property
    def images_differences_filepath(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.images_fluxes_path, "differences.dat")

    # -----------------------------------------------------------------

    @property
    def has_images_differences(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.images_differences_filepath)

    # -----------------------------------------------------------------

    def write_image_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing image flux differences ...")

        # Save
        self.images_differences.saveto(self.images_differences_filepath)

    # -----------------------------------------------------------------

    @property
    def fluxes_differences_filepath(self):

        """
        This function ...
        :return:
        """

        #return fs.join(self.images_fluxes_path, "simulated_differences.dat")
        return fs.join(self.images_fluxes_simulated_path, "differences.dat")

    # -----------------------------------------------------------------

    @property
    def has_fluxes_differences(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.fluxes_differences_filepath)

    # -----------------------------------------------------------------

    def write_fluxes_differences(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing fluxes differences ...")

        # Write
        self.fluxes_differences.saveto(self.fluxes_differences_filepath)

    # -----------------------------------------------------------------

    @property
    def sed_differences_filepath(self):

        """
        Thisnf unction ...
        :return:
        """

        #return fs.join(self.images_fluxes_path, "observed_differences.dat")
        return fs.join(self.images_fluxes_observed_path, "differences.dat")

    # -----------------------------------------------------------------

    @property
    def has_sed_differences(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sed_differences_filepath)

    # -----------------------------------------------------------------

    def write_sed_differences(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing SED differences ...")

        # Write
        self.sed_differences.saveto(self.sed_differences_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def images_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.evaluation_path, "images")

    # -----------------------------------------------------------------

    def get_image_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.images_path, str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_image_for_filter(self, fltr):

        """
        This funtion ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_image_filepath_for_filter(fltr))

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

            # Already present?
            if self.has_image_for_filter(frame.filter): continue

            # Determine the path
            path = self.get_image_filepath_for_filter(frame.filter)

            # Save the frame
            frame.saveto(path)

    # -----------------------------------------------------------------

    def get_errors_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.images_path, str(fltr) + "_errors.fits")

    # -----------------------------------------------------------------

    def has_errors_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_errors_filepath_for_filter(fltr))

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

            # Already present?
            if self.has_errors_for_filter(errors.filter): continue

            # Determine the path
            path = self.get_errors_filepath_for_filter(errors.filter)

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

    def get_residuals_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.residuals_path, str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_residuals_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_residuals_filepath_for_filter(fltr))

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

            # If already present
            if self.has_residuals_for_filter(residuals.filter): continue

            # Determine the path
            path = self.get_residuals_filepath_for_filter(residuals.filter)

            # Save
            residuals.saveto(path)

    # -----------------------------------------------------------------

    def get_weighed_filepath_for_filter(self, fltr):

        """
        Thisn function ...
        :param fltr:
        :return:
        """

        return fs.join(self.weighed_path, str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_weighed_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_weighed_filepath_for_filter(fltr))

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

            # If already present
            if self.has_weighed_for_filter(residuals.filter): continue

            # Determine the path
            path = self.get_weighed_filepath_for_filter(residuals.filter)

            # Save
            residuals.saveto(path)

    # -----------------------------------------------------------------

    def get_residuals_distribution_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.residuals_path, str(fltr) + "_distribution.dat")

    # -----------------------------------------------------------------

    def has_residuals_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_residuals_distribution_filepath_for_filter(fltr))

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

            # check if already present
            if self.has_residuals_distribution_for_filter(fltr): continue

            # Get the distribution
            distribution = self.residuals_distributions[fltr]

            # Determine the path
            path = self.get_residuals_distribution_filepath_for_filter(fltr)

            # Save
            distribution.saveto(path)

    # -----------------------------------------------------------------

    def get_weighed_distribution_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.weighed_path, str(fltr) + "_distribution.dat")

    # -----------------------------------------------------------------

    def has_weighed_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_weighed_distribution_filepath_for_filter(fltr))

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

            # Check if already present
            if self.has_weighed_distribution_for_filter(fltr): continue

            # Get the distribution
            distribution = self.weighed_distributions[fltr]

            # Determine the path
            path = self.get_weighed_distribution_filepath_for_filter(fltr)

            # Save
            distribution.saveto(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # 1. Plot the simulated SED together with the observed SED
        if not self.has_simulated_sed_plot: self.plot_simulated_sed()

        # Plot the simulated datacube SED together with the observed SED
        if not self.has_simulated_datacube_sed_plot: self.plot_simulated_datacube_sed()

        # 2. Plot the simulated fluxes together with the observed SED
        if not self.has_simulated_fluxes_plot: self.plot_simulated_fluxes()

        # Plot differences between observed and simulated fluxes
        if not self.has_differences_plot: self.plot_differences()

        # Plot the images
        if not self.has_all_image_plots: self.plot_images()

        # Plot the relative error maps
        if not self.has_all_relative_error_plots: self.plot_relative_errors()

        # Plot the comparison between the simulated fluxes (from SED and from images)
        if not self.has_fluxes_plot: self.plot_fluxes()

        # Plot the SED of the simulated fluxes from images, and the simulated SED
        if not self.has_fluxes_only_images_plot: self.plot_fluxes_only_images()

        # Proper
        if not self.has_fluxes_proper_plot: self.plot_fluxes_proper()

        # Plot differences between simulated fluxes (SED and from images)
        if not self.has_fluxes_differences_plot: self.plot_fluxes_differences()

        # Plot the comparison between the observed fluxes (from SED and from images)
        if not self.has_sed_plot: self.plot_seds()

        # Plot the SED of the observed fluxes derived from images
        if not self.has_sed_only_images_plot: self.plot_seds_only_images()

        # Plot differences between observed fluxes (SED and from images)
        if not self.has_sed_differences_plot: self.plot_sed_differences()

        # Plot the comparison between simulated fluxes and observed fluxes, both from IMAGES, and the simulated SED
        if not self.has_seds_images_plot: self.plot_seds_images()

        # Plot differences between simulated fluxes and observed fluxes, both from IMAGES
        if not self.has_images_differences_plot: self.plot_image_differences()

        # Plot residuals distributions
        self.plot_residuals_distributions()

        # Plot weighed distributions
        self.plot_weighed_distributions()

        # Plot residual maps
        self.plot_residuals()

        # Plot weighed residual maps
        self.plot_weighed_residuals()

    # -----------------------------------------------------------------

    @property
    def simulated_sed_plot_filepath(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.fluxes_path, "simulated_sed.pdf")
    
    # -----------------------------------------------------------------

    @property
    def has_simulated_sed_plot(self):
        
        """
        This function ...
        :return: 
        """
        
        return fs.is_file(self.simulated_sed_plot_filepath)
    
    # -----------------------------------------------------------------

    def plot_simulated_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the simulated SED and observed SED ...")

        # Initialize the plotter
        plotter = SEDPlotter()

        # Ignore filters
        plotter.config.ignore_filters = self.ignore_sed_plot_filters

        #print(self.observed_sed)
        #print(self.simulated_sed)

        # Add the SEDs
        plotter.add_sed(self.observed_sed, "Observation (SED)")
        plotter.add_sed(self.simulated_sed, "Simulation")
        plotter.format = "pdf"

        # Plot
        plotter.run(output=self.simulated_sed_plot_filepath)

    # -----------------------------------------------------------------

    @property
    def simulated_datacube_sed_plot_filepath(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.fluxes_path, "simulated_datacube_sed.pdf")

    # -----------------------------------------------------------------

    @property
    def has_simulated_datacube_sed_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.simulated_datacube_sed_plot_filepath)

    # -----------------------------------------------------------------

    def plot_simulated_datacube_sed(self):

        """
        Thisjfnction ...
        :return:
        """

        # Inform the user
        log.info("Plotting the simulated SED based on the datacube, and the observed SED ...")

        # Initialize the plotter
        plotter = SEDPlotter()

        # Ignore filters
        plotter.config.ignore_filters = self.ignore_sed_plot_filters

        # Add the SEDs
        plotter.add_sed(self.observed_sed, "Observation (SED)")
        plotter.add_sed(self.simulated_datacube_sed, "Simulation")
        plotter.format = "pdf"

        # Plot
        plotter.run(output=self.simulated_datacube_sed_plot_filepath)

    # -----------------------------------------------------------------

    @property
    def simulated_fluxes_plot_filepath(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.fluxes_path, "simulated_fluxes.pdf")

    # -----------------------------------------------------------------

    @property
    def has_simulated_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.simulated_fluxes_plot_filepath)

    # -----------------------------------------------------------------

    def plot_simulated_fluxes(self):

        """
        This functio n...
        :return:
        """

        # Inform the user
        log.info("Plotting the mock observed fluxes and observed SED ...")

        # Initialize the plotter
        plotter = SEDPlotter()

        # Ignore filters
        plotter.config.ignore_filters = self.ignore_sed_plot_filters

        # Add the SEDs
        plotter.add_sed(self.observed_sed, "Observation (SED)") # ObservedSED
        plotter.add_sed(self.simulated_fluxes, "Simulation (mock observations)") # ObservedSED
        plotter.add_sed(self.simulated_sed, "Simulation") # SED
        plotter.format = "pdf"

        # Plot
        plotter.run(output=self.simulated_fluxes_plot_filepath)

    # -----------------------------------------------------------------

    @property
    def differences_plot_filepath(self):

        """
        This fucntion ...
        :return:
        """

        return fs.join(self.fluxes_path, "differences.pdf")

    # -----------------------------------------------------------------

    @property
    def has_differences_plot(self):

        """
        Thisj fucntion ...
        :return:
        """

        return fs.is_file(self.differences_plot_filepath)

    # -----------------------------------------------------------------

    def plot_differences(self):

        """
        Thins function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the differences between observed and mock observed fluxes ...")

        # Plot the differences
        plot_differences(self.differences, self.differences_plot_filepath)

    # -----------------------------------------------------------------

    def get_image_plot_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.images_path, str(fltr) + ".png")

    # -----------------------------------------------------------------

    def has_image_plot_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_image_plot_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    @property
    def has_all_image_plots(self):

        """
        This function ...
        :return:
        """

        #print("FILTERS:", self.simulated_filters_no_iras_planck)

        for fltr in self.simulated_filters_no_iras_planck:
            #print(fltr, self.get_image_plot_filepath_for_filter(fltr), "has image plot:", self.has_image_plot_for_filter(fltr))
            if not self.has_image_plot_for_filter(fltr): return False
        return True

    # -----------------------------------------------------------------

    def plot_images(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the mock observed images ...")

        # Loop over the filters
        for fltr in self.simulated_filters_no_iras_planck:

            # Check whether present
            if self.has_image_plot_for_filter(fltr): continue

            # Get the image
            frame = self.images[fltr]

            # Determine the path
            path = self.get_image_plot_filepath_for_filter(fltr)

            # Make plot
            vmin, vmax = frame.saveto_png(path, colours=self.config.colours,
                                              interval=self.config.interval,
                                              scale=self.config.scale, alpha=self.config.alpha_method,
                                              peak_alpha=self.config.peak_alpha)

    # -----------------------------------------------------------------

    def get_relative_errors_plot_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.images_path, str(fltr) + "_relerrors.png")

    # -----------------------------------------------------------------

    def has_relative_errors_plot_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_relative_errors_plot_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    @property
    def has_all_relative_error_plots(self):

        """
        This function ...
        :return:
        """

        #print("FILTERS:", self.simulated_filters_no_iras_planck)

        for fltr in self.simulated_filters_no_iras_planck:
            #print(fltr, self.get_relative_errors_plot_filepath_for_filter(fltr), "has relative errors plot:", self.has_relative_errors_plot_for_filter(fltr))
            if not self.has_relative_errors_plot_for_filter(fltr): return False
        return True

    # -----------------------------------------------------------------

    def plot_relative_errors(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the relative error maps of the mock observed images ...")

        # Loop over the filters
        for fltr in self.simulated_filters_no_iras_planck:

            # Check whether present
            if self.has_relative_errors_plot_for_filter(fltr): continue

            # Get the image and error map
            frame = self.images[fltr]
            errors = self.errors[fltr]
            relerrors = errors / frame

            # Determine the path
            path = self.get_relative_errors_plot_filepath_for_filter(fltr)

            # Plot
            plotting.plot_box(relerrors, path=path, colorbar=True, scale="linear")

    # -----------------------------------------------------------------

    @property
    def fluxes_plot_filepath(self):

        """
        This function ...
        :return:
        """

        #return fs.join(self.images_fluxes_path, "fluxes.pdf")
        return fs.join(self.images_fluxes_simulated_path, "fluxes.pdf")

    # -----------------------------------------------------------------

    @property
    def has_fluxes_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.fluxes_plot_filepath)

    # -----------------------------------------------------------------

    def plot_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the fluxes ...")

        # Initialize the plotter
        plotter = SEDPlotter()

        # Ignore filters
        plotter.config.ignore_filters = self.ignore_sed_plot_filters

        # Add the SEDs
        plotter.add_sed(self.images_fluxes, "Simulation (from mock observed images)") # ObservedSED
        plotter.add_sed(self.simulated_fluxes, "Simulation (mock observations)") # ObservedSED
        plotter.add_sed(self.simulated_sed, "Simulation") # SED
        plotter.format = "pdf"

        # Plot
        plotter.run(output=self.fluxes_plot_filepath)

    # -----------------------------------------------------------------

    @property
    def fluxes_only_images_plot_filepath(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.images_fluxes_simulated_path, "fluxes_only_images.pdf")

    # -----------------------------------------------------------------

    @property
    def has_fluxes_only_images_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.fluxes_only_images_plot_filepath)

    # -----------------------------------------------------------------

    def plot_fluxes_only_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the fluxes only from the mock observed images ...")

        # Initialize the plotter
        plotter = SEDPlotter()

        # Ignore filters
        plotter.config.ignore_filters = self.ignore_sed_plot_filters

        # Add the SEDs
        plotter.add_sed(self.images_fluxes, "Simulation (from mock observed images)")  # ObservedSED
        plotter.add_sed(self.simulated_sed, "Simulation")  # SED
        plotter.format = "pdf"

        # Plot
        plotter.run(output=self.fluxes_only_images_plot_filepath)

    # -----------------------------------------------------------------

    @property
    def fluxes_proper_plot_filepath(self):

        """
        Thisf untion ...
        :return:
        """

        return fs.join(self.images_fluxes_simulated_path, "fluxes_proper.pdf")

    # -----------------------------------------------------------------

    @property
    def has_fluxes_proper_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.fluxes_proper_plot_filepath)

    # -----------------------------------------------------------------

    def plot_fluxes_proper(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the fluxes based on the proper mock observed images ...")

        # Initialize the plotter
        plotter = SEDPlotter()

        # Ignore filters
        plotter.config.ignore_filters = self.ignore_sed_plot_filters

        # Add the SEDs
        plotter.add_sed(self.proper_images_fluxes, "Simulation (from proper mock observed images)")  # ObservedSED
        plotter.add_sed(self.simulated_sed, "Simulation")  # SED
        plotter.format = "pdf"

        # Plot
        plotter.run(output=self.fluxes_proper_plot_filepath)

    # -----------------------------------------------------------------

    @property
    def fluxes_differences_plot_filepath(self):

        """
        Thisf ucntion ...
        :return:
        """

        #return fs.join(self.images_fluxes_path, "fluxes_differences.pdf")
        return fs.join(self.images_fluxes_simulated_path, "differences.pdf")

    # -----------------------------------------------------------------

    @property
    def has_fluxes_differences_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.fluxes_differences_plot_filepath)

    # -----------------------------------------------------------------

    def plot_fluxes_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the differences between simulated fluxes (from mock observved images and from mock observations) ...")

        # Plot
        plot_differences(self.fluxes_differences, self.fluxes_differences_plot_filepath)

    # -----------------------------------------------------------------

    @property
    def sed_plot_filepath(self):

        """
        This fucntion ...
        :return:
        """

        #return fs.join(self.images_fluxes_path, "seds.pdf")
        return fs.join(self.images_fluxes_observed_path, "seds.pdf")

    # -----------------------------------------------------------------

    @property
    def has_sed_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sed_plot_filepath)

    # -----------------------------------------------------------------

    def plot_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SEDs ...")

        # Initialize the plotter
        plotter = SEDPlotter()

        # Ignore filters
        plotter.config.ignore_filters = self.ignore_sed_plot_filters

        # Add the SEDs
        plotter.add_sed(self.images_sed, "Observation (from images)")
        plotter.add_sed(self.observed_sed, "Observation (SED)")
        plotter.format = "pdf"

        # Plot
        plotter.run(output=self.sed_plot_filepath)

    # -----------------------------------------------------------------

    @property
    def sed_only_images_plot_filepath(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.images_fluxes_observed_path, "seds_only_images.pdf")

    # -----------------------------------------------------------------

    @property
    def has_sed_only_images_plot(self):

        """
        Thisnf unction ...
        :return:
        """

        return fs.is_file(self.sed_only_images_plot_filepath)

    # -----------------------------------------------------------------

    def plot_seds_only_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SED only from images ...")

        # Initialize the plotter
        plotter = SEDPlotter()

        # Ignore filters
        plotter.config.ignore_filters = self.ignore_sed_plot_filters

        # Add the SEDs
        plotter.add_sed(self.images_sed, "Observation (from images)")
        plotter.format = "pdf"

        # Plot
        plotter.run(output=self.sed_only_images_plot_filepath)

    # -----------------------------------------------------------------

    @property
    def sed_differences_plot_filepath(self):

        """
        This function ...
        :return:
        """

        #return fs.join(self.images_fluxes_path, "sed_differences.pdf")
        return fs.join(self.images_fluxes_observed_path, "differences.pdf")

    # -----------------------------------------------------------------

    @property
    def has_sed_differences_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.sed_differences_plot_filepath)

    # -----------------------------------------------------------------

    def plot_sed_differences(self):

        """
        Thisj function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the differences between the observed fluxes (from images and from observed SED) ...")

        # Plot
        plot_differences(self.sed_differences, self.sed_differences_plot_filepath)

    # -----------------------------------------------------------------

    @property
    def seds_images_plot_filepath(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.images_fluxes_path, "images.pdf")

    # -----------------------------------------------------------------

    @property
    def has_seds_images_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.seds_images_plot_filepath)

    # -----------------------------------------------------------------

    def plot_seds_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SEDs based on images (and with the simulated SED as reference) ...")

        # Initialize the plotter
        plotter = SEDPlotter()

        # Ignore filters
        plotter.config.ignore_filters = self.ignore_sed_plot_filters

        # Add the SEDs
        plotter.add_sed(self.images_sed, "Observation (from images)")
        plotter.add_sed(self.images_fluxes, "Simulation (from mock observed images)")
        plotter.add_sed(self.simulated_sed, "Simulation")
        plotter.format = "pdf"

        # Plot
        plotter.run(output=self.seds_images_plot_filepath)

    # -----------------------------------------------------------------

    @property
    def images_differences_plot_filepath(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.join(self.images_fluxes_path, "images_differences.pdf")

    # -----------------------------------------------------------------

    @property
    def has_images_differences_plot(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.images_differences_plot_filepath)

    # -----------------------------------------------------------------

    def plot_image_differences(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the flux differences based on images ...")

        # Plot
        plot_differences(self.images_differences, self.images_differences_plot_filepath)

    # -----------------------------------------------------------------

    def get_residuals_distribution_plot_filepath_for_filter(self, fltr):

        """
        Thijs function ...
        :param fltr:
        :return:
        """

        return fs.join(self.residuals_path, str(fltr) + "_distribution.pdf")

    # -----------------------------------------------------------------

    def has_residuals_distribution_plot_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_residuals_distribution_plot_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    @property
    def distribution_x_limits(self):

        """
        This function ...
        :return:
        """

        return (self.distribution_x_min, self.distribution_x_max) # because should be less than an order of magnitude (10)

    # -----------------------------------------------------------------

    def plot_residuals_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the residuals distributions ...")

        # Loop over the distributions
        for fltr in self.simulated_filters_no_iras_planck:

            # Check if already present
            if self.has_residuals_distribution_plot_for_filter(fltr): continue

            # Get the distribution
            distribution = self.residuals_distributions[fltr]

            # Determine the plot path
            path = self.get_residuals_distribution_plot_filepath_for_filter(fltr)

            # Plot
            distribution.plot(path=path, x_limits=self.distribution_x_limits)

    # -----------------------------------------------------------------

    def get_weighed_distribution_plot_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.weighed_path, str(fltr) + "_distribution.pdf")

    # -----------------------------------------------------------------

    def has_weighed_distribution_plot_for_filter(self, fltr):

        """
        Thisj function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_weighed_distribution_plot_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    @property
    def weighed_distribution_x_limits(self):

        """
        This function ...
        :return:
        """

        return None

    # -----------------------------------------------------------------

    def plot_weighed_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the weighed residuals distributions ...")

        # Loop over the weighed distributions
        for fltr in self.simulated_filters_no_iras_planck:

            # Check if already present
            if self.has_weighed_distribution_plot_for_filter(fltr): continue

            # Get the distribution
            distribution = self.weighed_distributions[fltr]

            # Determine the plot path
            path = self.get_weighed_distribution_plot_filepath_for_filter(fltr)

            # Plot
            distribution.plot(path=path, x_limits=self.weighed_distribution_x_limits)

    # -----------------------------------------------------------------

    def get_residuals_plot_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.residuals_path, str(fltr) + ".png")

    # -----------------------------------------------------------------

    def has_residuals_plot_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_residuals_plot_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    @property
    def residuals_plot_interval(self):

        """
        This functio n...
        :return:
        """

        return [-100., 100.]

    # -----------------------------------------------------------------

    def plot_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the residuals ...")

        # Loop over the filters
        for fltr in self.simulated_filters_no_iras_planck:

            # Check whether present
            if self.has_residuals_plot_for_filter(fltr): continue

            # Get the residual map
            residuals = self.residuals[fltr]

            # Get residuals in percentage
            residuals_percentage = residuals * 100

            # Determine the path
            path = self.get_residuals_plot_filepath_for_filter(fltr)

            # # Make plot
            # vmin, vmax = residuals.saveto_png(path, colours=self.config.colours,
            #                                   interval=self.config.interval,
            #                                   scale=self.config.scale, alpha=self.config.alpha_method,
            #                                   peak_alpha=self.config.peak_alpha)

            # Plot
            plotting.plot_box(residuals_percentage, interval=self.residuals_plot_interval, path=path, colorbar=True, around_zero=True, scale="linear")

    # -----------------------------------------------------------------

    def get_weighed_plot_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.weighed_path, str(fltr) + ".png")

    # -----------------------------------------------------------------

    def has_weighed_plot_for_filter(self, fltr):

        """
        Thins function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_weighed_plot_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def plot_weighed_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the weighed residuals ...")

        # Loop over the filters
        for fltr in self.simulated_filters_no_iras_planck:

            # Ceck whether present
            if self.has_weighed_plot_for_filter(fltr): continue

            # Get the weighed residual map
            residuals = self.weighed[fltr]

            # Get residuals in percentage
            #residuals_percentage = residuals * 100

            # Determine the path
            path = self.get_weighed_plot_filepath_for_filter(fltr)

            # # Make plot
            # vmin, vmax = residuals.saveto_png(path, colours=self.config.colours,
            #                                   interval=self.config.interval,
            #                                   scale=self.config.scale, alpha=self.config.alpha_method,
            #                                   peak_alpha=self.config.peak_alpha)

            # Plot
            plotting.plot_box(residuals, path=path, colorbar=True, around_zero=True, scale="linear")

# -----------------------------------------------------------------

def plot_differences(table, path):

    """
    Thisfunction ...
    :param table:
    :param path:
    :return:
    """

    # Debugging
    log.debug("Plotting flux differences ...")

    import numpy as np
    import matplotlib.pyplot as plt

    # Initialize figure
    fig, ax1 = plt.subplots()

    # Set logarithmic
    ax1.set_xscale("log")

    # Set x label
    ax1.set_xlabel('Wavelength [micron]')

    # Get the filters, sorted on wavelength
    filters = table.filters(sort=True)
    wavelengths = [fltr.wavelength.to("micron").value for fltr in filters]

    # Get the arrays of differences and relative differences
    #differences = [table.difference_for_filter(fltr) for fltr in filters]
    rel_differences = [table.relative_difference_for_filter(fltr) * 100 for fltr in filters]

    # # Plot
    # ax1.plot(wavelengths, differences, 'b-')
    #
    # # Make the y-axis label, ticks and tick labels match the line color.
    # ax1.set_ylabel('Differences', color='b')
    # ax1.tick_params('y', colors='b')

    # Create twin axis
    #ax2 = ax1.twinx()

    # Plot relative differences
    ax1.plot(wavelengths, rel_differences, 'b-')

    # Set label
    ax1.set_ylabel('Relative differences (\%)')
    #ax1.tick_params('y', colors='r')

    # Save the figure
    fig.tight_layout()
    plt.savefig(path)

    # Close
    plt.close()

# -----------------------------------------------------------------
