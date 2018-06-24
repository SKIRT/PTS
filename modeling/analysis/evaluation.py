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
from ..fitting.modelanalyser import FluxDifferencesTable
from ...core.tools.utils import lazyproperty
from ...magic.core.list import FrameList
from ...core.tools import sequences
from ...core.basics.distribution import Distribution
from ...core.basics.containers import FilterBasedList
from ...core.data.sed import SED, ObservedSED
from ...magic.core.frame import Frame
from ...core.plot.sed import SEDPlotter
from ...magic.core.list import convert_to_same_unit, uniformize
from ...magic.tools import plotting
from ...magic.core.list import rebin_to_highest_pixelscale
from ...core.misc.fluxes import ObservedFluxCalculator
from ...core.misc.images import ObservedImageMaker
from ...core.units.parsing import parse_unit as u
from ...core.plot.distribution import plot_distribution

# -----------------------------------------------------------------

earth_name = "earth"

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

        # The reference SED
        self.reference_sed = None

        # The mock fluxes
        self.fluxes = None

        # Initialize the differences table
        self.differences = FluxDifferencesTable()

        # The simulated SED derived from the datacube
        self.simulated_datacube_sed = None

        # The mock observed images and their errors
        self.images = FrameList()
        #self.errors = FrameList()

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

        # The proper residual frames and weighed
        self.proper_residuals = FrameList()
        self.proper_weighed = FrameList()

        # The distributions
        self.residuals_distributions = FilterBasedList()
        self.weighed_distributions = FilterBasedList()

        # Proper distributions
        self.proper_residuals_distributions = FilterBasedList()
        self.proper_weighed_distributions = FilterBasedList()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # FLUXES

        # Get the observed SED
        self.get_reference_sed()

        # Get the mock fluxes
        self.get_fluxes()

        # 2. Get the differences
        self.get_differences()

        # 4. Get the simulated SED
        self.get_simulated_datacube_sed()

        # IMAGES

        # 5. Get images
        self.get_images()

        # 6. Get images created with spatial/spectral convolution etc.
        self.get_proper_images()

        # 7. Load observed images
        self.load_observed_images()

        # 8. Rebin the images to the same pixelscale
        self.rebin_images()

        # FLUXES FROM IMAGES

        # 9. Get image fluxes
        self.get_image_fluxes()

        # 10. Get proper image fluxes
        self.get_proper_image_fluxes()

        # 11. Get images SED
        self.get_images_sed()

        # 12. Get image differences
        self.get_image_differences()

        # Get fluxes differences
        self.get_fluxes_differences()

        # Get SED differences
        self.get_sed_differences()

        # RESIDUAL IMAGES

        # Get the residuals with the images
        self.get_residuals()

        # Get the residuals with the proper images
        self.get_proper_residuals()

        # 17. Calculate the weighed residual images
        self.get_weighed()

        # Calculate the weighed proper residual simages
        self.get_proper_weighed()

        # 18. Create distributions of the residual values
        self.get_residuals_distributions()

        # Create distributions of proper residual values
        self.get_proper_residuals_distributions()

        # 19. Create distributions of the weighed residual values
        self.get_weighed_distributions()

        # Create distributions of the proper weighed residual values
        self.get_proper_weighed_distributions()

        # 20. Writing
        self.write()

        # 21. Plotting
        if self.config.plot: self.plot()

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
    def evaluation_path(self):

        """
        Thisf unction ...
        :return:
        """

        return self.analysis_run.evaluation_path

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
        Thisnfunction ...
        :return:
        """

        return self.fluxes

    # -----------------------------------------------------------------

    @property
    def simulated_flux_filter_names(self):

        """
        This function ...
        :return:
        """

        return self.observed_filter_names_in_range_without_iras

    # -----------------------------------------------------------------

    @property
    def simulated_flux_filters(self):

        """
        This function ...
        :return:
        """

        return self.observed_filters_in_range_without_iras

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

    def get_reference_sed(self):
        
        """
        Thisf unction ...
        :return: 
        """
        
        if self.has_reference_sed: self.load_reference_sed()
        else: self.create_reference_sed()

    # -----------------------------------------------------------------

    def load_reference_sed(self):

        """
        This function ...
        :return:
        """

        # Load
        self.reference_sed = ObservedSED.from_file(self.reference_sed_path)

    # -----------------------------------------------------------------

    def create_reference_sed(self):

        """
        Thisf unction ...
        :return: 
        """

        # Inform the user
        log.info("Getting the reference observed SED ...")

        # Load the appropriate observed SED
        if self.config.not_clipped: self.reference_sed = self.truncated_sed.copy()
        else: self.reference_sed = self.observed_sed.copy()

        # Add additional relative error?
        if self.config.additional_error is not None: self.reference_sed.add_relative_error(self.config.additional_error)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filters_in_range(self):

        """
        This function ...
        :return:
        """

        filters = []
        for fltr in self.observed_filters:
            if not self.analysis_run.wavelength_grid.covers(fltr.wavelength):
                log.warning("The '" + str(fltr) + "' filter is not covered by the wavelength range: not making observations for this filter")
                continue
            filters.append(fltr)
        return filters

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filters_in_range_without_iras(self):

        """
        This function ...
        :return:
        """

        return [fltr for fltr in self.observed_filters_in_range if fltr not in self.iras_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filter_names_in_range(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.observed_filters_in_range]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filter_names_in_range_without_iras(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.observed_filters_in_range_without_iras]

    # -----------------------------------------------------------------

    def get_fluxes(self):

        """
        This function ...
        :return:
        """

        if self.has_fluxes: self.load_fluxes()
        else: self.create_fluxes()

    # -----------------------------------------------------------------

    def load_fluxes(self):

        """
        This function ...
        :return:
        """

        self.fluxes = ObservedSED.from_file(self.fluxes_filepath)

    # -----------------------------------------------------------------

    def create_fluxes(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Calculating the observed fluxes ...")

        # Create a ObservedFluxCalculator object
        flux_calculator = ObservedFluxCalculator()

        # Set spectral convolution flag
        flux_calculator.config.spectral_convolution = True #self.misc_options.fluxes_spectral_convolution

        # Set plot flags
        flux_calculator.config.plot = False
        flux_calculator.config.plot_seds = False
        flux_calculator.config.plot_images = False

        # EXTRA OPTIONS
        flux_calculator.config.check_wavelengths = True #self.config.check_wavelengths
        flux_calculator.config.ignore_bad = True #self.config.ignore_bad
        flux_calculator.config.skip_ignored_bad_convolution = False #self.config.skip_ignored_bad_convolution
        flux_calculator.config.skip_ignored_bad_closest = False #self.config.skip_ignored_bad_closest

        # DEPLOYMENT
        #self.flux_calculator.config.deploy_pts = self.config.deploy_pts
        #self.flux_calculator.config.update_dependencies = self.config.update_dependencies

        # Set input
        input_dict = dict()
        #input_dict["simulation"] = self.simulation
        input_dict["simulation_output_path"] = self.analysis_run.total_output_path
        #input_dict["output_path"] = self.fluxes_output_path
        input_dict["filter_names"] = self.simulated_flux_filters #self.filters_for_fluxes
        input_dict["instrument_names"] = [earth_name]
        #input_dict["instrument_names"] = self.misc_options.observation_instruments
        #input_dict["errors"] = self.misc_options.flux_errors
        #input_dict["no_spectral_convolution_filters"] = self.misc_options.no_fluxes_spectral_convolution_filters
        #input_dict["reference_seds"] = self.fluxes_reference_seds

        # Run
        flux_calculator.run(**input_dict)

        # Get the SED
        self.fluxes = flux_calculator.mock_seds[earth_name]

    # -----------------------------------------------------------------

    def get_differences(self):

        """
        This function ...
        :return:
        """

        # 3. Calculate flux differences
        if not self.has_differences: self.calculate_differences()
        else: self.load_differences()

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
            #observed_fluxdensity = self.observed_sed.photometry_for_filter(fltr, unit="Jy", add_unit=False)
            observed_fluxdensity = self.reference_sed.photometry_for_filter(fltr, unit="Jy", add_unit=False)

            # Find the corresponding flux error in the SED derived from observation
            #observed_fluxdensity_error = self.observed_sed.error_for_band(instrument, band, unit="Jy").average.to("Jy").value
            #observed_fluxdensity_error = self.observed_sed.error_for_filter(fltr, unit="Jy", add_unit=False).average

            # If no match with (instrument, band) is found in the observed SED
            if observed_fluxdensity is None:
                log.warning("The observed flux density could not be found for the " + str(fltr) + " filter")
                continue

            # Calculate the difference
            difference = fluxdensity - observed_fluxdensity
            relative_difference = difference / observed_fluxdensity

            # Find the index of the current band in the weights table
            #index = self.weights.index_for_filter(fltr, return_none=True)
            #if index is None:
            #    log.warning("A weight is not calculated for the " + str(fltr) + " filter")
            #    continue # Skip this band if a weight is not found

            # Get the weight
            #weight = self.weights["Weight"][index]

            # Calculate the chi squared term
            #chi_squared_term = weight * difference ** 2 / observed_fluxdensity_error ** 2

            # Add an entry to the differences table
            self.differences.add_entry(fltr.instrument, fltr.band, difference, relative_difference) #chi_squared_term)

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

    def get_simulated_datacube_sed(self):

        """
        This function ...
        :return:
        """

        # 5. Create simulated SED derived from the datacubes
        if not self.has_simulated_datacube_sed: self.create_simulated_datacube_sed()
        else: self.load_simulated_datacube_sed()

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

    def get_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the mock observed images ...")

        # Loop over the filters for which we want to create images
        for fltr in self.simulated_flux_filters:

            # Make image
            frame = self.get_image_for_filter(fltr)
            self.images.append(frame)

            # Make error map
            #errors = self.get_errors_for_filter(fltr, frame)
            #self.errors.append(errors)

            # Save the image
            if not self.has_image_for_filter(fltr):

                # Determine the path
                path = self.get_image_filepath_for_filter(frame.filter)

                # Save the frame
                frame.saveto(path)

            # Save the error map
            #if not self.has_errors_for_filter(fltr):
                # Determine the path
                #path = self.get_errors_filepath_for_filter(errors.filter)
                # Save the error map
                #errors.saveto(path)

    # -----------------------------------------------------------------

    def get_image_for_filter(self, fltr):

        """
        Thisn function ...
        :param fltr:
        :return:
        """

        # Present?
        if self.has_image_for_filter(fltr):

            # Needed?
            #if self.has_residuals_for_filter(fltr) and self.has_weighed_for_filter(fltr) and self.has_image_fluxes and self.has_all_image_plots and self.has_all_relative_error_plots: continue

            # Success
            log.success("The '" + str(fltr) + "' image is already present: loading ...")

            # Load the frame
            frame = Frame.from_file(self.get_image_filepath_for_filter(fltr))

        # Not yet present
        else:

            # Debugging
            log.debug("Making the mock observed '" + str(fltr) + "' image ...")

            # Get the appropriate simulated frame
            frame = self.get_frame_for_filter(fltr, convolve=False)

            # Print wavelengths to check
            # print(fltr.wavelength, fltr.center, fltr.pivot, frame._wavelength)

            # Give warning
            rel_difference = (frame.wavelength - fltr.pivot) / fltr.pivot
            if rel_difference > 0.05: log.warning("Wavelength of filter: " + str(fltr.pivot) + ", wavelength of frame: " + str(frame._wavelength))

            # Print unit BEFORE
            # print("UNIT BEFORE:", frame.unit)

            # Convert to non-brightness
            #conversion_factor = frame.convert_to_corresponding_non_brightness_unit()
            conversion_factor = frame.convert_to("Jy")

            # print("PIXELSCALE:", frame.average_pixelscale)
            # print("CONVERSION FACTOR", conversion_factor)

            # Print unit AFTER
            # print("UNIT AFTER:", frame.unit)

        # Return
        return frame

    # -----------------------------------------------------------------

    def get_errors_for_filter(self, fltr, frame):

        """
        This function ...
        :param fltr:
        :param frame:
        :return:
        """

        # Present?
        if self.has_errors_for_filter(fltr):

            # Needed?
            #if self.has_residuals_for_filter(fltr) and self.has_weighed_for_filter(fltr) and self.has_all_relative_error_plots: continue

            # Success
            log.success("The '" + str(fltr) + "' error map is already present: loading ...")

            # Load the error frame
            errors = Frame.from_file(self.get_errors_filepath_for_filter(fltr))

        # Not yet present
        else:

            # Debugging
            log.debug("Making an approximate error map for the '" + str(fltr) + "' image ...")

            # Create an approximate error frame
            errors = Frame(np.sqrt(frame.data), wcs=frame.wcs, filter=frame.filter, unit=frame.unit)

        # Return the error map
        return errors

    # -----------------------------------------------------------------

    def get_proper_images(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Getting the proper mock observed images ...")

        # Keep list of filters for which the image still has to be made
        filters = []

        # Loop over the filters for which we want to create images
        for fltr in self.simulated_flux_filters:

            # Has image?
            if self.has_proper_image_for_filter(fltr):
                frame = self.load_proper_image_for_filter(fltr)
                self.proper_images.append(frame, fltr=fltr)
            else: filters.append(fltr)

        # Any filters yet to create images for?
        nfilters = len(filters)
        if nfilters > 0:

            # Make the proper images
            images = self.make_proper_images(filters)

            # Add the images
            for fltr in images: self.proper_images.append(images[fltr], fltr=fltr)

    # -----------------------------------------------------------------

    def load_proper_image_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return Frame.from_file(self.get_proper_image_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def make_proper_images(self, filters):

        """
        This function ...
        :param filters:
        :return:
        """

        # Create the maker
        maker = ObservedImageMaker()

        # Set options
        maker.config.spectral_convolution = True

        # Write intermediate results
        maker.config.write_intermediate = True # write_intermediate_images
        maker.config.write_kernels = True # write_convolution_kernels

        # Set number of processes to one
        maker.config.nprocesses_local = 1

        # Settings for convolution
        maker.config.check_wavelengths = True
        maker.config.ignore_bad = True
        maker.config.skip_ignored_bad_convolution = False
        maker.config.skip_ignored_bad_closest = False

        # No plotting
        maker.config.plot = False

        # Group, but write earth instrument output into
        maker.config.group = True

        # Set input
        input_dict = dict()

        # The output path of the simulation
        input_dict["simulation_output_path"] = self.analysis_run.total_output_path

        # The output path for the images
        input_dict["output_path"] = self.proper_images_path
        input_dict["output_paths_instruments"] = {earth_name: self.proper_images_path}

        # Filters and instruments
        input_dict["filters"] = filters
        input_dict["instrument_names"] = [earth_name]

        # Set coordinate system of the datacube
        input_dict["wcs_path"] = self.analysis_run.reference_map_path
        input_dict["wcs_instrument"] = earth_name

        # Unit conversion
        input_dict["unit"] = u("Jy") #self.misc_options.images_unit

        # Convolution
        input_dict["auto_psfs"] = True
        #input_dict["kernel_paths"] = self.misc_options.images_kernels
        input_dict["fwhms_dataset"] = self.photometry_dataset #self.misc_options.fwhms_dataset # path or dataset is possible

        # Set dataset for rebinning
        input_dict["rebin_dataset"] = self.photometry_dataset  # path or dataset is possible
        input_dict["rebin_instrument"] = earth_name

        # NO SPECTRAL CONVOLUTION FOR CERTAIN IMAGES?
        input_dict["no_spectral_convolution_filters"] = self.planck_filters

        # Run
        maker.run(**input_dict)

        # Return the images
        return maker.images[earth_name] if earth_name in maker.images else {} # if all images were already made

    # -----------------------------------------------------------------

    def load_proper_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the proper mock observed images ...")

        # Loop over the filters for which we want to create images
        for fltr in self.simulated_flux_filters:
            
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
        image_names_for_filters = self.dataset.get_names_for_filters(self.simulated_flux_filters)
        for fltr, image_name in zip(self.simulated_flux_filters, image_names_for_filters):

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
        for fltr in self.simulated_flux_filters:

            # All already present
            if self.has_residuals_for_filter(fltr) and self.has_weighed_for_filter(fltr) and self.has_images_sed and self.has_image_fluxes: continue

            # Debugging
            log.debug("Rebinning the '" + str(fltr) + "' images if necessary ...")

            # Get the frames
            simulated = self.images[fltr]
            #simulated_errors = self.errors[fltr]
            observed = self.observed_images[fltr]
            observed_errors = self.observed_errors[fltr]

            # Rebin in-place
            names = ["simulated", "simulated_errors", "observed", "observed_errors"]
            #rebin_to_highest_pixelscale(simulated, simulated_errors, observed, observed_errors, names=names, in_place=True)
            rebin_to_highest_pixelscale(simulated, observed, observed_errors, names=names, in_place=True)

    # -----------------------------------------------------------------

    def get_image_fluxes(self):

        """
        This function ...
        :return:
        """

        # 10. Calculate fluxes from the images
        if not self.has_image_fluxes: self.calculate_image_fluxes()
        else: self.load_image_fluxes()

    # -----------------------------------------------------------------

    def calculate_image_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the fluxes based on the mock observed images ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

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

    def get_proper_image_fluxes(self):

        """
        This function ...
        :return:
        """

        # Calculate
        if not self.has_proper_image_fluxes: self.calculate_proper_image_fluxes()

        # Load
        else: self.load_proper_image_fluxes()

    # -----------------------------------------------------------------

    def calculate_proper_image_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the fluxes based on the proper mock observed images ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            #print(self.proper_images.filters)

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

    def get_images_sed(self):

        """
        This function ...
        :return:
        """

        # 12. Calculate SED
        if not self.has_images_sed: self.calculate_images_sed()
        else: self.load_images_sed()

    # -----------------------------------------------------------------

    def calculate_images_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the fluxes based on the observed images ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

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

    def get_image_differences(self):

        """
        This function ...
        :return:
        """

        # 13. Calcualte differences between fluxes calculated from the images and the observed fluxes
        if not self.has_images_differences: self.calculate_image_differences()
        else: self.load_image_differences()

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

    def get_fluxes_differences(self):

        """
        This function ...
        :return:
        """

        # 14. Calculate the differences between the simulated fluxes and the fluxes from the observed images
        if not self.has_fluxes_differences: self.calculate_fluxes_differences()
        else: self.load_fluxes_differences()

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

    def get_sed_differences(self):

        """
        This function ...
        :return:
        """

        # 15. Calculate the differences between the observed fluxes and the observed fluxes from the images
        if not self.has_sed_differences: self.calculate_sed_differences()
        else: self.load_sed_differences()

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
            #observed_flux = self.observed_sed.photometry_for_filter(fltr, unit="Jy", add_unit=False)
            observed_flux = self.reference_sed.photometry_for_filter(fltr, unit="Jy", add_unit=False)

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

    def get_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the residual images ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Load
            if self.has_residuals_for_filter(fltr): residuals = self.load_residuals_for_filter(fltr)

            # Calculate
            else: residuals = self.calculate_residuals_for_filter(fltr)

            # Add the residual image
            self.residuals.append(residuals)

    # -----------------------------------------------------------------

    def load_residuals_for_filter(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

        return Frame.from_file(self.get_residuals_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def calculate_residuals_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

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
        significance_mask = self.get_significance_mask(observed, errors, min_npixels=self.config.min_npixels,
                                                       connectivity=self.config.connectivity)

        # MASK
        residual[truncation_mask] = 0.0
        residual[significance_mask] = 0.0

        # Return
        return residual

    # -----------------------------------------------------------------

    def get_proper_residuals(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Getting the proper residual images ... ")
        
        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Load
            if self.has_proper_residuals_for_filter(fltr): residuals = self.load_proper_residuals_for_filter(fltr)

            # Calculate
            else: residuals = self.calculate_proper_residuals_for_filter(fltr)

            # Add the residuals image
            self.proper_residuals.append(residuals)

    # -----------------------------------------------------------------

    def load_proper_residuals_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return Frame.from_file(self.get_proper_residuals_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def calculate_proper_residuals_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Debugging
        log.debug("Creating the proper residual frame for the '" + str(fltr) + "' filter ...")

        # Get the images in the same units
        #simulated, observed, errors = convert_to_same_unit(self.proper_images[fltr], self.observed_images[fltr], self.observed_errors[fltr])
        simulated, observed, errors = uniformize(self.proper_images[fltr], self.observed_images[fltr], self.observed_errors[fltr], convolve=False)

        # Calculate the residual image
        residual = (simulated - observed) / observed

        # Set the filter
        residual.filter = fltr

        # Replace infs
        residual.replace_infs(0.0)

        # Get the truncation mask
        truncation_mask = self.get_truncation_mask(observed.wcs)

        # Get the significance mask
        significance_mask = self.get_significance_mask(observed, errors, min_npixels=self.config.min_npixels,
                                                       connectivity=self.config.connectivity)

        # MASK
        residual[truncation_mask] = 0.0
        residual[significance_mask] = 0.0

        # Return
        return residual

    # -----------------------------------------------------------------

    def get_weighed(self):

        """
        This ufnction ...
        :return:
        """

        # Inform the user
        log.info("Getting the weighed residual images ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Load
            if self.has_weighed_for_filter(fltr): residuals = self.load_weighed_residuals_for_filter(fltr)

            # Calculate
            else: residuals = self.calculate_weighed_residuals_for_filter(fltr)

            # Add the weighed residual image
            self.weighed.append(residuals)

    # -----------------------------------------------------------------

    def load_weighed_residuals_for_filter(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

        return Frame.from_file(self.get_weighed_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def calculate_weighed_residuals_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Debugging
        log.debug("Creating the weighed residual frame for the '" + str(fltr) + "' filter ...")

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
        significance_mask = self.get_significance_mask(observed, errors, min_npixels=self.config.min_npixels,
                                                       connectivity=self.config.connectivity)

        # MASK
        residual[truncation_mask] = 0.0
        residual[significance_mask] = 0.0

        # Return
        return residual

    # -----------------------------------------------------------------

    def get_proper_weighed(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Getting the proper weighed residual images ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Load
            if self.has_proper_weighed_for_filter(fltr): residuals = self.load_proper_weighed_residuals_for_filter(fltr)

            # Calculate
            else: residuals = self.calculate_proper_weighed_residuals_for_filter(fltr)

            # Add the weighed residual image
            self.proper_weighed.append(residuals)

    # -----------------------------------------------------------------

    def load_proper_weighed_residuals_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return Frame.from_file(self.get_proper_weighed_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def calculate_proper_weighed_residuals_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Debugging
        log.debug("Creating the proper weighed residual frame for the '" + str(fltr) + "' filter ...")

        # Get the images in the same units
        #simulated, observed, errors = convert_to_same_unit(self.proper_images[fltr], self.observed_images[fltr], self.observed_errors[fltr])
        simulated, observed, errors = uniformize(self.proper_images[fltr], self.observed_images[fltr], self.observed_errors[fltr], convolve=False)

        # Calculate the weighed residual image
        residual = (simulated - observed) / errors

        # Set the filter
        residual.filter = fltr

        # Replace infs
        residual.replace_infs(0.0)

        # Get the truncation mask
        truncation_mask = self.get_truncation_mask(observed.wcs)

        # Get the significance mask
        significance_mask = self.get_significance_mask(observed, errors, min_npixels=self.config.min_npixels,
                                                       connectivity=self.config.connectivity)

        # MASK
        residual[truncation_mask] = 0.0
        residual[significance_mask] = 0.0

        # Return
        return residual

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

    def get_residuals_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the residuals distributions ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Load
            if self.has_residuals_distribution_for_filter(fltr): distribution = self.load_residuals_distribution_for_filter(fltr)

            # Create
            else: distribution = self.create_residuals_distribution_for_filter(fltr)

            # Add the distribution
            if distribution is not None: self.residuals_distributions.append(fltr, distribution)

    # -----------------------------------------------------------------

    def load_residuals_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return Distribution.from_file(self.get_residuals_distribution_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def create_residuals_distribution_for_filter(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

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
            return None

        # Create distribution
        distribution = Distribution.from_values("Residual", values, nbins=self.config.nbins)

        # Reutrn
        return distribution

    # -----------------------------------------------------------------

    def get_proper_residuals_distributions(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Getting the proper residuals distributions ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Load
            if self.has_proper_residuals_distribution_for_filter(fltr): distribution = self.load_proper_residuals_distribution_for_filter(fltr)

            # Create
            else: distribution = self.create_proper_residuals_distribution_for_filter(fltr)

            # Add the distribution
            if distribution is not None: self.proper_residuals_distributions.append(fltr, distribution)

    # -----------------------------------------------------------------

    def load_proper_residuals_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return Distribution.from_file(self.get_proper_residuals_distribution_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def create_proper_residuals_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Debugging
        log.debug("Creating the proper residuals distribution for the '" + str(fltr) + "' filter ...")

        # Get the values within the truncation ellipse
        values = self.proper_residuals[fltr].values_in(self.truncation_ellipse)

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
            log.error("Cannot create proper distribution for the '" + str(fltr) + "' filter")
            return None

        # Create distribution
        distribution = Distribution.from_values("Residual", values, nbins=self.config.nbins)

        # Reutrn
        return distribution

    # -----------------------------------------------------------------

    def get_weighed_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the weighed residuals distributions ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Load
            if self.has_weighed_distribution_for_filter(fltr): distribution = self.load_weighed_distribution_for_filter(fltr)

            # Create
            else: distribution = self.create_weighed_distribution_for_filter(fltr)

            # Add the distribution
            if distribution is not None: self.weighed_distributions.append(fltr, distribution)

    # -----------------------------------------------------------------

    def load_weighed_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return Distribution.from_file(self.get_weighed_distribution_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def create_weighed_distribution_for_filter(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

        # Debugging
        log.debug("Creating the weighed residuals distribution for the '" + str(fltr) + "' filter ...")

        # Get the values within the truncation ellipse
        values = self.weighed[fltr].values_in(self.truncation_ellipse)

        # REMOVE EXACT ZEROES
        indices = np.argwhere(values == 0)
        values = np.delete(values, indices)

        # # REMOVE TOO LOW OR TOO HIGH (PROBABLY NOISE)
        # indices = np.argwhere(values < self.distribution_x_min)
        # values = np.delete(values, indices)
        # indices = np.argwhere(values > self.distribution_x_max)
        # values = np.delete(values, indices)

        # Check
        if len(values) == 0 or sequences.all_equal(values):
            log.error("Cannot create distribution for the '" + str(fltr) + "' filter")
            return None

        # Create distribution
        distribution = Distribution.from_values("Residual", values, nbins=self.config.nbins)

        # Retunr the distribution
        return distribution

    # -----------------------------------------------------------------

    def get_proper_weighed_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the proper weighed residuals distributions ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Load
            if self.has_proper_weighed_distribution_for_filter(fltr): distribution = self.load_proper_weighed_distribution_for_filter(fltr)

            # Create
            else: distribution = self.create_proper_weighed_distribution_for_filter(fltr)

            # Add the distribution
            if distribution is not None: self.proper_weighed_distributions.append(fltr, distribution)

    # -----------------------------------------------------------------

    def load_proper_weighed_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return Distribution.from_file(self.get_proper_weighed_distribution_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def create_proper_weighed_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Debugging
        log.debug("Creating the proper weighed residuals distribution for the '" + str(fltr) + "' filter ...")

        # Get the values within the truncation ellipse
        values = self.proper_weighed[fltr].values_in(self.truncation_ellipse)

        # REMOVE EXACT ZEROES
        indices = np.argwhere(values == 0)
        values = np.delete(values, indices)

        # Check
        if len(values) == 0 or sequences.all_equal(values):
            log.error("Cannot create distribution for the '" + str(fltr) + "' filter")
            return None

        # Create distribution
        distribution = Distribution.from_values("Residual", values, nbins=self.config.nbins)

        # Retunr the distribution
        return distribution

    # -----------------------------------------------------------------

    @property
    def do_write_reference_sed(self):

        """
        This function ...
        :return:
        """

        return not self.has_reference_sed

    # -----------------------------------------------------------------

    @property
    def do_write_fluxes(self):

        """
        Thisf unction ...
        :return:
        """

        return not self.has_fluxes

    # -----------------------------------------------------------------

    @property
    def do_write_differences(self):

        """
        Thisf unction ...
        :return:
        """

        return not self.has_differences

    # -----------------------------------------------------------------

    @property
    def do_write_simulated_datacube_sed(self):

        """
        This function ...
        :return:
        """

        return not self.has_simulated_datacube_sed

    # -----------------------------------------------------------------

    @property
    def do_write_image_fluxes(self):

        """
        This function ...
        :return:
        """

        return not self.has_image_fluxes

    # -----------------------------------------------------------------

    @property
    def do_write_proper_image_fluxes(self):

        """
        This function ...
        :return:
        """

        return not self.has_proper_image_fluxes

    # -----------------------------------------------------------------

    @property
    def do_write_images_sed(self):

        """
        This function ...
        :return:
        """

        return not self.has_images_sed

    # -----------------------------------------------------------------

    @property
    def do_write_images_differences(self):

        """
        This function ...
        :return:
        """

        return not self.has_images_differences

    # -----------------------------------------------------------------

    @property
    def do_write_fluxes_differences(self):

        """
        This function ...
        :return:
        """

        return not self.has_fluxes_differences

    # -----------------------------------------------------------------

    @property
    def do_write_sed_differences(self):

        """
        This fucntion ...
        :return:
        """

        return not self.has_sed_differences

    # -----------------------------------------------------------------

    @property
    def do_write_images(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_write_residuals(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_write_proper_residuals(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_write_weighed_residuals(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_write_proper_weighed_residuals(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_write_residuals_distributions(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_write_proper_residuals_distributions(self):

        """
        Thisf unction ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_write_weighed_distributions(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_write_proper_weighed_distributions(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the reference SED
        if self.do_write_reference_sed: self.write_reference_sed()

        # Write the mock fluxes
        if self.do_write_fluxes: self.write_fluxes()

        # 2. Write the flux differences
        if self.do_write_differences: self.write_differences()

        # Write the simulated datacube SED
        if self.do_write_simulated_datacube_sed: self.write_simulated_datacube_sed()

        # 3. Write the image fluxes
        if self.do_write_image_fluxes: self.write_image_fluxes()

        # Write proper image fluxes
        if self.do_write_proper_image_fluxes: self.write_proper_image_fluxes()

        # 4. Write the image SED
        if self.do_write_images_sed: self.write_image_sed()

        # 5. Write the image flux differences
        if self.do_write_images_differences: self.write_image_differences()

        # 6. Write fluxes differences
        if self.do_write_fluxes_differences: self.write_fluxes_differences()

        # 7. Write SED differences
        if self.do_write_sed_differences: self.write_sed_differences()

        # 9. Write the residual frames
        if self.do_write_residuals: self.write_residuals()

        # Write the proper residuals frames
        if self.do_write_proper_residuals: self.write_proper_residuals()

        # 10. Write the weighed residual frames
        if self.do_write_weighed_residuals: self.write_weighed_residuals()

        # Write the proper weighed residual frames
        if self.do_write_proper_weighed_residuals: self.write_proper_weighed_residuals()

        # 11. Write the residual distributions
        if self.do_write_residuals_distributions: self.write_residuals_distributions()

        # Write the proper residual distributions
        if self.do_write_proper_residuals_distributions: self.write_proper_residuals_distributions()

        # 12. Write the weighed residual distributions
        if self.do_write_weighed_distributions: self.write_weighed_distributions()

        # Write the proper weighed residual distributions
        if self.do_write_proper_weighed_distributions: self.write_proper_weighed_distributions()

    # -----------------------------------------------------------------

    @property
    def reference_sed_path(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.join(self.evaluation_path, "sed.dat")

    # -----------------------------------------------------------------

    @property
    def has_reference_sed(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.reference_sed_path)

    # -----------------------------------------------------------------

    def remove_reference_sed(self):

        """
        This function ...
        :return:
        """

        fs.remove_file(self.reference_sed_path)

    # -----------------------------------------------------------------

    def write_reference_sed(self):

        """
        This function ...
        :return:
        """

        self.reference_sed.saveto(self.reference_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def fluxes_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.evaluation_path, "fluxes")

    # -----------------------------------------------------------------

    @property
    def fluxes_filepath(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.fluxes_path, "fluxes.dat")

    # -----------------------------------------------------------------

    @property
    def has_fluxes(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.fluxes_filepath)

    # -----------------------------------------------------------------

    def write_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the mock fluxes ...")

        # Save
        self.fluxes.saveto(self.fluxes_filepath)

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

    def remove_image_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        fs.remove_file(self.get_image_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    @lazyproperty
    def proper_images_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.evaluation_path, "proper_images")

    # -----------------------------------------------------------------

    def get_proper_image_filepath_for_filter(self, fltr):

        """
        Thisf unction ....
        :param fltr:
        :return:
        """

        return fs.join(self.proper_images_path, str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_proper_image_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_proper_image_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def remove_proper_image_for_filter(self, fltr):

        """
        Thisf unction ..
        :param fltr:
        :return:
        """

        fs.remove_file(self.get_proper_image_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    # def get_errors_filepath_for_filter(self, fltr):
    #
    #     """
    #     This function ...
    #     :param fltr:
    #     :return:
    #     """
    #
    #     return fs.join(self.images_path, str(fltr) + "_errors.fits")
    #
    # # -----------------------------------------------------------------
    #
    # def has_errors_for_filter(self, fltr):
    #
    #     """
    #     This function ...
    #     :param fltr:
    #     :return:
    #     """
    #
    #     return fs.is_file(self.get_errors_filepath_for_filter(fltr))
    #
    # # -----------------------------------------------------------------
    #
    # def remove_errors_for_filter(self, fltr):
    #
    #     """
    #     This function ...
    #     :param fltr:
    #     :return:
    #     """
    #
    #     fs.remove_file(self.get_errors_filepath_for_filter(fltr))

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
    def proper_residuals_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.evaluation_path, "proper_residuals")

    # -----------------------------------------------------------------

    @lazyproperty
    def weighed_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.evaluation_path, "weighed")

    # -----------------------------------------------------------------

    @lazyproperty
    def proper_weighed_path(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.evaluation_path, "proper_weighed")

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

    def get_proper_residuals_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.proper_residuals_path, str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_proper_residuals_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_proper_residuals_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def write_proper_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the proper residuals ...")

        # Loop over the proper residual maps
        for residuals in self.proper_residuals:

            # If already present
            if self.has_proper_residuals_for_filter(residuals.filter): continue

            # Determine the path
            path = self.get_proper_residuals_filepath_for_filter(residuals.filter)

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

    def write_weighed_residuals(self):

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

    def get_proper_weighed_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.proper_weighed_path, str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_proper_weighed_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_proper_weighed_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def write_proper_weighed_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the proper weighed residuals ...")

        # Loop over the proper weighed residuals maps
        for residuals in self.proper_weighed:

            # If already present
            if self.has_proper_weighed_for_filter(residuals.filter): continue

            # Determine the path
            path = self.get_proper_weighed_filepath_for_filter(residuals.filter)

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
        for fltr in self.simulated_flux_filters:

            # Check if already present
            if self.has_residuals_distribution_for_filter(fltr): continue

            # Get the distribution
            distribution = self.residuals_distributions[fltr]

            # Determine the path
            path = self.get_residuals_distribution_filepath_for_filter(fltr)

            # Save
            distribution.saveto(path)

    # -----------------------------------------------------------------

    def get_proper_residuals_distribution_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.proper_residuals_path, str(fltr) + "_distribution.dat")

    # -----------------------------------------------------------------

    def has_proper_residuals_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_proper_residuals_distribution_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def write_proper_residuals_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the proper residuals distributions ...")

        # Loop over the distributions
        for fltr in self.simulated_flux_filters:

            # Check if already present
            if self.has_proper_residuals_distribution_for_filter(fltr): continue

            # Get the distribution
            distribution = self.proper_residuals_distributions[fltr]

            # Determine the path
            path = self.get_proper_residuals_distribution_filepath_for_filter(fltr)

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
        for fltr in self.simulated_flux_filters:

            # Check if already present
            if self.has_weighed_distribution_for_filter(fltr): continue

            # Get the distribution
            distribution = self.weighed_distributions[fltr]

            # Determine the path
            path = self.get_weighed_distribution_filepath_for_filter(fltr)

            # Save
            distribution.saveto(path)

    # -----------------------------------------------------------------

    def get_proper_weighed_distribution_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.proper_weighed_path, str(fltr) + "_distribution.dat")

    # -----------------------------------------------------------------

    def has_proper_weighed_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_proper_weighed_distribution_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def write_proper_weighed_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the proper weighed distributions ...")

        # Loop over the distributions
        for fltr in self.simulated_flux_filters:

            # Check if already present
            if self.has_proper_weighed_distribution_for_filter(fltr): continue

            # Get the distribution
            distribution = self.proper_weighed_distributions[fltr]

            # Determine the path
            path = self.get_proper_weighed_distribution_filepath_for_filter(fltr)

            # Save
            distribution.saveto(path)

    # -----------------------------------------------------------------

    @property
    def do_plot_simulated_sed(self):

        """
        This function ...
        :return:
        """

        return not self.has_simulated_sed_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_simulated_datacube_sed(self):

        """
        This function ...
        :return:
        """

        return not self.has_simulated_datacube_sed_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_simulated_fluxes(self):

        """
        This function ...
        :return:
        """

        return not self.has_simulated_fluxes_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_differences(self):

        """
        This function ...
        :return:
        """

        return not self.has_differences_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_images(self):

        """
        Thisf unction ...
        :return:
        """

        return not self.has_all_image_plots

    # -----------------------------------------------------------------

    @property
    def do_plot_proper_images(self):

        """
        This function ...
        :return:
        """

        return not self.has_all_proper_image_plots

    # -----------------------------------------------------------------

    @property
    def do_plot_relative_errors(self):

        """
        This function ...
        :return:
        """

        return not self.has_all_relative_error_plots

    # -----------------------------------------------------------------

    @property
    def do_plot_fluxes(self):

        """
        Thisf unction ...
        :return:
        """

        return not self.has_fluxes_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_fluxes_only_images(self):

        """
        This function ...
        :return:
        """

        return not self.has_fluxes_only_images_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_fluxes_proper(self):

        """
        This function ...
        :return:
        """

        return not self.has_fluxes_proper_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_fluxes_differences(self):

        """
        This function ...
        :return:
        """

        return not self.has_fluxes_differences_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_seds(self):

        """
        This function ...
        :return:
        """

        return not self.has_sed_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_seds_only_images(self):

        """
        Thisnfunction ...
        :return:
        """

        return not self.has_sed_only_images_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_sed_differences(self):

        """
        This function ...
        :return:
        """

        return not self.has_sed_differences_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_seds_images(self):

        """
        This function ...
        :return:
        """

        return not self.has_seds_images_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_image_differences(self):

        """
        Thisf unction ...
        :return:
        """

        return not self.has_images_differences_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_residuals_distributions(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_plot_proper_residuals_distributions(self):

        """
        Thisf unction ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_plot_weighed_distributions(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_plot_proper_weighed_distributions(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_plot_residuals(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_plot_proper_residuals(self):

        """
        Thisnf unction ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_plot_residuals_absolute(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_plot_proper_residuals_absolute(self):

        """
        This function ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_plot_weighed_residuals(self):

        """
        Thisf unction ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    @property
    def do_plot_proper_weighed_residuals(self):

        """
        Thisf unction ...
        :return:
        """

        return True

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # 1. Plot the simulated SED together with the observed SED
        if self.do_plot_simulated_sed: self.plot_simulated_sed()

        # Plot the simulated datacube SED together with the observed SED
        if self.do_plot_simulated_datacube_sed: self.plot_simulated_datacube_sed()

        # 2. Plot the simulated fluxes together with the observed SED
        if self.do_plot_simulated_fluxes: self.plot_simulated_fluxes()

        # Plot differences between observed and simulated fluxes
        if self.do_plot_differences: self.plot_differences()

        # Plot the images
        if self.do_plot_images: self.plot_images()

        # Plot the proper images
        if self.do_plot_proper_images: self.plot_proper_images()

        # Plot the relative error maps
        #if self.do_plot_relative_errors: self.plot_relative_errors()

        # Plot the comparison between the simulated fluxes (from SED and from images)
        if self.do_plot_fluxes: self.plot_fluxes()

        # Plot the SED of the simulated fluxes from images, and the simulated SED
        if self.do_plot_fluxes_only_images: self.plot_fluxes_only_images()

        # Proper
        if self.do_plot_fluxes_proper: self.plot_fluxes_proper()

        # Plot differences between simulated fluxes (SED and from images)
        if self.do_plot_fluxes_differences: self.plot_fluxes_differences()

        # Plot the comparison between the observed fluxes (from SED and from images)
        if self.do_plot_seds: self.plot_seds()

        # Plot the SED of the observed fluxes derived from images
        if self.do_plot_seds_only_images: self.plot_seds_only_images()

        # Plot differences between observed fluxes (SED and from images)
        if self.do_plot_sed_differences: self.plot_sed_differences()

        # Plot the comparison between simulated fluxes and observed fluxes, both from IMAGES, and the simulated SED
        if self.do_plot_seds_images: self.plot_seds_images()

        # Plot differences between simulated fluxes and observed fluxes, both from IMAGES
        if self.do_plot_image_differences: self.plot_image_differences()

        # Plot residual maps
        if self.do_plot_residuals: self.plot_residuals()

        # Plot proper residual maps
        if self.do_plot_proper_residuals: self.plot_proper_residuals()

        # Plot absolute residual maps
        if self.do_plot_residuals_absolute: self.plot_residuals_absolute()

        # Plot proper absolute residual maps
        if self.do_plot_proper_residuals_absolute: self.plot_proper_residuals_absolute()

        # Plot weighed residual maps
        if self.do_plot_weighed_residuals: self.plot_weighed_residuals()

        # Plot proper weighed residual maps
        if self.do_plot_proper_weighed_residuals: self.plot_proper_weighed_residuals()

        # Plot residual distributions
        if self.do_plot_residuals_distributions: self.plot_residuals_distributions()

        # Plot proper residual distributions
        if self.do_plot_proper_residuals_distributions: self.plot_proper_residuals_distributions()

        # Plot weighed distributions
        if self.do_plot_weighed_distributions: self.plot_weighed_distributions()

        # Plot proper weighed distributions
        if self.do_plot_proper_weighed_distributions: self.plot_proper_weighed_distributions()

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
        #plotter.add_sed(self.observed_sed, "Observation (SED)")
        plotter.add_sed(self.reference_sed, "Observation")
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
        #plotter.add_sed(self.observed_sed, "Observation (SED)")
        plotter.add_sed(self.reference_sed, "Observation")
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
        #plotter.add_sed(self.observed_sed, "Observation (SED)") # ObservedSED
        plotter.add_sed(self.reference_sed, "Observation (SED)")
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

        for fltr in self.simulated_flux_filters:
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
        for fltr in self.simulated_flux_filters:

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

    def get_proper_image_plot_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.proper_images_path, str(fltr) + ".png")

    # -----------------------------------------------------------------

    def has_proper_image_plot_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_proper_image_plot_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    @property
    def has_all_proper_image_plots(self):

        """
        Thisf unction ...
        :return:
        """

        for fltr in self.simulated_flux_filters:
            if not self.has_proper_image_plot_for_filter(fltr): return False
        return True

    # -----------------------------------------------------------------

    def plot_proper_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the proper mock observed images ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Check whether present
            if self.has_proper_image_plot_for_filter(fltr): continue

            # Get the image
            frame = self.proper_images[fltr]

            # Determine the path
            path = self.get_proper_image_plot_filepath_for_filter(fltr)

            # Make plot
            vmin, vmax = frame.saveto_png(path, colours=self.config.colours, interval=self.config.interval,
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

        for fltr in self.simulated_flux_filters:
            #print(fltr, self.get_relative_errors_plot_filepath_for_filter(fltr), "has relative errors plot:", self.has_relative_errors_plot_for_filter(fltr))
            if not self.has_relative_errors_plot_for_filter(fltr): return False
        return True

    # -----------------------------------------------------------------

    # def plot_relative_errors(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Plotting the relative error maps of the mock observed images ...")
    #
    #     # Loop over the filters
    #     for fltr in self.simulated_flux_filters:
    #
    #         # Check whether present
    #         if self.has_relative_errors_plot_for_filter(fltr): continue
    #
    #         # Get the image and error map
    #         frame = self.images[fltr]
    #         errors = self.errors[fltr]
    #         relerrors = errors / frame
    #
    #         # Determine the path
    #         path = self.get_relative_errors_plot_filepath_for_filter(fltr)
    #
    #         # Plot
    #         plotting.plot_box(relerrors, path=path, colorbar=True, scale="linear")

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
        #plotter.add_sed(self.observed_sed, "Observation (SED)")
        plotter.add_sed(self.reference_sed, "Observation (SED)")
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
        for fltr in self.simulated_flux_filters:

            # Check if already present
            if self.has_residuals_distribution_plot_for_filter(fltr): continue

            # Get the distribution
            distribution = self.residuals_distributions[fltr]

            # Determine the plot path
            path = self.get_residuals_distribution_plot_filepath_for_filter(fltr)

            # Plot
            plot_distribution(distribution, path=path, x_limits=self.distribution_x_limits)

    # -----------------------------------------------------------------

    def get_proper_residuals_distribution_plot_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.proper_residuals_path, str(fltr) + "_distribution.pdf")

    # -----------------------------------------------------------------

    def has_proper_residuals_distribution_plot_for_filter(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_proper_residuals_distribution_plot_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def plot_proper_residuals_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the proper residuals distributions ...")

        # Loop over the distributions
        for fltr in self.simulated_flux_filters:

            # Check if already present
            if self.has_proper_residuals_distribution_plot_for_filter(fltr): continue

            # Get the distribution
            distribution = self.proper_residuals_distributions[fltr]

            # Determine the plot path
            path = self.get_proper_residuals_distribution_plot_filepath_for_filter(fltr)

            # Plot
            plot_distribution(distribution, path=path, x_limits=self.distribution_x_limits)

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
        for fltr in self.simulated_flux_filters:

            # Check if already present
            if self.has_weighed_distribution_plot_for_filter(fltr): continue

            # Get the distribution
            distribution = self.weighed_distributions[fltr]

            # Determine the plot path
            path = self.get_weighed_distribution_plot_filepath_for_filter(fltr)

            # Plot
            plot_distribution(distribution, path=path, x_limits=self.weighed_distribution_x_limits)

    # -----------------------------------------------------------------

    def get_proper_weighed_distribution_plot_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.proper_weighed_path, str(fltr) + "_distribution.pdf")

    # -----------------------------------------------------------------

    def has_proper_weighed_distribution_plot_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_proper_weighed_distribution_plot_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def plot_proper_weighed_distributions(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Plotting the proper weighed residuals distributions ...")

        # Loop over the weighed distributions
        for fltr in self.simulated_flux_filters:

            # Check if already present
            if self.has_proper_weighed_distribution_plot_for_filter(fltr): continue

            # Get the distribution
            distribution = self.proper_weighed_distributions[fltr]

            # Determine the plot path
            path = self.get_proper_weighed_distribution_plot_filepath_for_filter(fltr)

            # Plot
            plot_distribution(distribution, path=path, x_limits=self.weighed_distribution_x_limits)

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

        return (-100., 100.)

    # -----------------------------------------------------------------

    @property
    def residuals_plot_cmap(self):

        """
        This function ...
        :return:
        """

        return "RdBu_r"

    # -----------------------------------------------------------------

    def plot_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the residuals ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Check whether present
            if self.has_residuals_plot_for_filter(fltr): continue

            # Debugging
            log.debug("Plotting the residuals for the '" + str(fltr) + "' filter ...")

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
            plotting.plot_box(residuals_percentage, interval=self.residuals_plot_interval, path=path, colorbar=True,
                              around_zero=True, scale="linear", cmap=self.residuals_plot_cmap, check_around_zero=False)

    # -----------------------------------------------------------------

    def get_proper_residuals_plot_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.proper_residuals_path, str(fltr) + ".png")

    # -----------------------------------------------------------------

    def has_proper_residuals_plot_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_proper_residuals_plot_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def plot_proper_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the proper residuals ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Check whether present
            if self.has_proper_residuals_plot_for_filter(fltr): continue

            # Debugging
            log.debug("Plotting the proper residuals for the '" + str(fltr) + "' filter ...")

            # Get the residual map
            residuals = self.proper_residuals[fltr]

            # Get residuals in percentage
            residuals_percentage = residuals * 100

            # Determine the path
            path = self.get_proper_residuals_plot_filepath_for_filter(fltr)

            # Plot
            plotting.plot_box(residuals_percentage, interval=self.residuals_plot_interval, path=path, colorbar=True,
                              around_zero=True, scale="linear", cmap=self.residuals_plot_cmap, check_around_zero=False)

    # -----------------------------------------------------------------

    def get_residuals_absolute_plot_filepath_for_filter(self, fltr):

        """
        Thisn fucntion ...
        :param fltr:
        :return:
        """

        return fs.join(self.residuals_path, str(fltr) + "_absolute.png")

    # -----------------------------------------------------------------

    def has_residuals_absolute_plot_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_residuals_absolute_plot_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    @property
    def residuals_absolute_plot_interval(self):

        """
        This function ...
        :return:
        """

        return (0, 100.)

    # -----------------------------------------------------------------

    @property
    def residuals_absolute_plot_cmap(self):

        """
        This function ...
        :return:
        """

        return "viridis"

    # -----------------------------------------------------------------

    def plot_residuals_absolute(self):

        """
        Thins function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the residuals in absolute values ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Check whether present
            if self.has_residuals_absolute_plot_for_filter(fltr): continue

            # Debugging
            log.debug("Plotting the residuals in absolute values for the '" + str(fltr) + "' filter ...")

            # Get the absolute residual map
            residuals = self.residuals[fltr].absolute

            # Get residuals in percentage
            residuals_percentage = residuals * 100

            # Determine the path
            path = self.get_residuals_absolute_plot_filepath_for_filter(fltr)

            # Plot
            plotting.plot_box(residuals_percentage, interval=self.residuals_absolute_plot_interval, path=path, colorbar=True, scale="linear", cmap=self.residuals_absolute_plot_cmap)

    # -----------------------------------------------------------------

    def get_proper_residuals_absolute_plot_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.proper_residuals_path, str(fltr) + "_absolute.png")

    # -----------------------------------------------------------------

    def has_proper_residuals_absolute_plot_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_proper_residuals_absolute_plot_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def plot_proper_residuals_absolute(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the proper residuals in absolute values ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Check whether present
            if self.has_proper_residuals_absolute_plot_for_filter(fltr): continue

            # Debugging
            log.debug("Plotting the proper residuals in absolute values for the '" + str(fltr) + "' filter ...")

            # Get the absolute residual map
            residuals = self.proper_residuals[fltr].absolute

            # Get residuals in percentage
            residuals_percentage = residuals * 100

            # Determine the path
            path = self.get_proper_residuals_absolute_plot_filepath_for_filter(fltr)

            # Plot
            plotting.plot_box(residuals_percentage, interval=self.residuals_absolute_plot_interval, path=path,
                              colorbar=True, scale="linear", cmap=self.residuals_absolute_plot_cmap)

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
        for fltr in self.simulated_flux_filters:

            # Ceck whether present
            if self.has_weighed_plot_for_filter(fltr): continue

            # Debugging
            log.debug("Plotting the weighed residuals for the '" + str(fltr) + "' filter ...")

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
            plotting.plot_box(residuals, path=path, colorbar=True, around_zero=True, scale="linear", check_around_zero=False)

    # -----------------------------------------------------------------

    def get_proper_weighed_plot_filepath_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.proper_weighed_path, str(fltr) + ".png")

    # -----------------------------------------------------------------

    def has_proper_weighed_plot_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_proper_weighed_plot_filepath_for_filter(fltr))

    # -----------------------------------------------------------------

    def plot_proper_weighed_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the proper weighed residuals ...")

        # Loop over the filters
        for fltr in self.simulated_flux_filters:

            # Ceck whether present
            if self.has_proper_weighed_plot_for_filter(fltr): continue

            # Debugging
            log.debug("Plotting the proper weighed residuals for the '" + str(fltr) + "' filter ...")

            # Get the weighed residual map
            residuals = self.proper_weighed[fltr]

            # Determine the path
            path = self.get_proper_weighed_plot_filepath_for_filter(fltr)

            # Plot
            plotting.plot_box(residuals, path=path, colorbar=True, around_zero=True, scale="linear", check_around_zero=False)

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
