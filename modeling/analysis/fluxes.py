#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.fluxes Contains the FluxesAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .component import AnalysisRunComponent
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.tools.utils import lazyproperty
from ...core.tools import sequences
from ...core.data.sed import ObservedSED
from ...magic.core.dataset import StaticDataSet
from ...core.tools.stringify import tostr
from ...core.plot.sed import plot_seds
from ...core.filter.filter import parse_filter_from_instrument_and_band
from ...core.tools import tables
from ..fitting.modelanalyser import FluxDifferencesTable
from ...magic.core.mask import Mask

# -----------------------------------------------------------------

class FluxesAnalyser(AnalysisRunComponent):

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
        super(AnalysisRunComponent, self).__init__(*args, **kwargs)

        # The simulated images
        self.simulated = dict()

        # The observed images
        self.observed = dict()

        # The masks
        self.masks = dict()

        # The SED
        self.sed = None

        # The flux differences table
        self.differences = None

    # -----------------------------------------------------------------

    @property
    def needs_images(self):
        return not self.has_sed

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get images
        if self.needs_images: self.get_images()

        # Calculate
        self.get_fluxes()

        # Calculate differences
        self.calculate_differences()

        # Writing
        self.write()

        # Plotting
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisRunComponent, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @property
    def fluxes_path(self):
        return self.analysis_run.fluxes_path

    # -----------------------------------------------------------------
    # OBSERVED IMAGES
    # -----------------------------------------------------------------

    @property
    def observed_images_dataset(self):
        return self.static_photometry_dataset

    # -----------------------------------------------------------------
    # SIMULATED IMAGES
    # -----------------------------------------------------------------

    @property
    def simulated_images_path(self):
        return self.analysis_run.images_path

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_images_dataset(self):
        return StaticDataSet.from_directory(self.simulated_images_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):
        return sequences.sorted_by_attribute(self.simulated_images_dataset.filters, "wavelength")

    # -----------------------------------------------------------------

    @lazyproperty
    def filter_names(self):
        return [str(fltr) for fltr in self.filters]

    # -----------------------------------------------------------------

    def get_images(self):

        """
        This function ...
        :return:
        """

        # Load the images
        self.load_images()

        # Mask the images
        self.mask_images()

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

        # Load masks
        self.load_masks()

    # -----------------------------------------------------------------

    def load_observed_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the observed images ...")

        # Get the frames
        self.observed = self.observed_images_dataset.get_frames_for_filters(self.filters, as_dict=True)

    # -----------------------------------------------------------------

    def load_simulated_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated images ...")

        # Get the frames
        self.simulated = self.simulated_images_dataset.get_frames_for_filters(self.filters, as_dict=True)

    # -----------------------------------------------------------------

    def load_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the image masks ...")

        # Get the frames
        for fltr in self.filters:

            # Get mask path
            path = fs.join(self.analysis_run.residuals_path, "masks", str(fltr) + ".fits")

            # Load
            self.masks[fltr] = Mask.from_file(path)

    # -----------------------------------------------------------------

    def mask_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning and masking the simulated images ...")

        # Loop over the filters
        for fltr in self.filters:

            # Get observation and model image
            observation = self.observed[fltr]
            model = self.simulated[fltr]

            # Convert to surface brightness units
            observation.convert_to("mJy/sr", distance=self.galaxy_distance)
            model.convert_to("mJy/sr", distance=self.galaxy_distance)

            # Mask
            #model.rebin(observation.wcs)
            #model.apply_mask_nans(observation.nans)

            # Mask
            model.rebin(self.masks[fltr].wcs)
            model.apply_mask_nans(self.masks[fltr])

            # Convert back to unit
            observation.convert_to(self.config.unit, distance=self.galaxy_distance)
            model.convert_to(self.config.unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    def get_fluxes(self):

        """
        This function ...
        :return:
        """

        # Load from file
        if self.has_sed: self.load_fluxes()

        # Calculate
        else: self.calculate_fluxes()

    # -----------------------------------------------------------------

    def load_fluxes(self):

        """
        This function ...
        :return:
        """

        # Load
        self.sed = ObservedSED.from_file(self.sed_path)

    # -----------------------------------------------------------------

    def calculate_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating fluxes from the mock images ...")

        # Create an observed SED for the mock fluxes
        self.sed = ObservedSED(photometry_unit=self.config.unit)

        # Loop over the filters
        for fltr in self.filters:

            # Set filter name
            filter_name = str(fltr)

            # Debugging
            log.debug("Calculating flux for the '" + filter_name + "' filter ...")

            # Get the image frame
            frame = self.simulated[fltr]

            # Calculate proper flux
            flux = frame.sum(add_unit=True)

            # Add this entry to the SED
            self.sed.add_point(fltr, flux)

    # -----------------------------------------------------------------

    @lazyproperty
    def fit_sed(self):

        """
        This function ...
        :return:
        """

        # Load the clipped observed SED
        sed = self.observed_sed.copy()

        # Add additional relative error?
        if self.config.additional_error is not None: sed.add_relative_error(self.config.additional_error)

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def has_observed_flux(self, fltr):

        """
        Thisn function ...
        :param fltr:
        :return:
        """

        observed_fluxdensity = self.fit_sed.photometry_for_filter(fltr)
        return observed_fluxdensity is not None

    # -----------------------------------------------------------------

    def get_observed_flux(self, fltr, unit="Jy"):

        """
        This function ...
        :param fltr:
        :param unit:
        :return:
        """

        # Find the corresponding flux in the SED derived from observation
        observed_fluxdensity = self.fit_sed.photometry_for_filter(fltr, unit=unit).value

        # Find the corresponding flux error in the SED derived from observation
        observed_fluxdensity_error = self.fit_sed.error_for_filter(fltr, unit=unit).average.to(unit).value

        # Return the values
        return observed_fluxdensity, observed_fluxdensity_error

    # -----------------------------------------------------------------

    @property
    def has_weights(self):
        return self.analysis_run.from_fitting

    # -----------------------------------------------------------------

    @property
    def fitting_run(self):
        return self.analysis_run.fitting_run

    # -----------------------------------------------------------------

    @property
    def weights(self):
        return self.fitting_run.weights

    # -----------------------------------------------------------------

    def has_weight(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        if not self.has_weights: return False

        # Find the index of the current band in the weights table
        index = tables.find_index(self.weights, key=[fltr.instrument, fltr.band], column_name=["Instrument", "Band"])
        return index is not None

    # -----------------------------------------------------------------

    def get_weight(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Get index
        index = tables.find_index(self.weights, key=[fltr.instrument, fltr.band], column_name=["Instrument", "Band"])
        weight = self.weights["Weight"][index] # apparently, this is a string, so parsing the table went wrong ...
        return float(weight)

    # -----------------------------------------------------------------

    def calculate_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the flux differences with the observed fluxes ...")

        # Initialize the differences table
        self.differences = FluxDifferencesTable()

        # Loop over the entries in the fluxdensity table (SED) derived from the simulation
        for i in range(len(self.sed)):

            # Get instrument, band and flux density
            instrument = self.sed["Instrument"][i]
            band = self.sed["Band"][i]
            fluxdensity = self.sed["Photometry"][i]
            fltr = parse_filter_from_instrument_and_band(instrument, band)

            # No observed flux?
            if not self.has_observed_flux(fltr):
                log.warning("The observed flux density could not be found for the " + instrument + " '" + str(fltr) + "' filter")
                continue

            # Get fluxdensity and error
            observed_fluxdensity, observed_fluxdensity_error = self.get_observed_flux(fltr)

            # Debugging
            log.debug("Calculating relative difference and chi squared term for the '" + str(fltr) + "' filter ...")
            log.debug("Observed flux: " + tostr(observed_fluxdensity) + " ± " + tostr(observed_fluxdensity_error))

            # Calculate difference and relative difference
            difference = fluxdensity - observed_fluxdensity
            relative_difference = difference / observed_fluxdensity

            # No weight?
            if self.has_weight(fltr):

                # Get the weight
                weight = self.get_weight(fltr)

                # Show the weight
                log.debug("The weight for the '" + str(fltr) + "' filter is " + tostr(weight))

                # Calculate the chi squared term
                chi_squared_term = weight * difference ** 2 / observed_fluxdensity_error ** 2

            else:

                # Give warnings
                log.warning("A weight is not found for the '" + str(fltr) + "' filter: skipping ...")

                # Set chi squared term to None
                chi_squared_term = None

            # Add an entry to the differences table
            self.differences.add_entry(instrument, band, difference, relative_difference, chi_squared_term)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the fluxes
        self.write_fluxes()

        # Write the differences table
        self.write_differences()

    # -----------------------------------------------------------------

    @property
    def sed_path(self):
        return fs.join(self.fluxes_path, "fluxes.dat")

    # -----------------------------------------------------------------

    @property
    def has_sed(self):
        return fs.is_file(self.sed_path)

    # -----------------------------------------------------------------

    def write_fluxes(self):

        """
        This function ....
        :return:
        """

        # Inform the user
        log.info("Writing the fluxes ...")

        # Write
        self.sed.saveto(self.sed_path)

    # -----------------------------------------------------------------

    @property
    def differences_path(self):
        return fs.join(self.fluxes_path, "differences.dat")

    # -----------------------------------------------------------------

    @property
    def has_differences(self):
        return fs.is_file(self.differences_path)

    # -----------------------------------------------------------------

    def write_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the flux differences ...")

        # Write
        self.differences.saveto(self.differences_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the observed and
        self.plot_fluxes()

    # -----------------------------------------------------------------

    @property
    def fluxes_plot_path(self):
        return fs.join(self.fluxes_path, "fluxes.pdf")

    # -----------------------------------------------------------------

    @property
    def has_fluxes_plot(self):
        return fs.is_file(self.fluxes_plot_path)

    # -----------------------------------------------------------------

    def plot_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the observed and mock fluxes ...")

        # Set SEDs
        seds = OrderedDict()
        seds["observed"] = self.fit_sed
        seds["mock"] = self.sed

        # Plot
        plot_seds(seds, path=self.fluxes_plot_path)

# -----------------------------------------------------------------
