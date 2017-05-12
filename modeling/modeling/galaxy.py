#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.modeler Contains the GalaxyModeler class.
#  Perform radiative transfer modeling for a certain galaxy

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ..data.properties import PropertyFetcher
from ..data.images import ImageFetcher
from ..data.seds import SEDFetcher
from ..preparation.initialization import PreparationInitializer
from ..preparation.preparer import DataPreparer
from ..decomposition.decomposition import GalaxyDecomposer
from ..truncation.truncation import Truncator
from ..photometry.photometry import PhotoMeter
from ..maps.old import OldStellarMapMaker
from ..maps.young import YoungStellarMapMaker
from ..maps.ionizing import IonizingStellarMapMaker
from ..maps.dust import DustMapMaker
from ..maps.colour import ColourMapMaker
from ..maps.attenuation import AttenuationMapMaker
from ..maps.tir import TIRMapMaker
from ..fitting.configuration import FittingConfigurer
from ..fitting.initialization.galaxy import GalaxyFittingInitializer
from ...core.basics.range import QuantityRange
from .base import ModelerBase
from ..config.parameters import units as parameter_units
from ..config.parameters import default_ranges, types, parameter_descriptions
from ...core.units.parsing import parse_unit as u
from ..build.model import ModelBuilder
from ..build.representation import RepresentationBuilder
from ..component.galaxy import get_galaxy_properties_path, get_data_seds_path, get_data_images_path, get_dustpedia_sed
from ...core.tools import filesystem as fs
from ...core.filter.filter import parse_filter
from ...magic.core.image import Image
from ..core.environment import GalaxyModelingEnvironment
from ...core.remote.utils import DetachedCalculation
from ...core.tools.utils import UserIntervention
from ..maps.ssfr import SSFRMapMaker
from ...core.tools import types
from ..maps.significance import SignificanceMaskCreator
from ..preparation.inspector import PreparationInspector
from ..component.galaxy import get_observed_sed_file_path, get_reference_seds
from ...core.plot.sed import SEDPlotter
from ...core.tools import parsing
from ...core.tools.logging import set_log_file, unset_log_file

# -----------------------------------------------------------------

# Define the different modeling methods
modeling_methods = dict()
modeling_methods["DL14"] = "High-resolution radiative transfer modeling of face-on galaxies based on De Looze et al. 2014"
modeling_methods["FitSKIRT"] = "FitSKIRT method for fitting edge-on galaxies"
modeling_methods["MGE"] = "Panchromatic radiative transfer modeling using Multi-Gaussian expansion"

# -----------------------------------------------------------------

# Set the free parameters for different modeling methods
free_parameters = dict()
free_parameters["DL14"] = ["fuv_young", "dust_mass", "fuv_ionizing"]

# -----------------------------------------------------------------

free_parameter_ranges = dict()
free_parameter_ranges["DL14"] = {"fuv_young": QuantityRange(0.0, 1e37, unit="W/micron"),
                                 "dust_mass": QuantityRange(0.5e7, 3.e7, unit="Msun"),
                                 "fuv_ionizing": QuantityRange(0.0, 1e34, unit="W/micron")}

# -----------------------------------------------------------------

parameter_grid_scales = dict()
parameter_grid_scales["DL14"] = {"fuv_young": "logarithmic",
                                 "dust_mass": "logarithmic",
                                 "fuv_ionizing": "logarithmic"}

# -----------------------------------------------------------------

parameter_grid_weights = dict()
weight_young, weight_mass, weight_ionizing = parsing.weights("7,10,5")
parameter_grid_weights["DL14"] = {"fuv_young": weight_young,
                                  "dust_mass": weight_mass,
                                  "fuv_ionizing": weight_ionizing}

# -----------------------------------------------------------------

# Set the filters for which the data can be used for fitting (data that can be trusted well enough)
fitting_filter_names = ["GALEX FUV", "GALEX NUV", "SDSS u", "SDSS g", "SDSS r", "SDSS i", "SDSS z", "WISE W1",
                        "IRAC I1", "IRAC I2", "WISE W2", "IRAC I3", "IRAC I4", "WISE W3", "WISE W4", "MIPS 24mu",
                        "Pacs blue", "Pacs red", "SPIRE PSW", "SPIRE PMW", "SPIRE PLW"]

# -----------------------------------------------------------------

# URLS for Halpha data for different galaxies (keys are the HYPERLEDA names)
halpha_urls = {"NGC3031": "https://ned.ipac.caltech.edu/img/2001ApJ...559..878H/MESSIER_081:I:Ha:hwb2001.fits.gz"}
halpha_fluxes = {"NGC3031": 7.8e40 * u("erg/s")}

# -----------------------------------------------------------------

# Set the filters that are not necessary for making maps and therefore shouldn't be included in the preparation steps
# that brings all data to the same resolution and pixelscale
# We want to exclude the SPIRE images from the procedures that bring all images to the same resolution
lower_resolution_filters = ["SPIRE PSW", "SPIRE PMW", "SPIRE PLW", "MIPS 70mu", "MIPS 160mu", "WISE W4"] # WISE W4 because it has just a slightly higher FWHM than Pacs 160mu

# -----------------------------------------------------------------

class GalaxyModeler(ModelerBase):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(GalaxyModeler, self).__init__(*args, **kwargs)

        # Attributes
        self.properties = None
        self.seds = None
        self.images = None

        # Models
        self.disk_model = None
        self.bulge_model = None

    # -----------------------------------------------------------------

    @property
    def galaxy_name(self):

        """
        This function ...
        :return: 
        """

        return self.modeling_config.name

    # -----------------------------------------------------------------

    @property
    def ngc_name(self):

        """
        This function ...
        :return:
        """

        return self.modeling_config.ngc_name

    # -----------------------------------------------------------------

    @property
    def hyperleda_name(self):

        """
        This function ...
        :return:
        """

        return self.modeling_config.hyperleda_name

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Get the data
        self.get_data()

        # 3. Data preparation
        self.prepare_data()

        # 4. Decomposition
        if "decompose" not in self.history: self.decompose()

        # 5. Truncation
        if "truncate" not in self.history: self.truncate()

        # 6. Do the photometry
        if "photometry" not in self.history: self.photometry()

        # 7. Make the maps
        self.make_maps()

        # 8. Build model and its representations
        self.build()

        # 9. Do the fitting
        self.fit(scales=parameter_grid_scales, sampling_weights=parameter_grid_weights)

        # 10. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(GalaxyModeler, self).setup(**kwargs)

        # Load the modeling environment
        self.environment = GalaxyModelingEnvironment(self.modeling_path)

        # Get arguments
        self.properties = kwargs.pop("properties", None)
        self.seds = kwargs.pop("seds", None)
        self.images = kwargs.pop("images", None)

        # Bulge and disk model
        self.bulge_model = kwargs.pop("bulge_model", None)
        self.disk_model = kwargs.pop("disk_model", None)

        # Check whether a remote is available for the heavy computations, if one was configured
        if not self.other_local and self.moderator.host_id_for_single("other") is None: raise RuntimeError("The desired remote(s) for heavy computations are currently unavailable")

    # -----------------------------------------------------------------

    def get_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the galaxy data ...")

        # Get the galaxy properties
        if "fetch_properties" not in self.history:
            if self.properties is not None: self.set_properties()
            else: self.get_properties()

        # Get the galaxy SEDs
        if "fetch_seds" not in self.history:
            if self.seds is not None: self.set_seds()
            else: self.get_seds()

        # Get the galaxy images
        if "fetch_images" not in self.history:
            if self.images is not None: self.set_images()
            else: self.get_images()

    # -----------------------------------------------------------------

    def get_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the general galaxy properties ...")

        # Create the fetcher
        fetcher = PropertyFetcher()

        # Add an entry to the history and save
        self.history.add_entry(PropertyFetcher.command_name())

        # Set the working directory
        fetcher.config.path = self.modeling_path

        # Run the fetcher
        fetcher.run()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def set_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the provided galaxy properties ...")

        # Add an entry to the history and save
        self.history.add_entry(PropertyFetcher.command_name())

        # Write the properties to the correct path
        self.properties.saveto(get_galaxy_properties_path(self.modeling_path))

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def get_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the galaxy SEDs ...")

        # Create the SED fetcher
        fetcher = SEDFetcher()

        # Add an entry to the history
        self.history.add_entry(SEDFetcher.command_name())

        # Set the working directory
        fetcher.config.path = self.modeling_path

        # Run the SED fetcher
        fetcher.run()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def set_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the provided SEDs ...")

        # Loop over the different SEDs
        for label in self.seds:

            # Debugging info
            log.debug("Writing " + label + " SED ...")

            # Determine the path to the new SED file
            sed_path = fs.join(get_data_seds_path(self.modeling_path), label + ".dat")

            # Save the SED at the specified location
            self.seds[label].saveto(sed_path)

    # -----------------------------------------------------------------

    def get_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the galaxy images (this can take a while) ...")

        # Create the configuration
        config = dict()
        config["remote"] = self.moderator.host_id_for_single("other")
        config["attached"] = self.config.attached

        # Get Halpha URL and flux
        if self.hyperleda_name in halpha_urls: halpha_url = halpha_urls[self.hyperleda_name]
        else: halpha_url = None
        if self.hyperleda_name in halpha_fluxes: halpha_flux = halpha_fluxes[self.hyperleda_name]
        else: halpha_flux = None
        config["halpha_url"] = halpha_url
        config["halpha_flux"] = halpha_flux

        # Advanced
        config["max_nobservations_mosaic"] = self.config.max_nobservations_mosaic

        # Set flag to create Poisson error maps
        config["make_poisson"] = self.config.make_poisson

        # Create the image fetcher
        fetcher = ImageFetcher(config)

        # Add an entry to the history
        self.history.add_entry(ImageFetcher.command_name())

        # Set the working directory
        fetcher.config.path = self.modeling_path

        # Run the image fetcher
        fetcher.run()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

        # If not running in attached remote mode, exit now
        if fetcher.detached:
            log.warning("The procedure that calculates the Poisson error maps for GALEX and SDSS is now running. "
                        "Wait for it to finished and resume the modeling afterwards")
            raise DetachedCalculation(self.__class__, "get_data")

    # -----------------------------------------------------------------

    def set_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the provided galaxy images ...")

        # Loop over the images
        for filter_name in self.images:

            # Open the image if necessary
            if types.is_string_type(self.images[filter_name]): image = Image.from_file(self.images[filter_name])
            else: image = self.images[filter_name]

            # Determine the path
            fltr = parse_filter(filter_name)
            origin = fltr.observatory
            images_path = get_data_images_path(self.modeling_path)
            fs.create_directory(images_path)

            # Create directory for origin
            origin_path = fs.create_directory_in(images_path, origin)
            image_path = fs.join(origin_path, filter_name)

            # Write the image
            image.saveto(image_path)

    # -----------------------------------------------------------------

    def prepare_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the galaxy data ...")

        # Initialize the preparation
        if "initialize_preparation" not in self.history: self.initialize_preparation()

        # Run the preparation
        if "prepare_data" not in self.history: self.prepare()

        # Inspect the preparation
        if "inspect_preparation" not in self.history: self.inspect_preparation()

    # -----------------------------------------------------------------

    def initialize_preparation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the data preparation ...")

        # Create the configuration
        config = dict()
        config["remote"] = self.moderator.host_id_for_single("other")
        config["attached"] = self.config.attached
        config["sources"] = dict()
        config["sources"]["weak"] = self.config.sources_weak
        config["sources"]["nprocesses"] = self.config.nprocesses

        # Create the initializer
        initializer = PreparationInitializer(config)

        # Add an entry to the history
        command_name = PreparationInitializer.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        initializer.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the initializer
        initializer.run()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

        # Give warning and exit
        message = "Check the result of the source detection, make adjustments where necessary, and resume the modeling afterwards"
        raise UserIntervention(message, self.__class__, "initialize_preparation")

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the data preparation ...")

        # Create the configuration
        config = dict()
        config["remote"] = self.moderator.host_id_for_single("other")
        config["attached"] = self.config.attached
        config["exclude_filters"] = lower_resolution_filters
        config["nprocesses"] = self.config.nprocesses

        # Create the data preparer
        preparer = DataPreparer()

        # Add an entry to the history
        command_name = DataPreparer.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        preparer.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the preparer
        preparer.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def inspect_preparation(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Inspecting the preparation ...")

        # Create the configuration
        config = dict()

        # Create the prepration inspector
        inspector = PreparationInspector(config)

        # Add an entry to the history
        PreparationInspector.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        inspector.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the inspector
        inspector.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def decompose(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Running decomposition on the galaxy ...")

        # Create the decomposer
        decomposer = GalaxyDecomposer()

        # Add an entry to the history
        command_name = GalaxyDecomposer.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        decomposer.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the decomposer
        decomposer.run(disk=self.disk_model, bulge=self.bulge_model)

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def truncate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making truncation masks for the galaxy images ...")

        # Create the truncator
        truncator = Truncator()

        # Add an entry to the history
        command_name = Truncator.command_name()
        self.history.add_entry()

        # Set the working directory
        truncator.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the truncator
        truncator.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def photometry(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the photometry ...")

        # Perform the photometry
        if self.config.perform_photometry: self.perform_photometry()

        # Or just get the DustPedia photometry
        else: self.get_photometry()

    # -----------------------------------------------------------------

    def perform_photometry(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Performing photometry on the galaxy images ...")

        # Create the photometer
        photometer = PhotoMeter()

        # Add an entry to the history
        command_name = PhotoMeter.command_name()
        self.history.add_entry()

        # Set the working directory
        photometer.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the photometer
        photometer.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def get_photometry(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Getting the DustPedia photometry for the galaxy ...")

        # Add an entry to the history
        self.history.add_entry("photometry")

        # Get the DustPedia SED
        sed = get_dustpedia_sed(self.modeling_path)

        # Save
        sed.saveto(get_observed_sed_file_path(self.modeling_path))

        # Plot
        plotter = SEDPlotter()
        plotter.add_sed(sed, "DustPedia")
        path = fs.join(self.environment.phot_path, "sed.pdf")
        plotter.run(output=path, title=self.galaxy_name)

        # Plot with references
        plotter = SEDPlotter()
        seds = get_reference_seds(self.modeling_path)
        for label in seds:
            if label == "DustPedia": continue
            plotter.add_sed(seds[label], label)
        path = fs.join(self.environment.phot_path, "sed_with_references.pdf")
        plotter.run(ouput=path, title=self.galaxy_name)

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the maps describing the model geometries ...")

        # Create colour maps
        if "make_colour_maps" not in self.history: self.make_colour_maps()

        # Create sSFR maps
        if "make_ssfr_maps" not in self.history: self.make_ssfr_maps()

        # Create the TIR map
        if "make_tir_maps" not in self.history: self.make_tir_map()

        # Create the attenuation map(s)
        if "make_attenuation_maps" not in self.history: self.make_attenuation_maps()

        # Create the dust map
        if "make_dust_map" not in self.history: self.make_dust_maps()

        # Create the map of the old stellar disk
        if "make_old_stars_map" not in self.history: self.make_old_stellar_maps()

        # Create the map of the young stellar population
        if "make_young_stars_map" not in self.history: self.make_young_stellar_maps()

        # Create the map of the ionizing stellar population
        if "make_ionizing_stars_map" not in self.history: self.make_ionizing_stellar_maps()

        # Calculate the significance masks
        if "create_significance_masks" not in self.history: self.create_significance_masks()

    # -----------------------------------------------------------------

    def make_colour_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the colour maps ...")

        # Create the colour maps
        maker = ColourMapMaker()

        # Add an entry to the history
        command_name = maker.command_name()
        self.history.add_entry(command_name)

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Set the working directory
        maker.config.path = self.modeling_path

        # Run maker
        maker.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def make_ssfr_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the sSFR maps ...")

        # Create the sSFR maps
        maker = SSFRMapMaker()

        # Add an entry to the history
        command_name = maker.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the maker
        maker.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def make_tir_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the TIR map ...")

        # Create the TIR map maker
        maker = TIRMapMaker()

        # Add an entry to the history
        command_name = maker.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the maker
        maker.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def make_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the attenuation maps ...")

        # Create the attenuation map maker
        maker = AttenuationMapMaker()

        # Add an entry to the history
        self.history.add_entry(maker.command_name())

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the maker
        maker.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def make_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the dust maps ...")

        # Create the configuration
        config = dict()
        config["black_body"] = dict()
        config["black_body"]["remote"] = self.moderator.host_id_for_single("other")

        # Create the dust map maker
        maker = DustMapMaker(config)

        # Add an entry to the history
        command_name = DustMapMaker.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the dust map maker
        maker.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def make_old_stellar_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the maps of old stars ...")

        # Create the old stellar map maker
        maker = OldStellarMapMaker()

        # Add an entry to the history
        command_name = maker.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the maker
        maker.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def make_young_stellar_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the maps of young stars ...")

        # Create the young stellar map maker
        maker = YoungStellarMapMaker()

        # Add an entry to the history
        command_name = YoungStellarMapMaker.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the maker
        maker.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def make_ionizing_stellar_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the maps of ionizing stars ...")

        # Create the ionizing stellar map maker
        maker = IonizingStellarMapMaker()

        # Add an entry to the history
        command_name = IonizingStellarMapMaker.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the ionizing stellar map maker
        maker.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def create_significance_masks(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating significance masks ...")

        # Create the configuration
        config = dict()

        # Create the significance mask creator
        creator = SignificanceMaskCreator(config)

        # Add an entry to the history
        command_name = SignificanceMaskCreator.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        creator.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the creator
        creator.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def build(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the model and its representations ...")

        # Build model
        self.build_model()

        # Build representations
        self.build_representations()

    # -----------------------------------------------------------------

    def build_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the model ...")

        # Create configuration
        config = dict()

        # Set model name
        config["name"] = self.model_name

        # Create the builder
        builder = ModelBuilder(config)

        # Add an entry to the history
        command_name = ModelBuilder.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        builder.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the model builder
        builder.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def build_representations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the representations ...")

        # Create configuration
        config = dict()

        # Set name for representation
        config["name"] = self.representation_name

        # Set model name
        config["model_name"] = self.model_name

        # Create the builder
        builder = RepresentationBuilder(config)

        # Add an entry to the history
        command_name = RepresentationBuilder.command_name()
        self.history.add_entry(command_name)

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the builder
        builder.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def configure_fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the fitting ...")

        # Create configuration
        config = dict()

        # Set the name for the fitting run
        config["name"] = self.fitting_run_name

        # Set the model name
        config["model_name"] = self.model_name

        # Set the representation name
        config["representation_name"] = self.representation_name

        # Set free parameters
        config["parameters"] = free_parameters[self.modeling_config.method]

        # Set parameter descriptions
        config["descriptions"] = parameter_descriptions

        # Set parameter types
        config["types"].name = types

        # Set parameter units
        config["units"] = parameter_units

        # Set ranges
        config["ranges"] = free_parameter_ranges[self.modeling_config.method]

        # Set fitting filters
        config["filters"] = fitting_filter_names

        # Ask for and set the fitting method
        config["fitting_method"] = self.fitting_method

        # Create the fitting configurer
        configurer = FittingConfigurer(config)

        # Add an entry to the history
        command_name = FittingConfigurer.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        configurer.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the fitting configurer
        configurer.run(default_ranges=default_ranges, settings=self.config.fitting_settings)

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def initialize_fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the fitting ...")

        # Create configuration
        config = dict()

        # Set the name for the fitting run
        config["name"] = self.fitting_run_name

        # Create the fitting initializer
        initializer = GalaxyFittingInitializer(config)

        # Add an entry to the history
        command_name = GalaxyFittingInitializer.command_name()
        self.history.add_entry(command_name)

        # Set the working directory
        initializer.config.path = self.modeling_path

        # Set log path
        log_path = fs.join(self.environment.log_path, command_name + "_" + self.timestamp + ".txt")
        set_log_file(log_path)

        # Run the fitting initializer
        initializer.run()

        # Unset log path
        unset_log_file()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
