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
from ..maps.stars.old import OldStellarMapMaker
from ..maps.stars.young import YoungStellarMapMaker
from ..maps.stars.ionizing import IonizingStellarMapMaker
from ..maps.dust.dust import DustMapMaker
from ..fitting.configuration import FittingConfigurer
from ..fitting.initialization.galaxy import GalaxyFittingInitializer
from ...core.basics.range import QuantityRange
from .base import ModelerBase
from ..config.parameters import units as parameter_units
from ..config.parameters import default_ranges, types, parameter_descriptions
from ...core.basics.unit import parse_unit as u
from ..build.model import ModelBuilder
from ..build.representation import RepresentationBuilder

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

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(GalaxyModeler, self).__init__(config, interactive)

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
        self.fit()

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

        # Check whether a remote is available for the heavy computations
        if self.moderator.host_id_for_single("other") is None: raise RuntimeError("The desired remote(s) for heavy computations are currently unavailable")

    # -----------------------------------------------------------------

    def get_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the galaxy data ...")

        # Get the galaxy properties
        if "fetch_properties" not in self.history: self.get_properties()

        # Get the galaxy SEDs
        if "fetch_seds" not in self.history: self.get_seds()

        # Get the galaxy images
        if "fetch_images" not in self.history: self.get_images_and_exit()

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

    def get_images_and_exit(self):

        """
        This function ...
        :return:
        """

        # Get the images
        self.get_images()

        # If not running in attached remote mode, exit now
        if not self.config.attached:
            log.warning("The procedure that calculates the Poisson error maps for GALEX and SDSS is now running. "
                        "Wait for it to finished and resume the modeling afterwards")
            exit()

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
        config["halpha_url"] = halpha_urls[self.hyperleda_name]
        config["halpha_flux"] = halpha_fluxes[self.hyperleda_name]
        config["max_nobservations_mosaic"] = self.config.max_nobservations_mosaic

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

    # -----------------------------------------------------------------

    def prepare_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the galaxy data ...")

        # Initialize the preparation
        if "initialize_preparation" not in self.history: self.initialize_preparation_and_exit()

        # Run the preparation
        if "prepare_data" not in self.history: self.prepare()

    # -----------------------------------------------------------------

    def initialize_preparation_and_exit(self):

        """
        This function ...
        :return:
        """

        self.initialize_preparation()

        # Give warning and exit
        log.warning("Check the result of the source detection, make adjustments where necessary, and resume the modeling afterwards")
        exit()

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

        # Create the initializer
        initializer = PreparationInitializer(config)

        # Add an entry to the history
        self.history.add_entry(PreparationInitializer.command_name())

        # Set the working directory
        initializer.config.path = self.modeling_path

        # Run the initializer
        initializer.run()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

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

        # Create the data preparer
        preparer = DataPreparer()

        # Add an entry to the history
        self.history.add_entry(DataPreparer.command_name())

        # Set the working directory
        preparer.config.path = self.modeling_path

        # Run the preparer
        preparer.run()

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
        self.history.add_entry(GalaxyDecomposer.command_name())

        # Set the working directory
        decomposer.config.path = self.modeling_path

        # Run the decomposer
        decomposer.run()

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
        self.history.add_entry(Truncator.command_name())

        # Set the working directory
        truncator.config.path = self.modeling_path

        # Run the truncator
        truncator.run()

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
        log.info("Running photometry on the galaxy images ...")

        # Create the photometer
        photometer = PhotoMeter()

        # Add an entry to the history
        self.history.add_entry(PhotoMeter.command_name())

        # Set the working directory
        photometer.config.path = self.modeling_path

        # Run the photometer
        photometer.run()

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

        # Create the map of the old stellar disk
        if "make_old_map" not in self.history: self.make_old_stellar_map()

        # Create the map of the young stellar population
        if "make_young_map" not in self.history: self.make_young_stellar_map()

        # Create the map of the ionizing stellar population
        if "make_ionizing_map" not in self.history: self.make_ionizing_stellar_map()

        # Create the dust map
        if "make_dust_map" not in self.history: self.make_dust_map()

    # -----------------------------------------------------------------

    def make_old_stellar_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of old stars ...")

        # Create the old stellar map maker
        maker = OldStellarMapMaker()

        # Add an entry to the history
        self.history.add_entry(maker.command_name())

        # Set the working directory
        maker.config.path = self.modeling_path

        # Run the maker
        maker.run()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def make_young_stellar_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of young stars ...")

        # Create the young stellar map maker
        maker = YoungStellarMapMaker()

        # Add an entry to the history
        self.history.add_entry(YoungStellarMapMaker.command_name())

        # Set the working directory
        maker.config.path = self.modeling_path

        # Run the maker
        maker.run()

        # Mark the end and save the history
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def make_ionizing_stellar_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of ionizing stars ...")

        # Create the ionizing stellar map maker
        maker = IonizingStellarMapMaker()

        # Add an entry to the history
        self.history.add_entry(IonizingStellarMapMaker.command_name())

        # Set the working directory
        maker.config.path = self.modeling_path

        # Run the ionizing stellar map maker
        maker.run()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def make_dust_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the dust map ...")

        # Create the configuration
        config = dict()
        config["black_body"] = dict()
        config["black_body"]["remote"] = self.moderator.host_id_for_single("other")

        # Create the dust map maker
        maker = DustMapMaker(config)

        # Add an entry to the history
        self.history.add_entry(DustMapMaker.command_name())

        # Set the working directory
        maker.config.path = self.modeling_path

        # Run the dust map maker
        maker.run()

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
        self.history.add_entry(ModelBuilder.command_name())

        # Set the working directory
        builder.config.path = self.modeling_path

        # Run the model builder
        builder.run()

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
        self.history.add_entry(RepresentationBuilder.command_name())

        # Run the builder
        builder.run()

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

        # Create the fitting configurer
        configurer = FittingConfigurer(config)

        # Add an entry to the history
        self.history.add_entry(FittingConfigurer.command_name())

        # Set the working directory
        configurer.config.path = self.modeling_path

        # Run the fitting configurer
        configurer.run(default_ranges=default_ranges, settings=self.config.fitting_settings)

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
        self.history.add_entry(GalaxyFittingInitializer.command_name())

        # Set the working directory
        initializer.config.path = self.modeling_path

        # Run the fitting initializer
        initializer.run()

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
