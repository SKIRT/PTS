#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.modeler Contains the GalaxyModeler class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ..core.basics.configurable import Configurable
from ..core.tools.logging import log
from ..core.tools import filesystem as fs
from .setup import ModelingSetupTool
from .data.properties import PropertyFetcher
from .data.images import ImageFetcher
from .data.datasetcreator import DataSetCreator
from .data.seds import SEDFetcher
from .preparation.initialization import PreparationInitializer
from .preparation.preparer import DataPreparer
from .decomposition.decomposition import GalaxyDecomposer
from .truncation.truncation import Truncator
from .photometry.photometry import PhotoMeter
from .maps.stars.old import OldStellarMapMaker
from .maps.stars.young import YoungStellarMapMaker
from .maps.stars.ionizing import IonizingStellarMapMaker
from .maps.dust.dust import DustMapMaker
from .fitting.configuration import FittingConfigurer
from .fitting.initialization import FittingInitializer
from .fitting.explorer import ParameterExplorer
from .core.component import load_modeling_history, get_meta_file_path
from ..magic.tools import catalogs
from ..core.basics.range import QuantityRange

# -----------------------------------------------------------------

# Define the different modeling methods
modeling_methods = ["DL14", "FitSKIRT"]
modeling_methods_descriptions = dict()
modeling_methods_descriptions["DL14"] = "High-resolution radiative transfer modeling of face-on galaxies based on De Looze et al. 2014"
modeling_methods_descriptions["FitSKIRT"] = "FitSKIRT method for fitting edge-on galaxies"

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
halpha_fluxes = {"NGC3031": 7.8e40 * Unit("erg/s")}

# -----------------------------------------------------------------

# Set the filters that are not necessary for making maps and therefore shouldn't be included in the preparation steps
# that brings all data to the same resolution and pixelscale
# We want to exclude the SPIRE images from the procedures that bring all images to the same resolution
lower_resolution_filters = ["SPIRE PSW", "SPIRE PMW", "SPIRE PLW"]

# -----------------------------------------------------------------

class GalaxyModeler(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(GalaxyModeler, self).__init__(config)

        # The HYPERLEDA name of the galaxy
        self.hyperleda_name = None

        # The path to the modeling directory
        self.modeling_path = None

        # The modeling history
        self.history = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Get the data
        self.get_data()

        # 3. Data preparation
        self.prepare_data()

        # 4. Decomposition
        self.decompose()

        # 5. Truncation
        self.truncate()

        # 6. Do the photometry
        self.photometry()

        # 7. Make the maps
        self.make_maps()

        # 8. Do the fitting
        self.fit()

        # 9. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(GalaxyModeler, self).setup()

        # Set the hyperleda name
        self.hyperleda_name = catalogs.get_hyperleda_name(self.config.galaxy_name)

        # Set the path to the modeling directory
        self.modeling_path = fs.join(self.config.path, self.config.galaxy_name)

        # Check whether the meta file is present
        if not fs.is_file(get_meta_file_path(self.modeling_path)):

            # Create the configuration
            config = dict()
            config["galaxy_name"] = self.config.galaxy_name

            # Create the modeling setup tool
            tool = ModelingSetupTool(config)

            # Run the setup tool
            tool.run()

        # Load the modeling history
        self.history = load_modeling_history(self.modeling_path)

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
        if "fetch_images" not in self.history: self.get_images()

        # Create the dataset
        if "create_dataset" not in self.history: self.create_dataset()

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

    def get_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the galaxy images ...")

        # Create the configuration
        config = dict()
        config["remote"] = self.config.host_id
        config["halpha_url"] = halpha_urls[self.hyperleda_name]
        config["halpha_flux"] = halpha_fluxes[self.hyperleda_name]

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

    def create_dataset(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dataset ...")

        # Create the dataset creator
        creator = DataSetCreator()

        # Add an entry to the history
        self.history.add_entry(DataSetCreator.command_name())

        # Set the working directory
        creator.config.path = self.modeling_path

        # Run the dataset creator
        creator.run()

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
        self.initialize_preparation()

        # Run the preparation
        self.prepare()

    # -----------------------------------------------------------------

    def initialize_preparation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Initializing the data preparation ...")

        # Create the initializer
        initializer = PreparationInitializer()

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
        self.make_old_stellar_map()

        # Create the map of the young stellar population
        self.make_young_stellar_map()

        # Create the map of the ionizing stellar population
        self.make_ionizing_stellar_map()

        # Create the dust map
        self.make_dust_map()

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

        # Create the dust map maker
        maker = DustMapMaker()

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

    def fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting radiative transfer models to the data ...")

        # Configure the fitting
        self.configure_fit()

        # Initialize the fitting
        self.initialize_fit()

        # Explore the parameter space
        self.explore()

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

        # Set free parameters
        config["parameters"] = free_parameters[self.config.method]
        config["ranges"] = free_parameter_ranges[self.config.method]
        config["filters"] = fitting_filter_names

        # Create the fitting configurer
        configurer = FittingConfigurer(config)

        # Add an entry to the history
        self.history.add_entry(FittingConfigurer.command_name())

        # Set the working directory
        configurer.config.path = self.modeling_path

        # Run the fitting configurer
        configurer.run()

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

        # Create the fitting initializer
        initializer = FittingInitializer()

        # Add an entry to the history
        self.history.add_entry(FittingInitializer.command_name())

        # Set the working directory
        initializer.config.path = self.modeling_path

        # Run the fitting initializer
        initializer.run()

        # Mark the end and save the history file
        self.history.mark_end()
        self.history.save()

    # -----------------------------------------------------------------

    def explore(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Exploring the parameter space ...")

        # Create the parameter explorer
        explorer = ParameterExplorer()

        # Add an entry to the history
        self.history.add_entry(ParameterExplorer.command_name())

        # Set the working directory
        explorer.config.path = self.modeling_path

        # Run the parameter explorer
        explorer.run()

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
