#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.modeler Contains the GalaxyModeler class, which runs the radiative transfer modelling procedure
#  for a certain galaxy.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
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
from ..fitting.initialization import FittingInitializer
from ..fitting.explorer import ParameterExplorer
from ..fitting.sedfitting import SEDFitter
from ..core.component import load_modeling_history, get_config_file_path, load_modeling_configuration
from ...core.basics.range import QuantityRange
from ..fitting.component import get_generations_table
from ...core.launch.synchronizer import RemoteSynchronizer
from ...core.remote.remote import is_available
from ...core.prep.deploy import Deployer

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
halpha_fluxes = {"NGC3031": 7.8e40 * Unit("erg/s")}

# -----------------------------------------------------------------

# Set the filters that are not necessary for making maps and therefore shouldn't be included in the preparation steps
# that brings all data to the same resolution and pixelscale
# We want to exclude the SPIRE images from the procedures that bring all images to the same resolution
lower_resolution_filters = ["SPIRE PSW", "SPIRE PMW", "SPIRE PLW", "MIPS 70mu", "MIPS 160mu", "WISE W4"] # WISE W4 because it has just a slightly higher FWHM than Pacs 160mu

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

        # The path to the modeling directory
        self.modeling_path = None

        # The modeling configuration
        self.modeling_config = None

        # Host ids of the available hosts
        self.available_host_ids = set()

        # The modeling history
        self.history = None

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

    @property
    def used_host_ids(self):

        """
        This function ...
        :return:
        """

        host_ids = set()

        # Add main host ID
        host_ids.add(self.host_id)

        # Add fitting host ids, if they are available
        for host_id in self.modeling_config.fitting_host_ids:
            if host_id in self.available_host_ids: host_ids.add(host_id)

        # Return the list of host IDs
        return list(host_ids)

    # -----------------------------------------------------------------

    @property
    def host_id(self):

        """
        This function ...
        :return:
        """

        # Loop over the preferred hosts
        for host_id in self.modeling_config.host_ids:
            if host_id in self.available_host_ids: return host_id

        # No host avilable
        return None

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
        if "decompose" not in self.history: self.decompose()

        # 5. Truncation
        if "truncate" not in self.history: self.truncate()

        # 6. Do the photometry
        if "photometry" not in self.history: self.photometry()

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

        # Set the path to the modeling directory
        self.modeling_path = self.config.path

        # Check for the presence of the configuration file
        if not fs.is_file(get_config_file_path(self.modeling_path)): raise ValueError("The current working directory is not a radiative transfer modeling directory (the configuration file is missing)")
        else: self.modeling_config = load_modeling_configuration(self.modeling_path)

        # Find which hosts are available
        if self.config.check_hosts: self.find_available_hosts()

        # Load the modeling history
        self.history = load_modeling_history(self.modeling_path)

        # Deploy SKIRT and PTS
        if self.config.deploy: self.deploy()

    # -----------------------------------------------------------------

    def find_available_hosts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding available hosts ...")

        # Find available hosts from host_ids list
        for host_id in self.modeling_config.host_ids:
            if is_available(host_id):
                log.debug("Host '" + host_id + "' is available")
                self.available_host_ids.add(host_id)
            else: log.debug("Host '" + host_id + "' is not available")

        # Find available hosts from fitting.host_ids list
        for host_id in self.modeling_config.fitting_host_ids:
            if host_id in self.modeling_config.host_ids: continue
            if is_available(host_id):
                log.debug("Host '" + host_id + "' is available")
                self.available_host_ids.add(host_id)
            else: log.debug("Host '" + host_id + "' is not available")

        # No available host in the list of preferred host ids
        if len(self.available_host_ids) == 0: raise RuntimeError("None of the preferred hosts are available at this moment")

    # -----------------------------------------------------------------

    def deploy(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deploying SKIRT and PTS ...")

        # Create the deployer
        deployer = Deployer()

        # Set the host ids
        deployer.config.host_ids = self.used_host_ids

        # Set the host id on which PTS should be installed
        deployer.config.pts_on = [self.host_id]

        # Set
        deployer.config.check = self.config.check_versions

        # Run the deployer
        deployer.run()

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
        if "fetch_images" not in self.history:

            self.get_images()
            log.warning("The procedure that calculates the Poisson error maps for GALEX and SDSS is now running. "
                        "Wait for it to finished and resume the modeling afterwards")
            exit()

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
        config["remote"] = self.host_id
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

    def prepare_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the galaxy data ...")

        # Initialize the preparation
        if "initialize_preparation" not in self.history:

            self.initialize_preparation()
            log.warning("Check the result of the source detection, make adjustments where necessary, and resume the modeling afterwards")
            exit()

        # Run the preparation
        if "prepare_data" not in self.history: self.prepare()

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
        config["remote"] = self.host_id

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
        config["remote"] = self.host_id
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
        config["black_body"]["remote"] = self.host_id

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

    def fit(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting radiative transfer models to the data ...")

        # Configure the fitting
        if "configure_fit" not in self.history: self.configure_fit()

        # Initialize the fitting
        if "initialize_fit" not in self.history: self.initialize_fit()

        # Load the generations table
        generations = get_generations_table(self.modeling_path)

        # If some generations have not finished, check the status of and retrieve simulations
        if generations.has_unfinished: self.synchronize()

        # If some generations have finished, fit the SED
        if generations.has_finished: self.fit_sed()

        # IF all generations have finished, explore new generation of models
        if generations.all_finished: self.explore()

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
        config["parameters"] = free_parameters[self.modeling_config.method]
        config["ranges"] = free_parameter_ranges[self.modeling_config.method]
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

    def synchronize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Synchronizing with the remotes (retrieving and analysing finished models) ...")

        # Create the remote synchronizer
        synchronizer = RemoteSynchronizer()

        # Set the host IDs
        synchronizer.config.host_ids = self.modeling_config.fitting_host_ids

        # Run the remote synchronizer
        synchronizer.run()

    # -----------------------------------------------------------------

    def fit_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fitting the SED to the finished generations ...")

        # Create the SED fitter
        fitter = SEDFitter()

        # Add an entry to the history
        self.history.add_entry(SEDFitter.command_name())

        # Run the fitter
        fitter.run()

        # Mark the end and save the history file
        self.history.mark_end()

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
