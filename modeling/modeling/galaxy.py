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
from ...core.basics.log import log
from ..data.properties import PropertyFetcher
from ..data.images import ImageFetcher
from ..data.seds import SEDFetcher
from ..data.inspector import DataInspector
from ..preparation.initialization import PreparationInitializer
from ..preparation.preparer import DataPreparer
from ..decomposition.decomposition import GalaxyDecomposer
from ..truncation.truncator import Truncator
from ..truncation.levels import SignificanceLevelsSetter
from ..photometry.photometry import PhotoMeter
from ..maps.old import OldStellarMapMaker
from ..maps.young import YoungStellarMapMaker
from ..maps.ionizing import IonizingStellarMapMaker
from ..maps.dust import DustMapMaker
from ..maps.colour import ColoursMapMaker
from ..maps.attenuation import AttenuationMapMaker
from ..maps.tir import TIRMapMaker
from ..maps.selector import ComponentMapsSelector
from ..maps.components import ComponentMapsMaker
from ..fitting.configuration import FittingConfigurer
from ..fitting.initialization.galaxy import GalaxyFittingInitializer
from ...core.basics.range import QuantityRange
from .base import ModelerBase
from ..config.parameters import units as parameter_units
from ..config.parameters import default_ranges, types, parameter_descriptions
from ...core.units.parsing import parse_unit as u
from ..build.models.galaxy import GalaxyModelBuilder
from ..build.representations.galaxy import GalaxyRepresentationBuilder
from ..build.representations.generator import RepresentationGenerator
from ..component.galaxy import get_galaxy_properties_path, get_data_seds_path, get_data_images_path, get_dustpedia_sed
from ...core.tools import filesystem as fs
from ...core.filter.filter import parse_filter
from ...magic.core.image import Image
from ...magic.core.frame import get_filter
from ..core.environment import GalaxyModelingEnvironment
from ...core.remote.utils import DetachedCalculation
from ...core.tools.utils import UserIntervention
from ..maps.ssfr import SSFRMapMaker
from ...core.tools import types
from ..preparation.inspector import PreparationInspector
from ..component.galaxy import get_observed_sed_file_path, get_reference_seds
from ...core.plot.sed import SEDPlotter
from ...core.tools import parsing
from ...dustpedia.core.database import get_mbb_dust_mass
from ..build.component import get_model_definition
from ...core.units.parsing import parse_quantity

# -----------------------------------------------------------------

# Define the different modeling methods
modeling_methods = dict()
modeling_methods["DL14"] = "High-resolution radiative transfer modeling of face-on galaxies based on De Looze et al. 2014"
modeling_methods["MGE"] = "Panchromatic radiative transfer modeling using Multi-Gaussian expansion"

# -----------------------------------------------------------------

# Set the free parameters for different modeling methods
free_parameters = dict()
free_parameters["DL14"] = ["fuv_young", "dust_mass", "fuv_ionizing"]

# -----------------------------------------------------------------

# Define maximum values for the standard free parameters
max_fuv_young = parse_quantity("1e37 W/micron")
max_fuv_ionizing = parse_quantity("1e35 W/micron")
max_dust_mass = parse_quantity("5.e7 Msun")

# -----------------------------------------------------------------

# Magnitude of fitting ranges
dust_mass_magnitudes = 1
fuv_young_magnitudes = 2
fuv_ionizing_magnitudes = 4

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

# URLS for Halpha data for different galaxies (keys are the HYPERLEDA [DustPedia] names)
halpha_urls = {"NGC3031": "https://ned.ipac.caltech.edu/img/2001ApJ...559..878H/MESSIER_081:I:Ha:hwb2001.fits.gz"}
halpha_fluxes = {"NGC3031": 7.8e40 * u("erg/s")}

# -----------------------------------------------------------------

# SF rates for different galaxies
SFRs = {"NGC3031": 0.8} # The average star formation rate # see Perez-Gonzalez 2006 (mentions Devereux et al 1995)

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
        return self.modeling_config.name

    # -----------------------------------------------------------------

    @property
    def ngc_name(self):
        return self.modeling_config.ngc_name

    # -----------------------------------------------------------------

    @property
    def hyperleda_name(self):
        return self.modeling_config.hyperleda_name

    # -----------------------------------------------------------------

    @property
    def needs_decomposition(self):
        return self.check_needs_step("decompose")

    # -----------------------------------------------------------------

    @property
    def needs_photometry(self):
        return self.check_needs_step("photometry")

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get the data
        self.get_data()

        # 3. Data preparation
        self.prepare_data()

        # 4. Decomposition
        if self.needs_decomposition: self.decompose()

        # 5. Truncation
        self.truncate()

        # 6. Do the photometry
        if self.needs_photometry: self.photometry()

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

        # CHECK IF RERUN IS DEFINED, IF SO, REMOVE COMMANDS FROM THE MODELLING HISTORY
        if self.config.rerun is not None: self.set_rerun()

    # -----------------------------------------------------------------

    @property
    def needs_properties(self):
        return self.check_needs_step("fetch_properties")

    # -----------------------------------------------------------------

    @property
    def needs_seds(self):
        return self.check_needs_step("fetch_seds")

    # -----------------------------------------------------------------

    @property
    def needs_images(self):
        return self.check_needs_step("fetch_images")

    # -----------------------------------------------------------------

    @property
    def needs_data_inspection(self):
        return self.check_needs_step("inspect_data")

    # -----------------------------------------------------------------

    def get_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the galaxy data ...")

        # Get the galaxy properties
        if self.needs_properties:
            if self.properties is not None: self.set_properties()
            else: self.get_properties()

        # Get the galaxy SEDs
        if self.needs_seds:
            if self.seds is not None: self.set_seds()
            else: self.get_seds()

        # Get the galaxy images
        if self.needs_images:
            if self.images is not None: self.set_images()
            else: self.get_images()

        # Inspect the data
        if self.needs_data_inspection: self.inspect_data()

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

        # Set the working directory
        fetcher.config.path = self.modeling_path

        # Run the fetcher
        with self.write_log(fetcher), self.register(fetcher), self.write_config(fetcher): fetcher.run()

    # -----------------------------------------------------------------

    def set_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the provided galaxy properties ...")

        # Add an entry to the history and save
        self.history.add_entry_and_save(PropertyFetcher.command_name())

        # Write the properties to the correct path
        self.properties.saveto(get_galaxy_properties_path(self.modeling_path))

        # Mark the end and save the history file
        self.history.mark_end_and_save()

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

        # Set the working directory
        fetcher.config.path = self.modeling_path

        # Set log path
        with self.write_log(fetcher), self.register(fetcher), self.write_config(fetcher): fetcher.run()

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
        config["nprocesses"] = self.config.nprocesses

        # Create the image fetcher
        fetcher = ImageFetcher(config)

        # Set the working directory
        fetcher.config.path = self.modeling_path

        # Run
        with self.write_log(fetcher), self.register(fetcher), self.write_config(fetcher): fetcher.run()

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

    def inspect_data(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Inspecting the data ...")

        # Create configuration
        config = dict()

        # Create the data inspector
        inspector = DataInspector(config)

        # Set the working directory
        inspector.config.path = self.modeling_path

        # Set log path
        with self.write_log(inspector), self.register(inspector), self.write_config(inspector): inspector.run()

    # -----------------------------------------------------------------

    @property
    def needs_preparation_initialization(self):
        return self.check_needs_step("initialize_preparation")

    # -----------------------------------------------------------------

    @property
    def needs_initialization_inspection(self):
        return self.check_needs_step("inspect_initialization")

    # -----------------------------------------------------------------

    @property
    def needs_preparation(self):
        return self.check_needs_step("prepare_data")

    # -----------------------------------------------------------------

    @property
    def needs_preparation_inspection(self):
        return self.check_needs_step("inspect_preparation")

    # -----------------------------------------------------------------

    def prepare_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the galaxy data ...")

        # Initialize the preparation
        if self.needs_preparation_initialization: self.initialize_preparation()

        # Inspect the initialization
        if self.needs_initialization_inspection: self.inspect_initialization()

        # Run the preparation
        if self.needs_preparation: self.prepare()

        # Inspect the preparation
        if self.needs_preparation_inspection: self.inspect_preparation()

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
        config["manual"] = self.config.sources_manual

        # CACHE DATA/IMAGES
        config["cache"] = self.config.cache

        # Create the initializer
        initializer = PreparationInitializer(config)

        # Set the working directory
        initializer.config.path = self.modeling_path

        # Set log path
        with self.write_log(initializer), self.register(initializer), self.write_config(initializer): initializer.run()

        # Give warning and exit
        message = "Check the result of the source detection, make adjustments where necessary, and resume the modeling afterwards"
        raise UserIntervention(message, self.__class__, "initialize_preparation")

    # -----------------------------------------------------------------

    def inspect_initialization(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Inspecting the initialization ...")

        # Add an entry to the history
        command_name = "inspect_initialization"
        self.history.add_entry_and_save(command_name)

        # Loop over the directories in the data/images path
        images_path = fs.join(self.environment.data_path, "images")
        for image_path in fs.files_in_path(images_path, returns="path", extension="fits", not_contains="poisson", recursive=True):

            # Load frame to get filter
            fltr = get_filter(image_path)

            # Get prep name
            prep_name = str(fltr)

            # Check whether directory present in prep path
            prep_path = fs.join(self.environment.prep_path, prep_name)

            # Preparation directory not found
            if not fs.is_directory(prep_path):
                self.history.remove_entry("initialize_preparation") # ? good idea?
                self.history.save()
                raise RuntimeError("Preparation directory was not found for the " + prep_name + " image. Run initialize_preparation again to solve this.")

            # Determine initialized file path
            initialized_path = fs.join(prep_path, "initialized.fits")

            # Initialized file not found
            if not fs.is_file(initialized_path):
                self.history.remove_entry("initialize_preparation") # ? good idea?
                self.history.save()
                raise RuntimeError("Initialized image was not found for the " + prep_name + " image. Run initialize_preparation again to fix this.")

            # Determine sources directory path
            sources_path = fs.join(prep_path, "sources")

            # Sources directory not found
            if not fs.is_directory(sources_path):
                self.history.remove_entry("initialize_preparation") # ? good idea?
                self.history.save()
                raise RuntimeError("Sources directory was not found for the " + prep_name + " image. Run initialize_preparation again to fix this.")

            # Empty sources directory
            if fs.is_empty(sources_path):
                self.history.remove_entry("initialize_preparation") # ? good idea?
                self.history.save()
                raise RuntimeError("Sources directory for the " + prep_name + " image is empty. Run initialize_preparation again to fix this.")

        # Mark the end and save the history file
        self.history.mark_end_and_save()

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

        config["rerun"] = self.config.rerun_preparation_step

        config["cache"] = self.config.cache

        # Create the data preparer
        preparer = DataPreparer(config)

        # Set the working directory
        preparer.config.path = self.modeling_path

        # Run
        with self.write_log(preparer), self.register(preparer), self.write_config(preparer): preparer.run()

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

        # Set the working directory
        inspector.config.path = self.modeling_path

        # Set log path
        with self.write_log(inspector), self.register(inspector), self.write_config(inspector): inspector.run()

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

        # Set the working directory
        decomposer.config.path = self.modeling_path

        # Set input dict
        input_dict = dict()
        input_dict["disk"] = self.disk_model
        input_dict["bulge"] = self.bulge_model

        # Set log path
        with self.write_log(decomposer), self.register(decomposer), self.write_config(decomposer), self.write_input(decomposer, **input_dict): decomposer.run(**input_dict) #decomposer.run(disk=self.disk_model, bulge=self.bulge_model)

    # -----------------------------------------------------------------

    @property
    def needs_truncation(self):
        return self.check_needs_step("truncate")

    # -----------------------------------------------------------------

    @property
    def needs_significance_levels(self):
        return self.check_needs_step("set_significance_levels")

    # -----------------------------------------------------------------

    def truncate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating ...")

        # Truncate at boundary
        if self.needs_truncation: self.truncate_boundary()

        # Set significance levels
        if self.needs_significance_levels: self.set_significance_levels()

    # -----------------------------------------------------------------

    def truncate_boundary(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting truncation boundary for the galaxy images ...")

        # Create configuration
        config = dict()
        #config["cache"] = self.config.cache

        # CHECK WHETHER TRUNCATION FACTOR IS DEFINED
        if self.config.truncation_factor is None: raise RuntimeError("Truncation factor has to be defined in order for the galaxy modeling to continue")

        # Set the truncation factor
        config["factor"] = self.config.truncation_factor

        # Create the truncator
        truncator = Truncator(config)

        # Set the working directory
        truncator.config.path = self.modeling_path

        # Set log path
        with self.write_log(truncator), self.register(truncator), self.write_config(truncator): truncator.run()

    # -----------------------------------------------------------------

    def set_significance_levels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting significance levels for the galaxy images ...")

        # Create configuration
        config = dict()

        # Check
        if self.config.significance_levels is None: raise RuntimeError("Significance levels have to be defined in order for the galaxy modeling to continue")

        # Set the levels
        config["levels"] = self.config.significance_levels

        # Create the levels setter
        setter = SignificanceLevelsSetter(config)

        # Set the working directory
        setter.config.path = self.modeling_path

        # Run
        with self.write_log(setter), self.register(setter), self.write_config(setter): setter.run()

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

        # Create configuration
        config = dict()
        config["remote"] = self.moderator.host_id_for_single("other")

        # Create the photometer
        photometer = PhotoMeter(config)

        # Set the working directory
        photometer.config.path = self.modeling_path

        # Set log path
        with self.write_log(photometer), self.register(photometer), self.write_config(photometer): photometer.run()

    # -----------------------------------------------------------------

    def get_photometry(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Getting the DustPedia photometry for the galaxy ...")

        # Add an entry to the history
        self.history.add_entry_and_save("photometry")

        # Get the DustPedia SED
        sed = get_dustpedia_sed(self.modeling_path)

        # Create phot path if not yet present
        if not fs.is_directory(self.environment.phot_path): fs.create_directory(self.environment.phot_path)

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
        plotter.run(output=path, title=self.galaxy_name)

        # Mark the end and save the history file
        self.history.mark_end_and_save()

    # -----------------------------------------------------------------

    @property
    def needs_colour_maps(self):
        return self.check_needs_step("make_colours_maps")

    # -----------------------------------------------------------------

    @property
    def needs_ssfr_maps(self):
        return self.check_needs_step("make_ssfr_maps")

    # -----------------------------------------------------------------

    @property
    def needs_tir_maps(self):
        return self.check_needs_step("make_tir_maps")

    # -----------------------------------------------------------------

    @property
    def needs_attenuation_maps(self):
        return self.check_needs_step("make_attenuation_maps")

    # -----------------------------------------------------------------

    @property
    def needs_old_stellar_maps(self):
        return self.check_needs_step("make_old_stellar_maps")

    # -----------------------------------------------------------------

    @property
    def needs_dust_maps(self):
        return self.check_needs_step("make_dust_map")

    # -----------------------------------------------------------------

    @property
    def needs_young_stellar_maps(self):
        return self.check_needs_step("make_young_stellar_maps")

    # -----------------------------------------------------------------

    @property
    def needs_ionizing_stellar_maps(self):
        return self.check_needs_step("make_ionizing_stellar_maps")

    # -----------------------------------------------------------------

    @property
    def needs_maps_selection(self):
        return self.check_needs_step("select_component_maps")

    # -----------------------------------------------------------------

    @property
    def needs_component_maps(self):
        return self.check_needs_step("make_component_maps")

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the maps describing the model geometries ...")

        # Create colour maps
        if self.needs_colour_maps: self.make_colour_maps()

        # Create sSFR maps
        if self.needs_ssfr_maps: self.make_ssfr_maps()

        # Create the TIR map
        if self.needs_tir_maps: self.make_tir_map()

        # Create the attenuation map(s)
        if self.needs_attenuation_maps: self.make_attenuation_maps()

        # Create the map of the old stellar disk
        if self.needs_old_stellar_maps: self.make_old_stellar_maps()

        # Create the dust map
        if self.needs_dust_maps: self.make_dust_maps()

        # Create the map of the young stellar population
        if self.needs_young_stellar_maps: self.make_young_stellar_maps()

        # Create the map of the ionizing stellar population
        if self.needs_ionizing_stellar_maps: self.make_ionizing_stellar_maps()

        # Make the component maps selection
        if self.needs_maps_selection: self.select_component_maps()

        # Create the maps for the different RT model components
        if self.needs_component_maps: self.make_component_maps()

    # -----------------------------------------------------------------

    def make_colour_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the colour maps ...")

        # Create the colour maps
        maker = ColoursMapMaker()

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        with self.write_log(maker), self.register(maker), self.write_config(maker): maker.run()

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

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        with self.write_log(maker), self.register(maker), self.write_config(maker): maker.run()

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

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        with self.write_log(maker), self.register(maker), self.write_config(maker): maker.run()

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

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        with self.write_log(maker), self.register(maker), self.write_config(maker): maker.run()

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

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        with self.write_log(maker), self.register(maker), self.write_config(maker): maker.run()

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

        # make_black_body
        # make_emission
        config["make_black_body"] = False
        config["make_emission"] = False

        # Create the dust map maker
        maker = DustMapMaker(config)

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        with self.write_log(maker), self.register(maker), self.write_config(maker): maker.run()

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

        # Set the working directory
        maker.config.path = self.modeling_path

        # Set log path
        with self.write_log(maker), self.register(maker), self.write_config(maker): maker.run()

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

        # Set the working directory
        maker.config.path = self.modeling_path

        # Run the ionizing stellar map maker
        with self.write_log(maker), self.register(maker), self.write_config(maker): maker.run()

    # -----------------------------------------------------------------

    def select_component_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Selecting the component maps ...")

        # Create configuration
        config = dict()
        config["auto"] = True  # automatic selection

        # Create the maps selector
        selector = ComponentMapsSelector(config)

        # Set the working directory
        selector.config.path = self.modeling_path

        # Run
        with self.write_log(selector), self.register(selector), self.write_config(selector): selector.run()

    # -----------------------------------------------------------------

    def make_component_maps(self):

        """
        Thins function ...
        :return:
        """

        # Inform the user
        log.info("Making the component maps ...")

        # Create the configuration
        config = dict()
        config["remote"] = self.moderator.host_id_for_single("other")
        #config["attached"] = self.config.attached

        # Create the map maker
        maker = ComponentMapsMaker(config)

        # Set the working directory
        maker.config.path = self.modeling_path

        # Run
        with self.write_log(maker), self.register(maker), self.write_config(maker): maker.run()

    # -----------------------------------------------------------------

    @property
    def needs_build_model(self):
        if self.is_galaxy_modeler: return self.check_needs_step("build_model_galaxy")
        elif self.is_sed_modeler: return self.check_needs_step("build_model_sed")
        elif self.is_images_modeler: return self.check_needs_step("build_model_images")
        else: raise RuntimeError("Something horribly went wrong")

    # -----------------------------------------------------------------

    @property
    def needs_representations(self):
        return self.check_needs_step("generate_representations")

    # -----------------------------------------------------------------

    def build(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the model and its representations ...")

        # 1. Build model
        if self.needs_build_model: self.build_model()

        # Generate the representations
        if self.needs_representations: self.generate_representations()

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

        # Get the estimated SFR
        config["sfr"] = SFRs[self.hyperleda_name] if self.hyperleda_name in SFRs else 1.0

        # Get the dust mass estimated by black body fitting (DustPedia)
        config["dust_mass"] = get_mbb_dust_mass(self.hyperleda_name)

        # Create the builder
        builder = GalaxyModelBuilder(config)

        # Set the working directory
        builder.config.path = self.modeling_path

        # Run the model builder
        with self.write_log(builder), self.register(builder), self.write_config(builder): builder.run()

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
        builder = GalaxyRepresentationBuilder(config)

        # Run the builder
        with self.write_log(builder), self.register(builder), self.write_config(builder): builder.run()

    # -----------------------------------------------------------------

    def generate_representations(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Generating the representations ...")

        # Create configuration
        config = dict()

        # Set model name
        config["model_name"] = self.model_name

        # Set number of representations
        config["nrepresentations"] = self.config.ndust_grids

        # Create the builder
        generator = RepresentationGenerator(config)

        # Run the generator
        with self.write_log(generator), self.register(generator), self.write_config(generator): generator.run()

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
        parameter_labels = free_parameters[self.modeling_config.method]
        config["parameters"] = parameter_labels

        # Set parameter descriptions
        config["descriptions"] = parameter_descriptions

        # Set parameter types
        config["types"].name = types

        # Set parameter units
        config["units"] = parameter_units

        # Calculate ranges
        free_parameter_ranges = dict()

        # Get the model definition
        definition = get_model_definition(self.modeling_path, self.model_name)

        # Get FUV young, FUV ionizing and dust mass
        fuv_young = definition.young_stars_luminosity
        fuv_ionizing = definition.ionizing_stars_luminosity
        dust_mass = definition.dust_mass

        # Define ranges
        for label in parameter_labels:
            if label == "fuv_young":
                #young_range = QuantityRange.around_magnitude(fuv_young, 3)
                #young_range.max = max_fuv_young
                young_range = QuantityRange.within_magnitude(fuv_young, fuv_young_magnitudes)
                free_parameter_ranges[label] = young_range
            elif label == "fuv_ionizing":
                #ionizing_range = QuantityRange.around_magnitude(fuv_ionizing, 3)
                #ionizing_range.max = max_fuv_ionizing
                ionizing_range = QuantityRange.within_magnitude(fuv_ionizing, fuv_ionizing_magnitudes)
                free_parameter_ranges[label] = ionizing_range
            elif label == "dust_mass":
                #mass_range = QuantityRange.around_magnitude(dust_mass, 2)
                #mass_range.max = max_dust_mass
                mass_range = QuantityRange.within_magnitude(dust_mass, dust_mass_magnitudes)
                free_parameter_ranges[label] = mass_range
            else: raise ValueError("Free parameter label not recognized: '" + label + "'")

        # Set ranges
        config["ranges"] = free_parameter_ranges  #free_parameter_ranges[self.modeling_config.method]

        # Set fitting filters
        config["filters"] = fitting_filter_names

        # Ask for and set the fitting method
        config["fitting_method"] = self.fitting_method

        # Create the fitting configurer
        configurer = FittingConfigurer(config)

        # Set the working directory
        configurer.config.path = self.modeling_path

        # Set the input dict
        input_dict = dict()
        input_dict["default_ranges"] = default_ranges
        input_dict["settings"] = self.config.fitting_settings

        # Run the fitting configurer
        with self.write_log(configurer), self.register(configurer), self.write_config(configurer), self.write_input(configurer, **input_dict): configurer.run(**input_dict)

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

        # Set the number of wavelength grids
        config["wg"] = dict()
        config["wg"]["ngrids"] = self.config.nwavelength_grids

        # Create the fitting initializer
        initializer = GalaxyFittingInitializer(config)

        # Set the working directory
        initializer.config.path = self.modeling_path

        # Run the fitting initializer
        with self.write_log(initializer), self.register(initializer), self.write_config(initializer): initializer.run()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
